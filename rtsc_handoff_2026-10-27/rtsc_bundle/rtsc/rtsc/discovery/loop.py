"""The discovery round loop (DISCOVERY.md section 2).

Each round has phases, each committed atomically to the ledger together with its results:

    proposed     the policy's generation requests (deterministic from seed + ledger)
    generated    generator outputs -> validation -> fingerprint dedupe against the WHOLE ledger
                 -> new materials (DEFERRED per goal, awaiting selection) or duplicate links
    prescreened  heuristic ML pre-screen (S2: DEFERRED with a recorded reason, or ordering only)
    selected     the top candidates within the per-round DFT budget get a dft_jobs row
    executed     the existing pipeline per goal (rtsc.pipeline.run: DFT existence filters,
                 DFPT, certify.classify, validity flags) via the executor; results -> ledger
    done         policy snapshot, reports

A crash loses at most the phase in flight: ``continue`` resumes the unfinished round at its phase;
the pipeline's own resume (record reuse only when material id, config hash AND input digest match)
means finished materials are never recomputed, and the ledger never appends the same pipeline
verdict twice.

Every material is evaluated for every enabled goal whose pressure range contains its pressure:
one DFT stage (in the pipeline config of its first in-scope goal), then ``certify.classify`` at each
further goal's T* on the SAME alternative spectra (``pipeline.compute_material`` on a record seeded
with those spectra and DFT stages: the pipeline's own certify-only resume path).  Materials outside
a goal's pressure range are DEFERRED (stage "scope", in_scope = 0) for that goal and can never be
certified for it.
"""
from __future__ import annotations

import copy
import hashlib
import json
import logging
import math
import time
from pathlib import Path
from typing import Any

from .. import pipeline
from ..report import write_reports as _pipeline_reports
from ..sources.base import Candidate, candidate_from_structure_arrays, safe_id
from ..sources.synthetic import gaussian_peak_spectrum
from ..status import PROVEN_STATUSES, Status, Verdict
from .executors import make_executor
from .fingerprint import (InvalidStructure, StructureMatcherBackend, coarse_key, composition_of,
                          fingerprint, invalid_fingerprint, reduced_formula, same_structure,
                          validate_structure)
from .generators import make_generator, make_prescreener
from .goals import DiscoveryConfig, Goal, counts_for
from .interfaces import GeneratedStructure, GenerationRequest, PrescreenResult
from .ledger import Ledger
from .policy import ArmStats, BanditPolicy, build_arms, reward

log = logging.getLogger("rtsc.discovery")

EXECUTOR_ORIGIN = "discovery.executor"
STOP_FILE = "STOP"


def _seed(*parts: Any) -> int:
    return int(hashlib.sha256(":".join(map(str, parts)).encode()).hexdigest()[:8], 16)


def _finite(x: Any) -> float | None:
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


def synthetic_spectrum(fp: str, chemsys: tuple[str, ...], pressure_GPa: float,
                       scfg: dict[str, Any]) -> dict[str, Any]:
    """Deterministic SYNTHETIC alpha^2F for demo runs (invented; NOT data for any material).
    Strength depends on the chemical system (so the bandit has something to learn) and on the
    fingerprint (variety within a system)."""
    import random
    a_lo, a_hi = scfg.get("A_range_K", [150.0, 3200.0])
    c_lo, c_hi = scfg.get("center_range_K", [600.0, 2200.0])
    chem = random.Random("chem:" + "-".join(chemsys)).uniform(0.25, 1.0)
    r = random.Random("fp:" + fp)
    A = a_lo + (a_hi - a_lo) * chem * r.uniform(0.5, 1.0)
    center = r.uniform(c_lo, c_hi)
    d = gaussian_peak_spectrum(center, center * float(scfg.get("width_frac", 0.06)), A,
                               label=f"SYNTHETIC discovery-demo {fp[:14]}")
    d["provenance"]["purpose"] = ("SYNTHETIC spectrum attached to a hypothetical structure by the "
                                  "discovery demo; invented, not computed, not data for this formula")
    d["provenance"]["pressure_GPa"] = pressure_GPa
    return d


class DiscoveryLoop:
    def __init__(self, dc: DiscoveryConfig, *, generators: dict[str, Any] | None = None,
                 prescreener: Any = None, executor: Any = None):
        self.dc = dc
        self.d = dc.d
        self.workdir = dc.workdir
        self.workdir.mkdir(parents=True, exist_ok=True)
        self.ledger = Ledger(self.workdir / "ledger.sqlite")
        self.generators = generators if generators is not None else {
            n: make_generator(n, g) for n, g in sorted(self.d["generators"].items())}
        self.prescreener = prescreener if prescreener is not None else make_prescreener(self.d["prescreen"])
        self.executor = executor if executor is not None else make_executor(self.d)
        self.goals: list[Goal] = dc.enabled_goals
        self.arms = build_arms(self.goals, sorted(self.generators))
        s = self.d["search"]
        self.policy = BanditPolicy(self.arms, int(self.d["seed"]), s["policy"],
                                   float(s["exploration_floor"]), float(s["ucb_c"]))
        dd = self.d["dedupe"]
        self.confirm_tol = float(dd.get("confirm_tol_A", 0.02))
        self.matcher = StructureMatcherBackend(dd["ltol"], dd["stol"], dd["angle_tol"],
                                               bool(dd["use_structure_matcher"]))
        self.dft_hash = dc.dft_hash()
        self.synthetic = bool(self.d["synthetic_spectra"].get("enabled", False))
        stored = self.ledger.get_meta("config")
        if stored is None:
            self.ledger.set_meta("config", dc.to_json())
            self.ledger.set_meta("goal_ids", [g.id for g in self.goals])
        else:
            old_goals = self.ledger.get_meta("goal_ids", [])
            missing = [g for g in old_goals if g not in [x.id for x in self.goals]]
            if missing:
                raise pipeline.ConfigError(f"goals {missing} exist in the ledger but are not enabled now")
            new = [g for g in self.goals if g.id not in old_goals]
            if new:
                self._add_goal_statuses(new)
                self.ledger.set_meta("goal_ids", old_goals + [g.id for g in new])

    # ================================================================== public API
    def run(self, command: str = "run", **overrides: Any) -> dict[str, Any]:
        limits = dict(self.d["stop"])
        limits.update({k: v for k, v in overrides.items() if v is not None})
        sid = self.ledger.start_session(command, limits)
        self.ledger.set_meta("config", self.dc.to_json())
        self.ledger.set_meta(f"config_session_{sid}", self.dc.to_json())
        self._reconcile(sid)
        t0 = time.monotonic()
        cert0 = {g.id: self.certified_ids(g) for g in self.goals}
        jobs0 = len(self.ledger.jobs())
        rounds_done = 0
        reason = None
        self.write_reports()
        while True:
            last = self.ledger.last_round()
            in_flight = last is not None and last["phase"] != "done"
            new_cert = {g.id: len(self.certified_ids(g) - cert0[g.id]) for g in self.goals}
            dft_used = len(self.ledger.jobs()) - jobs0
            reason = self._stop_reason(limits, rounds_done, t0, new_cert, dft_used, in_flight)
            if reason:
                break
            k = last["round"] if in_flight else (0 if last is None else last["round"]) + 1
            if not in_flight:
                self.ledger.start_round(k, sid, self._propose(k))
            budget_left = None if limits.get("max_dft_jobs") is None else \
                max(0, int(limits["max_dft_jobs"]) - dft_used)
            self._run_round(k, budget_left)
            rounds_done += 1
            r = self.ledger.round_row(k)
            c = json.loads(r["counts"])
            if c.get("new", 0) == 0 and c.get("selected", 0) == 0 and not self._pool():
                reason = ("exhausted: the last round produced no new structure and no material is "
                          "waiting for DFT (every in-scope goal of every computed material is evaluated)")
                break
        self.write_reports()
        self.ledger.finish_session(sid, reason)
        return {"stop_reason": reason, "rounds_this_session": rounds_done,
                "new_certified": {g.id: sorted(self.certified_ids(g) - cert0[g.id]) for g in self.goals},
                "certified_total": {g.id: len(self.certified_ids(g)) for g in self.goals},
                "accounting": self.ledger.accounting([g.id for g in self.goals])}

    def counts(self, goal: Goal, row: dict) -> tuple[bool, str]:
        """``goals.counts_for`` with the CURRENT config hashes (the single counting predicate)."""
        return counts_for(goal, row, certify_hash=self.dc.certify_hash(goal), dft_hash=self.dft_hash)

    def certified_ids(self, goal: Goal) -> set[str]:
        """Materials that satisfy ``goal`` (``counts_for`` on their current row)."""
        return {r["material_id"] for r in self.ledger.current(goal.id) if self.counts(goal, r)[0]}

    def _final_under_current_config(self, cur: dict | None, g: Goal) -> bool:
        """Reconcile leaves a (material, goal) alone iff its current row is final under the
        CURRENT hashes: a pipeline verdict (classify-origin ones also need the current
        certify_hash), or a non-retryable executor ERROR (e.g. DFT_PRESSURE_MISMATCH) made under
        the current dft_hash and certify_hash.  This makes reconcile idempotent."""
        if cur is None:
            return False
        if cur["stage"] == "pipeline":
            return cur["origin"] != pipeline.CLASSIFY_ORIGIN or cur["certify_hash"] == self.dc.certify_hash(g)
        if cur["stage"] == "executor" and cur["status"] == Status.ERROR.value and not cur["retryable"]:
            return cur["dft_hash"] == self.dft_hash and cur["certify_hash"] == self.dc.certify_hash(g)
        return False

    # ================================================================== config reconciliation
    def _reconcile(self, sid: int) -> None:
        """Bring the ledger in line with the CURRENT config before any new round (audit F4-F6):
        scope changes -> scope / awaiting statuses; spectra computed under another dft_hash ->
        DEFERRED (stage config) and re-queued for DFT; an in-scope goal of a computed material
        whose current row is not a pipeline verdict under the current certify_hash (new goal,
        changed T* or [certify]) -> certify.classify only, on the stored spectra."""
        last = self.ledger.last_round()
        if last is None:
            return
        k = int(last["round"])
        done = {j["material_id"] for j in self.ledger.jobs(None, ("DONE",))}
        touched: dict[str, tuple[dict, list[str]]] = {}
        for m in self.ledger.materials():
            if not m["valid"] or m["dup_of"] is not None:
                continue
            mid, P = m["material_id"], float(m["pressure_GPa"])
            with self.ledger.tx():
                for g in self.goals:
                    cur = self.ledger.current_one(mid, g.id)
                    if not g.in_scope(P) and (cur is None or cur["stage"] != "scope"):
                        self._initial_status_for(g, mid, k, P, "")
                    elif g.in_scope(P) and (cur is None or cur["stage"] == "scope"):
                        self._initial_status_for(g, mid, k, P, f"goal {g.id} now covers {P:g} GPa "
                                                 "(config changed); awaiting evaluation")
            if mid not in done:
                continue
            pipe = [self.ledger.current_one(mid, g.id) for g in self.goals]
            pipe = [r for r in pipe if r and r["stage"] == "pipeline"]
            if any(r["dft_hash"] != self.dft_hash for r in pipe):
                with self.ledger.tx():
                    for g in self.goals:
                        if g.in_scope(P):
                            self.ledger.append_verdict(
                                material_id=mid, goal_id=g.id, round_=k, stage="config",
                                origin="discovery.config", status=Status.DEFERRED, T_star_K=g.T_star_K,
                                pressure_GPa=P, flags=["HYPOTHETICAL_STRUCTURE"],
                                reasons=["DFT / spectra / filter configuration changed since this "
                                         "material was computed; re-queued for the DFT stage"])
                    for j in self.ledger.jobs(None, ("DONE",)):
                        if j["material_id"] == mid:
                            self.ledger.set_job_state(mid, j["round"], "STALE")
                continue
            src = next((r for r in pipe if r["origin"] == pipeline.CLASSIFY_ORIGIN), None) or \
                next((r for r in pipe if not r["retryable"]), None)
            if src is None or not src["record_path"]:
                continue
            todo = []
            for g in self.goals:
                if not g.in_scope(P):
                    continue
                if self._final_under_current_config(self.ledger.current_one(mid, g.id), g):
                    continue
                todo.append(g)
            if not todo:
                continue
            sp = Path(src["record_path"])
            src_rec = pipeline.read_record(sp.parent.parent, mid)
            c = self.candidate(mid)
            for g in todo:
                cfg = self.dc.goal_pipeline_config(g, self.workdir / "reconcile" / f"s{sid:03d}" / g.id)
                touched.setdefault(g.id, (cfg, []))[1].append(mid)
                rec, copied = self._other_goal_record(mid, c, src_rec, cfg, sp)
                with self.ledger.tx():
                    self._record(k, mid, g, rec, cfg, c, src_record_path=str(sp) if copied else None,
                                 note=f"re-evaluated at session {sid} start: goal added, T* or "
                                      "[certify] changed (certify.classify on the stored spectra)")
        for cfg, mids in touched.values():
            self._rewrite_goal_reports(cfg, mids)

    def write_reports(self) -> None:
        from .report import write_discovery_reports
        write_discovery_reports(self.ledger, self.dc, self.goals, self.workdir)

    # ================================================================== stop criteria
    def _stop_reason(self, limits, rounds_done, t0, new_cert, dft_used, in_flight) -> str | None:
        if (self.workdir / STOP_FILE).exists():
            return f"STOP file present ({self.workdir / STOP_FILE})"
        n = limits.get("certified_per_goal")
        if n is not None and all(new_cert[g.id] >= int(n) for g in self.goals):
            return f"certified_per_goal={n} reached for every enabled goal ({new_cert})"
        if limits.get("max_rounds") is not None and rounds_done >= int(limits["max_rounds"]):
            return f"max_rounds={limits['max_rounds']} reached"
        if limits.get("max_wall_hours") is not None and \
                time.monotonic() - t0 >= float(limits["max_wall_hours"]) * 3600:
            return f"max_wall_hours={limits['max_wall_hours']} reached"
        if limits.get("max_dft_jobs") is not None and dft_used >= int(limits["max_dft_jobs"]) \
                and not in_flight:
            return f"max_dft_jobs={limits['max_dft_jobs']} reached"
        return None

    # ================================================================== policy
    def arm_stats(self) -> dict[str, ArmStats]:
        st = {a.key: ArmStats() for a in self.arms}
        for r in self.ledger.rounds():
            per_new = r["counts"].get("per_request_new")
            for i, rq in enumerate(r["requests"]):
                s = st.setdefault(rq["arm"], ArmStats())
                s.pulls += 1
                if per_new is not None and i < len(per_new) and per_new[i] == 0:
                    s.F += 0.5
        goal_T = {g.id: g.T_star_K for g in self.goals}
        for g in self.goals:
            for r in self.ledger.current(g.id):
                arm_goal = r["arm"].split("|", 1)[0]
                if arm_goal != g.id or r["stage"] != "pipeline":
                    continue
                if r["status"] == Status.ERROR.value and r["retryable"]:
                    continue
                x = reward(r["status"], r["T_low_K"], r["T_up_K"], goal_T[g.id])
                s = st.setdefault(r["arm"], ArmStats())
                s.S += x
                s.F += 1.0 - x
                s.outcomes += 1
        for r in self.ledger.q("SELECT m.arm FROM ml_screen s JOIN materials m "
                               "ON m.material_id=s.material_id WHERE s.ml_deferred=1 ORDER BY m.seq"):
            st.setdefault(r["arm"], ArmStats()).F += 0.5
        return st

    def _propose(self, k: int) -> list[dict[str, Any]]:
        stats = self.arm_stats()
        arms = self.policy.propose(k, int(self.d["search"]["requests_per_round"]), stats)
        seen: dict[str, int] = {}
        out = []
        for i, a in enumerate(arms):
            pulls = stats.get(a.key, ArmStats()).pulls + seen.get(a.key, 0)
            seen[a.key] = seen.get(a.key, 0) + 1
            out.append({"arm": a.key, "generator": a.generator, "chemical_system": list(a.chemsys),
                        "space_group": a.space_group, "pressure_GPa": a.pressure_GPa,
                        "goal_id": a.goal_id, "n": int(self.d["search"]["structures_per_request"]),
                        "seed": _seed(self.d["seed"], k, i), "extra": {"arm": a.key, "arm_pulls": pulls}})
        return out

    @staticmethod
    def request_of(rq: dict[str, Any]) -> GenerationRequest:
        return GenerationRequest(tuple(rq["chemical_system"]), rq["space_group"], float(rq["pressure_GPa"]),
                                 rq["goal_id"], int(rq["n"]), int(rq["seed"]), dict(rq["extra"]))

    # ================================================================== round phases
    def _run_round(self, k: int, budget_left: int | None) -> None:
        phase = self.ledger.round_row(k)["phase"]
        if phase == "proposed":
            self._generate(k)
            phase = "generated"
        if phase == "generated":
            self._prescreen(k)
            phase = "prescreened"
        if phase == "prescreened":
            self._select(k, budget_left)
            phase = "selected"
        if phase == "selected":
            self._execute(k)
            phase = "executed"
        if phase == "executed":
            self.ledger.save_policy_state(k, {"arms": {a: s.to_dict() for a, s in self.arm_stats().items()},
                                              "seed": self.d["seed"], "kind": self.policy.kind})
            self.write_reports()
            self.ledger.set_phase(k, "done")

    # ---------------------------------------------------------------- generate + dedupe
    def _same(self, a: dict, b: dict) -> bool:
        if self.matcher.available:
            return self.matcher.fit(a, b)
        try:
            return same_structure(a, b, self.confirm_tol)
        except Exception:  # noqa: BLE001 - unconfirmable => not a duplicate
            return False

    def _match(self, s: dict, fp: str, ckey: str) -> tuple[str | None, str, list[str]]:
        """(existing material, how, unconfirmed fingerprint hits).  A fingerprint hit is only a
        candidate; it is a duplicate only when confirmed (StructureMatcher or same_structure)."""
        hits = self.ledger.materials_by_fingerprint(fp)
        unconfirmed = []
        for row in hits:
            if self._same(json.loads(row["structure"]), s):
                return row["material_id"], "fingerprint+confirmed", []
            unconfirmed.append(row["material_id"])
        if self.matcher.available:
            for m in self.ledger.materials_with_coarse_key(ckey):
                if m["material_id"] not in unconfirmed and self.matcher.fit(json.loads(m["structure"]), s):
                    return m["material_id"], "structure_matcher", []
        return None, "", unconfirmed

    def _initial_statuses(self, mid: str, k: int, P: float, why: str, *, error: str | None = None) -> None:
        for g in self.goals:
            self._initial_status_for(g, mid, k, P, why, error=error)

    def _initial_status_for(self, g: Goal, mid: str, k: int, P: float, why: str, *,
                            error: str | None = None) -> None:
        if error is not None:
            self.ledger.append_verdict(material_id=mid, goal_id=g.id, round_=k, stage="generate",
                                       origin="discovery.generate", status=Status.ERROR,
                                       T_star_K=g.T_star_K, pressure_GPa=P, reasons=[error],
                                       flags=["HYPOTHETICAL_STRUCTURE"], in_scope=g.in_scope(P))
        elif g.in_scope(P):
            self.ledger.append_verdict(material_id=mid, goal_id=g.id, round_=k, stage="generate",
                                       origin="discovery.generate", status=Status.DEFERRED,
                                       T_star_K=g.T_star_K, pressure_GPa=P, reasons=[why],
                                       flags=["HYPOTHETICAL_STRUCTURE"])
        else:
            self.ledger.append_verdict(
                material_id=mid, goal_id=g.id, round_=k, stage="scope", origin="discovery.scope",
                status=Status.DEFERRED, T_star_K=g.T_star_K, pressure_GPa=P, in_scope=False,
                flags=["HYPOTHETICAL_STRUCTURE"],
                reasons=[f"outside goal {g.id}: pressure {P:g} GPa > P_max {g.P_max_GPa:g} GPa; "
                         "never evaluated or counted for this goal"])

    def _add_goal_statuses(self, goals: list[Goal]) -> None:
        with self.ledger.tx():
            for m in self.ledger.materials():
                for g in goals:
                    self._initial_status_for(g, m["material_id"], m["round"], m["pressure_GPa"],
                                             f"goal {g.id} added after generation; awaiting selection")

    def _generate(self, k: int) -> None:
        row = self.ledger.round_row(k)
        requests = json.loads(row["requests"])
        produced: list[tuple[int, dict, GeneratedStructure | Any]] = []
        notes = []
        for i, rq in enumerate(requests):
            gen = self.generators.get(rq["generator"])
            try:
                if gen is None:
                    raise KeyError(f"generator {rq['generator']!r} not configured")
                outs = list(gen.generate(self.request_of(rq)))
            except Exception as exc:  # noqa: BLE001 - a failing generator only loses its proposals
                notes.append(f"request {i} ({rq['arm']}): generator failed: {type(exc).__name__}: {exc}")
                outs = []
            produced.extend((i, rq, o) for o in outs)
        prefix = ("SYNTH-" if self.synthetic else "") + str(self.d["id_prefix"])
        per_new = [0] * len(requests)
        n_new = n_dup = n_bad = 0
        with self.ledger.tx():
            seq = self.ledger.next_seq()
            for i, rq, gs in produced:
                P = float(rq["pressure_GPa"])
                gname = str(getattr(gs, "generator", rq["generator"]) or rq["generator"])
                prov = {"generator_provenance": getattr(gs, "provenance", {}),
                        "generator": gname, "request_id": getattr(gs, "request_id", ""),
                        "flags": ["HYPOTHETICAL_STRUCTURE"]}
                raw = getattr(gs, "structure", gs)
                try:
                    s = validate_structure(raw)
                    fp = fingerprint(s, P)
                    ckey = coarse_key(s, P)
                    err = None
                except (InvalidStructure, Exception) as exc:  # noqa: BLE001
                    s, fp, ckey = {"raw": raw}, invalid_fingerprint([raw, P]), "invalid"
                    err = f"invalid generated structure: {type(exc).__name__}: {exc}"
                if err is None:
                    existing, how, possible = self._match(s, fp, ckey)
                else:
                    existing = (self.ledger.material_by_fingerprint(fp) or {"material_id": None})["material_id"]
                    how, possible = "invalid (identical raw output)", []
                if possible:
                    prov["possible_duplicate_of"] = possible
                    prov["dedupe_note"] = ("equal fingerprint but NOT confirmed as the same structure; "
                                           "kept as a separate material")
                if existing is not None:
                    self.ledger.add_duplicate(material_id=existing, fingerprint=fp, generator=gname,
                                              request_id=prov["request_id"], arm=rq["arm"], round_=k,
                                              how=how, provenance=prov)
                    n_dup += 1
                    continue
                formula = reduced_formula(s) if err is None else "invalid"
                chemsys = tuple(sorted(composition_of(s))) if err is None else tuple(rq["chemical_system"])
                mid = f"{prefix}{seq:06d}-{formula}"
                self.ledger.add_material(
                    material_id=mid, seq=seq, fingerprint=fp, coarse_key=ckey, formula=formula,
                    chemsys=chemsys, pressure_GPa=P, structure=s, generator=gname,
                    request_id=prov["request_id"], request=rq, arm=rq["arm"], round_=k,
                    provenance=prov, valid=err is None)
                seq += 1
                if err is None:
                    n_new += 1
                    per_new[i] += 1
                else:
                    n_bad += 1
                self._initial_statuses(mid, k, P, f"generated in round {k} by {gname} "
                                       f"(request {rq['arm']}); not computed yet, awaiting selection",
                                       error=err)
            self.ledger.set_phase(k, "generated", counts={
                "requests": len(requests), "generated": len(produced), "new": n_new,
                "duplicates": n_dup, "invalid": n_bad, "per_request_new": per_new},
                note="; ".join(notes) or None)

    # ---------------------------------------------------------------- ML pre-screen (S2)
    def _ml_decision(self, res: PrescreenResult | None) -> tuple[bool, list[str]]:
        pc = self.d["prescreen"]
        if res is None or not getattr(res, "ok", False):
            return False, []
        why = []
        if pc.get("defer_if_ml_unstable", True) and getattr(res, "dynamically_stable_ml", None) is False:
            mf = _finite(getattr(res, "min_freq_cm1", None))
            why.append("ML pre-screen (HEURISTIC, S2): dynamically unstable"
                       + (f" (min ML frequency {mf:g} cm^-1)" if mf is not None else ""))
        e = _finite(getattr(res, "e_above_hull_ml_eV", None))
        thr = pc.get("max_e_above_hull_ml_eV")
        if e is not None and thr is not None and e > float(thr):
            why.append(f"ML pre-screen (HEURISTIC, S2): ML energy above hull {e:.4g} eV/atom > {thr} "
                       "(ordering/deferral only, never an exclusion)")
        return bool(why), why

    @staticmethod
    def _result_dict(res: PrescreenResult | None) -> dict[str, Any]:
        if res is None:
            return {"ok": False, "notes": ["no pre-screen result returned"]}
        out = {}
        for k in ("fingerprint", "relaxed_structure", "e_above_hull_ml_eV", "dynamically_stable_ml",
                  "min_freq_cm1", "phonon_dos", "ok", "notes"):
            out[k] = getattr(res, k, None)
        out["heuristic"] = True
        return json.loads(json.dumps(out, default=str))

    def _prescreen(self, k: int) -> None:
        mats = [m for m in self.ledger.materials(k) if m["valid"] and m["dup_of"] is None
                and self.ledger.ml_screen(m["material_id"]) is None]
        name = str(getattr(self.prescreener, "name", type(self.prescreener).__name__))
        by_p: dict[float, list] = {}
        for m in mats:
            by_p.setdefault(m["pressure_GPa"], []).append(m)
        for P, group in sorted(by_p.items()):
            gss = [GeneratedStructure(json.loads(m["structure"]), m["generator"], m["request_id"],
                                      {"fingerprint": m["fingerprint"], "material_id": m["material_id"]})
                   for m in group]
            try:
                results = list(self.prescreener.prescreen(gss, float(P)))
            except Exception as exc:  # noqa: BLE001 - heuristic: a failure only loses ordering info
                results = []
                self.ledger.add_round_note(k, f"pre-screener failed at {P} GPa: {type(exc).__name__}: {exc}")
            by_fp = {}
            for r in results:
                fpk = getattr(r, "fingerprint", None)
                if isinstance(fpk, str) and fpk:
                    by_fp.setdefault(fpk, r)
            for j, m in enumerate(group):
                res = by_fp.get(m["fingerprint"])
                if res is None and len(results) == len(group) and not by_fp:
                    res = results[j]
                self._record_prescreen(k, m, res, name)
        self.ledger.set_phase(k, "prescreened", counts={"prescreened": len(mats)})

    def _record_prescreen(self, k: int, m, res, name: str) -> None:
        mid, P = m["material_id"], m["pressure_GPa"]
        rd = self._result_dict(res)
        deferred, why = self._ml_decision(res)
        dup_of = fpr = None
        rel = rd.get("relaxed_structure")
        if rel is not None and rd.get("ok"):
            try:
                relv = validate_structure(rel)
                fpr = fingerprint(relv, P)
                for other in self.ledger.materials_by_fingerprint(fpr):
                    if other["material_id"] == mid:
                        continue
                    ml_o = self.ledger.ml_screen(other["material_id"])
                    o_s = (ml_o or {}).get("result", {}).get("relaxed_structure") \
                        if other["fingerprint_relaxed"] == fpr else None
                    if self._same(o_s or json.loads(other["structure"]), relv):
                        dup_of = other["material_id"]
                        break
            except Exception:  # noqa: BLE001
                rd["relaxed_structure"] = None
                rd.setdefault("notes", []).append("relaxed structure invalid; ignored")
        with self.ledger.tx():
            self.ledger.add_ml_screen(mid, k, name, bool(rd.get("ok")), rd, deferred)
            if fpr is not None:
                self.ledger.set_relaxed_fingerprint(mid, fpr, dup_of)
            reasons = why
            if dup_of:
                reasons = [f"duplicate of {dup_of} after ML relaxation (HEURISTIC); not computed separately"]
            if reasons:
                for g in self.goals:
                    if g.in_scope(P):
                        self.ledger.append_verdict(material_id=mid, goal_id=g.id, round_=k,
                                                   stage="prescreen", origin=f"discovery.prescreen:{name}",
                                                   status=Status.DEFERRED, T_star_K=g.T_star_K,
                                                   pressure_GPa=P, reasons=reasons,
                                                   flags=["HYPOTHETICAL_STRUCTURE"])

    # ---------------------------------------------------------------- select (budget, S2)
    def _pool(self) -> list[dict]:
        done = {r["material_id"] for r in self.ledger.q("SELECT DISTINCT material_id FROM dft_jobs "
                                                        "WHERE state IN ('DONE', 'PENDING', 'SUBMITTED')")}
        out = []
        for m in self.ledger.materials():
            if not m["valid"] or m["dup_of"] is not None or m["material_id"] in done:
                continue
            if not any(g.in_scope(m["pressure_GPa"]) for g in self.goals):
                continue
            out.append(m)
        return out

    def _select(self, k: int, budget_left: int | None) -> None:
        stats = self.arm_stats()
        pool = self._pool()

        retry_ids = {j["material_id"] for j in self.ledger.jobs(None, ("RETRY",))}

        def key(m):
            ml = self.ledger.ml_screen(m["material_id"])
            mldef = int(bool(ml and ml["ml_deferred"]))
            e = _finite(ml["result"].get("e_above_hull_ml_eV")) if ml else None
            return (int(m["material_id"] not in retry_ids), mldef,
                    -stats.get(m["arm"], ArmStats()).mean, e if e is not None else 0.0, m["seq"])

        ranked = sorted(pool, key=key)
        n = int(self.d["budget"]["dft_per_round"])
        if budget_left is not None:
            n = min(n, budget_left)
        chosen = ranked[:n]
        chosen_ids = {m["material_id"] for m in chosen}
        with self.ledger.tx():
            for m in chosen:
                self.ledger.add_dft_job(m["material_id"], k, self.executor.name)
            for rank, m in enumerate(ranked):
                if m["material_id"] in chosen_ids or m["round"] != k:
                    continue
                for g in self.goals:
                    if g.in_scope(m["pressure_GPa"]):
                        self.ledger.append_verdict(
                            material_id=m["material_id"], goal_id=g.id, round_=k, stage="select",
                            origin="discovery.select", status=Status.DEFERRED, T_star_K=g.T_star_K,
                            pressure_GPa=m["pressure_GPa"], flags=["HYPOTHETICAL_STRUCTURE"],
                            reasons=[f"not selected for DFT in round {k} (per-round budget "
                                     f"dft_per_round={self.d['budget']['dft_per_round']}, heuristic "
                                     f"priority rank {rank}); stays in the pool for later rounds"])
            self.ledger.set_phase(k, "selected", counts={"selected": len(chosen), "pool": len(pool)})

    # ---------------------------------------------------------------- execute (existing pipeline)
    def candidate(self, mid: str) -> Candidate:
        m = self.ledger.material(mid)
        s = json.loads(m["structure"])
        used = "generated"
        ml = self.ledger.ml_screen(mid)
        mlsum = None
        if ml is not None:
            r = ml["result"]
            mlsum = {k: r.get(k) for k in ("ok", "e_above_hull_ml_eV", "dynamically_stable_ml",
                                           "min_freq_cm1", "notes")}
            mlsum["prescreener"] = ml["prescreener"]
            if self.d["prescreen"].get("use_ml_relaxed", True) and r.get("ok") and r.get("relaxed_structure"):
                s = validate_structure(r["relaxed_structure"])
                used = "ml_relaxed (heuristic starting geometry; the DFT stage relaxes again)"
        P = float(m["pressure_GPa"])
        prov = json.loads(m["provenance"])
        c = candidate_from_structure_arrays(
            mid, s["lattice"], s["species"], s["frac_coords"], source=f"discovery:{m['generator']}",
            pressure_GPa=P,
            provenance={"source": "discovery", "generator": m["generator"], "request_id": m["request_id"],
                        "arm": m["arm"], "round": m["round"], "fingerprint": m["fingerprint"],
                        "flags": ["HYPOTHETICAL_STRUCTURE"], "structure_used": used,
                        "generator_provenance": prov.get("generator_provenance", {}),
                        # HEURISTIC ML results live in provenance only: no filter reads them (S2)
                        "ml_prescreen_heuristic": mlsum})
        c.properties = {"target_pressure_GPa": P}
        if self.synthetic:
            c.synthetic = True
            c.spectra = [synthetic_spectrum(m["fingerprint"], tuple(m["chemsys"].split("-")), P,
                                            self.d["synthetic_spectra"])]
        return c

    def _goal_dir(self, k: int, g: Goal) -> Path:
        return self.workdir / "rounds" / f"{k:04d}" / g.id

    def _first_goal(self, P: float) -> Goal | None:
        for g in self.goals:
            if g.in_scope(P):
                return g
        return None

    def _execute(self, k: int) -> None:
        jobs = self.ledger.jobs(k, ("PENDING", "SUBMITTED", "RUNNING"))
        mids = [j["material_id"] for j in jobs]
        ext = {j["material_id"]: j["external_id"] for j in jobs}
        cands = {mid: self.candidate(mid) for mid in mids}
        first: dict[str, Goal] = {}
        for mid in mids:
            g = self._first_goal(cands[mid].pressure_GPa)
            if g is not None:
                first[mid] = g
        records: dict[str, dict[str, dict | None]] = {}    # goal -> mid -> record
        for g in self.goals:
            group = [mid for mid in mids if first.get(mid) is g]
            if not group:
                continue
            cfg = self.dc.goal_pipeline_config(g, self._goal_dir(k, g))
            existing = next((ext[m] for m in group if ext.get(m)), None)

            def on_submit(jid: str, group=group) -> None:
                with self.ledger.tx():
                    for m in group:
                        self.ledger.set_job_state(m, k, "SUBMITTED", external_id=jid, attempt=True)

            recs = self.executor.run_pipeline(cfg, [cands[m] for m in group], on_submit=on_submit,
                                              existing_job=existing)
            records[g.id] = {m: recs.get(m) for m in group}
            with self.ledger.tx():
                for m in group:
                    self._record(k, m, g, records[g.id][m], cfg, cands[m])
        # further in-scope goals: certify the SAME spectra at their T* (certify-only resume path)
        for g in self.goals:
            todo = [m for m in mids if m in first and first[m] is not g and g.in_scope(cands[m].pressure_GPa)]
            if not todo:
                continue
            cfg = self.dc.goal_pipeline_config(g, self._goal_dir(k, g))
            for m in todo:
                src = records.get(first[m].id, {}).get(m)
                src_path = pipeline.record_path(self._goal_dir(k, first[m]), m)
                rec, copied = self._other_goal_record(m, cands[m], src, cfg, src_path)
                with self.ledger.tx():
                    self._record(k, m, g, rec, cfg, cands[m], src_record_path=str(src_path) if copied else None)
            self._rewrite_goal_reports(cfg, [j["material_id"] for j in self.ledger.jobs(k)])
        with self.ledger.tx():
            for mid in mids:
                retry = any(r["retryable"] and r["status"] == Status.ERROR.value
                            for r in self.ledger.current(None, mid))
                self.ledger.set_job_state(mid, k, "RETRY" if retry else "DONE")
            self.ledger.set_phase(k, "executed", counts={"dft_jobs": len(mids)})

    def _other_goal_record(self, mid: str, c: Candidate, src: dict | None, cfg: dict,
                           src_path: Path) -> tuple[dict | None, str | None]:
        if src is None or not src.get("verdict"):
            return None, None
        wd = cfg["run"]["workdir"]
        chash = pipeline.config_hash(cfg)
        digest = pipeline.input_digest(c)
        old = pipeline.read_record(wd, mid)
        if pipeline.record_reusable(old, c, chash, digest) and old.get("verdict"):
            return old, None                                  # resumed: never recomputed
        if src.get("verdict_origin") == pipeline.CLASSIFY_ORIGIN:
            if not (src.get("spectra") and src.get("input_digest") == digest):
                return None, None                             # cannot seed: retry (never copy a proof)
            rec = pipeline.new_record(c, chash, digest)
            rec["spectra"] = copy.deepcopy(src["spectra"])
            rec["stages"] = {k: copy.deepcopy(v) for k, v in src.get("stages", {}).items()
                             if k not in ("certify", "validity", "error")}
            rec["stages"]["seeded_from"] = {
                "record": str(src_path), "config_hash": src.get("config_hash"), "input_digest": src.get("input_digest"),
                "note": "alternative spectra and DFT stages of the same material in the same round; "
                        "only certify (at this goal's T*) and validity flags run here"}
            return pipeline.compute_material(c.to_dict(), rec, cfg, wd), None
        return src, src.get("config_hash")                    # non-classify outcome: reused as is

    def _rewrite_goal_reports(self, cfg: dict, mids) -> None:
        wd = Path(cfg["run"]["workdir"])
        recs = []
        for mid in mids:
            r = pipeline.read_record(wd, mid)
            if r and r.get("verdict"):
                recs.append(r)
        if not recs:
            return
        verdicts = {r["material_id"]: Verdict.from_dict(r["verdict"]) for r in recs}
        origins = {r["material_id"]: r.get("verdict_origin") for r in recs}
        summ = pipeline.RunSummary(list(verdicts), verdicts, origins, {}, pipeline.config_hash(cfg),
                                   "", pipeline._now(), None)
        pipeline.account(list(verdicts), verdicts, origins)
        _pipeline_reports(cfg, summ, recs)

    @staticmethod
    def dft_pressure(rec: dict, c: Candidate) -> float | None:
        """Pressure at which the record's spectra were actually computed (None = unknown):
        DFPT provenance (pressure_GPa; None if not relaxed), else the spectra's own provenance
        pressure, else (attached spectra without pressure) the candidate's target pressure."""
        d = (rec.get("stages") or {}).get("dfpt")
        if d is not None:
            prov = d.get("provenance") or {}
            if prov.get("relaxed") is False or prov.get("pressure_GPa") is None:
                return None
            return float(prov["pressure_GPa"])
        ps = {(sd.get("provenance") or {}).get("pressure_GPa") for sd in rec.get("spectra") or []}
        ps.discard(None)
        if len(ps) > 1:
            return None
        if ps:
            return float(ps.pop())
        return float(c.properties.get("target_pressure_GPa", c.pressure_GPa))

    def _record(self, k: int, mid: str, g: Goal, rec: dict | None, cfg: dict, c: Candidate,
                src_record_path: str | None = None, note: str | None = None) -> None:
        P = float(c.pressure_GPa)
        wd = Path(cfg["run"]["workdir"])
        if rec is None or not rec.get("verdict"):
            fn = getattr(self.executor, "failure_note", None)
            extra = fn() if callable(fn) else ""
            self.ledger.append_verdict(
                material_id=mid, goal_id=g.id, round_=k, stage="executor", origin=EXECUTOR_ORIGIN,
                status=Status.ERROR, T_star_K=g.T_star_K, pressure_GPa=P, retryable=True,
                flags=["HYPOTHETICAL_STRUCTURE"],
                reasons=["no finished pipeline record (job failed, was killed or is still running); "
                         "retried on continue" + (f"; {extra}" if extra else "")])
            return
        v = Verdict.from_dict(rec["verdict"])
        origin = str(rec.get("verdict_origin"))
        dp = self.dft_pressure(rec, c)
        if v.status in PROVEN_STATUSES:
            T_used = float(rec.get("stages", {}).get("certify", {}).get("params", {}).get("T_target_K", float("nan")))
            if T_used != g.T_star_K:
                raise AssertionError(f"{mid}: proven verdict made at T*={T_used}, goal {g.id} needs {g.T_star_K}")
            if dp is None or abs(dp - P) > 1e-9:
                self.ledger.append_verdict(
                    material_id=mid, goal_id=g.id, round_=k, stage="executor", origin=EXECUTOR_ORIGIN,
                    status=Status.ERROR, T_star_K=g.T_star_K, pressure_GPa=P, dft_pressure_GPa=dp,
                    dft_hash=self.dft_hash, certify_hash=self.dc.certify_hash(g),
                    flags=["HYPOTHETICAL_STRUCTURE", "DFT_PRESSURE_MISMATCH"],
                    record_path=str(pipeline.record_path(wd, mid)),
                    reasons=[f"certify.classify returned {v.status.value}, but the spectra were "
                             f"computed at {dp} GPa (unknown if None), not at the material's {P:g} "
                             "GPa: the verdict is not about this material and is not recorded"])
                return
        rp = src_record_path or str(pipeline.record_path(wd, mid))
        cur = self.ledger.current_one(mid, g.id)
        if cur and cur["stage"] == "pipeline" and cur["record_path"] == rp and \
                cur["input_digest"] == rec.get("input_digest") and cur["config_hash"] == rec.get("config_hash") \
                and cur["status"] == v.status.value and cur["certify_hash"] == self.dc.certify_hash(g) \
                and cur["dft_hash"] == self.dft_hash:
            return                                            # already recorded (resumed round)
        flags = list(v.flags)
        for f in ("HYPOTHETICAL_STRUCTURE",) + (("SYNTHETIC_DATA",) if c.synthetic else ()):
            if f not in flags:
                flags.append(f)
        reasons = list(v.reasons)
        if src_record_path is not None:
            reasons.append(f"outcome of the goal-independent stages (existence filters / DFT) of the "
                           f"same structure and pressure ({src_record_path}); certify.classify was not reached")
        if note:
            reasons.append(note)
        retry = v.status is Status.ERROR and (bool(rec.get("retryable")) or
                                              origin in ("pipeline.worker", "pipeline.collect"))
        cert_dir = Path(cfg.get("report", {}).get("outdir") or wd) if src_record_path is None else \
            Path(src_record_path).parent.parent
        cert = cert_dir / "certificates" / f"{safe_id(mid)}.json"
        self.ledger.append_verdict(
            material_id=mid, goal_id=g.id, round_=k, stage="pipeline", origin=origin, status=v.status,
            T_star_K=g.T_star_K, pressure_GPa=P, T_low_K=v.T_low_K, T_up_K=v.T_up_K,
            rule_ids=v.rule_ids, assumptions=v.assumptions, flags=flags, reasons=reasons,
            certificate_path=str(cert), record_path=rp, input_digest=rec.get("input_digest"),
            config_hash=rec.get("config_hash"), dft_hash=self.dft_hash,
            certify_hash=self.dc.certify_hash(g), dft_pressure_GPa=dp, retryable=retry,
            in_scope=g.in_scope(P))
