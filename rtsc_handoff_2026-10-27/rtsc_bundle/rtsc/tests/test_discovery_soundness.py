"""Soundness tests of the discovery loop (DISCOVERY.md + ARCH S1-S6).

(1) S2 extended: generators, pre-screeners, policy and budget never produce REJECTED_PROVEN,
    CERTIFIED_CANDIDATE or EXCLUDED_ASSUMPTION (property test with random adversarial mocks).
(2) Continuation after a kill mid-round: nothing lost, nothing duplicated, no proven verdict
    recomputed or recorded with a mismatched digest; same final ledger as an uninterrupted run.
(3) Ledger accounting: every generated structure = one material or one duplicate link; exactly one
    current status per material x goal.
(4) Goal semantics: CERTIFIED at T* = 77 K never counts for G-RT; ambient needs P <= 1e-4 GPa.
"""
from __future__ import annotations

import sqlite3
import sys
from collections import Counter
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from discovery_mocks import (FixedGenerator, Kill, KillingExecutor, RandomGenerator,  # noqa: E402
                             RandomPrescreener, make_dc)

import rtsc.discovery.loop as L  # noqa: E402
import rtsc.pipeline as P  # noqa: E402
from rtsc.discovery.generators import NullPrescreener  # noqa: E402
from rtsc.discovery.goals import normalize_discovery  # noqa: E402
from rtsc.discovery.ledger import Ledger, LedgerSoundnessError  # noqa: E402
from rtsc.discovery.loop import DiscoveryLoop  # noqa: E402
from rtsc.pipeline import CLASSIFY_ORIGIN, ConfigError  # noqa: E402
from rtsc.sources.synthetic import gaussian_peak_spectrum  # noqa: E402
from rtsc.status import Status  # noqa: E402

PROVEN_OR_EXCLUDED = {"REJECTED_PROVEN", "CERTIFIED_CANDIDATE", "EXCLUDED_ASSUMPTION"}


def _all_verdicts(loop):
    return loop.ledger.history()


# ------------------------------------------------------------------------------- (1) S2 extended
@pytest.mark.parametrize("seed", [0, 1, 2, 3, 4, 5])
def test_s2_heuristics_never_prove_or_exclude(tmp_path, seed):
    """No spectra (synthetic off, DFPT off): with random adversarial generators / ML / policy
    seeds / budgets, the ledger may contain only DEFERRED and ERROR."""
    import random
    r = random.Random(seed)
    dc = make_dc(tmp_path / "w", seed=seed, synthetic_spectra={"enabled": False},
                 budget={"dft_per_round": r.randint(0, 4)},
                 search={"policy": r.choice(["thompson", "ucb"]), "exploration_floor": r.random(),
                         "requests_per_round": r.randint(1, 4), "structures_per_request": r.randint(1, 4)},
                 prescreen={"defer_if_ml_unstable": r.random() < 0.5,
                            "max_e_above_hull_ml_eV": r.choice([0.0, 0.5]),
                            "use_ml_relaxed": r.random() < 0.5},
                 stop={"max_rounds": 3, "certified_per_goal": 100})
    gen = RandomGenerator(seed)
    loop = DiscoveryLoop(dc, generators={"randgen": gen}, prescreener=RandomPrescreener(seed))
    loop.run()
    vs = _all_verdicts(loop)
    assert vs, "the loop recorded nothing"
    bad = [v for v in vs if v["status"] in PROVEN_OR_EXCLUDED]
    assert not bad, bad[:3]
    for v in vs:
        if v["stage"] != "pipeline":
            assert v["status"] in ("DEFERRED", "ERROR")
            if v["status"] == "ERROR":
                assert v["stage"] in ("generate", "executor")
    acc = loop.ledger.accounting(["G-RT", "G-AMB"])
    assert acc["generated"] == gen.count          # (3) every generated structure accounted for


def test_s2_ml_outputs_do_not_change_dft_verdicts(tmp_path):
    """With SYNTHETIC spectra: proven statuses come only from certify.classify, no EXCLUDED, and the
    verdict of every material computed under both a random adversarial ML pre-screener and no ML is
    the same (ML may only reorder / defer)."""
    res = {}
    for name, pre in (("null", NullPrescreener()), ("random", RandomPrescreener(11))):
        dc = make_dc(tmp_path / name, budget={"dft_per_round": 20},
                     search={"chemical_systems": [["H", "La"], ["B", "Mg"]], "requests_per_round": 2,
                             "structures_per_request": 3},
                     stop={"max_rounds": 1, "certified_per_goal": 100})
        loop = DiscoveryLoop(dc, prescreener=pre)
        loop.run()
        for v in _all_verdicts(loop):
            if v["status"] in ("REJECTED_PROVEN", "CERTIFIED_CANDIDATE"):
                assert v["stage"] == "pipeline" and v["origin"] == CLASSIFY_ORIGIN
            assert v["status"] != "EXCLUDED_ASSUMPTION"
        fps = {m["material_id"]: m["fingerprint"] for m in loop.ledger.materials()}
        res[name] = {(fps[r["material_id"]], r["goal_id"]): (r["status"], r["T_low_K"], r["T_up_K"])
                     for g in ("G-RT", "G-AMB") for r in loop.ledger.current(g) if r["stage"] == "pipeline"}
    common = set(res["null"]) & set(res["random"])
    assert common
    for k in common:
        assert res["null"][k] == res["random"][k], k


def test_ledger_guard_and_append_only(tmp_path):
    led = Ledger(tmp_path / "l.sqlite")
    led.add_material(material_id="M1", seq=1, fingerprint="fp", coarse_key="H|P0", formula="H",
                     chemsys=("H",), pressure_GPa=0.0, structure={}, generator="g", request_id="r",
                     request={}, arm="a", round_=1, provenance={})
    for stage in ("generate", "scope", "prescreen", "select", "policy", "executor"):
        for st in (Status.REJECTED_PROVEN, Status.CERTIFIED_CANDIDATE, Status.EXCLUDED_ASSUMPTION):
            with pytest.raises(LedgerSoundnessError):
                led.append_verdict(material_id="M1", goal_id="G-RT", round_=1, stage=stage,
                                   origin=CLASSIFY_ORIGIN, status=st, T_star_K=300, pressure_GPa=0)
    with pytest.raises(LedgerSoundnessError):     # pipeline stage, wrong origin
        led.append_verdict(material_id="M1", goal_id="G-RT", round_=1, stage="pipeline",
                           origin="discovery.ml", status=Status.CERTIFIED_CANDIDATE, T_star_K=300,
                           pressure_GPa=0)
    with pytest.raises(LedgerSoundnessError):
        led.append_verdict(material_id="M1", goal_id="G-RT", round_=1, stage="pipeline",
                           origin="discovery.ml", status=Status.EXCLUDED_ASSUMPTION, T_star_K=300,
                           pressure_GPa=0)
    vid = led.append_verdict(material_id="M1", goal_id="G-RT", round_=1, stage="prescreen",
                             origin="x", status=Status.DEFERRED, T_star_K=300, pressure_GPa=0)
    with pytest.raises(sqlite3.DatabaseError):
        led.con.execute("UPDATE verdicts SET status='CERTIFIED_CANDIDATE' WHERE id=?", (vid,))
    with pytest.raises(sqlite3.DatabaseError):
        led.con.execute("DELETE FROM verdicts WHERE id=?", (vid,))


# ------------------------------------------------------------------------------- (2) continuation
def _snapshot(loop):
    mats = {m["material_id"]: m["fingerprint"] for m in loop.ledger.materials()}
    cur = {(r["material_id"], r["goal_id"]): (r["status"], r["T_low_K"], r["T_up_K"], r["stage"])
           for g in ("G-RT", "G-AMB") for r in loop.ledger.current(g)}
    return mats, cur


@pytest.mark.parametrize("after", [0, 2])
def test_continue_after_kill(tmp_path, monkeypatch, after):
    calls: Counter = Counter()
    orig = P.certify_spectra

    def counting(mid, spectra, ccfg, **kw):
        calls[(mid, float(ccfg["T_target_K"]))] += 1
        return orig(mid, spectra, ccfg, **kw)

    monkeypatch.setattr(P, "certify_spectra", counting)
    stop = {"max_rounds": 3, "certified_per_goal": 100}
    # reference: uninterrupted
    ref = DiscoveryLoop(make_dc(tmp_path / "ref", stop=stop))
    ref.run()
    ref_snap = _snapshot(ref)
    calls.clear()
    # killed in round 2 after `after` certify calls, then continued
    dc = make_dc(tmp_path / "run", stop=stop)
    with pytest.raises(Kill):
        DiscoveryLoop(dc, executor=KillingExecutor("/0002/", after)).run()
    led = Ledger(dc.workdir / "ledger.sqlite")
    assert led.last_round()["round"] == 2 and led.last_round()["phase"] == "selected"
    n_before = len(led.materials())
    led.close()
    loop = DiscoveryLoop(make_dc(tmp_path / "run", stop={"max_rounds": 2, "certified_per_goal": 100}))
    res = loop.run("continue")
    assert res["rounds_this_session"] == 2
    # nothing lost / duplicated; identical to the uninterrupted run
    assert len(loop.ledger.materials()) >= n_before
    assert _snapshot(loop) == ref_snap
    # no material x T* was certified twice (finished work is never recomputed)
    assert calls and max(calls.values()) == 1, [k for k, v in calls.items() if v > 1]
    # each pipeline verdict recorded once, with digest/config hash equal to the record on disk
    seen = Counter()
    for v in loop.ledger.history():
        if v["stage"] != "pipeline":
            continue
        seen[(v["material_id"], v["goal_id"], v["input_digest"], v["config_hash"])] += 1
        rec = P.read_record(Path(v["record_path"]).parent.parent, v["material_id"])
        assert rec["input_digest"] == v["input_digest"] and rec["config_hash"] == v["config_hash"]
        assert rec["input_digest"] == P.input_digest(loop.candidate(v["material_id"]))
        assert rec["verdict"]["status"] == v["status"]
    assert max(seen.values()) == 1
    loop.ledger.accounting(["G-RT", "G-AMB"])


def test_retryable_executor_error_is_retried(tmp_path):
    """A job whose record never appeared (crashed node) -> retryable ERROR -> recomputed later."""
    from rtsc.discovery.executors import LocalExecutor

    class Lossy(LocalExecutor):
        def run_pipeline(self, cfg, candidates, **kw):
            recs = super().run_pipeline(cfg, candidates, **kw)
            if "/0001/" in cfg["run"]["workdir"]:
                recs[candidates[0].material_id] = None
            return recs

    dc = make_dc(tmp_path / "w", stop={"max_rounds": 1, "certified_per_goal": 100})
    loop = DiscoveryLoop(dc, executor=Lossy())
    loop.run()
    lost = [r for r in loop.ledger.current("G-RT") if r["status"] == "ERROR"]
    assert len(lost) == 1 and lost[0]["retryable"] and lost[0]["stage"] == "executor"
    mid = lost[0]["material_id"]
    loop2 = DiscoveryLoop(make_dc(tmp_path / "w", stop={"max_rounds": 1, "certified_per_goal": 100}))
    loop2.run("continue")
    assert loop2.ledger.current_one(mid, "G-RT")["stage"] == "pipeline"
    assert loop2.ledger.current_one(mid, "G-RT")["round"] == 2


# ------------------------------------------------------------------------------- (3) accounting
def test_accounting_every_structure_one_status_per_goal(tmp_path):
    gen = RandomGenerator(42)
    dc = make_dc(tmp_path / "w", synthetic_spectra={"enabled": False}, stop={"max_rounds": 4},
                 search={"requests_per_round": 4, "structures_per_request": 4})
    loop = DiscoveryLoop(dc, generators={"randgen": gen}, prescreener=RandomPrescreener(1))
    loop.run()
    acc = loop.ledger.accounting(["G-RT", "G-AMB"])
    assert acc["generated"] == gen.count and acc["duplicates"] > 0
    assert acc["current_pairs"] == 2 * acc["materials"]
    for m in loop.ledger.materials():
        for g in ("G-RT", "G-AMB"):
            assert loop.ledger.current_one(m["material_id"], g) is not None
    dup_rows = loop.ledger.q("SELECT d.material_id, d.fingerprint, d.how, m.fingerprint AS mfp, "
                             "m.fingerprint_relaxed AS mfr FROM duplicates d "
                             "JOIN materials m ON m.material_id=d.material_id")
    assert {r["how"] for r in dup_rows} <= {"fingerprint+confirmed", "structure_matcher",
                                            "invalid (identical raw output)"}
    # a fingerprint link must be a real fingerprint hit; a StructureMatcher-confirmed duplicate may have a
    # different fingerprint (the fingerprint is only a candidate key)
    assert all(r["fingerprint"] in (r["mfp"], r["mfr"]) for r in dup_rows if r["how"] != "structure_matcher")


# ------------------------------------------------------------------------------- (4) goal semantics
_CUBE = {"lattice": [[3.0, 0, 0], [0, 3.0, 0], [0, 0, 3.0]], "species": ["La", "H"],
         "frac_coords": [[0, 0, 0], [0.5, 0.5, 0.5]]}


def _spec(center, A):
    def f(fp, chemsys, p, scfg):
        return gaussian_peak_spectrum(center, center * 0.06, A, label="SYNTHETIC goal-test")
    return f


def test_77K_certificate_never_counts_for_G_RT(tmp_path, monkeypatch):
    monkeypatch.setattr(L, "synthetic_spectrum", _spec(800.0, 900.0))   # A/a0 ~ 208 K < 300 K
    dc = make_dc(tmp_path / "w", goals={"G-RT": {"pressures_GPa": [0.0]}},
                 search={"chemical_systems": [["H", "La"]], "requests_per_round": 1,
                         "structures_per_request": 1},
                 stop={"max_rounds": 1, "certified_per_goal": 100})
    loop = DiscoveryLoop(dc, generators={"fixed": FixedGenerator([_CUBE])})
    loop.run()
    (mid,) = [m["material_id"] for m in loop.ledger.materials()]
    amb = loop.ledger.current_one(mid, "G-AMB")
    rt = loop.ledger.current_one(mid, "G-RT")
    assert amb["status"] == "CERTIFIED_CANDIDATE" and amb["T_star_K"] == 77.0
    assert rt["status"] == "REJECTED_PROVEN" and rt["T_star_K"] == 300.0
    assert loop.certified_ids(dc.goal("G-AMB")) == {mid}
    assert loop.certified_ids(dc.goal("G-RT")) == set()
    # a (forged) classify verdict at 77 K filed under G-RT is still not counted for G-RT
    loop.ledger.append_verdict(material_id=mid, goal_id="G-RT", round_=1, stage="pipeline",
                               origin=CLASSIFY_ORIGIN, status=Status.CERTIFIED_CANDIDATE,
                               T_star_K=77.0, pressure_GPa=0.0, T_low_K=100.0)
    assert loop.certified_ids(dc.goal("G-RT")) == set()


def test_ambient_goal_requires_ambient_pressure(tmp_path, monkeypatch):
    monkeypatch.setattr(L, "synthetic_spectrum", _spec(2000.0, 3000.0))  # strong: certified at 300 K
    dc = make_dc(tmp_path / "w", goals={"G-RT": {"pressures_GPa": [5.0]}},
                 search={"chemical_systems": [["H", "La"]], "requests_per_round": 1,
                         "structures_per_request": 1},
                 stop={"max_rounds": 1, "certified_per_goal": 100})
    loop = DiscoveryLoop(dc, generators={"fixed": FixedGenerator([_CUBE])})
    loop.run()
    rows = loop.ledger.materials()
    hp = [m["material_id"] for m in rows if m["pressure_GPa"] == 5.0]
    assert hp
    for mid in hp:
        assert loop.ledger.current_one(mid, "G-RT")["status"] == "CERTIFIED_CANDIDATE"
        a = loop.ledger.current_one(mid, "G-AMB")
        assert a["status"] == "DEFERRED" and a["stage"] == "scope" and not a["in_scope"]
        assert not any(v["goal_id"] == "G-AMB" and v["stage"] == "pipeline"
                       for v in loop.ledger.history(mid))
        assert mid not in loop.certified_ids(dc.goal("G-AMB"))
        with pytest.raises(LedgerSoundnessError):
            loop.ledger.append_verdict(material_id=mid, goal_id="G-AMB", round_=1, stage="pipeline",
                                       origin=CLASSIFY_ORIGIN, status=Status.CERTIFIED_CANDIDATE,
                                       T_star_K=77.0, pressure_GPa=5.0, in_scope=False)
    for bad in ({"G-AMB": {"P_max_GPa": 1e-3}}, {"G-AMB": {"T_star_K": 350.0}},
                {"G-AMB": {"pressures_GPa": [1.0]}}):
        with pytest.raises(ConfigError):
            make_dc(tmp_path / "x", goals=bad)


# ------------------------------------------------------------------------------- DFT pressure
def test_dft_pressure_policy_and_mismatch(tmp_path, monkeypatch):
    import rtsc.dfpt as D
    from discovery_mocks import base_raw
    from rtsc.spectrum import Spectrum

    for settings in ({"pressure_GPa": 150.0}, {"press_kbar": 10.0}, {"do_relax": False},
                     {"cell_dofree": "shape"}):
        raw = base_raw(tmp_path / "x", synthetic_spectra={"enabled": False})
        raw["dfpt"] = {"enabled": True, "pseudos": {"dir": str(tmp_path)}, "settings": settings}
        with pytest.raises(ConfigError):
            normalize_discovery(raw, Path(__file__).resolve().parents[1])

    class Res:                                 # a DFT stage that (wrongly) ran at 150 GPa
        spectra = [Spectrum.from_dict(gaussian_peak_spectrum(800.0, 48.0, 700.0, label="fake DFPT"))]

        def to_dict(self):
            return {"provenance": {"pressure_GPa": 150.0, "relaxed": True}}

    class WF:
        def __init__(self, *a, **k):
            pass

        def run(self):
            return Res()

    monkeypatch.setattr(D, "EPhWorkflow", WF)
    raw = base_raw(tmp_path / "w", goals={"G-RT": {"pressures_GPa": [0.0]}},
                   search={"chemical_systems": [["H", "La"]], "requests_per_round": 1,
                           "structures_per_request": 1},
                   stop={"max_rounds": 1}, synthetic_spectra={"enabled": False})
    raw["dfpt"] = {"enabled": True, "pseudos": {"dir": str(tmp_path)}}
    dc = normalize_discovery(raw, Path(__file__).resolve().parents[1])
    loop = DiscoveryLoop(dc, generators={"fixed": FixedGenerator([_CUBE])})
    loop.run()
    (m,) = loop.ledger.materials()
    for g in ("G-RT", "G-AMB"):
        cur = loop.ledger.current_one(m["material_id"], g)
        assert cur["status"] == "ERROR" and "DFT_PRESSURE_MISMATCH" in cur["flags"], cur
        assert cur["dft_pressure_GPa"] == 150.0
    assert loop.certified_ids(dc.goal("G-AMB")) == set()


# ------------------------------------------------------------------------------- D5 / D3 (Lean)
def _ledger_counts(loop):
    return (loop.ledger.q("SELECT COUNT(*) AS n FROM verdicts")[0]["n"],
            sorted(p.name for p in (loop.workdir / "reconcile").glob("*")) if (loop.workdir / "reconcile").exists() else [])


def test_reconcile_idempotent_with_pressure_mismatch_error(tmp_path, monkeypatch):
    """One goal has a pipeline row (UNDECIDED at 300 K), the other the non-retryable
    DFT_PRESSURE_MISMATCH ERROR (proven at 77 K, spectra at 150 GPa).  Continuing three times with
    an unchanged config appends no verdict and creates no reconcile directory."""
    import rtsc.dfpt as D
    from discovery_mocks import base_raw
    from rtsc.spectrum import Spectrum

    class Res:
        spectra = [Spectrum.from_dict(gaussian_peak_spectrum(4000.0, 200.0, 1700.0, label="fake DFPT"))]

        def to_dict(self):
            return {"provenance": {"pressure_GPa": 150.0, "relaxed": True}}

    class WF:
        def __init__(self, *a, **k):
            pass

        def run(self):
            return Res()

    monkeypatch.setattr(D, "EPhWorkflow", WF)
    raw = base_raw(tmp_path / "w", goals={"G-RT": {"pressures_GPa": [0.0]}},
                   search={"chemical_systems": [["H", "La"]], "requests_per_round": 1,
                           "structures_per_request": 1},
                   stop={"max_rounds": 1}, synthetic_spectra={"enabled": False})
    raw["dfpt"] = {"enabled": True, "pseudos": {"dir": str(tmp_path)}}
    root = Path(__file__).resolve().parents[1]
    loop = DiscoveryLoop(normalize_discovery(raw, root), generators={"fixed": FixedGenerator([_CUBE])})
    loop.run()
    (m,) = loop.ledger.materials()
    rt = loop.ledger.current_one(m["material_id"], "G-RT")
    amb = loop.ledger.current_one(m["material_id"], "G-AMB")
    assert rt["stage"] == "pipeline" and rt["status"] == "UNDECIDED", rt
    assert amb["status"] == "ERROR" and "DFT_PRESSURE_MISMATCH" in amb["flags"] and not amb["retryable"]
    before = _ledger_counts(loop)
    loop.ledger.close()
    for _ in range(3):
        lp = DiscoveryLoop(normalize_discovery(raw, root), generators={"fixed": FixedGenerator([_CUBE])})
        lp.run("continue", max_rounds=0)
        assert _ledger_counts(lp) == before
        lp.ledger.close()


def test_reconcile_idempotent_synthetic_run(tmp_path):
    dc = make_dc(tmp_path / "w", stop={"max_rounds": 2, "certified_per_goal": 100})
    loop = DiscoveryLoop(dc)
    loop.run()
    before = _ledger_counts(loop)
    for _ in range(3):
        lp = DiscoveryLoop(make_dc(tmp_path / "w", stop={"max_rounds": 0}))
        lp.run("continue")
        assert _ledger_counts(lp) == before
    assert before[1] == []


def test_accounting_checks_generator_outputs_per_round(tmp_path):
    gen = RandomGenerator(7)
    dc = make_dc(tmp_path / "w", synthetic_spectra={"enabled": False}, stop={"max_rounds": 2})
    loop = DiscoveryLoop(dc, generators={"randgen": gen}, prescreener=RandomPrescreener(2))
    loop.run()
    acc = loop.ledger.accounting(["G-RT", "G-AMB"])
    assert acc["generated"] == gen.count
    # a duplicate row that no generator output produced breaks accounting
    m = loop.ledger.materials()[0]
    loop.ledger.add_duplicate(material_id=m["material_id"], fingerprint=m["fingerprint"],
                              generator="forged", request_id="r", arm="a", round_=1, how="forged",
                              provenance={})
    with pytest.raises(AssertionError, match="generator outputs"):
        loop.ledger.accounting(["G-RT", "G-AMB"])
