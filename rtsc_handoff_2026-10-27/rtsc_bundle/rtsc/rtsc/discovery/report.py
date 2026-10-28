"""Discovery outputs (DISCOVERY.md section 6), rewritten after every round:

* ``discovery_report.md``: assumption banner, goals, progress per round, top candidates per goal,
  open (UNDECIDED) list ranked by certified T_up, rejected list with certificates, accounting;
* ``candidates_<goal>.csv``: every material with its current status for that goal;
* ``export`` (JSON dump of the ledger) and ``export_finetune`` (DFT-level labels, see below).
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any

from ..status import Status

ORDER = [Status.CERTIFIED_CANDIDATE, Status.UNDECIDED, Status.REJECTED_PROVEN,
         Status.EXCLUDED_ASSUMPTION, Status.DEFERRED, Status.ERROR]
_RANK = {s.value: i for i, s in enumerate(ORDER)}

NOT_COUNTING = "CERTIFIED_CANDIDATE (does not count)"

CSV_COLUMNS = ["material_id", "formula", "chemsys", "pressure_GPa", "dft_pressure_GPa", "status",
               "counts_for_goal", "not_counting_reason", "T_low_K", "T_up_K", "T_star_K", "in_scope", "stage", "origin", "generator", "arm", "round_generated",
               "round_verdict", "synthetic", "flags", "rule_ids", "assumptions", "reasons",
               "certificate_path", "record_path", "input_digest", "config_hash", "dft_hash",
               "certify_hash", "fingerprint", "dup_of"]

BANNER = """> **Model assumptions (read first).** Every REJECTED / CERTIFIED / UNDECIDED verdict is a
> statement about the MODEL: isotropic Migdal-Eliashberg theory, linearised gap equation, block-form
> Coulomb cutoff, for the given alpha^2F (the DFPT spectra, or attached spectra), within the stated
> envelope and mu* range (see each verdict's assumptions).  Out-of-model risks (non-adiabatic,
> anharmonic, multiband, ...) are FLAGS and never upgrade a verdict.  All structures here are
> HYPOTHETICAL (generated); ML pre-screen results and the search policy are HEURISTIC and only
> order candidates or mark them DEFERRED.  A material satisfies a goal only if certify.classify
> returned CERTIFIED_CANDIDATE at that goal's T* and a pressure inside the goal's range."""

SYNTH_BANNER = """> **SYNTHETIC DATA.** This run attaches INVENTED alpha^2F spectra (labelled SYNTHETIC) to
> generated structures to exercise the loop without DFT.  No number in this report is data for
> any real compound, whatever its formula."""


def _sort_key(r: dict) -> tuple:
    st = r["status"]
    if st == Status.CERTIFIED_CANDIDATE.value:
        t = -(r["T_low_K"] or 0.0)
    elif st == Status.UNDECIDED.value:
        t = -(r["T_up_K"] if r["T_up_K"] is not None else math.inf)
    else:
        t = -(r["T_up_K"] or 0.0)
    return (_RANK.get(st, 99), t, r["seq"])


def goal_rows(ledger, goal_id: str, dc=None) -> list[dict]:
    """Current rows of a goal, each with ``counts`` / ``not_counting_reason`` from the single
    predicate ``goals.counts_for`` (CURRENT config hashes).  A CERTIFIED row that does not count
    gets ``status_shown`` = "CERTIFIED_CANDIDATE (does not count)"."""
    from .goals import counts_for
    rows = ledger.current(goal_id)
    goal = dc.goal(goal_id) if dc is not None else None
    ch = dc.certify_hash(goal) if goal is not None else None
    dh = dc.dft_hash() if dc is not None else None
    for r in rows:
        if goal is not None:
            ok, why = counts_for(goal, r, certify_hash=ch, dft_hash=dh)
        else:
            ok, why = False, "no configuration given"
        r["counts_for_goal"] = ok
        r["not_counting_reason"] = "" if ok or r["status"] != "CERTIFIED_CANDIDATE" else why
        r["status_shown"] = NOT_COUNTING if (r["status"] == "CERTIFIED_CANDIDATE" and not ok) else r["status"]
        r["synthetic"] = "SYNTHETIC_DATA" in r["flags"] or r["material_id"].startswith("SYNTH-")
        r["round_generated"] = r["gen_round"]
        r["round_verdict"] = r["round"]
    return sorted(rows, key=_sort_key)


def _cell(v: Any) -> Any:
    if isinstance(v, (list, tuple)):
        return "; ".join(map(str, v))
    return "" if v is None else v


def _fmt(x: Any) -> str:
    if x is None:
        return "-"
    if isinstance(x, float):
        return f"{x:.4g}"
    return str(x)


def write_discovery_reports(ledger, dc, goals, workdir: Path) -> None:
    workdir = Path(workdir)
    rows_by_goal = {g.id: goal_rows(ledger, g.id, dc) for g in goals}
    for g in goals:
        with open(workdir / f"candidates_{g.id}.csv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=CSV_COLUMNS, extrasaction="ignore")
            w.writeheader()
            for r in rows_by_goal[g.id]:
                w.writerow({k: _cell(r.get("status_shown" if k == "status" else k)) for k in CSV_COLUMNS})
    synthetic = any(r["synthetic"] for rows in rows_by_goal.values() for r in rows) or \
        bool(dc.d.get("synthetic_spectra", {}).get("enabled"))
    try:
        acc = ledger.accounting([g.id for g in goals])
        acc_line = (f"{acc['generated']} generated structures = {acc['materials']} unique materials + "
                    f"{acc['duplicates']} duplicates (linked, never re-run); every material has exactly "
                    f"one current status per goal ({acc['current_pairs']} pairs). S4 checked.")
    except AssertionError as exc:
        acc_line = f"ACCOUNTING FAILED: {exc}"
    L = [f"# Discovery report: {dc.pipeline['run'].get('name', 'rtsc')}", "", BANNER, ""]
    if synthetic:
        L += [SYNTH_BANNER, ""]
    L += ["## Goals", "",
          "| goal | T* (K) | P range (GPa) | certified (counts) | certified, does not count | undecided | rejected | excluded | deferred | error | out of scope |",
          "|---|---|---|---|---|---|---|---|---|---|---|"]
    for g in goals:
        rows = rows_by_goal[g.id]
        inn = [r for r in rows if r["in_scope"]]
        c = {s.value: sum(1 for r in inn if r["status"] == s.value) for s in Status}
        ncount = sum(1 for r in rows if r["counts_for_goal"])
        L.append(f"| {g.id} | {g.T_star_K:g} | 0 .. {g.P_max_GPa:g} | {ncount} | "
                 f"{sum(1 for r in rows if r['status_shown'] == NOT_COUNTING)} | "
                 f"{c['UNDECIDED']} | {c['REJECTED_PROVEN']} | {c['EXCLUDED_ASSUMPTION']} | "
                 f"{c['DEFERRED']} | {c['ERROR']} | {len(rows) - len(inn)} |")
    L += ["", f"Accounting: {acc_line}", "", "## Progress", "",
          "| round | phase | requests | generated | new | duplicates | invalid | selected for DFT | notes |",
          "|---|---|---|---|---|---|---|---|---|"]
    for r in ledger.rounds():
        c = r["counts"]
        L.append(f"| {r['round']} | {r['phase']} | {len(r['requests'])} | {c.get('generated', '-')} | "
                 f"{c.get('new', '-')} | {c.get('duplicates', '-')} | {c.get('invalid', '-')} | "
                 f"{c.get('selected', '-')} | {'; '.join(r['notes'])[:200]} |")
    for g in goals:
        rows = rows_by_goal[g.id]
        cert = [r for r in rows if r["counts_for_goal"]]
        nc = [r for r in rows if r["status_shown"] == NOT_COUNTING]
        und = [r for r in rows if r["status"] == Status.UNDECIDED.value]
        rej = [r for r in rows if r["status"] == Status.REJECTED_PROVEN.value]
        L += ["", f"## {g.id}: certified candidates (Tc_model >= {g.T_star_K:g} K proved), by certified T_low", ""]
        L += _table(cert[:25], g) if cert else ["(none yet)"]
        if nc:
            L += ["", f"## {g.id}: CERTIFIED verdicts that do NOT count for this goal (stale config, other "
                  "T*, pressure)", ""]
            L += [f"- {r['material_id']} ({r['formula']}, {r['pressure_GPa']:g} GPa, T*={r['T_star_K']:g} K): "
                  f"does not count: {r['not_counting_reason']}" for r in nc[:50]]
        L += ["", f"## {g.id}: open (UNDECIDED), ranked by certified T_up", ""]
        L += _table(und[:25], g) if und else ["(none)"]
        L += ["", f"## {g.id}: rejected (Tc_model < {g.T_star_K:g} K proved), with certificates", ""]
        L += _table(rej[:25], g) if rej else ["(none)"]
    L += ["", f"Per-goal CSV: candidates_<goal>.csv.  Ledger: {workdir / 'ledger.sqlite'}.  "
          "Per-round pipeline outputs: rounds/<k>/<goal>/ (summary.md, candidates.csv, certificates/)."]
    (workdir / "discovery_report.md").write_text("\n".join(L) + "\n")


def _table(rows: list[dict], g) -> list[str]:
    out = ["| material | formula | P (GPa) | T_low (K) | T_up (K) | rules | flags | certificate |",
           "|---|---|---|---|---|---|---|---|"]
    for r in rows:
        out.append(f"| {r['material_id']} | {r['formula']} | {r['pressure_GPa']:g} | {_fmt(r['T_low_K'])} | "
                   f"{_fmt(r['T_up_K'])} | {','.join(r['rule_ids'])} | {','.join(r['flags'])} | "
                   f"{r['certificate_path'] or '-'} |")
    return out


# ---------------------------------------------------------------------------------- export
def export_json(ledger, goals, out: Path) -> Path:
    data = {"materials": [], "duplicates": [dict(r) for r in ledger.q("SELECT * FROM duplicates ORDER BY id")],
            "rounds": ledger.rounds(), "accounting": ledger.accounting([g.id for g in goals])}
    for m in ledger.materials():
        d = dict(m)
        d["structure"] = json.loads(d["structure"])
        d["request"] = json.loads(d["request"])
        d["provenance"] = json.loads(d["provenance"])
        d["current"] = {g.id: ledger.current_one(m["material_id"], g.id) for g in goals}
        d["ml_screen"] = ledger.ml_screen(m["material_id"])
        data["materials"].append(d)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(data, indent=1, default=str))
    return out


def _cif(s: dict, name: str) -> str:
    import numpy as np
    L = np.asarray(s["lattice"], float)
    a, b, c = (float(np.linalg.norm(v)) for v in L)
    ang = lambda u, v: math.degrees(math.acos(float(np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v)))))
    al, be, ga = ang(L[1], L[2]), ang(L[0], L[2]), ang(L[0], L[1])
    lines = [f"data_{name}", "_symmetry_space_group_name_H-M   'P 1'", "_symmetry_Int_Tables_number   1",
             f"_cell_length_a   {a:.8f}", f"_cell_length_b   {b:.8f}", f"_cell_length_c   {c:.8f}",
             f"_cell_angle_alpha   {al:.8f}", f"_cell_angle_beta   {be:.8f}", f"_cell_angle_gamma   {ga:.8f}",
             "loop_", " _symmetry_equiv_pos_as_xyz", "  'x, y, z'", "loop_", " _atom_site_type_symbol",
             " _atom_site_label", " _atom_site_occupancy", " _atom_site_fract_x", " _atom_site_fract_y",
             " _atom_site_fract_z"]
    for i, (sp, f) in enumerate(zip(s["species"], s["frac_coords"])):
        lines.append(f"  {sp}  {sp}{i}  1  {f[0]:.8f}  {f[1]:.8f}  {f[2]:.8f}")
    return "\n".join(lines) + "\n"


def export_finetune(ledger, goals, outdir: Path, *, val_fraction: float = 0.1,
                    include_synthetic: bool = False) -> dict[str, Any]:
    """DFT-level labels for generator fine-tuning (DISCOVERY.md section 4).  Documented, not run.

    Writes ``train.csv`` / ``val.csv`` (deterministic split by fingerprint hash) with a ``cif`` column
    (P1 CIF of the structure the pipeline computed) and property columns: pressure_GPa, A_K
    (min/max over alternative spectra), lambda, omega_log_K, and per goal the current status and
    certified bracket.  This is the CSV layout (``material_id``, ``cif``, property columns) used by
    MatterGen's ``datasets/*/{train,val}.csv``; property names must be registered in your MatterGen
    installation before conditioning on them (verify against the installed version).  Only
    materials whose spectra came from the pipeline are exported; SYNTHETIC spectra are excluded
    unless ``include_synthetic``.  ``labels.jsonl`` carries the full per-goal verdict records."""
    from ..pipeline import read_record
    from ..spectrum import Spectrum
    outdir.mkdir(parents=True, exist_ok=True)
    rows, labels = [], []
    for m in ledger.materials():
        mid = m["material_id"]
        cur = {g.id: ledger.current_one(mid, g.id) for g in goals}
        rec = None
        for g in goals:
            c = cur[g.id]
            if c and c["stage"] == "pipeline" and c["record_path"]:
                p = Path(c["record_path"])
                rec = read_record(p.parent.parent, mid)
                if rec and rec.get("spectra"):
                    break
        if not rec or not rec.get("spectra"):
            continue
        synth = bool(rec.get("candidate", {}).get("synthetic"))
        if synth and not include_synthetic:
            continue
        sp = [Spectrum.from_dict(d) for d in rec["spectra"]]
        As = [s.A() for s in sp]
        lams = [s.lambda_() for s in sp]
        wl = [s.omega_log() for s in sp]
        cand = rec["candidate"]
        s = {"lattice": cand["lattice"], "species": cand["species"], "frac_coords": cand["frac_coords"]}
        row = {"material_id": mid, "cif": _cif(s, mid.replace("-", "_")), "chemical_system": m["chemsys"],
               "pressure_GPa": m["pressure_GPa"], "A_K_min": min(As), "A_K_max": max(As),
               "lambda_min": min(lams), "lambda_max": max(lams), "omega_log_K_min": min(wl),
               "omega_log_K_max": max(wl), "synthetic": synth}
        for g in goals:
            c = cur[g.id] or {}
            row[f"status_{g.id}"] = c.get("status")
            row[f"T_low_K_{g.id}"] = c.get("T_low_K")
            row[f"T_up_K_{g.id}"] = c.get("T_up_K")
        rows.append(row)
        labels.append({"material_id": mid, "structure": s, "pressure_GPa": m["pressure_GPa"],
                       "A_K": As, "lambda": lams, "omega_log_K": wl, "synthetic": synth,
                       "verdicts": cur})
    cols = list(rows[0]) if rows else ["material_id", "cif"]
    split = {"train": [], "val": []}
    for r in rows:
        h = int(hashlib.sha256(r["material_id"].encode()).hexdigest()[:8], 16) / 0xFFFFFFFF
        split["val" if h < val_fraction else "train"].append(r)
    for name, rs in split.items():
        with open(outdir / f"{name}.csv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols)
            w.writeheader()
            w.writerows(rs)
    with open(outdir / "labels.jsonl", "w") as fh:
        for x in labels:
            fh.write(json.dumps(x, default=str) + "\n")
    (outdir / "README.txt").write_text(
        "DFT-level labels exported by `rtsc discover export --finetune`.\n"
        "train.csv / val.csv: material_id, cif (P1), property columns (see rtsc.discovery.report."
        "export_finetune).\nStatuses are model statements (isotropic Migdal-Eliashberg); T_low/T_up are "
        "certified brackets on Tc_model.\nTo fine-tune MatterGen, register the property columns you "
        "condition on in your MatterGen installation and convert the CSVs with its dataset tooling; "
        "this export is not run or validated against MatterGen here.\n")
    return {"n": len(rows), "train": len(split["train"]), "val": len(split["val"]), "outdir": str(outdir)}
