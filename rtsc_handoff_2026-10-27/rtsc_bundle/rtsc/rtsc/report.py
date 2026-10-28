"""Reports (ARCH section 6, stage 9).

Writes into ``[report] outdir`` (default: the run's workdir):

* ``candidates.csv`` / ``candidates.json``: one row per input, sorted by status
  (CERTIFIED_CANDIDATE, UNDECIDED, REJECTED_PROVEN, EXCLUDED_ASSUMPTION, DEFERRED, ERROR), then by
  certified T_low (descending, missing last), then id;
* ``certificates/<id>.json``: the Verdict, its origin, the spectra and the pipeline context
  (filters, priority, validity details, config hash); ``<workdir>/<id>/verdict.json`` holds the
  bare Verdict for ``rtsc recheck``;
* ``summary.md``: counts (accounting for every input, S4), the model assumptions in plain words
  (S6), the top candidates, exclusion criteria, deferrals, errors and flags;
* ``run_summary.json``: machine-readable counts and per-status id lists.
"""
from __future__ import annotations

import csv
import json
import math
import re
from collections import Counter
from pathlib import Path
from typing import Any, Sequence

from .sources.base import safe_id
from .status import Status

STATUS_ORDER = [Status.CERTIFIED_CANDIDATE, Status.UNDECIDED, Status.REJECTED_PROVEN,
                Status.EXCLUDED_ASSUMPTION, Status.DEFERRED, Status.ERROR]
_RANK = {s.value: i for i, s in enumerate(STATUS_ORDER)}

CSV_COLUMNS = ["material_id", "formula", "source", "synthetic", "status", "T_low_K", "T_up_K",
               "bracket_model",
               "T_target_K", "mu_hi", "rel_envelope", "cutoff", "n_spectra", "A_max_K",
               "lambda_min", "lambda_max", "dropped_nonpos_mass_K", "pressure_GPa", "priority", "proof_basis", "rule_ids",
               "assumptions", "flags", "reasons", "verdict_origin", "record"]

#: assumption markers that mean "this step has a hand proof only / is pending audit"
NON_LEAN_MARKERS = ("NOT in Lean", "hand proof", "pending audit")

MODEL_STATEMENT = """\
**What a verdict is about.** Every REJECTED_PROVEN / CERTIFIED_CANDIDATE / UNDECIDED verdict is a
statement about a MODEL, not directly about the material:

* isotropic Migdal-Eliashberg theory, linearised gap equation (Tc_model := sup of unstable T,
  which is >= the usual highest-solution Tc).  With the default rejection_model = "isotropic"
  every rejection uses the DOS-averaged single-band alpha^2F (assumption M-ISO-SINGLEBAND, flag
  MULTIBAND_RISK): in a multiband material the averaged A can underestimate Tc;
* the given alpha^2F (from DFPT or an attached file), within the stated relative weight
  envelope e (true atom weights in [(1-e)w, (1+e)w]); frequency uncertainty enters only through
  the alternative spectra (e.g. DFPT smearings), all of which must agree (S5);
* Coulomb pseudopotential mu* in [0, mu_hi] with the block-form cutoff
  N_c(T) = #{n >= 0 : (2n+1) pi T <= omega_c} (CERTIFIED depends on this rule and omega_c;
  REJECTED holds for every cutoff and every mu* >= 0).

**What a verdict is not.** It is not a prediction of the measured Tc.  Out-of-model physics -
anharmonicity and quantum nuclear motion (important for hydrogen), non-adiabatic effects,
anisotropic or multiband gaps, spin fluctuations, vertex corrections, the accuracy of DFT itself,
thermodynamic and kinetic stability - is NOT covered.  Such risks appear as flags; flags never
change a status.  EXCLUDED_ASSUMPTION verdicts rest on physical, non-proved criteria (hull,
metallicity, magnetism, imaginary modes) and DEFERRED means "not computed".  REJECTED_PROVEN and
CERTIFIED_CANDIDATE are produced only by rtsc.certify.classify, from rigorous ball-arithmetic
checks of rules proved in Lean (library PipelineLogic; docs/CORRESPONDENCE_LEAN.md) -- except
where a verdict's assumptions say otherwise (e.g. row_check = "cform": hand proof, not in Lean);
such verdicts are listed separately below and marked in the proof_basis column.
"""


def _f(x: Any) -> float | None:
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


def _spectra_stats(rec: dict[str, Any]) -> dict[str, Any]:
    sp = rec.get("spectra") or []
    if not sp:
        return {"n_spectra": 0, "A_max_K": None, "lambda_min": None, "lambda_max": None,
                "rel_envelope": None}
    As, lams = [], []
    for d in sp:
        w = [float(x) for x in d["weight_K"]]
        o = [float(x) for x in d["omega_K"]]
        As.append(math.fsum(w))
        lams.append(2.0 * math.fsum(wi / oi for wi, oi in zip(w, o)))
    return {"n_spectra": len(sp), "A_max_K": max(As), "lambda_min": min(lams),
            "lambda_max": max(lams), "rel_envelope": max(float(d.get("rel_envelope", 0)) for d in sp)}


def _is_synthetic(rec: dict[str, Any], v: dict[str, Any]) -> bool:
    cand = rec.get("candidate", {})
    return bool(cand.get("synthetic") or "SYNTHETIC_DATA" in v.get("flags", [])
                or any((d.get("provenance") or {}).get("synthetic") for d in rec.get("spectra") or []))


def proof_basis(v: dict[str, Any]) -> str:
    if v["status"] not in ("REJECTED_PROVEN", "CERTIFIED_CANDIDATE", "UNDECIDED"):
        return "-" if v["status"] != "EXCLUDED_ASSUMPTION" else "physical criterion (not proved)"
    hand = [a for a in v.get("assumptions", []) if any(m in a for m in NON_LEAN_MARKERS)]
    return "Arb + hand proof (NOT all in Lean)" if hand else "Arb + Lean"


def _dropped_mass(rec: dict[str, Any]) -> float | None:
    """Largest alpha^2F weight (K) at omega <= 0 dropped from any alternative (team DFT
    provenance 'dropped_nonpositive_omega' or the A-DROP-NONPOS tag)."""
    vals = []
    for d in rec.get("spectra") or []:
        pv = d.get("provenance") or {}
        dn = pv.get("dropped_nonpositive_omega")
        if isinstance(dn, dict) and dn.get("signed_mass_K") is not None:
            vals.append(float(dn["signed_mass_K"]))
        for a in pv.get("assumptions", []) or []:
            m = re.match(r"^A-DROP-NONPOS\(([-+0-9.eE]+)K\)", str(a))
            if m:
                vals.append(float(m.group(1)))
    return max(vals) if vals else None


def bracket_model(v: dict[str, Any]) -> str:
    """Which model the [T_low, T_up] bracket refers to (Lean note MB-N1)."""
    if v.get("T_low_K") is None and v.get("T_up_K") is None:
        return "-"
    if "R-MULTI" in v.get("rule_ids", []):
        return "multiband (R-MULTI upper bound)"
    rm = ((v.get("certificate") or {}).get("rejection_model")
          or ((v.get("certificate") or {}).get("params") or {}).get("rejection_model"))
    return "single-band (DOS-averaged)" + (" -- multiband mode" if rm == "multiband" else "")


def row_for(rec: dict[str, Any], cfg: dict[str, Any], workdir: Path) -> dict[str, Any]:
    v = rec["verdict"]
    cand = rec.get("candidate", {})
    st = _spectra_stats(rec)
    cert = cfg.get("certify", {})
    flags = list(v.get("flags", []))
    # T* and mu_hi come from what certify actually used (audit F2); disagreement is flagged
    used = (v.get("certificate") or {}).get("params") or {}
    T_used, mu_used = used.get("T_star"), used.get("mu_hi")
    T_cfg, mu_cfg = cert.get("T_target_K"), cert.get("mu_hi")
    if (T_used is not None and T_cfg is not None and float(T_used) != float(T_cfg)) or \
            (mu_used is not None and mu_cfg is not None and float(mu_used) != float(mu_cfg)):
        flags.append("PARAMS_MISMATCH")
    cutoff = next((a for a in v.get("assumptions", []) if "omega_c" in a), "")
    prio = rec.get("stages", {}).get("priority", {}).get("total")
    return {
        "material_id": rec["material_id"], "formula": cand.get("formula", ""),
        "source": cand.get("source", ""), "synthetic": _is_synthetic(rec, v),
        "status": v["status"], "T_low_K": _f(v.get("T_low_K")), "T_up_K": _f(v.get("T_up_K")),
        "bracket_model": bracket_model(v),
        "T_target_K": T_used if T_used is not None else T_cfg,
        "mu_hi": mu_used if mu_used is not None else mu_cfg,
        "rel_envelope": st["rel_envelope"], "cutoff": cutoff, "n_spectra": st["n_spectra"],
        "A_max_K": st["A_max_K"], "lambda_min": st["lambda_min"], "lambda_max": st["lambda_max"],
        "dropped_nonpos_mass_K": _dropped_mass(rec),
        "pressure_GPa": cand.get("properties", {}).get("target_pressure_GPa",
                                                       cand.get("pressure_GPa")),
        "priority": prio, "proof_basis": proof_basis(v), "rule_ids": list(v.get("rule_ids", [])),
        "assumptions": list(v.get("assumptions", [])), "flags": flags,
        "reasons": list(v.get("reasons", [])), "verdict_origin": rec.get("verdict_origin"),
        "record": str(workdir / safe_id(rec["material_id"]) / "record.json"),
    }


def sort_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return sorted(rows, key=lambda r: (_RANK.get(r["status"], 99),
                                       -(r["T_low_K"] if r["T_low_K"] is not None else -math.inf),
                                       r["material_id"]))


def _csv_value(v: Any) -> Any:
    if isinstance(v, list):
        return "; ".join(str(x) for x in v)
    return "" if v is None else v


def write_reports(cfg: dict[str, Any], summary, records: Sequence[dict[str, Any]]) -> Path:
    workdir = Path(cfg["run"]["workdir"])
    outdir = Path(cfg.get("report", {}).get("outdir") or workdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rows = sort_rows([row_for(r, cfg, workdir) for r in records if r.get("verdict")])
    with open(outdir / "candidates.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=CSV_COLUMNS)
        w.writeheader()
        for r in rows:
            w.writerow({k: _csv_value(r[k]) for k in CSV_COLUMNS})
    (outdir / "candidates.json").write_text(json.dumps(rows, indent=1, default=str))
    cdir = outdir / "certificates"
    cdir.mkdir(exist_ok=True)
    for rec in records:
        if not rec.get("verdict"):
            continue
        sid = safe_id(rec["material_id"])
        cert = {"material_id": rec["material_id"], "verdict": rec["verdict"],
                "verdict_origin": rec.get("verdict_origin"),
                "synthetic": bool(rec.get("candidate", {}).get("synthetic", False)),
                "candidate": {k: rec.get("candidate", {}).get(k) for k in
                              ("formula", "composition", "source", "pressure_GPa", "spacegroup",
                               "provenance", "band_coupling_files")},
                "spectra": rec.get("spectra", []),
                "pipeline_context": {"config_hash": rec.get("config_hash"),
                                     "stages": {k: v for k, v in rec.get("stages", {}).items()
                                                if k != "error"}},
                "recheck": f"rtsc recheck {workdir / sid / 'verdict.json'}"}
        (cdir / f"{sid}.json").write_text(json.dumps(cert, indent=1, default=str))
        vpath = workdir / sid / "verdict.json"
        vpath.parent.mkdir(parents=True, exist_ok=True)
        vpath.write_text(json.dumps(rec["verdict"], indent=1, default=str))
    (outdir / "run_summary.json").write_text(json.dumps(summary.to_dict(), indent=1, default=str))
    (outdir / "summary.md").write_text(summary_markdown(cfg, summary, rows, records))
    return outdir


def _fmt(x: Any) -> str:
    """Full precision for bounds: 299.9999999997 must not print as 300.0."""
    return "-" if x is None else (repr(x) if isinstance(x, float) else str(x))


def summary_markdown(cfg: dict[str, Any], summary, rows: list[dict[str, Any]],
                     records: Sequence[dict[str, Any]]) -> str:
    cert = cfg.get("certify", {})
    counts = summary.counts()
    n = len(summary.inputs)
    L: list[str] = []
    L.append(f"# rtsc run summary: {cfg.get('run', {}).get('name', '')}\n")
    n_syn = sum(1 for r in rows if r["synthetic"])
    if n_syn:
        L.append(f"> **SYNTHETIC DATA**: {n_syn} of {n} inputs use invented (synthetic) spectra. "
                 "Their verdicts exercise the software and say nothing about any real material.\n")
    L.append(f"- config hash: `{summary.config_hash}`; finished {summary.finished_utc}")
    L.append(f"- target T* = {cert.get('T_target_K')} K; mu* in [0, {cert.get('mu_hi')}]; "
             f"weight envelope e >= {cfg.get('spectra', {}).get('rel_envelope')}")
    from .pipeline import cutoff_description
    L.append(f"- Coulomb cutoff: {cutoff_description(cert)}")
    L.append(f"- rejection model: {cert.get('rejection_model', 'isotropic')} "
             "(isotropic: rejections assume M-ISO-SINGLEBAND)")
    L.append(f"- inputs: **{n}**, each with exactly one status (checked):\n")
    L.append("| status | count |\n|---|---|")
    for s in STATUS_ORDER:
        L.append(f"| {s.value} | {counts.get(s.value, 0)} |")
    L.append(f"| **total** | **{sum(counts.values())}** |\n")
    L.append("## What the verdicts mean (and do not mean)\n")
    L.append(MODEL_STATEMENT)
    L.append("| status | meaning |\n|---|---|")
    L.append("| CERTIFIED_CANDIDATE | proved: Tc_model >= T* for every spectrum in the envelope, "
             "every alternative spectrum and every mu* <= mu_hi (with the stated cutoff) |")
    L.append("| REJECTED_PROVEN | proved: Tc_model < T* for every spectrum in the envelope, every "
             "alternative and every mu* >= 0 |")
    L.append("| UNDECIDED | the certified bracket [T_low, T_up] contains T*, or alternatives disagree |")
    L.append("| EXCLUDED_ASSUMPTION | removed by a physical, non-proved criterion (recorded) |")
    L.append("| DEFERRED | not computed (budget, ordering, duplicate, no spectrum) |")
    L.append("| ERROR | a stage failed; see the record |\n")
    top_n = int(cfg.get("report", {}).get("top_n", 25))
    good = [r for r in rows if r["status"] in ("CERTIFIED_CANDIDATE", "UNDECIDED")][:top_n]
    if good:
        L.append("## Top candidates (certified first, then undecided; by certified T_low)\n")
        L.append("Brackets are single-band (DOS-averaged) unless the bracket column says otherwise.\n")
        L.append("| material | formula | status | T_low [K] | T_up [K] | bracket | A_max [K] | lambda | flags |")
        L.append("|---|---|---|---|---|---|---|---|---|")
        for r in good:
            lam = (f"{r['lambda_min']:.2f}-{r['lambda_max']:.2f}" if r["lambda_min"] is not None
                   else "-")
            L.append(f"| {r['material_id']} | {r['formula']} | {r['status']} | {_fmt(r['T_low_K'])} "
                     f"| {_fmt(r['T_up_K'])} | {r['bracket_model']} | {r['A_max_K']:.6g} | {lam} | "
                     f"{', '.join(r['flags']) or '-'} |")
        L.append("")
    rej = Counter(rid for r in rows if r["status"] == "REJECTED_PROVEN" for rid in r["rule_ids"])
    if rej:
        L.append("## Rejections by rule\n")
        L.extend(f"- {k}: {v}" for k, v in rej.most_common())
        n_sb = sum(1 for r in rows if r["status"] == "REJECTED_PROVEN"
                   and any(a.startswith("M-ISO-SINGLEBAND") for a in r["assumptions"]))
        if n_sb:
            L.append(f"\n{n_sb} rejection(s) hold in the single-band (DOS-averaged) model only "
                     "(M-ISO-SINGLEBAND); see rejection_model = \"multiband\".")
        L.append("")
    hand = [r for r in rows if "hand proof" in r["proof_basis"]]
    if hand:
        L.append("## Verdicts NOT fully backed by Lean\n")
        L.extend(f"- {r['material_id']} ({r['status']}): "
                 + "; ".join(a for a in r["assumptions"] if any(m in a for m in NON_LEAN_MARKERS))
                 for r in hand)
        L.append("")
    mism = [r for r in rows if "PARAMS_MISMATCH" in r["flags"]]
    if mism:
        L.append("## Parameter mismatches\n")
        L.append("These verdicts were computed with T*/mu_hi different from the current config; "
                 "the table shows the values actually used.\n")
        L.extend(f"- {r['material_id']}: T*={r['T_target_K']}, mu_hi={r['mu_hi']}" for r in mism)
        L.append("")
    exc = Counter(rid for r in rows if r["status"] == "EXCLUDED_ASSUMPTION" for rid in r["rule_ids"])
    if exc:
        L.append("## Exclusions by criterion (physical assumptions, not proofs)\n")
        L.extend(f"- {k}: {v}" for k, v in exc.most_common())
        thr = cfg.get("filters", {})
        L.append(f"\nThresholds: `{json.dumps(thr, default=str)}`\n")
    dfr = Counter((r["reasons"][0].split(":")[0] if r["reasons"] else "?")
                  for r in rows if r["status"] == "DEFERRED")
    if dfr:
        L.append("## Deferred (not computed)\n")
        L.extend(f"- {k}: {v}" for k, v in dfr.most_common())
        L.append("")
    errs = [r for r in rows if r["status"] == "ERROR"]
    if errs:
        L.append("## Errors\n")
        L.extend(f"- {r['material_id']}: {'; '.join(r['reasons'])[:300]}" for r in errs[:50])
        L.append("")
    fl = Counter(f for r in rows for f in r["flags"])
    if fl:
        L.append("## Out-of-model flags (never change a status)\n")
        L.extend(f"- {k}: {v}" for k, v in fl.most_common())
        L.append("")
    L.append("Files: `candidates.csv`, `candidates.json`, `certificates/<id>.json`, "
             "`<workdir>/<id>/record.json`, `<workdir>/<id>/verdict.json` (input of `rtsc recheck`).")
    return "\n".join(L) + "\n"
