#!/usr/bin/env python3
"""analyze_benchmark.py -- retrodiction-benchmark analysis for one material (HANDOFF.md section 5, step 2).

Input (positional INPUT):
  * a finished ``rtsc run`` material directory (contains record.json; verdict.json optional), or
  * a directory with QE matdyn output ``a2F.dos1..N`` plus the ``lambda`` file.

Output: <out>/<key>_benchmark.json and <key>_benchmark.md, plus every produced verdict
(<key>_verdict_<tag>.json) and the exact spectra it was computed on (<key>_spectra.json) so each
verdict can be re-verified independently with
    rtsc recheck <key>_verdict_<tag>.json --spectra <key>_spectra.json

What is PROVED and what is not
------------------------------
* Every status / bracket labelled "PROVED" comes from ``rtsc.pipeline.certify_spectra`` -- the
  exact function ``rtsc certify`` calls, which calls ``rtsc.certify.classify`` (ARCH S1) -- and is
  re-verified with ``rtsc.certify.recheck`` (what ``rtsc recheck -c`` calls).  This script contains
  NO decision logic: per-smearing brackets are read back from recheck's re-derived evidence.
* Allen-Dynes and the float Eliashberg solution below are FLOAT REFERENCES, NOT VERDICTS.  They are
  only compared with the proved bracket as a sanity check (a float Tc outside the certified bracket
  would be a serious finding and is reported loudly).
* Literature values are quoted from memory with the stated confidence; they are marked LITERATURE.

This script never modifies the rtsc bundle.
"""
from __future__ import annotations

import os

# at most the requested cores: one BLAS thread per process (set before numpy is imported)
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "FLINT_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import dataclasses
import datetime as _dt
import json
import math
import re
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import scipy.linalg as sla
from scipy.optimize import brentq

from rtsc import pipeline
from rtsc.certify import CertifyParams, recheck
from rtsc.spectrum import Spectrum
from rtsc.status import Verdict

SMOKE_BANNER = "SMOKE TEST, unconverged, no physical meaning"
FLOAT_LABEL = "float reference, not a verdict"

# ============================================================================ literature
# LITERATURE values, quoted from memory.  "confidence" says how sure we are of each number; nothing
# marked 'uncertain' should be quoted without checking the cited paper.
LITERATURE: dict[str, dict] = {
    "Al": {
        "Tc_exp_K": 1.175, "Tc_used_K": 1.2, "pressure_GPa": 0.0,
        "Tc_exp_ref": "standard tabulated value (e.g. Kittel; CRC Handbook): 1.175 K",
        "lambda": "0.38-0.44 (DFPT, e.g. Savrasov & Savrasov, PRB 54, 16487 (1996): ~0.44); "
                  "~0.38-0.43 from specific heat / tunnelling (Grimvall, 'The electron-phonon "
                  "interaction in metals', 1981)",
        "omega_log_K": "~270-320 K (uncertain; typical DFPT values)",
        "notes": "Weak coupling: Tc is exponentially sensitive to lambda and mu*; a 10% change in "
                 "lambda moves Tc by a factor ~2. mu* ~0.10-0.13 commonly fitted.",
        "confidence": "Tc: high. lambda range: medium. omega_log: low/uncertain.",
    },
    "Pb": {
        "Tc_exp_K": 7.19, "Tc_used_K": 7.2, "pressure_GPa": 0.0,
        "Tc_exp_ref": "standard tabulated value: 7.19-7.20 K",
        "lambda": "~1.55 from tunnelling inversion (McMillan & Rowell, PRL 14, 108 (1965); "
                  "Carbotte, RMP 62, 1027 (1990) table: 1.55); DFPT ~1.1-1.2 without spin-orbit "
                  "coupling and ~1.5-1.6 with SOC (Heid, Bohnen, Sklyadneva, Chulkov, "
                  "PRB 81, 174527 (2010))",
        "omega_log_K": "~56 K (Carbotte 1990 table: 4.83 meV; medium confidence)",
        "notes": "SOC matters for lambda: a scalar-relativistic DFPT run is expected to give the "
                 "smaller lambda and a lower Tc_model.",
        "confidence": "Tc: high. tunnelling lambda 1.55: high. DFPT no-SOC value: medium. "
                      "omega_log: medium.",
    },
    "Nb": {
        "Tc_exp_K": 9.25, "Tc_used_K": 9.3, "pressure_GPa": 0.0,
        "Tc_exp_ref": "standard tabulated value: 9.25-9.26 K",
        "lambda": "~1.0-1.3 (DFPT; e.g. Savrasov & Savrasov 1996 report ~1.26; tunnelling/"
                  "specific-heat estimates ~1.0)",
        "omega_log_K": "~150-200 K (uncertain)",
        "notes": "Spin fluctuations are thought to suppress Tc relative to the phonon-only "
                 "model (out of model; would appear as experiment < Tc_model).",
        "confidence": "Tc: high. lambda range: medium. omega_log: low/uncertain.",
    },
    "MgB2": {
        "Tc_exp_K": 39.0, "Tc_used_K": 39.0, "pressure_GPa": 0.0,
        "Tc_exp_ref": "Nagamatsu et al., Nature 410, 63 (2001): 39 K",
        "lambda": "isotropic (DOS-averaged) lambda ~0.7-0.9 (e.g. Kong, Dolgov, Jepsen, Andersen, "
                  "PRB 64, 020501(R) (2001): ~0.87; Choi, Roundy, Sun, Cohen, Louie, Nature 418, "
                  "758 (2002): ~0.73 with anharmonic phonons)",
        "omega_log_K": "~60-75 meV, i.e. ~700-870 K (uncertain)",
        "notes": "The canonical multiband case: isotropic single-band Eliashberg with mu*~0.1-0.15 "
                 "gives Tc roughly 15-25 K too low (order of magnitude from Choi et al. 2002 and "
                 "Golubov et al., J. Phys.: Condens. Matter 14, 1353 (2002)); two-band/anisotropic "
                 "treatments recover ~39 K. Experiment ABOVE the single-band T_up is therefore the "
                 "EXPECTED outcome here, not a pipeline failure. Expect MULTIBAND_RISK and "
                 "M-ISO-SINGLEBAND.",
        "confidence": "Tc: high. isotropic lambda range: medium. Size of the isotropic "
                      "underestimate: medium (order of magnitude). omega_log: low.",
    },
    "H3S": {
        "Tc_exp_K": 203.0, "Tc_used_K": 200.0, "pressure_GPa": 155.0,
        "Tc_exp_ref": "Drozdov, Eremets, Troyan, Ksenofontov, Shylin, Nature 525, 73 (2015): "
                      "~203 K at ~155 GPa",
        "lambda": "~2.0-2.2 harmonic DFPT for Im-3m H3S near 150-200 GPa (Duan et al., Sci. Rep. "
                  "4, 6968 (2014)); anharmonicity lowers it (~1.8; Errea et al., PRL 114, 157004 "
                  "(2015)) -- values uncertain",
        "omega_log_K": "~1300 K (harmonic; uncertain)",
        "notes": "Anharmonic and quantum-nuclear effects are large (flag ANHARMONIC_H expected); "
                 "the harmonic Tc_model may overestimate.",
        "confidence": "Tc/pressure: high. lambda, omega_log: medium/low.",
    },
}

ALIASES = {"al": "Al", "pb": "Pb", "nb": "Nb", "mgb2": "MgB2", "h3s": "H3S"}


# ============================================================================ helpers
def _now() -> str:
    return _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _fmt(x, nd=4):
    if x is None:
        return "-"
    if isinstance(x, float):
        if not math.isfinite(x):
            return str(x)
        return f"{x:.{nd}g}"
    return str(x)


def _pct(x):
    return "-" if x is None else f"{x:+.2%}"


def _smearing_index(path: Path) -> int:
    m = re.search(r"a2F\.dos(\d+)$", path.name)
    return int(m.group(1)) if m else 0


def find_qe_dir(root: Path) -> Path | None:
    """Directory holding a2F.dosN files and lambda (the dir itself or below it). Any subset of
    smearings is accepted, e.g. an explicit selection that excludes a2F.dos1."""
    if any(_smearing_index(p) > 0 for p in root.glob("a2F.dos*")):
        return root
    for p in sorted(root.rglob("a2F.dos*")):
        if _smearing_index(p) > 0:
            return p.parent
    return None


NEGATIVE_POLICY = "error"   # rtsc default; --negative-policy clip is an explicit, recorded opt-in


def load_qe_dir(qd: Path) -> list[Spectrum]:
    """Each a2F.dosN through rtsc.dfpt.load_attached_spectra (parse + verification), with the
    lambda file -- the same reader ``rtsc certify --format qe_a2f_dos`` and attached files use."""
    from rtsc.dfpt import load_attached_spectra
    files = sorted(qd.glob("a2F.dos*"), key=_smearing_index)
    files = [f for f in files if _smearing_index(f) > 0]
    if not files:
        raise FileNotFoundError(f"no a2F.dos* in {qd}")
    lf = qd / "lambda"
    out = []
    for f in files:
        kw = {"lambda_file": str(lf)} if lf.is_file() else {}
        kw["negative_policy"] = NEGATIVE_POLICY
        out.extend(load_attached_spectra(str(f), "qe_a2f_dos", None, label=f.name, **kw))
    return out


def qe_lambda_entries(qd: Path | None):
    if qd is None or not (qd / "lambda").is_file():
        return []
    from rtsc.dfpt.parsers import parse_qe_lambda_file
    return parse_qe_lambda_file(qd / "lambda")


def detect_smoke(root: Path, spectra) -> bool:
    if "smoke" in str(root).lower():
        return True
    for cand in (root / "PROVENANCE.txt", root.parent / "PROVENANCE.txt"):
        if cand.is_file() and "SMOKE TEST" in cand.read_text(errors="replace").upper():
            return True
    return False


def is_synthetic(spectra) -> bool:
    return any(s.provenance.get("synthetic") or "SYNTHETIC" in s.label.upper() for s in spectra)


# ============================================================================ float references
def allen_dynes(lam: float, wlog: float, w2: float, mu: float) -> dict:
    """Allen & Dynes, PRB 12, 905 (1975), with strong-coupling (f1) and shape (f2) factors.
    FLOAT REFERENCE, NOT A VERDICT."""
    den = lam - mu * (1 + 0.62 * lam)
    if den <= 0:
        return {"Tc_K": 0.0, "f1": None, "f2": None, "note": "lambda <= mu*(1+0.62 lambda)"}
    L1 = 2.46 * (1 + 3.8 * mu)
    L2 = 1.82 * (1 + 6.3 * mu) * (w2 / wlog)
    f1 = (1 + (lam / L1) ** 1.5) ** (1 / 3)
    f2 = 1 + (w2 / wlog - 1) * lam ** 2 / (lam ** 2 + L2 ** 2)
    tc_mcm = wlog / 1.2 * math.exp(-1.04 * (1 + lam) / den)
    return {"Tc_K": f1 * f2 * tc_mcm, "Tc_McMillan_wlog_K": tc_mcm, "f1": f1, "f2": f2}


def _lam_k(om: np.ndarray, p: np.ndarray, T: float, kmax: int) -> np.ndarray:
    """lambda(k; T) = sum_j p_j 2 w_j / (w_j^2 + (2 pi k T)^2), k = 0..kmax (float)."""
    k = np.arange(kmax + 1, dtype=float)
    x2 = (2 * math.pi * T * k) ** 2
    out = np.empty(kmax + 1)
    step = max(1, 400_000 // max(1, len(om)))
    for i0 in range(0, kmax + 1, step):
        i1 = min(kmax + 1, i0 + step)
        out[i0:i1] = (2 * p[None, :] * om[None, :] / (om[None, :] ** 2 + x2[i0:i1, None])).sum(1)
    return out


def eliashberg_eta(om, p, T, mu, omega_c, N=None, fac=12.0, cap=4000):
    """Largest eigenvalue eta(T) of the symmetrised isotropic linearised Eliashberg kernel
    (even gap, positive Matsubara frequencies w_n = pi T (2n+1), n < N):

        Z_n = 1 + [lambda(0) + 2 sum_{k=1..n} lambda(k)] / (2n+1)
        Z_n Delta_n = sum_m [lambda(n-m) + lambda(n+m+1) - 2 mu* theta(w_c - w_m)] Delta_m / (2m+1)

    symmetrised with c_n = (2n+1) Z_n:  S = K / sqrt(c c^T).  Tc is where eta = 1.
    Independent of rtsc.certify.floatcore (different formulation and code).  A truncated N gives a
    principal submatrix, so eta (and the resulting Tc) can only be underestimated.
    Returns (eta, N used, N_c, truncated?).  FLOAT REFERENCE, NOT A VERDICT."""
    wmax = float(om.max())
    nc = 0
    if mu > 0 and omega_c >= math.pi * T:
        nc = int(math.floor((omega_c / (math.pi * T) - 1) / 2)) + 1
    truncated = False
    if N is None:
        N = max(64, int(math.ceil(fac * wmax / (2 * math.pi * T))), nc + 32)
        if N > cap:
            N, truncated = cap, True
    lam = _lam_k(om, p, T, 2 * N + 1)
    n = np.arange(N)
    K = lam[np.abs(n[:, None] - n[None, :])] + lam[n[:, None] + n[None, :] + 1]
    if nc:
        c = min(nc, N)
        K[:c, :c] -= 2 * mu
    S1 = np.concatenate([[0.0], np.cumsum(lam[1:N])])     # S1[n] = sum_{k=1..n} lambda(k)
    cvec = (2 * n + 1) + lam[0] + 2 * S1[:N]               # (2n+1) Z_n
    d = 1.0 / np.sqrt(cvec)
    K *= d[:, None]
    K *= d[None, :]
    if N <= 900:
        w = sla.eigh(K, eigvals_only=True, subset_by_index=[N - 1, N - 1])
        return float(w[0]), N, nc, truncated
    from scipy.sparse.linalg import eigsh
    w = eigsh(K, k=1, which="LA", return_eigenvectors=False, tol=1e-10)
    return float(w[0]), N, nc, truncated


def eliashberg_tc(spec: Spectrum, mu: float, omega_c: float, T_min=0.01) -> dict:
    """Highest T with eta(T) = 1 (downward scan by 0.7, then Brent in log T); truncation checked
    by re-evaluating eta at Tc with 1.5x the Matsubara count.  FLOAT REFERENCE, NOT A VERDICT."""
    om = np.asarray(spec.omega_K, float)
    p = np.asarray(spec.weight_K, float)
    keep = p > 0
    om, p = om[keep], p[keep]
    t0 = time.time()
    trunc = {"any": False}

    def f(lt):
        e, _, _, tr = eliashberg_eta(om, p, math.exp(lt), mu, omega_c)
        trunc["any"] |= tr
        return e - 1.0
    T_hi = float(p.sum()) / 4.33148145 * 1.0001     # R-UNIV scale: stable above (Lean TcR1)
    lt = math.log(T_hi)
    if f(lt) >= 0:
        return {"Tc_K": T_hi, "note": "eta>=1 at the R-UNIV scale?! (float bug?)",
                "seconds": time.time() - t0}
    while True:
        lo = lt + math.log(0.7)
        if math.exp(lo) < T_min:
            return {"Tc_K": 0.0, "below_K": T_min, "truncated": trunc["any"],
                    "note": f"no eta>=1 found above {T_min} K"
                            + (" (Matsubara truncation hit: lower estimate only)" if trunc["any"] else ""),
                    "seconds": time.time() - t0}
        if f(lo) >= 0:
            break
        lt = lo
    root = brentq(f, lo, lt, xtol=1e-7, rtol=1e-9)
    Tc = math.exp(root)
    e1, N1, nc, tr1 = eliashberg_eta(om, p, Tc, mu, omega_c)
    out = {"Tc_K": Tc, "N_matsubara": N1, "N_c": nc, "truncated": tr1, "eta_at_Tc": e1,
           "seconds": None}
    if 1.5 * N1 <= 4500:
        e2, N2, _, _ = eliashberg_eta(om, p, Tc, mu, omega_c, N=int(1.5 * N1))
        out["eta_at_Tc_1p5N"] = e2
        out["truncation_shift_eta"] = e2 - e1
    out["seconds"] = time.time() - t0
    return out


def rtsc_float_tc(spec: Spectrum, mu: float, omega_c: float) -> float:
    """rtsc's own float proposal (rtsc.certify.floatcore.float_tc) -- shown for comparison only."""
    from rtsc.certify import FloatSpectrum, float_tc
    fs = FloatSpectrum(spec.omega_K, spec.weight_K, 0.0)
    return float(float_tc(fs, mu=mu, omega_c=omega_c))


# ============================================================================ proved part
def run_certify(job: dict) -> dict:
    """One call of the public certify entry point (the function ``rtsc certify`` uses), then
    ``recheck`` (the function ``rtsc recheck -c`` uses).  No decision logic here."""
    spectra = [Spectrum.from_dict(d) for d in job["spectra"]]
    ccfg_raw = {"T_target_K": job["T_star"], "mu_hi": job["mu_hi"], "cutoff_rule": "block",
                "mode": job["mode"], "rejection_model": "isotropic",
                "omega_c_factor": job["omega_c_factor"]}
    if job.get("time_budget_s") is not None:
        ccfg_raw["time_budget_s"] = job["time_budget_s"]
    cfg = pipeline.normalize_config({"certify": ccfg_raw,
                                     "spectra": {"rel_envelope": job["envelope"]}})
    ccfg = cfg["certify"]
    t0 = time.time()
    v = pipeline.certify_spectra(job["material_id"], spectra, ccfg)
    t_cert = time.time() - t0
    t1 = time.time()
    res = recheck(v, spectra, expected_params=pipeline._certify_params(ccfg))
    t_re = time.time() - t1
    per_alt = []
    for s, a, alt in zip(spectra, res.get("alternatives", []), v.certificate.get("alternatives", [])):
        per_alt.append({"label": s.label, "smearing_index": s.provenance.get("smearing_index"),
                        "T_low_K": a.get("T_low"), "T_up_K": a.get("T_up"),
                        "classifyAlt": a.get("classifyAlt"), "evidence": a.get("evidence"),
                        "certify_seconds": alt.get("seconds"),
                        "upper_search": alt.get("search_upper"),
                        "lower_search": alt.get("search_lower")})
    return {"tag": job["tag"], "T_star": job["T_star"], "mu_hi": job["mu_hi"], "mode": job["mode"],
            "envelope": job["envelope"], "time_budget_s": ccfg.get("time_budget_s", 120.0),
            "verdict": v.to_dict(), "status": v.status.value, "T_low_K": v.T_low_K,
            "T_up_K": v.T_up_K, "flags": list(v.flags), "assumptions": list(v.assumptions),
            "rule_ids": list(v.rule_ids), "per_alternative": per_alt,
            "omega_c_K": (v.certificate or {}).get("params", {}).get("omega_c"),
            "recheck": {"ok": bool(res.get("ok")), "status": res.get("status"),
                        "T_low_K": res.get("T_low_K"), "T_up_K": res.get("T_up_K"),
                        "errors": res.get("errors", []), "seconds": t_re},
            "wall_certify_s": t_cert, "wall_recheck_s": t_re}


# ============================================================================ main analysis
def load_input(root: Path):
    """-> (spectra, kind, qe_dir, run_record_or_None, run_verdict_or_None)"""
    rec = None
    run_verdict = None
    if (root / "record.json").is_file():
        rec = json.loads((root / "record.json").read_text())
        if (root / "verdict.json").is_file():
            run_verdict = json.loads((root / "verdict.json").read_text())
        if not rec.get("spectra"):
            raise SystemExit(f"{root}/record.json has no spectra (status "
                             f"{(rec.get('verdict') or {}).get('status')}); nothing to analyse")
        spectra = [Spectrum.from_dict(d) for d in rec["spectra"]]
        return spectra, "rtsc_run", find_qe_dir(root), rec, run_verdict
    qd = find_qe_dir(root)
    if qd is None:
        raise SystemExit(f"{root}: neither record.json nor a2F.dos1 found")
    return load_qe_dir(qd), "qe_dir", qd, None, None


def analyse(args) -> dict:
    root = Path(args.input).resolve()
    key = ALIASES.get(args.material.lower(), args.material)
    lit = LITERATURE.get(key)
    if lit is None and args.tc_lit is None:
        raise SystemExit(f"unknown material key {args.material!r}; known {sorted(LITERATURE)} "
                         "(or pass --tc-lit)")
    tc_lit = float(args.tc_lit if args.tc_lit is not None else lit["Tc_used_K"])
    spectra_raw, kind, qd, rec, run_verdict = load_input(root)
    smoke = args.smoke or detect_smoke(root, spectra_raw)
    synth = is_synthetic(spectra_raw)
    data_label = (SMOKE_BANNER if smoke else
                  "SYNTHETIC DATA, no physical meaning" if synth else
                  "real DFPT data (convergence as stated by the producer of the input)")
    # envelope: max(input's own, default 0.10) -- conservative, never narrowed (pipeline rule)
    spectra = pipeline._with_envelope(spectra_raw, float(args.envelope))
    out: dict = {"tool": "bench_tools/analyze_benchmark.py", "created_utc": _now(),
                 "material_key": key, "material_id": args.material_id or f"{key}-bench",
                 "input": str(root), "input_kind": kind, "qe_dir": str(qd) if qd else None,
                 "data_label": data_label, "smoke": smoke, "synthetic": synth,
                 "T_lit_K": tc_lit, "envelope_applied": sorted({s.rel_envelope for s in spectra}),
                 "n_smearings": len(spectra)}
    try:
        import flint
        out["trust"] = f"python-flint {flint.__version__}"
    except Exception:  # pragma: no cover
        pass

    # ---------------- 1. per-smearing spectrum metrics + QE lambda-file cross-check
    entries = qe_lambda_entries(qd)
    rows = []
    for s in spectra:
        isig = s.provenance.get("smearing_index")
        r = {"label": s.label, "smearing_index": isig,
             "broadening_Ry": s.provenance.get("el_ph_broadening_Ry"),
             "n_atoms": s.n_atoms(), "lambda": s.lambda_(), "omega_log_K": s.omega_log(),
             "omega2_K": s.omega2(), "A_K": s.A(), "omega_max_K": s.omega_max(),
             "qe_footer_lambda": s.provenance.get("qe_reported_lambda_histogram")}
        if entries and isig and isig <= len(entries):
            e = entries[isig - 1]
            r["qe_lambda_file"] = {"broadening_Ry": e.broadening_Ry, "lambda": e.lambda_,
                                   "omega_ln_K": e.omega_ln_K, "dos_Ef": getattr(e, "dos_ef", None)}
            r["d_lambda_rel"] = (r["lambda"] - e.lambda_) / e.lambda_ if e.lambda_ else None
            if e.omega_ln_K:
                r["d_omega_log_rel"] = (r["omega_log_K"] - e.omega_ln_K) / e.omega_ln_K
        rows.append(r)
    out["smearings"] = rows

    def spread(key_):
        v = [r[key_] for r in rows if r.get(key_) is not None and math.isfinite(r[key_])]
        if not v:
            return None
        mid = 0.5 * (max(v) + min(v))
        return {"min": min(v), "max": max(v), "rel_spread": (max(v) - min(v)) / mid if mid else None}
    out["spread"] = {k: spread(k) for k in ("lambda", "omega_log_K", "A_K", "omega_max_K")}

    # ---------------- 2. PROVED results via the public certify entry point
    sdicts = [s.to_dict() for s in spectra]
    base = {"spectra": sdicts, "material_id": out["material_id"], "envelope": float(args.envelope),
            "omega_c_factor": float(args.omega_c_factor), "time_budget_s": args.time_budget}
    jobs = []
    if not args.skip_300:
        jobs.append(dict(base, tag="T300", T_star=300.0, mu_hi=0.15, mode=args.mode_300))
    for mu in args.mu_list:
        jobs.append(dict(base, tag=f"Tlit_mu{mu:.2f}", T_star=tc_lit, mu_hi=mu, mode="bracket"))
    t0 = time.time()
    if args.jobs > 1 and len(jobs) > 1:
        with ProcessPoolExecutor(max_workers=args.jobs) as ex:
            results = list(ex.map(run_certify, jobs))
    else:
        results = [run_certify(j) for j in jobs]
    out["certify_wall_total_s"] = time.time() - t0
    outdir = Path(args.out).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / f"{key}_spectra.json").write_text(json.dumps(sdicts, indent=1))
    proved = {}
    for r in results:
        vp = outdir / f"{key}_verdict_{r['tag']}.json"
        vp.write_text(json.dumps(r["verdict"], indent=1, default=str))
        r["verdict_file"] = str(vp)
        r["recheck_cli"] = f"rtsc recheck {vp} --spectra {outdir / (key + '_spectra.json')}"
        r.pop("verdict")
        proved[r["tag"]] = r
    out["proved"] = proved

    # the run's own verdict (if any): recheck it too (Arb only)
    if run_verdict is not None:
        try:
            v = Verdict.from_dict(run_verdict)
            rspec = [Spectrum.from_dict(d) for d in rec["spectra"]]
            t1 = time.time()
            res = recheck(v, rspec)
            out["run_verdict"] = {"status": v.status.value, "T_low_K": v.T_low_K,
                                  "T_up_K": v.T_up_K, "flags": list(v.flags),
                                  "T_star": (v.certificate or {}).get("params", {}).get("T_star"),
                                  "mu_hi": (v.certificate or {}).get("params", {}).get("mu_hi"),
                                  "recheck_ok": bool(res.get("ok")),
                                  "recheck_errors": res.get("errors", []),
                                  "recheck_s": time.time() - t1}
        except Exception as exc:
            out["run_verdict"] = {"error": repr(exc)}

    # ---------------- 3. FLOAT references (not verdicts)
    omega_c = None
    for r in results:
        if r.get("omega_c_K"):
            omega_c = float(r["omega_c_K"])
            break
    if omega_c is None:
        omega_c = float(args.omega_c_factor) * max(s.omega_max() for s in spectra)
    fl = {"label": FLOAT_LABEL, "omega_c_K": omega_c,
          "note": "Eliashberg float uses the same block cutoff N_c(T) as the certificate "
                  "(omega_c = factor x max omega over all smearings); Allen-Dynes uses its own "
                  "implicit mu* convention, so equal mu* does not mean the same model.",
          "per_smearing": []}
    mus = sorted(set([0.0] + list(args.mu_list)))
    for s, row in zip(spectra, rows):
        e = {"label": s.label, "smearing_index": row["smearing_index"], "allen_dynes": {},
             "eliashberg_float": {}, "rtsc_floatcore": {}}
        for mu in mus:
            e["allen_dynes"][f"{mu:.2f}"] = allen_dynes(row["lambda"], row["omega_log_K"],
                                                        row["omega2_K"], mu)
            e["eliashberg_float"][f"{mu:.2f}"] = eliashberg_tc(s, mu, omega_c)
            if args.rtsc_float:
                t1 = time.time()
                e["rtsc_floatcore"][f"{mu:.2f}"] = {"Tc_K": rtsc_float_tc(s, mu, omega_c),
                                                     "seconds": time.time() - t1}
        fl["per_smearing"].append(e)
    out["float_reference"] = fl

    # ---------------- 5. consistency checks
    out["checks"] = consistency(out, key, lit)
    return out


def consistency(out: dict, key: str, lit: dict | None) -> dict:
    chk: dict = {"items": [], "SERIOUS": []}

    def add(name, ok, detail, serious=False):
        chk["items"].append({"check": name, "result": ok, "detail": detail})
        if serious and ok is False:
            chk["SERIOUS"].append(f"{name}: {detail}")

    proved = out["proved"]
    fl = out["float_reference"]["per_smearing"]
    T_exp = lit["Tc_exp_K"] if lit else out["T_lit_K"]
    # recheck of every produced verdict
    for tag, r in proved.items():
        add(f"recheck[{tag}]", r["recheck"]["ok"],
            f"stored {r['status']} [{r['T_low_K']}, {r['T_up_K']}] -> rechecked "
            f"{r['recheck']['status']} [{r['recheck']['T_low_K']}, {r['recheck']['T_up_K']}]"
            + (f"; errors {r['recheck']['errors']}" if r["recheck"]["errors"] else ""),
            serious=True)
        # F7 bookkeeping (read-only): combined == [min T_low, max T_up] over alternatives
        lows = [a["T_low_K"] for a in r["per_alternative"]]
        ups = [a["T_up_K"] for a in r["per_alternative"]]
        f7_low = min(lows) if lows and all(x is not None for x in lows) else None
        f7_up = max(ups) if ups and all(x is not None for x in ups) else None
        add(f"F7 combination[{tag}]", (f7_low == r["T_low_K"] or r["T_low_K"] is None) and
            (f7_up == r["T_up_K"] or r["T_up_K"] is None),
            f"verdict [{r['T_low_K']}, {r['T_up_K']}] vs [min T_low, max T_up] = [{f7_low}, {f7_up}]",
            serious=True)
    brk = {t: r for t, r in proved.items() if t.startswith("Tlit_")}
    # experiment vs T_up (T_up holds for every mu* >= 0)
    for tag, r in brk.items():
        if r["T_up_K"] is not None:
            ok = T_exp <= r["T_up_K"]
            add(f"T_exp <= combined T_up [{tag}]", ok,
                f"T_exp={T_exp} K, T_up={_fmt(r['T_up_K'])} K"
                + ("" if ok else " -> experiment ABOVE the model's proved ceiling"
                   + (" (EXPECTED for MgB2: single-band underestimate)" if key == "MgB2" else
                      " (model/DFPT inaccuracy or unconverged input; not a soundness issue)")))
    for tag, r in brk.items():
        if abs(r["mu_hi"] - 0.15) < 1e-12:
            ok = None if r["T_low_K"] is None else r["T_low_K"] <= T_exp
            add("T_low(mu*=0.15) <= T_exp", ok,
                f"T_low={_fmt(r['T_low_K'])} K, T_exp={T_exp} K")
    # float Tc inside the certified bracket (SERIOUS if violated)
    for tag, r in brk.items():
        mu = r["mu_hi"]
        for a, f in zip(r["per_alternative"], fl):
            fe = f["eliashberg_float"][f"{mu:.2f}"]
            tc_mu = fe["Tc_K"]
            tc_0 = f["eliashberg_float"]["0.00"]["Tc_K"]
            if a["T_low_K"] is not None:
                ok = tc_mu >= a["T_low_K"]
                tr = bool(fe.get("truncated"))
                add(f"float Tc(mu*={mu:.2f}) >= T_low [{tag}, {a['label']}]",
                    ok if (ok or not tr) else None,
                    f"float {_fmt(tc_mu, 7)} K vs proved T_low {_fmt(a['T_low_K'], 7)} K"
                    + (" (float Matsubara-truncated: lower estimate)" if tr else ""),
                    serious=not tr)
            else:
                ls = a.get("lower_search") or {}
                add(f"no proved T_low [{tag}, {a['label']}]", None,
                    f"rtsc lower search found no C-RAYLEIGH certificate; its float proposal at the "
                    f"envelope minimum (weights x(1-e), mu*={mu:.2f}) was "
                    f"{_fmt(ls.get('tc_float_lo'), 5)} K; certify_lower only searches T >= "
                    f"max(1e-3*T_hi, 0.05 K) (search.py), so Tc_model below ~0.05 K cannot get a "
                    f"lower bound (sound, uninformative). Nominal-weight float Tc: {_fmt(tc_mu, 5)} K")
            if a["T_up_K"] is not None:
                add(f"float Tc(mu*=0) <= T_up [{tag}, {a['label']}]",
                    tc_0 <= a["T_up_K"], f"float {_fmt(tc_0, 7)} K vs proved T_up "
                    f"{_fmt(a['T_up_K'], 7)} K", serious=True)
    # Allen-Dynes vs Eliashberg float (sanity of the float code)
    # mu* = 0 is the clean comparison (no cutoff-convention ambiguity): expect ~10-20 %.
    # For mu* > 0 the two use different mu* conventions (AD's fitted mu* vs our block cutoff at
    # omega_c = factor x max omega), and in weak coupling Tc is exponentially sensitive to that,
    # so the mu* > 0 differences are reported as information only.
    for f in fl:
        for mu in [m for m in f["allen_dynes"]]:
            ad = f["allen_dynes"].get(mu, {}).get("Tc_K")
            el = f["eliashberg_float"].get(mu, {}).get("Tc_K")
            if ad and el:
                rel = (el - ad) / ad
                clean = mu == "0.00"
                add(f"AD vs Eliashberg float [{f['label']}, mu*={mu}]",
                    (abs(rel) <= 0.20) if clean else None,
                    f"AD={_fmt(ad)} K, Eliashberg={_fmt(el)} K, rel diff {rel:+.1%}"
                    + ("" if clean else " (info: different mu* conventions)"))
    # rtsc floatcore vs independent solver (both floats)
    for f in fl:
        for mu, d in f.get("rtsc_floatcore", {}).items():
            el = f["eliashberg_float"][mu]["Tc_K"]
            tr = bool(f["eliashberg_float"][mu].get("truncated"))
            if el > 0:
                ok = abs(d["Tc_K"] - el) / el <= 1e-3
                add(f"independent float vs rtsc floatcore [{f['label']}, mu*={mu}]",
                    ok if (ok or not tr) else None,
                    f"{_fmt(el, 7)} K vs {_fmt(d['Tc_K'], 7)} K"
                    + (" (Matsubara truncation in both: ours N<=4000, floatcore cap 3000; truncated "
                       "values are lower estimates, so the larger one is closer)" if tr else ""))
    # QE lambda-file cross-check
    for r in out["smearings"]:
        if r.get("d_lambda_rel") is not None:
            add(f"lambda vs QE lambda file [{r['label']}]", abs(r["d_lambda_rel"]) <= 0.05,
                f"parsed {r['lambda']:.5g} vs file {r['qe_lambda_file']['lambda']:.5g} "
                f"({r['d_lambda_rel']:+.2%}; histogram vs direct q-sum, rtsc tolerance 5%)")
    sp = out["spread"]
    add("smearing spread", None,
        "; ".join(f"{k}: {_fmt(v['min'])}..{_fmt(v['max'])} ({v['rel_spread']:.1%})"
                  for k, v in sp.items() if v))
    for tag, r in brk.items():
        w = [(a["T_low_K"], a["T_up_K"]) for a in r["per_alternative"]]
        add(f"per-smearing brackets [{tag}]", None,
            ", ".join(f"[{_fmt(lo)}, {_fmt(up)}]" for lo, up in w))
    # MgB2 / multiband bookkeeping (and generally: report what is present)
    all_flags = sorted({f for r in proved.values() for f in r["flags"]})
    sb = {t: any(a.startswith("M-ISO-SINGLEBAND") for a in r["assumptions"]) for t, r in proved.items()}
    mbr = {t: "MULTIBAND_RISK" in r["flags"] for t, r in proved.items()}
    add("M-ISO-SINGLEBAND present per verdict", all(sb.values()) if sb else None, json.dumps(sb),
        serious=(key == "MgB2"))
    rej = {t: r["status"] for t, r in proved.items()}
    exp_mbr = all(mbr[t] for t, s in rej.items() if s == "REJECTED_PROVEN")
    add("MULTIBAND_RISK on every REJECTED_PROVEN verdict", exp_mbr,
        f"statuses {rej}; MULTIBAND_RISK {mbr} (classify sets it only on single-band rejections)",
        serious=True)
    if key == "MgB2":
        add("MgB2: MULTIBAND_RISK present on at least one verdict", any(mbr.values()),
            "expected on the T*=300 K rejection; UNDECIDED/CERTIFIED brackets carry only "
            "M-ISO-SINGLEBAND", serious=False)
    rv = out.get("run_verdict")
    if rv and "flags" in rv:
        add("run verdict flags (incl. validity stage)", None, ", ".join(rv["flags"]) or "(none)")
        add("run verdict recheck", rv.get("recheck_ok"), "; ".join(rv.get("recheck_errors", [])),
            serious=True)
    chk["flags_seen"] = all_flags
    return chk


# ============================================================================ markdown
def to_markdown(out: dict) -> str:
    key = out["material_key"]
    L = []
    banner = out["data_label"]
    L.append(f"## Benchmark: {key} ({out['material_id']})")
    L.append("")
    if out["smoke"] or out["synthetic"]:
        L.append(f"> **{banner.upper()}.** Every number below is program output on data that "
                 "carries no physical meaning for this material. Do not quote.")
    else:
        L.append(f"> Data: {banner}.")
    L.append("")
    L.append(f"Input: `{out['input']}` ({out['input_kind']}); smearings: {out['n_smearings']}; "
             f"envelope e = {out['envelope_applied']}; literature T* = {out['T_lit_K']} K; "
             f"generated {out['created_utc']} ({out.get('trust', '')}).")
    L.append("")
    L.append("### 1. Spectra per smearing (parsed Spectrum; QE lambda file for comparison)")
    L.append("")
    L.append("| # | broadening (Ry) | λ | ω_log (K) | A (K) | ω_max (K) | QE λ (lambda file) | Δλ | QE ω_ln (K) | Δω_log |")
    L.append("|---|---|---|---|---|---|---|---|---|---|")
    for r in out["smearings"]:
        q = r.get("qe_lambda_file") or {}
        L.append(f"| {r['smearing_index']} | {_fmt(r['broadening_Ry'])} | {r['lambda']:.4f} | "
                 f"{r['omega_log_K']:.1f} | {r['A_K']:.2f} | {r['omega_max_K']:.1f} | "
                 f"{_fmt(q.get('lambda'))} | {_pct(r.get('d_lambda_rel'))} | "
                 f"{_fmt(q.get('omega_ln_K'))} | {_pct(r.get('d_omega_log_rel'))} |")
    sp = out["spread"]
    L.append("")
    L.append("Spread over smearings: " + "; ".join(
        f"{k} {_fmt(v['min'])}–{_fmt(v['max'])} ({v['rel_spread']:.1%})" for k, v in sp.items() if v))
    L.append("")
    L.append("### 2. PROVED results (rtsc.pipeline.certify_spectra → rtsc.certify.classify; "
             "re-verified with rtsc.certify.recheck)")
    L.append("")
    L.append("Statements about Tc_model (isotropic Migdal–Eliashberg, linearised gap equation, block "
             "Coulomb cutoff, given α²F, weight envelope e). T_up holds for every μ* ≥ 0; T_low holds "
             "for every μ* ≤ mu_hi under the stated cutoff. Combined over smearings per ERRATA F7.")
    L.append("")
    L.append("| run | T* (K) | mu_hi | mode | status | combined T_low (K) | combined T_up (K) | recheck | certify wall (s) | recheck wall (s) |")
    L.append("|---|---|---|---|---|---|---|---|---|---|")
    for tag, r in out["proved"].items():
        L.append(f"| {tag} | {r['T_star']} | {r['mu_hi']} | {r['mode']} | **{r['status']}** | "
                 f"{_fmt(r['T_low_K'], 6)} | {_fmt(r['T_up_K'], 6)} | "
                 f"{'OK' if r['recheck']['ok'] else 'FAILED'} | {r['wall_certify_s']:.1f} | "
                 f"{r['wall_recheck_s']:.1f} |")
    L.append("")
    for tag, r in out["proved"].items():
        if not tag.startswith("Tlit_"):
            continue
        L.append(f"Per-smearing proved brackets, {tag} (from recheck's re-derived evidence):")
        L.append("")
        L.append("| # | T_low (K) | T_up (K) | classifyAlt | certify s |")
        L.append("|---|---|---|---|---|")
        for a in r["per_alternative"]:
            L.append(f"| {a['smearing_index']} | {_fmt(a['T_low_K'], 6)} | {_fmt(a['T_up_K'], 6)} | "
                     f"{a['classifyAlt']} | {_fmt(a['certify_seconds'], 3)} |")
        L.append("")
    t300 = out["proved"].get("T300")
    if t300:
        L.append(f"T* = 300 K: **{t300['status']}**, flags {t300['flags'] or '(none)'}, rules "
                 f"{t300['rule_ids']}.")
        L.append("")
    L.append(f"omega_c used: {_fmt(out['float_reference']['omega_c_K'], 6)} K. Verdict files and "
             "re-check commands:")
    L.append("")
    for tag, r in out["proved"].items():
        L.append(f"- `{r['recheck_cli']}`")
    rv = out.get("run_verdict")
    if rv:
        L.append("")
        L.append(f"The `rtsc run` verdict itself: {rv.get('status')} [{rv.get('T_low_K')}, "
                 f"{rv.get('T_up_K')}] at T*={rv.get('T_star')} K, mu_hi={rv.get('mu_hi')}; flags "
                 f"{rv.get('flags')}; recheck {'OK' if rv.get('recheck_ok') else 'FAILED ' + str(rv.get('recheck_errors') or rv.get('error'))}.")
    L.append("")
    L.append(f"### 3. Float reference, not a verdict")
    L.append("")
    L.append("Allen–Dynes (with f1·f2) and an independent float solution of the isotropic linearised "
             "Eliashberg equations on the same spectrum (nominal weights, same block cutoff ω_c). "
             "**These numbers are NOT proved and are NOT verdicts.**")
    L.append("")
    fl = out["float_reference"]["per_smearing"]
    mus = list(fl[0]["eliashberg_float"].keys()) if fl else []
    hdr = "| # | " + " | ".join(f"AD μ*={m}" for m in mus if m != "0.00") + " | " + \
        " | ".join(f"Eliashberg μ*={m}" for m in mus) + " |"
    L.append(hdr)
    L.append("|" + "---|" * (1 + 2 * len(mus) - 1))
    for f in fl:
        L.append(f"| {f['smearing_index']} | " +
                 " | ".join(_fmt(f["allen_dynes"][m]["Tc_K"]) for m in mus if m != "0.00") + " | " +
                 " | ".join(_fmt(f["eliashberg_float"][m]["Tc_K"])
                            + (" †" if f["eliashberg_float"][m].get("truncated") else "")
                            for m in mus) + " |")
    L.append("")
    L.append("† Matsubara sum truncated at N = 4000 (needed more): a LOWER estimate of the float Tc.")
    L.append("")
    L.append("### 4. Literature (LITERATURE values, quoted from memory; see confidence)")
    L.append("")
    L.append("| material | Tc exp (K) | P (GPa) | λ (literature) | ω_log (literature) | confidence |")
    L.append("|---|---|---|---|---|---|")
    for k, d in LITERATURE.items():
        L.append(f"| {k} | {d['Tc_exp_K']} ({d['Tc_exp_ref']}) | {d['pressure_GPa']} | {d['lambda']} | "
                 f"{d['omega_log_K']} | {d['confidence']} |")
    L.append("")
    if key in LITERATURE:
        L.append(f"Note for {key}: {LITERATURE[key]['notes']}")
        L.append("")
    L.append("### 5. Consistency checks")
    L.append("")
    ch = out["checks"]
    if ch["SERIOUS"]:
        L.append("**!!! SERIOUS FINDINGS !!!**")
        for s in ch["SERIOUS"]:
            L.append(f"- **{s}**")
        L.append("")
    else:
        L.append("No serious finding (all rechecks OK, F7 bookkeeping consistent, every float Tc "
                 "inside its proved bracket).")
        L.append("")
    L.append("| check | result | detail |")
    L.append("|---|---|---|")
    for it in ch["items"]:
        res = {True: "yes", False: "**NO**", None: "info"}[it["result"]]
        L.append(f"| {it['check']} | {res} | {it['detail']} |")
    L.append("")
    return "\n".join(L)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", help="rtsc run material dir (record.json) or QE dir (a2F.dos* + lambda)")
    ap.add_argument("material", help=f"material key: {', '.join(LITERATURE)}")
    ap.add_argument("--out", default="bench_out", help="output directory")
    ap.add_argument("--material-id")
    ap.add_argument("--tc-lit", type=float, help="override the literature T* used for the bracket runs")
    ap.add_argument("--mu-list", default="0.10,0.15")
    ap.add_argument("--envelope", type=float, default=0.10,
                    help="minimum weight envelope (the input's own is kept if larger)")
    ap.add_argument("--omega-c-factor", type=float, default=10.0)
    ap.add_argument("--mode-300", default="decide", choices=["decide", "bracket"])
    ap.add_argument("--skip-300", action="store_true")
    ap.add_argument("--time-budget", type=float, default=None,
                    help="CertifyParams.time_budget_s (default: rtsc's 120 s; only refinement is cut)")
    ap.add_argument("--jobs", type=int, default=1, help="parallel certify processes (<= 2 advised)")
    ap.add_argument("--no-rtsc-float", dest="rtsc_float", action="store_false",
                    help="skip rtsc's own float_tc comparison")
    ap.add_argument("--smoke", action="store_true", help="force the SMOKE TEST label")
    ap.add_argument("--negative-policy", default="error", choices=["error", "clip"],
                    help="rtsc loader policy for negative alpha^2F bins (clip = explicit opt-in)")
    a = ap.parse_args(argv)
    a.mu_list = [float(x) for x in a.mu_list.split(",") if x.strip()]
    global NEGATIVE_POLICY
    NEGATIVE_POLICY = a.negative_policy
    try:
        os.nice(10)
    except Exception:
        pass
    t0 = time.time()
    out = analyse(a)
    out["wall_total_s"] = time.time() - t0
    od = Path(a.out).resolve()
    key = out["material_key"]
    (od / f"{key}_benchmark.json").write_text(json.dumps(out, indent=1, default=str))
    md = to_markdown(out)
    (od / f"{key}_benchmark.md").write_text(md)
    print(md)
    print(f"\nwritten: {od / (key + '_benchmark.json')}, {od / (key + '_benchmark.md')}")
    return 2 if out["checks"]["SERIOUS"] else 0


if __name__ == "__main__":
    sys.exit(main())
