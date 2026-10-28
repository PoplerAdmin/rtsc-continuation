"""Validation of rtsc.certify on SYNTHETIC spectra (not data for any real compound).

(a) Einstein spectrum: bracket contains the float Tc; Tc/A -> 0.23087 at Omega = 7.5749 Tc.
(b) ~30 random multi-peak spectra (incl. mu* > 0, e > 0): independent float reference Tc
    (N = 600) must lie in [T_low, T_up]; report relative bracket widths and timings.
(c) negative controls: tampered certificates must fail recheck.
(d) weak- and strong-coupling cases.

Usage:  nice python3 tests/validate_certify.py [out.json]
"""
from __future__ import annotations

import copy
import json
import math
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np  # noqa: E402

from synthetic import einstein, float_tc_reference, multipeak  # noqa: E402
from rtsc.certify import CertifyParams, classify, recheck  # noqa: E402
from rtsc.status import Verdict  # noqa: E402


def run_case(spec, mu, T_star, label, rel_tol=1e-4, refN=600):
    wc = 10 * max(spec.omega_K)
    t = time.time()
    v = classify([spec], CertifyParams(T_star=T_star, mu_hi=mu, rel_tol=rel_tol), material_id=label)
    tcls = time.time() - t
    t = time.time()
    rc = recheck(Verdict.from_json(v.to_json()), [spec])
    trc = time.time() - t
    ref = float_tc_reference(spec, mu=mu, omega_c=wc, N=refN)          # nominal, mu* = mu
    ref0 = float_tc_reference(spec, N=refN) if mu > 0 else ref         # nominal, mu* = 0
    inside = (v.T_low_K is not None and v.T_up_K is not None and v.T_low_K <= ref <= v.T_up_K
              and v.T_low_K <= ref0 <= v.T_up_K)
    a = v.certificate["alternatives"][0]
    row = {"label": label, "atoms": spec.n_atoms(), "lambda": round(spec.lambda_(), 4),
           "omega_log_K": round(spec.omega_log(), 2), "e": spec.rel_envelope, "mu": mu,
           "T_star": T_star, "status": v.status.value, "T_low": v.T_low_K, "T_up": v.T_up_K,
           "Tc_ref_mu": ref, "Tc_ref_mu0": ref0, "inside": inside,
           "rel_width": (v.T_up_K - v.T_low_K) / ref0 if inside else None,
           "up_gap": v.T_up_K / ref0 - 1 if v.T_up_K else None,
           "low_gap": 1 - v.T_low_K / ref if v.T_low_K and ref > 0 else None,
           "n_cells": len(a["cells"]["cells"]), "t_classify_s": round(tcls, 2),
           "t_recheck_s": round(trc, 2), "recheck_ok": rc["ok"]}
    return row, v


def tamper_tests(v, spec):
    """(c) negative controls: every tampered certificate must FAIL recheck."""
    out = []

    def attempt(name, mutate):
        w = Verdict.from_json(v.to_json())
        try:
            mutate(w)
        except Exception as exc:  # pragma: no cover
            out.append({"tamper": name, "applied": False, "err": repr(exc)})
            return
        r = recheck(w, [spec])
        out.append({"tamper": name, "recheck_ok": r["ok"], "errors": r["errors"][:2]})

    alt = lambda w: w.certificate["alternatives"][0]

    def scale_u(w):
        c = alt(w)["cells"]["cells"][0]
        i = len(c["u"]) // 3
        c["u"][i] *= 1.5
    attempt("head entry u_i *= 1.5 in first cell", scale_u)

    def eta_up(w):
        c = alt(w)["cells"]["cells"][0]
        c["eta"] *= 1e3
    attempt("eta *= 1000 in first cell", eta_up)

    def gap(w):
        cells = alt(w)["cells"]["cells"]
        if len(cells) >= 2:
            T = float.fromhex(cells[1]["T1"])
            cells[1]["T1"] = math.nextafter(T, math.inf).hex()
    attempt("1-ulp gap between cells 0 and 1", gap)

    def lower_Tup(w):
        cblk = alt(w)["cells"]
        cells = cblk["cells"]
        T1 = float.fromhex(cells[0]["T1"]) * 0.99
        cells[0]["T1"] = T1.hex()
        cblk["T_up"] = T1.hex()
        w.T_up_K = T1
    attempt("move T_up 1% lower (first cell stretched)", lower_Tup)

    def ray_T(w):
        r = alt(w)["rayleigh"]
        T = float.fromhex(r["T"]) * 1.02
        r["T"] = T.hex()
        w.T_low_K = T
    attempt("raise Rayleigh T by 2% (with matching stored bracket)", ray_T)

    def ray_x(w):
        r = alt(w)["rayleigh"]
        r["x"] = [x * (1 + 0.3 * ((i % 2) * 2 - 1)) for i, x in enumerate(r["x"])]
    attempt("perturb Rayleigh vector +/-30% alternating", ray_x)

    def status_flip(w):
        from rtsc.status import Status
        w.status = Status.CERTIFIED_CANDIDATE if w.status != Status.CERTIFIED_CANDIDATE else Status.REJECTED_PROVEN
    attempt("flip stored status", status_flip)

    def spec_change(w):
        pass
    # different spectrum (weights +1%) with the same certificate
    w = Verdict.from_json(v.to_json())
    from rtsc.spectrum import Spectrum
    sp2 = Spectrum(spec.omega_K, tuple(x * 1.01 for x in spec.weight_K), spec.rel_envelope, spec.label)
    r = recheck(w, [sp2])
    out.append({"tamper": "spectrum weights *1.01, same certificate", "recheck_ok": r["ok"],
                "errors": r["errors"][:2]})

    def mu_change(w):
        w.certificate["params"]["mu_hi"] = 0.0 if w.certificate["params"]["mu_hi"] > 0 else 0.1
    attempt("change mu_hi in stored params", mu_change)
    return out


def main(outpath="validation_certify.json"):
    results = {"a_einstein": [], "b_random": [], "c_tamper": [], "d_regimes": []}
    t_all = time.time()
    # ---------------- (a) Einstein
    A = 1000.0
    ratio_opt = 7.5749
    T_guess = A / 4.33148145
    for Om_over in [ratio_opt, 3.0, 20.0]:
        spec = einstein(Om_over * 0.2308679 * A, A)
        row, v = run_case(spec, 0.0, 300.0, f"einstein Omega/T~{Om_over}", rel_tol=1e-6)
        row["Tlow_over_A"] = v.T_low_K / A
        row["Tup_over_A"] = v.T_up_K / A
        results["a_einstein"].append(row)
        print("(a)", json.dumps(row), flush=True)
    # ---------------- (b) random multipeak
    rng = np.random.default_rng(12345)
    cases = []
    for i in range(30):
        lam = float(np.exp(rng.uniform(np.log(0.45), np.log(4.0))))
        om_hi = float(rng.uniform(300, 2500))
        n_pts = int(rng.choice([200, 600, 1500, 3000]))
        mu = float(rng.choice([0.0, 0.0, 0.1, 0.13]))
        e = float(rng.choice([0.0, 0.0, 0.0, 0.05, 0.1]))
        cases.append((i, lam, om_hi, n_pts, mu, e))
    for (i, lam, om_hi, n_pts, mu, e) in cases:
        spec = multipeak(1000 + i, lam, n_pts=n_pts, om_lo=0.05 * om_hi, om_hi=om_hi, e=e)
        T_star = float(rng.choice([20.0, 77.0, 300.0]))
        row, v = run_case(spec, mu, T_star, f"rand{i}")
        results["b_random"].append(row)
        print("(b)", json.dumps(row), flush=True)
        if i == 3:
            results["c_tamper"] = tamper_tests(v, spec)
            for t in results["c_tamper"]:
                print("(c)", json.dumps(t), flush=True)
    # ---------------- (d) regimes
    for (lam, om_hi, mu, lab) in [(0.3, 1500.0, 0.0, "weak lambda=0.3, Omega_max=1500K"),
                                  (0.5, 2000.0, 0.1, "weak lambda=0.5 mu*=0.1"),
                                  (3.5, 400.0, 0.13, "strong lambda=3.5 mu*=0.13"),
                                  (2.5, 2500.0, 0.1, "strong high-freq lambda=2.5 (Tc~300K)")]:
        spec = multipeak(77, lam, n_pts=2000, om_lo=0.05 * om_hi, om_hi=om_hi)
        row, v = run_case(spec, mu, 300.0, lab)
        results["d_regimes"].append(row)
        print("(d)", json.dumps(row), flush=True)
    results["total_seconds"] = time.time() - t_all
    Path(outpath).write_text(json.dumps(results, indent=1))
    print("TOTAL", results["total_seconds"])


if __name__ == "__main__":
    main(*(sys.argv[1:]))
