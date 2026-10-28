"""Cross-checks of the oracle (independent constructions + literature/limit checks).

    python crosscheck.py [--quick]      -> writes crosscheck_results.json

1. matrix: even-sector matrix vs brute-force full Matsubara line (float, direct sgn sums) folded
   by an explicit even embedding; same in mpmath (30 digits) with exact finite remainder; odd
   sector top eigenvalue < 0 (R-ODD); float top eigenvalue confirmed in mp by LDL inertia at
   +-1e-11 shifts.
2. Tc sign in mp: brute-force mp matrix (N=200) is unstable at Tc(1-1e-6), stable at Tc(1+1e-6).
3. eta formulation (top eig of D^-1/2 K D^-1/2 = 1) gives the same Tc for every suite case.
4. exact similarity: (omega, w, omega_c) -> s*(...) gives Tc -> s*Tc.
5. literature/limits: Einstein optimum Tc/A in PROOF.md's interval [0.2308678918, 0.2308678940];
   strong-coupling limit Tc/(omega sqrt(lambda)) -> 0.1827 (Allen-Dynes 1975); Allen-Dynes
   formula ratios; rho+(A/a0) < 0 for every case; Tc monotone non-increasing in mu*.
"""
from __future__ import annotations

import glob
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from oracle import (Model, even_matrix, find_tc, fold_even, fold_odd,  # noqa: E402
                    full_line_matrix_float, mp_even_matrix_fullline, mp_inertia_positive,
                    tc_eta, top_eig)


def matrix_checks():
    import mpmath as mp
    cfgs = [
        ("einstein", [1000.0], [500.0], 120.0, 0.0),
        ("einstein_mu", [1000.0], [500.0], 80.0, 0.1),
        ("twopeak_soft", [150.0, 900.0, 3.0], [60.0, 500.0, 2.0], 37.3, 0.13),
        ("hydride_like3", [300.0, 1700.0, 2300.0], [90.0, 700.0, 500.0], 250.0, 0.15),
        ("low_T_many_Nc", [100.0, 50.0], [22.0, 5.0], 2.0, 0.15),
    ]
    out = []
    for name, o, w, T, mu in cfgs:
        o, w = np.array(o), np.array(w)
        wc = 10 * o.max()
        L = 60
        ns, H = full_line_matrix_float(o, w, T, mu, wc, L)
        M = even_matrix(o, w, T, mu, wc, L + 1)
        d_float = float(np.max(np.abs(fold_even(ns, H, L) - M)))
        odd_top = top_eig(fold_odd(ns, H, L))
        Mp = mp_even_matrix_fullline(o, w, T, mu, wc, L)
        d_mp = max(abs(float(Mp[i, j] - mp.mpf(M[i, j]))) for i in range(L + 1)
                   for j in range(L + 1))
        r = top_eig(M)
        above = mp_inertia_positive(Mp, mp.mpf(r) - mp.mpf("1e-11"))
        below = mp_inertia_positive(Mp, mp.mpf(r) + mp.mpf("1e-11"))
        ok = d_mp < 1e-12 * max(1.0, float(np.max(np.abs(M)))) and above == 1 and below == 0 \
            and odd_top < 0
        out.append(dict(cfg=name, T=T, mu=mu, L=L, maxdiff_float_fullline=d_float,
                        maxdiff_mp_fullline=d_mp, top_even=r, top_odd=odd_top,
                        mp_count_above_r_minus=above, mp_count_above_r_plus=below, ok=ok))
        print("matrix", out[-1], flush=True)
    return out


def tc_sign_mp():
    import mpmath as mp
    out = []
    for name, o, w, mu in [("einstein_lam1_mu0.1", [1000.0], [500.0], 0.1),
                           ("twopeak_a_mu0", [200.0, 1200.0], [50.0, 420.0], 0.0)]:
        m = Model.from_spectrum(o, w, mu)
        r = find_tc(m)
        res = {"case": name, "Tc": r.Tc, "err": r.err}
        L = 199
        for tag, f in (("below", 1 - 1e-6), ("above", 1 + 1e-6)):
            T = r.Tc * f
            t = time.time()
            Mp = mp_even_matrix_fullline(o, w, T, mu, m.omega_c, L, dps=25)
            npos = mp_inertia_positive(Mp, mp.mpf(0), dps=25)
            res[f"npos_{tag}"] = npos
            res[f"sec_{tag}"] = round(time.time() - t, 1)
        res["ok"] = res["npos_below"] >= 1 and res["npos_above"] == 0
        print("tc_sign_mp", res, flush=True)
        out.append(res)
    return out


def similarity():
    out = []
    for o, w, mu in [([200.0, 1200.0], [50.0, 420.0], 0.1), ([1000.0], [1000.0], 0.15)]:
        t1 = find_tc(Model.from_spectrum(o, w, mu)).Tc
        s = 1.7
        t2 = find_tc(Model.from_spectrum([s * x for x in o], [s * x for x in w], mu)).Tc
        out.append(dict(o=o, w=w, mu=mu, Tc=t1, Tc_scaled_over_s=t2 / s,
                        rel=abs(t2 / s - t1) / t1))
        print("similarity", out[-1], flush=True)
    return out


def literature():
    out = {}
    # (a) Einstein optimum vs PROOF.md interval
    A = 1.0
    best = None
    for lam in (1.13, 1.14, 1.1436, 1.15, 1.16):
        tc = find_tc(Model.from_spectrum([2 * A / lam], [A], 0.0)).Tc
        best = max(best or 0, tc)
        out[f"einstein_TcA_lam{lam}"] = tc / A
    out["einstein_TcA_best"] = best
    out["proof_interval"] = [0.2308678918, 0.2308678940]
    out["einstein_opt_within_proof_interval_tol1e-10"] = \
        0.2308678918 - 1e-10 <= best <= 0.2308678940 + 1e-10
    # (b) strong coupling: Tc -> 0.1827 sqrt(lambda) omega
    for lam in (25.0, 100.0, 400.0):
        om = 100.0
        tc = find_tc(Model.from_spectrum([om], [lam * om / 2], 0.0)).Tc
        out[f"strong_Tc_over_omega_sqrtlam_lam{lam:g}"] = tc / (om * math.sqrt(lam))
    out["strong_limit_literature"] = 0.1827
    print("literature", out, flush=True)
    return out


def suite_checks():
    rows = []
    for f in sorted(glob.glob(os.path.join(HERE, "cases", "*.json"))):
        d = json.load(open(f))
        o, w = d["spectrum"]["omega_K"], d["spectrum"]["weight_K"]
        prev = None
        for mu, r in d["oracle"].items():
            row = dict(name=d["name"], mu=float(mu), Tc=r["Tc"], err=r["err"],
                       status=r["status"], rho_at_ceiling=r["rho_at_ceiling"],
                       reentrant=r["reentrant"], AD=d["allen_dynes"][mu])
            row["ceiling_negative"] = r["rho_at_ceiling"] < 0
            if r["Tc"] is not None:
                m = Model.from_spectrum(o, w, float(mu), d["omega_c_K"])
                te = tc_eta(m, r["Tc"], r["N_used"])
                row["Tc_eta"] = te
                row["rel_diff_eta"] = abs(te - r["Tc"]) / r["Tc"] if te == te else None
                row["AD_ratio"] = row["AD"] / r["Tc"]
            cur = r["Tc"] if r["Tc"] is not None else 0.0
            row["monotone_in_mu"] = prev is None or cur <= prev + 1e-9 * max(prev, 1)
            prev = cur
            rows.append(row)
            print("suite", {k: row[k] for k in ("name", "mu", "Tc", "rel_diff_eta",
                                                "ceiling_negative", "monotone_in_mu")
                            if k in row}, flush=True)
    return rows


if __name__ == "__main__":
    parts = [a for a in sys.argv[1:] if not a.startswith("-")] or \
        ["matrix", "similarity", "literature", "tc_sign_mp", "suite"]
    fn = os.path.join(HERE, "crosscheck_results.json")
    res = json.load(open(fn)) if os.path.exists(fn) else {}
    funcs = {"matrix": matrix_checks, "similarity": similarity, "literature": literature,
             "tc_sign_mp": tc_sign_mp, "suite": suite_checks}
    for p in parts:
        t = time.time()
        res[p] = funcs[p]()
        res[p + "_seconds"] = time.time() - t
        json.dump(res, open(fn, "w"), indent=1, default=float)
