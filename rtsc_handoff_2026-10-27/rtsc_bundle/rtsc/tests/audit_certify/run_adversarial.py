"""Adversarial synthetic spectra: engine verdict vs independent oracle (indep.py).

usage: nice python run_adversarial.py <group> [mode]
Writes results_<group>_<mode>.json.  Critical findings printed with 'CRITICAL'.
"""
from __future__ import annotations

import json
import math
import sys
import time

import mpmath as mp
import numpy as np

sys.path.insert(0, "/home/claude/rtsc")
sys.path.insert(0, "/home/claude/rtsc/tests/audit_certify")

from rtsc.spectrum import Spectrum  # noqa: E402
from rtsc.certify import classify, CertifyParams, recheck  # noqa: E402
from rtsc.certify.rayleigh import certified_n_c  # noqa: E402  (only for N_c comparison print)
import indep  # noqa: E402

A0 = 4.33148145


def einstein_A_for_tc(Om, Tc, N=None):
    """A such that the Einstein mode Om has oracle Tc (mu=0) exactly Tc (bisection on A)."""
    lo, hi = A0 * Tc * 0.999999, A0 * Tc * 50
    T = Tc
    NN = N or indep.auto_N([Om], T, fac=12)
    f = lambda A: indep.rho_N([Om], [A], T, NN)
    if f(lo) >= 0:
        return lo
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(mid) >= 0:
            hi = mid
        else:
            lo = mid
        if hi - lo < 1e-13 * hi:
            break
    return hi


def cases(group):
    rng = np.random.default_rng(12345)
    out = []
    if group in ("a", "all"):
        om = np.logspace(0, 5, 40)
        w = np.full(40, 1500.0 / 40)
        out.append(("wide_1K_1e5K", Spectrum(tuple(om), tuple(w), 0.0), {}))
        out.append(("tiny_low_atom", Spectrum((1e-3, 1000.0), (1e-6, 1250.0), 0.0), {}))
        out.append(("tiny_low_atom2", Spectrum((1e-6, 900.0), (1e-3, 1380.0), 0.0), {}))
        oms = rng.uniform(50, 2000, 3000)
        ws = rng.uniform(0, 1, 3000)
        ws *= 1400.0 / ws.sum()
        out.append(("many_atoms_3000", Spectrum(tuple(oms), tuple(ws), 0.05), {}))
    if group in ("b", "all"):
        out.append(("huge_atom_w301", Spectrum((301.0,), (1299.4,), 0.0), {}))
        out.append(("huge_atom_w301_A2000", Spectrum((301.0,), (2000.0,), 0.0), {}))
        out.append(("e_0.99", Spectrum((800.0, 2000.0), (400.0, 300.0), 0.99), {}))
        out.append(("e_0.999_einstein", Spectrum((1000.0,), (700.0,), 0.999), {}))
        out.append(("mu_large_10", Spectrum((1000.0,), (2500.0,), 0.0), {"mu_hi": 10.0}))
        out.append(("wc_tiny", Spectrum((1000.0,), (2500.0,), 0.0), {"mu_hi": 5.0, "omega_c": 1.0}))
        out.append(("wc_huge", Spectrum((1000.0,), (2500.0,), 0.0), {"mu_hi": 0.15, "omega_c": 1e7}))
    if group in ("c", "all"):
        for Om in (1000.0, 2271.0, 5000.0):
            A_c = einstein_A_for_tc(Om, 300.0)
            for rel in (-1e-4, -1e-6, 1e-6, 1e-4):
                A = A_c * (1 + rel)
                out.append((f"near_Tstar_Om{Om:g}_rel{rel:g}", Spectrum((Om,), (A,), 0.0),
                            {"mu_hi": 0.0, "_A_c": A_c}))
    if group in ("d", "all"):
        out.append(("two_peak_extreme", Spectrum((50.0, 50000.0), (600.0, 1500.0), 0.0), {}))
        out.append(("low_peak_strong", Spectrum((30.0,), (1400.0,), 0.0), {}))
        out.append(("tiny_T_scale", Spectrum((2.0,), (3.0,), 0.0), {"T_star": 0.5}))
        out.append(("high_T_scale", Spectrum((1e6,), (2e6,), 0.1), {"T_star": 3e5}))
    return out


def oracle_compare(name, spec, kw, v, prm):
    rec = {}
    om = list(spec.omega_K)
    p = np.asarray(spec.weight_K)
    e = spec.rel_envelope
    T_star = prm.T_star
    mu_hi = prm.mu_hi
    wc = prm.resolved_omega_c([spec])
    # oracle Tc at mu=0 for envelope vertices (upper side) -- REJECTED must be beyond all of them
    verts = {"hi": p * (1 + e), "nom": p}
    if len(p) > 1 and e > 0:
        r2 = np.random.default_rng(7)
        for i in range(3):
            s = r2.integers(0, 2, len(p))
            verts[f"rv{i}"] = p * np.where(s == 1, 1 + e, 1 - e)
        verts["lo"] = p * (1 - e)
    tcs = {}
    for k, pp in verts.items():
        t0 = time.time()
        tcs[k] = indep.tc_lower(om, pp, 0.0, None, cap=2000)
        rec.setdefault("oracle_s", 0)
        rec["oracle_s"] += time.time() - t0
    rec["oracle_tc_mu0"] = tcs
    tc_max = max(tcs.values())
    crit = []
    if v.status.value == "REJECTED_PROVEN" and tc_max >= T_star * (1 - 1e-9):
        crit.append(f"REJECTED but oracle Tc(mu=0) max over vertices = {tc_max} >= T*={T_star}")
    if v.T_up_K is not None and tc_max > v.T_up_K * (1 + 1e-9):
        crit.append(f"T_up={v.T_up_K} below oracle Tc {tc_max}")
    # lower side: oracle Tc at the (1-e) vertex with mu_hi (float lower bound on Tc for that member)
    if v.T_low_K is not None:
        tl = indep.tc_lower(om, p * (1 - e), mu_hi, wc, T_hi=v.T_low_K * 1.5, cap=2000)
        rec["oracle_tc_lo_muhi"] = tl
        # independent mpmath check of every stored Rayleigh certificate
        for alt in v.certificate["alternatives"]:
            rc = alt.get("rayleigh")
            if not rc:
                continue
            T = float.fromhex(rc["T"])
            x = rc["x"]
            Nc = indep.nc_count(T, rc["omega_c"])
            rec["Nc_indep"] = Nc
            rec["Nc_cert"] = rc["N_c"]
            if Nc != rc["N_c"]:
                crit.append(f"N_c mismatch indep {Nc} vs cert {rc['N_c']}")
            if len(x) <= 160:
                mp.mp.dps = 40
                Fs = [indep.mp_atomF(o, T, x) for o in om]
                Nx = len(x)
                diag = mp.fsum((2 * n + 1) * mp.mpf(x[n]) ** 2 for n in range(Nx))
                sx = mp.fsum(mp.mpf(t) for t in x[:min(Nc, Nx)])
                boxmin = mp.fsum(min((1 - mp.mpf(e)) * mp.mpf(pj) * F, (1 + mp.mpf(e)) * mp.mpf(pj) * F)
                                 for pj, F in zip(p, Fs)) - diag - 2 * mp.mpf(rc["mu_hi"]) * sx ** 2
                rec["mp_boxminQ"] = float(boxmin)
                if boxmin < 0:
                    crit.append(f"Rayleigh cert boxminQ (mpmath) = {boxmin} < 0")
            else:
                rec["mp_boxminQ"] = "skipped (Nx > 160)"
    rec["critical"] = crit
    return rec


def main():
    group = sys.argv[1]
    mode = sys.argv[2] if len(sys.argv) > 2 else "decide"
    res = []
    for name, spec, kw in cases(group):
        kw = dict(kw)
        A_c = kw.pop("_A_c", None)
        prm = CertifyParams(mode=mode, time_budget_s=float(kw.pop("budget", 40.0)), **kw)
        t0 = time.time()
        try:
            v = classify([spec], prm, material_id=name)
        except Exception as exc:
            print(f"{name}: EXCEPTION {exc!r}", flush=True)
            res.append({"name": name, "exception": repr(exc)})
            continue
        dt = time.time() - t0
        rc = recheck(v, [spec])
        rec = {"name": name, "mode": mode, "status": v.status.value, "T_low": v.T_low_K,
               "T_up": v.T_up_K, "seconds": dt, "recheck_ok": rc["ok"], "A_c": A_c,
               "params": {k: getattr(prm, k) for k in ("T_star", "mu_hi", "omega_c")}}
        try:
            rec.update(oracle_compare(name, spec, kw, v, prm))
        except Exception as exc:
            rec["oracle_exception"] = repr(exc)
        tag = "CRITICAL " if rec.get("critical") else ""
        print(f"{tag}{name}: {rec['status']} [{rec['T_low']}, {rec['T_up']}] {dt:.1f}s "
              f"recheck={rc['ok']} oracle={rec.get('oracle_tc_mu0')} lo_mu={rec.get('oracle_tc_lo_muhi')} "
              f"mpQ={rec.get('mp_boxminQ')} {rec.get('critical') or ''} {rec.get('oracle_exception') or ''}",
              flush=True)
        res.append(rec)
    with open(f"/home/claude/rtsc/tests/audit_certify/results_{group}_{mode}.json", "w") as fh:
        json.dump(res, fh, indent=1, default=str)


if __name__ == "__main__":
    main()
