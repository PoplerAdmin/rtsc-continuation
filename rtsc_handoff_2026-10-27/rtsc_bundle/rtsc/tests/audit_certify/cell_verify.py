"""Independent float verification of check_cell (Theorem C) and an eta-push attack.

For a cell certificate, compute in float64 (numpy, independent code) the Mhat(hi, lo) rows
    r_n = sum_m Mhat_nm v_m,  hi = (1+e) lam(.;T1), lo = (1-e) lam(.;T2),
for n < n_rows over columns m < M_big, plus an upper tail estimate, and return
    margin = max_n (r_n / v_n) + eta      (must be <= 0 if the certificate is true).
Attack: raise eta (and separately other knobs) to the largest value check_cell still accepts,
then confirm the float margin is still <= 0.  A positive margin at an accepted eta = bug.
"""
from __future__ import annotations

import copy
import math
import sys

import numpy as np

sys.path.insert(0, "/home/claude/rtsc")
sys.path.insert(0, "/home/claude/rtsc/tests/audit_certify")
import indep  # noqa: E402
from rtsc.certify.arbcore import ArbSpectrum  # noqa: E402
from rtsc.certify.cells import check_cell  # noqa: E402


def vec(cert, M):
    p = int(cert.get("p", 3))
    N = int(cert["N"])
    v = np.empty(M)
    v[:N] = cert["u"]
    m = np.arange(N, M)
    v[N:] = cert["C"] * (2 * m + 1.0) ** (-p)
    return v


def float_margin(om, w, e, cert, n_rows=None, M_big=None):
    T1 = float.fromhex(cert["T1"]); T2 = float.fromhex(cert["T2"])
    N_e = int(cert["N_e"])
    n_rows = n_rows or 3 * N_e
    M_big = M_big or max(20 * n_rows, 4000)
    K = n_rows + M_big + 2
    hi = indep.lam_vec(om, np.asarray(w) * (1 + e), T1, K)
    lo = indep.lam_vec(om, np.asarray(w) * (1 - e), T2, K)
    hi[0] = 0.0
    v = vec(cert, M_big)
    worst = -np.inf
    Slo = np.concatenate([[0.0], np.cumsum(lo[1:n_rows + 1])])
    m = np.arange(M_big)
    for n in range(n_rows):
        row = hi[np.abs(n - m)] @ v + hi[n + m + 1] @ v          # hi(0)=0 kills m = n Toeplitz
        row += -2 * Slo[n] * v[n] - (2 * n + 1) * v[n]
        # tail of columns m >= M_big: hi(k) <= b/k^2 bound, v_m <= v[M_big-1]
        b = (1 + e) * float(np.dot(w, om)) / (2 * math.pi ** 2 * T1 ** 2)
        row += 2 * v[-1] * b / (M_big - n - 1)
        worst = max(worst, row / v[n])
    return worst + float(cert["eta"])


def push(asp, cert, key, factor_hi=1e6):
    """largest multiplicative factor f on cert[key] (geometric bisection) that check_cell accepts"""
    lo, hi = 1.0, factor_hi
    c = copy.deepcopy(cert)
    c[key] = cert[key] * hi
    if check_cell(asp, c)["ok"]:
        return hi, c
    for _ in range(40):
        mid = math.sqrt(lo * hi)
        c = copy.deepcopy(cert)
        c[key] = cert[key] * mid
        if check_cell(asp, c)["ok"]:
            lo = mid
        else:
            hi = mid
        if hi / lo < 1 + 1e-6:
            break
    c = copy.deepcopy(cert)
    c[key] = cert[key] * lo
    return lo, c


if __name__ == "__main__":
    from rtsc.spectrum import Spectrum
    from rtsc.certify import classify, CertifyParams
    specs = [Spectrum((1000.0,), (1250.0,), 0.0, "E"),
             Spectrum((800.0, 2000.0), (400.0, 300.0), 0.5, "e05"),
             Spectrum((50.0, 3000.0), (300.0, 900.0), 0.1, "2pk")]
    for sp in specs:
        v = classify([sp], CertifyParams(mode="bracket", mu_hi=0.0, time_budget_s=30))
        asp = ArbSpectrum(sp.omega_K, sp.weight_K, sp.rel_envelope)
        cells = v.certificate["alternatives"][0]["cells"]["cells"]
        print(sp.label, v.status.value, v.T_up_K, len(cells), "cells", flush=True)
        for c in cells[:3]:
            m0 = float_margin(sp.omega_K, sp.weight_K, sp.rel_envelope, c)
            f, cpush = push(asp, c, "eta")
            m1 = float_margin(sp.omega_K, sp.weight_K, sp.rel_envelope, cpush)
            tag = "BUG" if m1 > 1e-9 * max(1.0, abs(cpush["eta"])) else "ok"
            print(f"  cell [{float.fromhex(c['T1']):.6g},{float.fromhex(c['T2']):.6g}] N={c['N']} "
                  f"N_e={c['N_e']} M_v={c['M_v']} p={c.get('p',3)}: margin(eta)={m0:.3e}; "
                  f"eta pushed x{f:.6g} -> eta={cpush['eta']:.4e}, float margin={m1:.3e} {tag}",
                  flush=True)
