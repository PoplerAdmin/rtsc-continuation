"""Independent referee oracle (written from SPEC.md only; does NOT import rtsc.certify.floatcore
or tests/oracle).  Float64 numpy for Tc scans, mpmath for high-precision spot checks.

Even sector, n,m = 0..N-1:
    M_nm = lam(|n-m|) + lam(n+m+1)            (n != m)
    M_nn = lam(2n+1) - (2n+1) - 2 sum_{k=1}^n lam(k)
    M_nm -= 2 mu   for n, m < N_c(T) = #{n >= 0 : (2n+1) pi T <= omega_c}
Truncation is a principal submatrix, so rho_N <= rho (a truncated rho_N >= 0 PROVES instability up
to float error).
"""
from __future__ import annotations

import math

import mpmath as mp
import numpy as np


def nc_count(T, wc):
    if wc is None or not math.isfinite(wc):
        return 10**9
    n = 0
    # direct definition (small loops only used with moderate counts)
    x = wc / (math.pi * T)
    if x < 1:
        return 0
    n = int((x - 1) // 2) + 1
    while n > 0 and (2 * (n - 1) + 1) * math.pi * T > wc:
        n -= 1
    while (2 * n + 1) * math.pi * T <= wc:
        n += 1
    return n


def lam_vec(om, p, T, kmax):
    om = np.asarray(om, float)
    p = np.asarray(p, float)
    k = np.arange(1, kmax + 1, dtype=float)
    x2 = (2 * math.pi * T * k) ** 2
    out = np.zeros(kmax + 1)
    step = max(1, 3_000_000 // max(1, len(om)))
    for i0 in range(0, kmax, step):
        xx = x2[i0:i0 + step]
        out[1 + i0:1 + i0 + len(xx)] = (2 * p[None, :] * om[None, :]
                                        / (om[None, :] ** 2 + xx[:, None])).sum(1)
    return out


def even_M(om, p, T, N, mu=0.0, wc=None):
    lam = lam_vec(om, p, T, 2 * N + 1)
    n = np.arange(N)
    M = lam[np.abs(n[:, None] - n[None, :])] + lam[n[:, None] + n[None, :] + 1]
    S = np.concatenate([[0.0], np.cumsum(lam[1:N])])
    M[n, n] = lam[2 * n + 1] - (2 * n + 1) - 2 * S
    if mu:
        c = min(nc_count(T, wc), N)
        if c > 0:
            M[:c, :c] -= 2 * mu
    return M


def rho_N(om, p, T, N, mu=0.0, wc=None):
    M = even_M(om, p, T, N, mu, wc)
    return float(np.linalg.eigvalsh(M)[-1])


def auto_N(om, T, mu=0.0, wc=None, fac=8.0, cap=2500):
    k0 = max(om) / (2 * math.pi * T)
    N = max(80, int(fac * k0) + 40)
    if mu and wc:
        N = max(N, min(nc_count(T, wc), 100000) + 40)
    return min(N, cap)


def tc_lower(om, p, mu=0.0, wc=None, T_hi=None, N=None, rtol=1e-9, T_min=1e-4, fac=8.0, cap=2500):
    """Largest T (found by a downward scan + bisection) with rho_N(T) >= 0.
    Because rho_N <= rho this is (up to float error) a LOWER bound on the largest unstable T
    near the scan; returns 0.0 if nothing unstable is found above T_min."""
    A = float(np.sum(p))
    if T_hi is None:
        T_hi = A / 4.33148145 * 1.0001
    f = lambda T: rho_N(om, p, T, N or auto_N(om, T, mu, wc, fac, cap), mu, wc)
    hi = T_hi
    if f(hi) >= 0:
        return hi
    lo = hi * 0.97
    while f(lo) < 0:
        hi = lo
        lo *= 0.97
        if lo < T_min:
            return 0.0
    while hi / lo - 1 > rtol:
        mid = math.sqrt(lo * hi)
        if f(mid) >= 0:
            lo = mid
        else:
            hi = mid
    return lo


# --------------------------------------------------------------------------- mpmath checks
def mp_lam(om, p, T, k):
    x2 = (2 * mp.pi * mp.mpf(T) * k) ** 2
    return mp.fsum(2 * mp.mpf(w) * mp.mpf(o) / (mp.mpf(o) ** 2 + x2) for o, w in zip(om, p))


def mp_Q(om, p, T, x, mu, Nc):
    """Q^+_T(x) straight from the SPEC definition (finite x), high precision."""
    N = len(x)
    xs = [mp.mpf(v) for v in x]
    lam = [mp.mpf(0)] + [mp_lam(om, p, T, k) for k in range(1, 2 * N + 1)]
    Q = mp.mpf(0)
    for n in range(N):
        for m in range(N):
            if n != m:
                Q += lam[abs(n - m)] * xs[n] * xs[m]
            Q += lam[n + m + 1] * xs[n] * xs[m]
    S = mp.mpf(0)
    for n in range(N):
        if n >= 1:
            S += lam[n]
        Q -= ((2 * n + 1) + 2 * S) * xs[n] ** 2
    sx = mp.fsum(xs[:min(Nc, N)])
    Q -= 2 * mp.mpf(mu) * sx ** 2
    return Q


def mp_atomF(o, T, x):
    """F_j(x): Q is sum_j p_j F_j - diag - mu term.  Computed from the definition with p = 1."""
    N = len(x)
    xs = [mp.mpf(v) for v in x]
    lam = [mp.mpf(0)] + [2 * mp.mpf(o) / (mp.mpf(o) ** 2 + (2 * mp.pi * mp.mpf(T) * k) ** 2)
                         for k in range(1, 2 * N + 1)]
    F = mp.mpf(0)
    for n in range(N):
        for m in range(N):
            if n != m:
                F += lam[abs(n - m)] * xs[n] * xs[m]
            F += lam[n + m + 1] * xs[n] * xs[m]
    S = mp.mpf(0)
    for n in range(N):
        if n >= 1:
            S += lam[n]
        F -= 2 * S * xs[n] ** 2
    return F
