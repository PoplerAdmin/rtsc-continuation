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
    if 1.5 * N1 <= 6000:
        e2, N2, _, _ = eliashberg_eta(om, p, Tc, mu, omega_c, N=int(1.5 * N1))
        out["eta_at_Tc_1p5N"] = e2
        out["truncation_shift_eta"] = e2 - e1
    out["seconds"] = time.time() - t0
    return out


