"""Rigorous (ball-arithmetic) building blocks for the Einstein-spectrum even-sector
problem in units T = 1:  lam(k) = 2 a w / (w^2 + (2 pi k)^2),  mu* = 0.

Matrix (n, m = 0..N-1), lam(0)-terms already cancelled:
   B_nm = lam(|n-m|) + lam(n+m+1)                    (n != m)
   B_nn = lam(2n+1) - (2n+1) - 2 sum_{k=1}^n lam(k)
All arithmetic is arb (python-flint) except the float Cholesky factor used as a
certificate, whose residual is then bounded rigorously in arb.
"""
from flint import arb, arb_mat, fmpq, ctx
import numpy as np

ctx.prec = 128
PI = arb.pi()
FOURPI2 = 4 * PI * PI


def lam(k, a, w):
    return 2 * a * w / (w * w + FOURPI2 * k * k)


def dlam(k, a, w):
    c2 = FOURPI2 * k * k
    d = w * w + c2
    return 2 * a * (c2 - w * w) / (d * d)


def d2lam(k, a, w):
    c2 = FOURPI2 * k * k
    d = w * w + c2
    return 4 * a * w * (w * w - 3 * c2) / (d * d * d)


def Ssum(a, w):
    """S = sum_{k>=1} lam(k) = a (coth(w/2)/2 - 1/w)  (closed form)."""
    return a * ((w / 2).coth() / 2 - 1 / w)


def build(f, a, w, N, with_diag_shift):
    """matrix with entries built from f(k) = lam/dlam/d2lam at w."""
    vals = [None] + [f(k, a, w) for k in range(1, 2 * N + 1)]
    cum = [arb(0)]
    for k in range(1, N):
        cum.append(cum[-1] + vals[k])
    rows = []
    for n in range(N):
        row = []
        for m in range(N):
            if n == m:
                e = vals[2 * n + 1] - 2 * cum[n]
                if with_diag_shift:
                    e = e - (2 * n + 1)
            else:
                e = vals[abs(n - m)] + vals[n + m + 1]
            row.append(e)
        rows.append(row)
    return arb_mat(rows)


def ghead(a, w, N, L):
    """G = sum_{m=N}^{N+L-1} c_m c_m^T / beta_m,  c_m(i) = lam(m-i)+lam(m+i+1),
    beta_m = 2m+1-2S(w).  Evaluated in ball arithmetic at (possibly ball) w."""
    vals = [None] + [lam(k, a, w) for k in range(1, 2 * N + L + 1)]
    S = Ssum(a, w)
    C = arb_mat([[vals[m - i] + vals[m + i + 1] for m in range(N, N + L)] for i in range(N)])
    CD = arb_mat([[C[i, j] / (2 * (N + j) + 1 - 2 * S) for j in range(L)] for i in range(N)])
    return CD * C.transpose(), S


def gtail(a, wmax, N, L, S):
    """trace bound for sum_{m>=N+L} c_m c_m^T/beta_m (see PROOF.md, Lemma T)."""
    q = a * wmax / (2 * PI * PI)
    return 4 * q * q / (6 * (L - 1) ** 2 * (2 * N + 2 * L + 1 - 2 * S))


def rowabs_upper(Mb):
    n = Mb.nrows()
    out = []
    for i in range(n):
        s = arb(0)
        for j in range(n):
            s += Mb[i, j].abs_upper()
        out.append(s.upper())
    return out


def certify_nsd(Mb, diag_extra, eps=1e-11):
    """Rigorously certify  x^T Mb' x + sum_i diag_extra[i] x_i^2 <= 0  for every matrix Mb'
    in the ball matrix Mb (symmetric).  Returns (ok, info)."""
    n = Mb.nrows()
    # For symmetric M' in the ball matrix Mb:
    #   x^T M' x <= x^T sym(mid) x + sum_i x_i^2 * (rowrad_i + colrad_i)/2 .
    # P = -sym(mid(Mb)) - diag((rowrad+colrad)/2 + diag_extra)  (exactly symmetric arb entries)
    rad = [[Mb[i, j].rad() for j in range(n)] for i in range(n)]
    mid = [[Mb[i, j].mid() for j in range(n)] for i in range(n)]
    P = [[None] * n for _ in range(n)]
    for i in range(n):
        rr = arb(0)
        for j in range(n):
            rr += rad[i][j] + rad[j][i]
        rr = rr / 2
        for j in range(i, n):
            e = -(mid[i][j] + mid[j][i]) / 2
            if i == j:
                e = e - rr - diag_extra[i]
            P[i][j] = e
            P[j][i] = e
    Pf = np.array([[float(P[i][j].mid()) for j in range(n)] for i in range(n)])
    Pf = 0.5 * (Pf + Pf.T)
    try:
        Lf = np.linalg.cholesky(Pf - eps * np.eye(n))
    except np.linalg.LinAlgError:
        return False, "cholesky failed; minimal eig approx %.3e" % np.linalg.eigvalsh(Pf)[0]
    La = arb_mat([[arb(float(Lf[i, j])) for j in range(n)] for i in range(n)])
    LL = La * La.transpose()
    worst = arb(0)
    for i in range(n):
        s = arb(0)
        for j in range(n):
            e = P[i][j] - LL[i, j]
            if i == j:
                e = e - eps
            s += e.abs_upper()
        if s.upper() > worst.upper():
            worst = s.upper()
    ok = worst < eps  # arb comparison: True only if certainly
    return bool(ok), "resid rowsum %s, approx min eig %.3e" % (worst.str(3), np.linalg.eigvalsh(Pf)[0])
