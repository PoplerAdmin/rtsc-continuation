"""Rigorous (ball-arithmetic) core for the Einstein-spectrum even-sector certificate.

Dimensionless: l(k) = 2 a w / (w^2 + (2 pi k)^2), k >= 1.
Even sector  M+_{nm} = l(|n-m|) + l(n+m+1) (n != m),
             M+_{nn} = l(2n+1) - (2n+1) - 2 sum_{k=1}^n l(k).

All bounds below are justified in CERT.md (lemma numbers refer to it).
Everything that enters a verdict is computed with python-flint `arb` balls;
floating point is used only to *propose* (a shift s and a Cholesky factor L),
and the proposal is then checked rigorously (Lemma 6).
"""
import numpy as np
from flint import arb, arb_mat, ctx

ctx.prec = 128
PI = arb.pi()
TWO_PI = 2 * PI


def ball_interval(w1, w2):
    """Ball containing the closed interval [w1, w2] (w1, w2 exact floats)."""
    return arb(w1).union(arb(w2))


def lvals(a, w, kmax):
    """[None, l(1), ..., l(kmax)] as balls (w may be a ball)."""
    out = [None]
    for k in range(1, kmax + 1):
        b2 = (TWO_PI * k) ** 2
        out.append(2 * a * w / (w * w + b2))
    return out


def dlvals(a, w, kmax):
    """first derivative d l(k)/dw = 2a (b^2 - w^2)/(w^2+b^2)^2."""
    out = [None]
    for k in range(1, kmax + 1):
        b2 = (TWO_PI * k) ** 2
        out.append(2 * a * (b2 - w * w) / (w * w + b2) ** 2)
    return out


def d2l_abs_upper(a, w, kmax):
    """upper bounds (arb, exact) of sup_{w in ball} |l''(k)|,
    l''(k) = -4 a w (3 b^2 - w^2)/(w^2 + b^2)^3."""
    out = [None]
    for k in range(1, kmax + 1):
        b2 = (TWO_PI * k) ** 2
        v = -4 * a * w * (3 * b2 - w * w) / (w * w + b2) ** 3
        out.append(arb(v.abs_upper()))
    return out


def L_total(a, w):
    """L = sum_{k>=1} l(k) = (a/2)(coth(w/2) - 2/w)   (Lemma 1)."""
    return (a / 2) * ((w / 2).coth() - 2 / w)


def T_upper(a, w, j, lv=None, S=None, L=None):
    """rigorous upper bound (exact arb) for T(j) = sum_{k>=j} l(k), j>=1, valid for all
    w in the ball: min of  (i) L - S_{j-1} (exact identity) and
    (ii) l(j) + (a/pi) arctan(w/(2 pi j))  (integral comparison, Lemma 1)."""
    if lv is not None and j < len(lv):
        lj = lv[j]
    else:
        lj = 2 * a * w / (w * w + (TWO_PI * j) ** 2)
    u2 = lj + (a / PI) * (w / (TWO_PI * j)).atan()
    best = u2.upper()
    if S is not None and L is not None and j - 1 < len(S):
        u1 = (L - S[j - 1]).upper()
        best = min(best, u1)
    return arb(best)


def cumsums(lv, nmax):
    """S[n] = sum_{k=1}^n l(k), n = 0..nmax"""
    S = [arb(0)]
    for k in range(1, nmax + 1):
        S.append(S[-1] + lv[k])
    return S


def even_block(lv, S, N, deriv=False):
    """N x N leading block of M+ (deriv=True: of dM+/dw, where lv holds l', S its cumsums)."""
    rows = []
    for n in range(N):
        row = []
        for m in range(N):
            if n == m:
                v = lv[2 * n + 1] - 2 * S[n]
                if not deriv:
                    v = v - (2 * n + 1)
            else:
                v = lv[abs(n - m)] + lv[n + m + 1]
            row.append(v)
        rows.append(row)
    return arb_mat(rows)


def rigorous_negdef(S_ball, tag=""):
    """Prove every symmetric matrix in the ball matrix S_ball is positive definite.
    Returns (ok, lower_bound_on_lambda_min or None).
    Method (Lemma 6): s>0, L float with S_mid - sI ~ L L^T; E := S - sI - L L^T (ball);
    if ||E||_inf (rigorous upper bound) < s, then S >= (s - ||E||) I > 0."""
    N = S_ball.nrows()
    Sm = np.array([[float(S_ball[i, j].mid()) for j in range(N)] for i in range(N)])
    # exact diagonal scaling by powers of two (S > 0  <=>  D S D > 0, D = diag(2^-e_i))
    dg = np.diag(Sm)
    if np.any(dg <= 0):
        return False, None, float(np.min(dg))
    e = np.round(0.5 * np.log2(dg)).astype(int)
    sc = [arb(2.0 ** (-int(t))) for t in e]
    S_ball = arb_mat([[S_ball[i, j] * sc[i] * sc[j] for j in range(N)] for i in range(N)])
    Sm = np.array([[float(S_ball[i, j].mid()) for j in range(N)] for i in range(N)])
    Sm = 0.5 * (Sm + Sm.T)
    ev = np.linalg.eigvalsh(Sm)
    if ev[0] <= 0:
        return False, None, ev[0]
    s = 0.5 * ev[0]
    try:
        Lf = np.linalg.cholesky(Sm - s * np.eye(N))
    except np.linalg.LinAlgError:
        return False, None, ev[0]
    Lb = arb_mat([[arb(float(Lf[i, j])) for j in range(N)] for i in range(N)])
    sb = arb(s)
    E = S_ball - Lb * Lb.transpose()
    for i in range(N):
        E[i, i] = E[i, i] - sb
    worst = arb(0)
    for i in range(N):
        r = arb(0)
        for j in range(N):
            r += arb(E[i, j].abs_upper())
        worst = worst.max(r)
    marg = sb - worst
    ok = bool(marg > 0)
    return ok, (float(marg.lower()) if ok else None), ev[0]
