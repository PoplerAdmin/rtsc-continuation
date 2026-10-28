"""Certify one w-cell [w1, w2] for the even sector:  Q+_w(x) <= -eta |x|^2  for all w in the
cell and all finitely supported x.   (Lemmas 2-6 of CERT.md.)"""
import time
from flint import arb, arb_mat
from core import (ball_interval, lvals, dlvals, d2l_abs_upper, L_total, T_upper,
                  cumsums, even_block, rigorous_negdef)


def certify_cell(a, eta, w1, w2, N=150, K=600, verbose=False):
    t0 = time.time()
    wc = (w1 + w2) / 2.0          # exact: w1, w2 dyadic with few bits
    h = (w2 - w1) / 2.0
    assert wc - h == w1 and wc + h == w2
    wcb = arb(wc)
    hb = arb(h)
    wb = ball_interval(w1, w2)

    # --- Taylor data of the N x N block at the centre (Lemma 5)
    lc = lvals(a, wcb, 2 * N)
    Sc = cumsums(lc, N)
    dlc = dlvals(a, wcb, 2 * N)
    dSc = cumsums(dlc, N)
    A0 = even_block(lc, Sc, N)
    A1 = even_block(dlc, dSc, N, deriv=True)
    mu = d2l_abs_upper(a, wb, 2 * N)
    Smu = cumsums(mu, N)
    half_h2 = hb * hb / 2
    rem = []
    for n in range(N):
        u = mu[2 * n + 1] + 2 * Smu[n]
        for m in range(N):
            if m != n:
                u += mu[abs(n - m)] + mu[n + m + 1]
        rem.append(arb((half_h2 * u).upper()))

    # --- tail and coupling, valid for every w in the cell (Lemmas 2-4)
    kmax = 2 * N + K + 2
    lb = lvals(a, wb, kmax)
    Sb = cumsums(lb, kmax)
    Lb = L_total(a, wb)

    def Tu(j):
        return T_upper(a, wb, j, lb, Sb, Lb)

    gam = []
    for m in range(N, N + K):
        g = arb((2 * m + 1 - 2 * Tu(m + 1) - eta).lower())
        if not g > 0:
            return False, {"why": "gamma<=0"}
        gam.append(g)
    TNK1 = Tu(N + K + 1)
    gmin = arb((2 * (N + K) + 1 - 2 * TNK1 - eta).lower())
    if not gmin > 0:
        return False, {"why": "gmin<=0"}
    # G1 = sum_{m=N}^{N+K-1} b_m b_m^T / gamma_m
    Brows, Bgrows = [], []
    for n in range(N):
        r, rg = [], []
        for j in range(K):
            m = N + j
            v = lb[m - n] + lb[n + m + 1]
            r.append(v)
            rg.append(v / gam[j])
        Brows.append(r)
        Bgrows.append(rg)
    B1 = arb_mat(Brows)
    G1 = arb_mat(Bgrows) * B1.transpose()
    # diagonal Schur bound for m >= N+K
    cmax = Tu(K + 1) + TNK1
    dfac = cmax / gmin
    d = [arb((dfac * (Tu(N + K - n) + Tu(N + K + n + 1))).upper()) for n in range(N)]

    info = {"N": N, "K": K, "w1": w1, "w2": w2, "rem0": float(rem[0]), "d0": float(d[0]),
            "G1_00": float(G1[0, 0].mid())}
    worst_margin = None
    for sgn in (+1, -1):
        P = A0 + (hb * sgn) * A1 + G1
        for n in range(N):
            P[n, n] = P[n, n] + rem[n] + d[n] + eta
        ok, marg, ev0 = rigorous_negdef(-P)
        info["ev_min_float_%+d" % sgn] = float(ev0)
        if not ok:
            info["time"] = time.time() - t0
            return False, info
        worst_margin = marg if worst_margin is None else min(worst_margin, marg)
    info["margin"] = worst_margin
    info["time"] = time.time() - t0
    return True, info


if __name__ == "__main__":
    import sys
    a = arb("4.33148")
    eta = arb("1e-7")
    for (w1, w2) in [(7.5625, 7.578125), (7.5, 7.625), (3.0, 3.5), (200.0, 232.0)]:
        print(w1, w2, certify_cell(a, eta, w1, w2, N=int(sys.argv[1]) if len(sys.argv) > 1 else 150))
