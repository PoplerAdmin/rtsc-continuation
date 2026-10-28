"""COMPUTER-ASSISTED LEMMA C (middle range).
Certifies: for a = A_LO and every w in [W_L, W_R], every finitely supported x,
    Q_{a,w}(x) <= 0      (even sector, Einstein spectrum, mu* = 0, T = 1).
Method per w-interval I = [w0-h, w0+h]  (see PROOF.md, Lemmas T and C):
  Q(x) <= u^T (B(w) + G(w)) u + tau |u|^2,   u = x restricted to {0..N-1}
  B(w) = B(w0) + t B'(w0) + t^2 R,  R in ball (1/2) B''(I)
  => Q(x) <= max_{s=+-1} u^T [B0 + s h B1 + G(I)] u + sum_i u_i^2 (h^2 rowabs(R)_i + tau)
and NSD of the two finite ball matrices is certified by float Cholesky + rigorous residual.
Usage: python3 cert_upper.py p q WL WR [N L]   (a = p/q, WL, WR rationals 'num/den')
"""
import sys, time
from flint import arb, fmpq
from rig_common import *


def check_interval(a, wl, wr, N, L):
    w0 = arb((wl + wr) / 2)
    h = arb((wr - wl) / 2)
    wI = arb(wl).union(arb(wr))
    B0 = build(lam, a, w0, N, True)
    B1 = build(dlam, a, w0, N, False)
    R = build(d2lam, a, wI, N, False) * arb(fmpq(1, 2))
    G, S = ghead(a, wI, N, L)
    tau = gtail(a, arb(wr), N, L, S)
    Rrow = rowabs_upper(R)
    extra = [h * h * Rrow[i] + tau for i in range(N)]
    infos = []
    for s in (1, -1):
        Mb = B0 + B1 * (s * h) + G
        ok, info = certify_nsd(Mb, extra)
        infos.append(info)
        if not ok:
            return False, infos
    return True, infos


def run(a, WL, WR, N=60, L=300, h0=fmpq(1, 4), hmin=fmpq(1, 2 ** 14), log=None):
    stack = []
    x = WL
    while x < WR:
        y = min(x + h0, WR)
        stack.append((x, y))
        x = y
    stack.reverse()
    done = []
    t0 = time.time()
    while stack:
        wl, wr = stack.pop()
        ok, infos = check_interval(a, wl, wr, N, L)
        if ok:
            done.append((wl, wr))
            msg = "OK   [%s, %s]  %s" % (wl, wr, " | ".join(infos))
        else:
            if wr - wl <= hmin:
                print("FAIL at minimal width", wl, wr, infos, flush=True)
                return False, done
            mid = (wl + wr) / 2
            stack.append((mid, wr))
            stack.append((wl, mid))
            msg = "split [%s, %s]  %s" % (wl, wr, infos[-1])
        print("%7.1fs %s" % (time.time() - t0, msg), flush=True)
    # coverage check
    done.sort()
    assert done[0][0] == WL and done[-1][1] == WR
    for (l1, r1), (l2, r2) in zip(done, done[1:]):
        assert r1 == l2
    print("CERTIFIED: [%s, %s] covered by %d intervals, a = %s, N = %d, L = %d"
          % (WL, WR, len(done), a, N, L), flush=True)
    return True, done


if __name__ == "__main__":
    p, q = int(sys.argv[1]), int(sys.argv[2])
    WL, WR = fmpq(sys.argv[3]), fmpq(sys.argv[4])
    N = int(sys.argv[5]) if len(sys.argv) > 5 else 60
    L = int(sys.argv[6]) if len(sys.argv) > 6 else 300
    a = arb(fmpq(p, q))
    run(a, WL, WR, N, L)
