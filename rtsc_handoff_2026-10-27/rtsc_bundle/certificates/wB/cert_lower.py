"""COMPUTER-ASSISTED LEMMA L (sharpness).
For the Einstein spectrum a2F = a * delta(omega - w)  at T = 1 (a = p/q, w = r/s rational),
mu* = 0, exhibit a finitely supported x (float vector, exactly representable) with
    Q^+_{T=1}(x) > 0     (computed exactly up to rigorous arb rounding).
Hence T = 1 is unstable, Tc >= 1 = A/a, so the best constant c >= 1/a.
Usage: python3 cert_lower.py p q r s [N]
"""
import sys
import numpy as np
from flint import arb, arb_mat, fmpq, ctx
from rig_common import lam, build

ctx.prec = 128


def float_vec(a, w, N):
    n = np.arange(N)
    lv = lambda k: 2 * a * w / (w ** 2 + (2 * np.pi * k) ** 2)
    K = lv(np.abs(np.subtract.outer(n, n))) + lv(np.add.outer(n, n) + 1)
    np.fill_diagonal(K, lv(2 * n + 1))
    cum = np.concatenate([[0], np.cumsum(lv(np.arange(1, N)))])
    Mf = K - np.diag(2 * n + 1 + 2 * cum)
    ev, V = np.linalg.eigh(Mf)
    x = V[:, -1]
    return x * np.sign(x[0]), ev[-1]


if __name__ == "__main__":
    p, q, r, s = map(int, sys.argv[1:5])
    N = int(sys.argv[5]) if len(sys.argv) > 5 else 300
    a, w = arb(fmpq(p, q)), arb(fmpq(r, s))
    x, rf = float_vec(p / q, r / s, N)
    print("float top eigenvalue of section:", rf)
    B = build(lam, a, w, N, True)          # exact finite section; x supported in 0..N-1
    X = arb_mat([[arb(float(v))] for v in x])
    Qx = (X.transpose() * B * X)[0, 0]
    nx = (X.transpose() * X)[0, 0]
    print("a =", fmpq(p, q), " w =", fmpq(r, s), " N =", N)
    print("Q(x)      =", Qx)
    print("Q(x)/|x|^2=", Qx / nx)
    print("CERTIFIED Q(x) > 0:", bool(Qx > 0))
