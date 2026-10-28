"""Sharpness (Lemma 10): for a1 > a*, exhibit w1 and finitely supported x (dyadic floats)
with x^T M+(a1, w1) x > 0, evaluated rigorously in arb.
usage: python3 sharp.py A1 [N]"""
import sys
import numpy as np
from scipy.optimize import minimize_scalar
from flint import arb, arb_mat
from core import lvals, cumsums, even_block


def lam(k, a, w):
    return 2 * a * w / (w ** 2 + (2 * np.pi * k) ** 2)


def spec_matrix(a, w, N):
    """float matrix built literally from SPEC.md's Q+ (incl. lam(0), which cancels)."""
    n = np.arange(N)
    Lm = lam(np.abs(np.subtract.outer(n, n)), a, w) + lam(np.add.outer(n, n) + 1, a, w)
    cum = np.concatenate([[0.0], np.cumsum(lam(np.arange(1, N), a, w))])
    return Lm - np.diag((2 * n + 1) + lam(0, a, w) + 2 * cum)


def main():
    A1 = sys.argv[1]
    N = int(sys.argv[2]) if len(sys.argv) > 2 else 200
    a1f = float(A1)
    top = lambda w: np.linalg.eigvalsh(spec_matrix(a1f, w, N))[-1]
    r = minimize_scalar(lambda w: -top(w), bounds=(7.0, 8.0), method="bounded",
                        options={"xatol": 1e-9})
    w1 = float(np.round(r.x * 2 ** 20) / 2 ** 20)       # dyadic w
    ev, V = np.linalg.eigh(spec_matrix(a1f, w1, N))
    x = V[:, -1]
    x = np.round(x * 2 ** 40) / 2 ** 40                 # dyadic x (exact in arb)
    a = arb(A1); wb = arb(w1)
    lv = lvals(a, wb, 2 * N)
    S = cumsums(lv, N)
    A = even_block(lv, S, N)
    # cross-check: arb block vs SPEC-literal float matrix
    diff = np.max(np.abs(np.array([[float(A[i, j].mid()) for j in range(N)] for i in range(N)])
                         - spec_matrix(a1f, w1, N)))
    X = arb_mat([[arb(float(t))] for t in x])
    Q = (X.transpose() * A * X)[0, 0]
    nx2 = sum(arb(float(t)) ** 2 for t in x)
    print("a1 = %s, w = %r (dyadic), N = %d, float top eig = %.3e" % (A1, w1, N, ev[-1]))
    print("max |arb block - SPEC float matrix| = %.2e" % diff)
    print("Q+(x) = %s ;  Q+(x)/|x|^2 = %s" % (Q, Q / nx2))
    ok = bool(Q > 0)
    print("VERDICT sharpness: %s" % ("CERTIFIED  (rho+(a1,w) > 0)" if ok else "NOT CERTIFIED"))
    np.savetxt("sharp_x_%s.txt" % A1, x, header="w=%r a1=%s ; x_n n=0..%d (dyadic, exact)" % (w1, A1, N - 1))
    return ok


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
