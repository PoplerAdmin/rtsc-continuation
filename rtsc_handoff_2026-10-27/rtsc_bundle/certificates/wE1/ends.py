"""Analytic end ranges (Lemmas 7, 8, 9 of CERT.md), evaluated with arb.
usage: python3 ends.py A0 ETA [WLO WHI]"""
import sys
from flint import arb


def main():
    A0 = sys.argv[1]; ETA = sys.argv[2]
    WLO = arb(sys.argv[3]) if len(sys.argv) > 3 else arb(1)
    WHI = arb(sys.argv[4]) if len(sys.argv) > 4 else arb(400)
    a = arb(A0); eta = arb(ETA)
    ok = True
    # Lemma 7 (small w): row sums of unit-weight Schur test <= 2L - 1 = a*Lang(w/2) - 1,
    # Lang(u) = coth u - 1/u increasing  =>  for 0 < w <= WLO bounded by value at WLO.
    small = a * ((WLO / 2).coth() - 2 / WLO) - 1
    print("small w <= %s : sup Q/|x|^2 <= %s" % (WLO, small))
    ok &= bool(small < -eta)
    # Lemma 8 (large w): weights h0=H0, h1=H1, h_n=1 (n>=2); l(k) <= 2a/W, T(j) <= L <= a/2.
    H0 = arb(8); H1 = arb(5) / 2
    lmax = 2 * a / WHI
    R0 = -1 + lmax * (1 + 2 * H1 / H0) + a / H0
    R1 = -3 + lmax * (1 + 2 * H0 / H1) + a / H1
    Rn = a - 5 + 2 * lmax * (H0 + H1 - 2)
    print("large w >= %s : R0 <= %s, R1 <= %s, Rn(n>=2) <= %s" % (WHI, R0, R1, Rn))
    for R in (R0, R1, Rn):
        ok &= bool(R < -eta)
    # also require a0 < 5 for Rn argument consistency (a - 5 < 0) -- implied by Rn < 0.
    # Lemma 9 (odd sector): sup Q-/|x|^2 <= -1 for all a, w  (no numerics needed).
    print("odd sector: Q-(x) <= -sum (2n+1) x_n^2 <= -|x|^2 for all a>0, w>0 (analytic)")
    print("VERDICT ends: %s" % ("CERTIFIED" if ok else "NOT CERTIFIED"))
    return ok


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
