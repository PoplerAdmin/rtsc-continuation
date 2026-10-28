"""COMPUTER-ASSISTED (scalar, rigorous arb) checks for the analytic end lemmas.
Lemma G (small w):  rho(a,w) <= -1 + a*Lgv(w/2),  Lgv(x) = coth x - 1/x increasing.
   Check: a_lo * Lgv(W_L/2) < 1.
Lemma H (large w):  HS^2 <= (4 a^2/pi^2) P(Lw) e^{-2 Lw},  Lw = ln(w/(2 pi)),
   P(L) = (7/16) L^2 + b L + c0,  b = c1 + (e-1/2)/2,  c0 = c1 (1/2 + e) + 1/32,
   c1 = 1 + ln(2)/2,  e = (2 pi / W)^2  (valid for all w >= W >= 2 pi).
   Check: 2b - 7/8 > 0, 2 c0 - b > 0  (=> bound decreasing in L >= 0)  and bound(W) < 1.
Usage: python3 cert_ends.py p q WLnum/WLden W
"""
import sys
from flint import arb, fmpq, ctx
ctx.prec = 128
p, q = int(sys.argv[1]), int(sys.argv[2])
a = arb(fmpq(p, q))
WL = arb(fmpq(sys.argv[3]))
W = arb(fmpq(sys.argv[4]))
pi = arb.pi()

x = WL / 2
gersh = a * (x.coth() - 1 / x)
print("Lemma G: a*Lgv(WL/2) =", gersh, " < 1 :", bool(gersh < 1))

c1 = 1 + arb(2).log() / 2
e = (2 * pi / W) ** 2
b = c1 + (e - arb(fmpq(1, 2))) / 2
c0 = c1 * (arb(fmpq(1, 2)) + e) + arb(fmpq(1, 32))
LW = (W / (2 * pi)).log()
print("W >= 2 pi:", bool(W > 2 * pi))
print("2b-7/8 =", 2 * b - arb(fmpq(7, 8)), bool(2 * b - arb(fmpq(7, 8)) > 0))
print("2c0-b  =", 2 * c0 - b, bool(2 * c0 - b > 0))
P = arb(fmpq(7, 16)) * LW * LW + b * LW + c0
HS2 = 4 * a * a / (pi * pi) * P * (-2 * LW).exp()
print("Lemma H: HS^2 bound at W =", HS2, " < 1 :", bool(HS2 < 1), "  (HS <=", HS2.sqrt(), ")")
