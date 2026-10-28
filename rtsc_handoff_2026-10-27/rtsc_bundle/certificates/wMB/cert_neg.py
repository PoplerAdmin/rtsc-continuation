"""Rigorous (Arb ball arithmetic) check of the scalar inequalities behind Proposition N
(negative-coupling single-band Einstein form is <= 0 at a0 for every w > 0):

  (P0)  a0 * g0(w) <= 1,  g0(w) = 2/w - (1 + (w/2pi)^2)/sinh(w/2)
  (P1)  a0 * (c(w) + 2 l1(w)) <= 3,  c(w) = 2/w - 1/sinh(w/2),  l1(w) = 2w/(w^2+4pi^2)
  (P2)  a0 * tanh(w/4) <= 5   (trivial: tanh < 1 < 5/a0)

Analytic facts used (proved in PROOF.md):
  * c(w) <= tanh(w/4) <= w/4 ;  c + 2 l1 <= tanh(w/4)   (since c + 2L = tanh(w/4), L >= l1)
  * h(w) = (1+(w/2pi)^2)/sinh(w/2) is decreasing; 2/w decreasing; 1/sinh(w/2) decreasing
  * g0 <= c ;  g0 <= 2/w ;  c + 2 l1 <= 6/w.
On a cell [al, be] monotonicity gives the upper enclosures
  g0 <= 2/al - h(be),   c + 2 l1 <= 2/al - 1/sinh(be/2) + 4 be/(al^2 + 4 pi^2).
"""
from flint import arb, ctx
ctx.prec = 128
a0 = arb("4.33148145")
pi = arb.pi()
def h(w): return (1 + (w/(2*pi))**2)/(w/2).sinh()
ok = True
# ---- (P0)
# w in (0, 0.9]: a0*g0 <= a0*w/4 <= a0*0.225
v = a0*arb("0.9")/4; print("P0 small-w bound a0*0.9/4 =", v.str(5)); ok &= bool(v < 1)
# w >= 8.7: a0*g0 <= a0*2/w <= a0*2/8.7
v = a0*2/arb("8.7"); print("P0 large-w bound a0*2/8.7 =", v.str(5)); ok &= bool(v < 1)
# cells on [0.9, 8.7]
lo = arb("0.9"); hi = arb("8.7"); ncell = 780; worst = arb(-10)
for k in range(ncell):
    al = lo + (hi-lo)*k/ncell; be = lo + (hi-lo)*(k+1)/ncell
    ub = a0*(2/al - h(be))
    if not bool(ub < 1): ok = False; print("P0 FAIL cell", al, be, ub)
    worst = ub if bool(ub > worst) else worst
print("P0 max cell upper bound on [0.9,8.7]:", worst.str(6))
# ---- (P1)
t = (arb("3.4")/4).tanh(); v = a0*t; print("P1 w<=3.4: a0*tanh(0.85) =", v.str(6)); ok &= bool(v < 3)
v = a0*6/arb(40); print("P1 w>=40: a0*6/40 =", v.str(5)); ok &= bool(v < 3)
lo = arb("3.4"); hi = arb(40); ncell = 2000; worst = arb(-10)
for k in range(ncell):
    al = lo + (hi-lo)*k/ncell; be = lo + (hi-lo)*(k+1)/ncell
    ub = a0*(2/al - 1/(be/2).sinh() + 4*be/(al*al + 4*pi*pi))
    if not bool(ub < 3): ok = False; print("P1 FAIL cell", al, be, ub)
    worst = ub if bool(ub > worst) else worst
print("P1 max cell upper bound on [3.4,40]:", worst.str(6))
# ---- (P2)
ok &= bool(a0 < 5); print("P2: a0 < 5:", bool(a0 < 5))
print("VERDICT Proposition N scalar inequalities:", "PASS" if ok else "FAIL")
