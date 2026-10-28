# Single-band Einstein form with NEGATIVE coupling: Q^E_{-a}(x) = -C(x) - a G_w(x).
# Find a_-^* = inf_w crit a (even and odd sectors).
import numpy as np
from scipy.linalg import eigh
TWO_PI = 2*np.pi
def G(w, nmax, sector=+1):
    k = np.arange(2*nmax+2); g = 2*w/(w*w+(TWO_PI*k)**2)
    n = np.arange(nmax)
    Bm = g[np.abs(n[:, None]-n[None, :])] + sector*g[n[:, None]+n[None, :]+1]
    D = g[0] + 2*np.concatenate([[0], np.cumsum(g[1:nmax])])
    return Bm - np.diag(D)
def top(a, w, nmax, sector):
    n = np.arange(nmax)
    M = -np.diag(2*n+1.0) + a*G(w, nmax, sector)
    return eigh(M, eigvals_only=True, subset_by_index=[nmax-1, nmax-1])[0]
def crit(sign, w, nmax, sector):
    lo, hi = 0.0, 1.0
    while top(sign*hi, w, nmax, sector) < 0:
        lo = hi; hi *= 2
        if hi > 1e4: return np.inf
    for _ in range(50):
        mid = (lo+hi)/2
        if top(sign*mid, w, nmax, sector) >= 0: hi = mid
        else: lo = mid
    return hi
ws = np.exp(np.linspace(np.log(0.05), np.log(2000), 90))
for sector in (+1, -1):
    vals = [crit(-1, w, 300 if w < 300 else 800, sector) for w in ws]
    i = int(np.argmin(vals))
    print('sector', sector, 'negative coupling: min crit a =', vals[i], 'at w =', ws[i])
    print('   samples:', [(round(w, 2), round(v, 3)) for w, v in list(zip(ws, vals))[::10]])
    print('   top eig at a=a0 max over w:', max(top(-4.33148145, w, 300, sector) for w in ws[ws < 300]))
