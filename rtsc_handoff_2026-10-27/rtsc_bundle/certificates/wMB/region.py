# Mismatched single-band form F(p,r;w) = -(2n+1) + p(B - s I) - r(D - s I), s = sigma_pi = 1/sinh(w/2)
# (unit-coupling kernel l(k) = 2w/(w^2+(2 pi k)^2), T=1). p*(rho) = inf_w crit p with r = rho p.
import numpy as np
from scipy.linalg import eigh
TWO_PI = 2*np.pi
def F(p, r, w, nmax=200):
    k = np.arange(2*nmax+2); g = 2*w/(w*w+(TWO_PI*k)**2)
    n = np.arange(nmax)
    Bm = g[np.abs(n[:, None]-n[None, :])] + g[n[:, None]+n[None, :]+1]
    D = g[0] + 2*np.concatenate([[0], np.cumsum(g[1:nmax])])
    s = 1/np.sinh(w/2)
    M = p*(Bm - s*np.eye(nmax)) - np.diag((2*n+1) + r*(D - s))
    return eigh(M, eigvals_only=True, subset_by_index=[nmax-1, nmax-1])[0]
def pcrit(rho, w, nmax=200):
    lo, hi = 0.0, 1.0
    while F(hi, rho*hi, w, nmax) < 0:
        lo = hi; hi *= 2
        if hi > 1e5: return np.inf
    for _ in range(45):
        mid = (lo+hi)/2
        if F(mid, rho*mid, w, nmax) >= 0: hi = mid
        else: lo = mid
    return hi
if __name__ == "__main__":
    ws = np.exp(np.linspace(np.log(0.2), np.log(300), 60))
    for rho in [0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0]:
        vals = [pcrit(rho, w, 200 if w < 100 else 400) for w in ws]
        i = int(np.argmin(vals))
        print(rho, 'p* =', round(vals[i], 5), 'at w =', round(ws[i], 3), ' p(w->large)=', round(vals[-1], 4))
