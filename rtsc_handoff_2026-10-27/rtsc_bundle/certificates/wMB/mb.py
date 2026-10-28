# Float exploration of multiband linearized Eliashberg (even sector), T = 1 units.
import numpy as np
from scipy.linalg import eigh
TWO_PI = 2*np.pi

def gk(w, kmax):
    k = np.arange(kmax+1)
    return 2*w/(w*w + (TWO_PI*k)**2)

def build(atoms, N, nmax, ext=0):
    """atoms: list of (w, Amat) with Amat[i,j] = int alpha^2F_ij over that atom (units of T).
    N: DOS weights. Returns symmetric matrix of even-sector form in x = sqrt(N) phi variables."""
    N = np.asarray(N, float); B = len(N)
    s = np.sqrt(np.outer(N, 1/N))           # s_ij = sqrt(N_i/N_j)
    n = np.arange(nmax)
    M = np.zeros((B*nmax, B*nmax))
    diag = np.zeros((B, nmax)); diag += (2*n+1)[None, :]
    for (w, A) in atoms:
        A = np.asarray(A, float)
        K = s*A
        r = A.sum(1)
        g = gk(w, 2*nmax+2)
        Bm = g[np.abs(n[:, None]-n[None, :])] + g[n[:, None]+n[None, :]+1]
        M += np.kron(K, Bm)    # index (i,n) -> i*nmax+n
        S = np.concatenate([[0], np.cumsum(g[1:nmax])])   # S(n)=sum_{k=1}^n g(k)
        D = g[0] + 2*S
        diag += r[:, None]*D[None, :]
    M -= np.diag(diag.ravel())
    return M

def topeig(atoms, N, nmax):
    M = build(atoms, N, nmax)
    return eigh(M, eigvals_only=True, subset_by_index=[M.shape[0]-1, M.shape[0]-1])[0]

def acrit(atoms_shape, N, nmax=200, lo=1e-3, hi=None, tol=1e-9):
    """smallest scale a with topeig(a*shape) >= 0 (bisection, assumes single crossing)."""
    f = lambda a: topeig([(w, a*np.asarray(A)) for (w, A) in atoms_shape], N, nmax)
    if hi is None:
        hi = 1.0
        while f(hi) < 0:
            hi *= 2
            if hi > 1e6: return np.inf
    while hi-lo > tol*hi:
        mid = 0.5*(lo+hi)
        if f(mid) >= 0: hi = mid
        else: lo = mid
    return hi

def funcs(atoms_shape, N):
    N = np.asarray(N, float); s = np.sqrt(np.outer(N, 1/N))
    Atot = sum(np.asarray(A, float) for _, A in atoms_shape)
    K = s*Atot
    lmax = np.linalg.eigvalsh(K)[-1]
    rowmax = Atot.sum(1).max()
    avg = (N[:, None]*Atot).sum()/N.sum()
    intl = sum(np.linalg.eigvalsh(s*np.asarray(A, float))[-1] for _, A in atoms_shape)
    return dict(lmax=lmax, rowmax=rowmax, avg=avg, int_lmax=intl)

if __name__ == "__main__":
    # sanity: single band, w = 7.5749, expect a* = 4.33148
    print(acrit([(7.5749, [[1.0]])], [1.0], nmax=150))
    print(acrit([(7.5749, np.ones((2, 2))/2)], [1.0, 1.0], nmax=150))  # uniform: rows sum to 1
