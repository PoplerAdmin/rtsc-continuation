# Independent check: build linearized multiband Eliashberg on the full Matsubara line
# n in [-Nh, Nh), with |w~_{n,i}| computed by direct (symmetric) summation, then compare the
# critical coupling scale with the reduced even-sector form of mb.py.
import numpy as np
from scipy.linalg import eigh
from mb import acrit
TWO_PI = 2*np.pi
def full_top(atoms, N, Nh):
    N = np.asarray(N, float); B = len(N)
    n = np.arange(-Nh, Nh); L = 2*Nh
    sgn = np.sign(2*n+1)
    lam = lambda d, w, A: A[..., None, None]*(2*w/(w*w+(TWO_PI*d)**2))  # T=1
    d = n[:, None]-n[None, :]
    Kfull = np.zeros((B, B, L, L))
    for w, A in atoms:
        Kfull += np.asarray(A, float)[:, :, None, None]*(2*w/(w*w+(TWO_PI*d)**2))[None, None]
    # |w~_{n,i}| = |(2n+1) + sum_{m,j} lam_ij(n-m) sgn(w_m)|  (pi T = pi; divide by pi T -> units of pi T)
    # Use a much longer m-range for the Z sum to avoid truncation error
    Mz = 20*Nh; mm = np.arange(-Mz, Mz); sg = np.sign(2*mm+1)
    wt = np.zeros((B, L))
    for w, A in atoms:
        A = np.asarray(A, float)
        g = 2*w/(w*w+(TWO_PI*(n[:, None]-mm[None, :]))**2)
        wt += A.sum(1)[:, None]*(g*sg[None, :]).sum(1)[None, :]
    wt = np.abs((2*n+1)[None, :] + wt)
    # equation: wt_{n,i} phi_{n,i} = sum_{m,j} lam_ij(n-m) phi_{m,j}; symmetric form with N_i:
    M = np.zeros((B*L, B*L))
    for i in range(B):
        for j in range(B):
            M[i*L:(i+1)*L, j*L:(j+1)*L] = np.sqrt(N[i]/N[j])*Kfull[i, j]
    M -= np.diag(wt.ravel())
    return eigh(M, eigvals_only=True, subset_by_index=[B*L-1, B*L-1])[0]
def acrit_full(atoms, N, Nh):
    f = lambda a: full_top([(w, a*np.asarray(A)) for w, A in atoms], N, Nh)
    lo, hi = 0.01, 1.0
    while f(hi) < 0: lo = hi; hi *= 2
    for _ in range(40):
        mid = (lo+hi)/2
        if f(mid) >= 0: hi = mid
        else: lo = mid
    return hi
if __name__ == "__main__":
    N = [0.3, 0.7]
    atoms = [(5.0, np.array([[1.0, 0.4], [0.4*0.3/0.7, 0.2]])), (12.0, np.array([[0.3, 0.5], [0.5*0.3/0.7, 0.6]]))]
    for Nh in [60, 120]:
        print('full line Nh=%d' % Nh, acrit_full(atoms, N, Nh))
    print('reduced even form   ', acrit(atoms, N, nmax=120, tol=1e-12))
