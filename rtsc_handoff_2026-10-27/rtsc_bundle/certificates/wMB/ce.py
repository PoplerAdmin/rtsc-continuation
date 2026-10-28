# Counterexamples and an illustrative two-band example (float; mu* = 0 unless stated).
import numpy as np
from scipy.linalg import eigh
from mb import build, topeig, acrit
a0 = 4.33148145; wstar = 7.574866
def Tc(atoms_phys, N, nmax=160):
    """atoms_phys: list of (Omega[K], Amat[K]). Bisection on T for top eigenvalue = 0."""
    f = lambda T: topeig([(W/T, A/T) for W, A in atoms_phys], N, nmax)
    lo, hi = 1e-3, 1.0
    while f(hi) >= 0: lo = hi; hi *= 2
    if f(lo) < 0: return 0.0
    for _ in range(50):
        mid = (lo+hi)/2
        if f(mid) >= 0: lo = mid
        else: hi = mid
    return lo
def funcs(atoms, N):
    N = np.asarray(N, float); s = np.sqrt(np.outer(N, 1/N)); A = sum(np.asarray(M, float) for _, M in atoms)
    K = s*A; R = np.diag(A.sum(1))
    return dict(avg=(N[:, None]*A).sum()/N.sum(), kappa=np.linalg.eigvalsh(K)[-1],
                LamA=np.linalg.eigvalsh((R+K)/2)[-1], rowmax=A.sum(1).max(), maxdiag=np.diag(A).max())
def report(name, atoms, N):
    T = Tc(atoms, N); f = funcs(atoms, N)
    print('%s: Tc = %.3f K | bounds Phi/a0:' % (name, T), ', '.join('%s %.3f' % (k, v/a0) for k, v in f.items()))
    return T, f
A = 1000.0; Om = wstar*A/a0     # extremal single-band Einstein: Tc = A/a0 = 230.87 K
print('--- CE1: decoupled bands, small-DOS band carries all coupling (average fails)')
for eps in [0.2, 0.05, 0.01]:
    report(' eps=%.2f' % eps, [(Om, np.array([[A, 0], [0, 0.0]]))], [eps, 1-eps])
print('--- CE1b: same but connected (interband b=50 K on band-1 row)')
eps = 0.05; b = 50.0
report(' eps=0.05', [(Om, np.array([[A, b], [b*eps/(1-eps), 0.0]]))], [eps, 1-eps])
print('--- CE2: pure interband, equal DOS (intraband-only functionals fail; kappa sharp)')
report(' ', [(Om, np.array([[0, A], [A, 0.0]]))], [0.5, 0.5])
print('--- CE3: equal row sums, unequal DOS: kappa = rowmax = LamA, bound attained')
N1 = 0.2; x = 300.0  # A = [[A-x, x],[x*N1/N2, A - x*N1/N2]]
N2 = 1-N1; report(' ', [(Om, np.array([[A-x, x], [x*N1/N2, A-x*N1/N2]]))], [N1, N2])
print('--- CE4: Coulomb matrix with nonneg entries but N_i mu_ij not PSD, NO phonons')
# pure interband repulsion mu, equal DOS, block cutoff omega_c: instability iff 2 mu sum_{n<Nc} 1/(2n+1) >= 1
mu = 0.2; wc = 1000.0
def unstable(T):
    Nc = int(np.floor(wc/(2*np.pi*T) - 0.5)) + 1   # number of n>=0 with (2n+1) pi T < wc
    return Nc >= 1 and 2*mu*np.sum(1/(2*np.arange(Nc)+1)) >= 1
Ts = np.linspace(0.5, 40, 4000); Tun = max([T for T in Ts if unstable(T)], default=0)
print('  mu_12=mu_21=%.2f, mu_ii=0, omega_c=%g K, alpha2F=0: unstable up to T ~ %.2f K (A = 0)' % (mu, wc, Tun))
print('--- Illustrative MgB2-like two-band Einstein model (NOT fitted literature values)')
Nm = np.array([0.43, 0.57]); lam = np.array([[1.0, 0.20], [0.20*0.43/0.57, 0.45]])
Omg = np.array([[700.0, 800.0], [800.0, 800.0]])
atoms = []
for (i, j) in [(0, 0), (1, 1)]:
    M = np.zeros((2, 2)); M[i, j] = lam[i, j]*Omg[i, j]/2; atoms.append((Omg[i, j], M))
M = np.zeros((2, 2)); M[0, 1] = lam[0, 1]*800/2; M[1, 0] = lam[1, 0]*800/2; atoms.append((800.0, M))
T2, f = report(' two-band', atoms, Nm)
# isotropic average: single band with N-weighted spectrum
avg_atoms = [(W, np.array([[(Nm[:, None]*Mx).sum()/Nm.sum()]])) for W, Mx in atoms]
T1 = Tc(avg_atoms, [1.0]); print('  isotropic-average single band Tc = %.3f K  (multiband/average = %.3f)' % (T1, T2/T1))
