# Numerical check of the exact identity (Theorem B):
# Q(x) = (1/kappa) sum_w sum_ij m_ij(w) [Q^E_{kappa,w}(psi_ij) + Q^E_{-kappa,w}(chi_ij)]
#        - 1/2 sum_w sum_ij W_ij(w) sum_n D_{w,n} (phi_ni - phi_nj)^2
import numpy as np
from mb import build
TP = 2*np.pi
def unitG(w, nmax):
    k = np.arange(2*nmax+2); g = 2*w/(w*w+(TP*k)**2); n = np.arange(nmax)
    Bm = g[np.abs(n[:, None]-n[None, :])] + g[n[:, None]+n[None, :]+1]
    D = g[0] + 2*np.concatenate([[0], np.cumsum(g[1:nmax])])
    return Bm - np.diag(D), D
rng = np.random.default_rng(7)
B, nat, nmax = 3, 3, 40
N = rng.dirichlet(np.ones(B)); atoms = []
for _ in range(nat):
    U = rng.exponential(1, (B, B)); W = (U+U.T)/2
    atoms.append((rng.uniform(1, 20), W/N[:, None]))
M = build(atoms, N, nmax)
Atot = sum(A for _, A in atoms)
ev, V = np.linalg.eig(Atot); i = np.argmax(ev.real); kappa = ev[i].real; v = np.abs(V[:, i].real)
Ks = np.sqrt(np.outer(N, 1/N))*Atot
print('kappa', kappa, 'lmax(K_sym)', np.linalg.eigvalsh(Ks)[-1])
n = np.arange(nmax); C = lambda y: ((2*n+1)*y*y).sum()
for trial in range(3):
    x = rng.normal(size=(B, nmax))           # x[i,n]
    Q = x.ravel() @ M @ x.ravel()
    phi = x/np.sqrt(N)[:, None]; th = phi/v[:, None]
    rhs = 0.0
    for w, A in atoms:
        G, D = unitG(w, nmax); W = N[:, None]*A; m = W*np.outer(v, v)
        for a in range(B):
            for b in range(B):
                psi = (th[a]+th[b])/2; chi = (th[a]-th[b])/2
                rhs += m[a, b]/kappa*((-C(psi)+kappa*psi@G@psi) + (-C(chi)-kappa*chi@G@chi))
                rhs -= 0.5*W[a, b]*(D*(phi[a]-phi[b])**2).sum()
    print('Q = %.12f   identity RHS = %.12f   diff = %.2e' % (Q, rhs, Q-rhs))
