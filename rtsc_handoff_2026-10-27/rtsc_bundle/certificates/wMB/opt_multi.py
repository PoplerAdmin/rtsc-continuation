import numpy as np, sys
from mb import acrit, funcs
from scipy.optimize import minimize
a0 = 4.33148145
key = sys.argv[1]; seed = int(sys.argv[2]); nat = int(sys.argv[3])
rng = np.random.default_rng(seed)
def shape(p):
    N1 = 1/(1+np.exp(-p[0])); N = np.array([N1, 1-N1]); atoms = []
    for k in range(nat):
        lw, u11, u12, u22 = p[1+4*k:5+4*k]
        W = np.array([[u11**2, u12**2], [u12**2, u22**2]])
        atoms.append((np.exp(lw), W/N[:, None]))
    return atoms, N
def ratio(p):
    at, N = shape(p)
    if sum(A.sum() for _, A in at) < 1e-12: return 10.
    return acrit(at, N, nmax=110, tol=1e-7)*funcs(at, N)[key]/a0
best = None
for s in range(25):
    p = [rng.uniform(-4, 4)] + sum([[rng.uniform(0.5, 3.5), *rng.uniform(0, 1, 3)] for _ in range(nat)], [])
    r = ratio(np.array(p))
    if best is None or r < best[0]: best = (r, np.array(p))
res = minimize(ratio, best[1], method='Nelder-Mead', options=dict(maxiter=600, fatol=1e-8, xatol=1e-5))
at, N = shape(res.x)
print(key, 'start', best[0], 'opt', res.fun, 'N', N)
for w, A in at: print('  w', w, 'A', A.tolist())
print(' ', funcs(at, N))
