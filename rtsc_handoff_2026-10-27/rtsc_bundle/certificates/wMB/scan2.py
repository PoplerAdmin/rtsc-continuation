import numpy as np, sys
from mb import acrit, funcs
from scipy.optimize import minimize
a0 = 4.33148145
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv)>1 else 0)
def shape(p):
    lw, lt, u11, u12, u22 = p
    w = np.exp(lw); N1 = 1/(1+np.exp(-lt)); N = np.array([N1, 1-N1])
    W = np.array([[u11**2, u12**2], [u12**2, u22**2]])
    A = W/N[:, None]
    return [(w, A)], N
def ratio(p, key='lmax', nmax=120):
    at, N = shape(p)
    if at[0][1].sum() < 1e-12: return 10.0
    a = acrit(at, N, nmax=nmax, tol=1e-7)
    return a*funcs(at, N)[key]/a0
best = []
for it in range(60):
    p = np.array([rng.uniform(0, 4), rng.uniform(-5, 5), *rng.uniform(0, 1, 3)])
    best.append((ratio(p), p))
best.sort(key=lambda t: t[0])
for r, p in best[:5]:
    print('rand', r, p)
r, p = best[0]
res = minimize(ratio, p, method='Nelder-Mead', options=dict(maxiter=400, xatol=1e-4, fatol=1e-7))
at, N = shape(res.x)
print('opt lmax ratio', res.fun, 'w=', at[0][0], 'N=', N, 'A=', at[0][1].tolist(), funcs(at, N))
print('rowmax ratio at that point', ratio(res.x, 'rowmax'))
