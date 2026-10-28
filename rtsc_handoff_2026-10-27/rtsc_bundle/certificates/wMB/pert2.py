import numpy as np
from mb import acrit, funcs
from scipy.optimize import minimize
a0 = 4.33148145
def ratio(p):
    N1 = 1/(1+np.exp(-p[0])); b = p[1]; c = p[2]; w = np.exp(p[3])
    N = np.array([N1, 1-N1])
    W = np.array([[N1, b*np.sqrt(N1*(1-N1))], [b*np.sqrt(N1*(1-N1)), c*c*(1-N1)]])
    A = W/N[:, None]
    a = acrit([(w, A)], N, nmax=130, tol=1e-10)
    return a*funcs([(w, A)], N)['lmax']/a0 - 1
best = None
for N1 in [0.8, 0.9, 0.95, 0.98]:
    for b in [0.2, 0.3, 0.5, 0.8]:
        r = ratio([np.log(N1/(1-N1)), b, 0.0, np.log(7.5749)])
        print(N1, b, '%.3e' % r)
        if best is None or r < best[0]: best = (r, [np.log(N1/(1-N1)), b, 0.0, np.log(7.5749)])
res = minimize(ratio, best[1], method='Nelder-Mead', options=dict(maxiter=300, fatol=1e-10, xatol=1e-6))
print(res.fun, res.x, 1/(1+np.exp(-res.x[0])), np.exp(res.x[3]))
