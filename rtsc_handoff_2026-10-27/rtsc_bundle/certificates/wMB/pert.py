import numpy as np
from mb import acrit, funcs
a0 = 4.33148145
w = 7.5749
for N1 in [0.5, 0.1, 0.9]:
    N = np.array([N1, 1-N1])
    for b in [0.3, 0.1, 0.03, 0.01]:
        # W = [[1, b],[b, 0]] (symmetric N_i A_ij); intraband only in band 1
        W = np.array([[N1*1.0, b*np.sqrt(N1*(1-N1))], [b*np.sqrt(N1*(1-N1)), 0.0]])
        A = W/N[:, None]
        a = acrit([(w, A)], N, nmax=150, tol=1e-11)
        f = funcs([(w, A)], N)
        print(N1, b, 'ratio_lmax-1 = %.3e' % (a*f['lmax']/a0-1), ' r2/kappa=%.3f' % (A.sum(1)[1]/f['lmax']))
