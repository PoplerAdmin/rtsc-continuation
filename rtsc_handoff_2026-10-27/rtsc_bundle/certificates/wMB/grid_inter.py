import numpy as np
from mb import acrit, funcs
a0 = 4.33148145
# pure interband: W12 = 1, N = (N1, 1-N1). lmax(K)=sqrt(A12*A21)
for N1 in [0.5, 0.2, 0.05, 0.01, 1e-3]:
    N = np.array([N1, 1-N1])
    row = []
    for w in [1, 3, 7.57, 15, 40]:
        A = np.array([[0, 1/N1], [1/(1-N1), 0]])
        a = acrit([(w, A)], N, nmax=150, tol=1e-6)
        f = funcs([(w, A)], N)
        row.append((w, round(a*f['lmax']/a0, 4), round(a*f['rowmax']/a0, 3)))
    print(N1, row)
