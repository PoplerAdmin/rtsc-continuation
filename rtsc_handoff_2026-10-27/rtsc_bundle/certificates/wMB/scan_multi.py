# random multi-band, multi-frequency tests of candidate functionals (float)
import numpy as np, sys
from mb import acrit, funcs
a0 = 4.33148145
rng = np.random.default_rng(int(sys.argv[1]))
B = int(sys.argv[2]); nat = int(sys.argv[3]); trials = int(sys.argv[4])
worst = {}
for t in range(trials):
    N = rng.dirichlet(np.ones(B)*0.5) + 1e-6
    atoms = []
    for k in range(nat):
        w = np.exp(rng.uniform(np.log(1.0), np.log(60)))
        U = rng.exponential(1, (B, B)) * (rng.random((B, B)) < 0.6)
        W = (U + U.T)/2
        atoms.append((w, W/N[:, None]))
    if sum(A.sum() for _, A in atoms) == 0: continue
    a = acrit(atoms, N, nmax=120, tol=1e-6)
    f = funcs(atoms, N)
    for key, v in f.items():
        r = a*v/a0
        if key not in worst or r < worst[key][0]:
            worst[key] = (r, N, atoms)
for key, (r, N, atoms) in worst.items():
    print(key, 'min ratio', round(r, 5))
r, N, atoms = worst['lmax']
print('worst lmax case N=', N); [print(w, A.tolist()) for w, A in atoms]
