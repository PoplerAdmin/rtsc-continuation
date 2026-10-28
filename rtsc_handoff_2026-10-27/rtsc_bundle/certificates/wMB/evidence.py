# Float evidence: ratios a* * Phi / a0 (>= 1 means "Tc <= Phi/a0" holds for that spectrum).
import numpy as np, sys, json
from mb import acrit
a0 = 4.33148145
def functionals(atoms, N):
    N = np.asarray(N, float); s = np.sqrt(np.outer(N, 1/N))
    A = sum(np.asarray(M, float) for _, M in atoms)
    K = s*A; R = np.diag(A.sum(1))
    return dict(avg=(N[:, None]*A).sum()/N.sum(), maxdiag=np.max(np.diag(A)),
                kappa=np.linalg.eigvalsh(K)[-1], thmA=np.linalg.eigvalsh((R+K)/2)[-1],
                rowmax=A.sum(1).max(), trace_half=np.trace((R+K)/2))
def ratios(atoms, N, nmax=110):
    a = acrit(atoms, N, nmax=nmax, tol=1e-8)
    return {k: a*v/a0 for k, v in functionals(atoms, N).items()}
def rand_case(rng, B, nat):
    N = rng.dirichlet(np.ones(B)*0.5) + 1e-4; N /= N.sum()
    atoms = []
    for _ in range(nat):
        w = np.exp(rng.uniform(np.log(2), np.log(40)))
        U = rng.exponential(1, (B, B))*(rng.random((B, B)) < 0.7)
        W = (U+U.T)/2
        if W.sum() == 0: W[0, 0] = 1
        atoms.append((w, W/N[:, None]))
    return atoms, N
if __name__ == "__main__":
    seed, B, nat, trials = map(int, sys.argv[1:5])
    rng = np.random.default_rng(seed)
    mins = {}
    for t in range(trials):
        atoms, N = rand_case(rng, B, nat)
        r = ratios(atoms, N)
        for k, v in r.items():
            mins[k] = min(mins.get(k, np.inf), v)
    print(json.dumps(dict(seed=seed, B=B, nat=nat, trials=trials, min_ratio={k: round(v, 5) for k, v in mins.items()})))
