import numpy as np
a0 = 4.33148145; TP = 2*np.pi
def worst(w, mode, nmax=4000):
    k = np.arange(1, nmax+1); l = 2*w/(w*w+(TP*k)**2)
    S = np.concatenate([[0], np.cumsum(l)])[:nmax]   # S(n), n=0..nmax-1
    n = np.arange(nmax)
    base = 2/w if mode == 'I' else 2/w - 1/np.sinh(w/2)
    return np.max(a0*(base + 2*S) - (2*n+1))   # need <= 0
for w in [0.5, 1, 2, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 8, 8.66, 9, 10, 12, 15, 20, 50, 200]:
    print(w, 'I: %+.4f' % worst(w, 'I'), ' II: %+.4f' % worst(w, 'II'))
