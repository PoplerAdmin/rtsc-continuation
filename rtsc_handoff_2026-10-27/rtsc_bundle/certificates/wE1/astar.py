# float estimate of a* (not part of the certificate)
import numpy as np
from scipy.optimize import brentq, minimize_scalar
from fexplore import M
def rho(a,w,N): return np.linalg.eigvalsh(M(a,w,N))[-1]
for N in [200,400,800]:
    ac=lambda w: brentq(lambda a: rho(a,w,N),4.3,4.4,xtol=1e-14)
    r=minimize_scalar(ac,bounds=(7.4,7.7),method='bounded',options={'xatol':1e-7})
    print(N, repr(r.fun), r.x)
