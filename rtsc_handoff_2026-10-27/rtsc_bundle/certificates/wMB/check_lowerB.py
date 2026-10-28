# float sanity: (i) Poisson formula for the symbol; (ii) B_w - sigma_pi I - gamma e0 e0^T >= 0 (even sector)
import numpy as np
TP = 2*np.pi
for w in [0.3, 2.0, 5.0, 7.57, 20.0]:
    k = np.arange(-200000, 200001); l = 2*w/(w*w+(TP*k)**2)
    th = np.array([0, 1.0, np.pi]); s_num = [(l*np.cos(k*t)).sum() for t in th]
    s_formula = np.cosh((w/TP)*(np.pi-th))/np.sinh(w/2)
    n = np.arange(400); g = 2*w/(w*w+(TP*np.arange(1000))**2)
    B = g[np.abs(n[:, None]-n[None, :])] + g[n[:, None]+n[None, :]+1]
    gam = (w/TP)**2/np.sinh(w/2)
    M = B - np.eye(400)/np.sinh(w/2); M[0, 0] -= gam
    print('w=%5.2f symbol err %.1e  min eig(B - s_pi - gam e0e0) = %.3e' % (w, np.max(np.abs(np.array(s_num)-s_formula)), np.linalg.eigvalsh(M)[0]))
