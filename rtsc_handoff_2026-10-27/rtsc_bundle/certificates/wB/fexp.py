# float exploration (guidance only)
import numpy as np
from scipy.optimize import brentq
def lamv(k,a,w): return 2*a*w/(w**2+(2*np.pi*k)**2)
def M(a,w,N):
    n=np.arange(N)
    K=lamv(np.abs(np.subtract.outer(n,n)),a,w)+lamv(np.add.outer(n,n)+1,a,w)
    np.fill_diagonal(K, lamv(2*n+1,a,w))
    cum=np.concatenate([[0],np.cumsum(lamv(np.arange(1,N),a,w))])
    return K-np.diag(2*n+1+2*cum)
def S(a,w): return a*(0.5/np.tanh(w/2)-1/w)
def gdiag(a,w,N,extra=4000):
    # g_i = sum_{m>=N} K_im * (2S)/beta_m
    i=np.arange(N)[:,None]; m=np.arange(N,N+extra)[None,:]
    K=lamv(m-i,a,w)+lamv(m+i+1,a,w)
    beta=2*m+1-2*S(a,w)
    return (K*2*S(a,w)/beta).sum(1)
def rho(a,w,N): return np.linalg.eigvalsh(M(a,w,N))[-1]
a0=4.3314814875; w0=7.5748661
for N in [30,50,80,120,200]:
    r=rho(a0,w0,N); g=gdiag(a0,w0,N); rg=np.linalg.eigvalsh(M(a0,w0,N)+np.diag(g))[-1]
    print(N, r, rg, g[-1])
eps=1e-5
print("drho/da", (rho(a0+eps,w0,300)-rho(a0-eps,w0,300))/(2*eps))
for w in [2,4,6,7,7.5,8,9,12,20,40,80]:
    print(w, rho(4.33145,w,200))
