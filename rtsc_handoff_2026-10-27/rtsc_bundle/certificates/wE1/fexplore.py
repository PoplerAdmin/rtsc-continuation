import numpy as np
from scipy.optimize import brentq, minimize_scalar
def lam(k,a,w): return 2*a*w/(w**2+(2*np.pi*k)**2)
def M(a,w,N=200,s=1):
    n=np.arange(N); K=lam(np.abs(np.subtract.outer(n,n)),a,w)+s*lam(np.add.outer(n,n)+1,a,w)
    np.fill_diagonal(K,0)
    cum=np.concatenate([[0],np.cumsum(lam(np.arange(1,N),a,w))])
    d=s*lam(2*n+1,a,w)-(2*n+1)-2*cum
    return K+np.diag(d)
def rho(a,w,N=200,s=1): return np.linalg.eigvalsh(M(a,w,N,s))[-1]
if __name__=="__main__":
  a=4.33148149
  r=minimize_scalar(lambda w:-rho(a,w),bounds=(5,10),method='bounded',options={'xatol':1e-10})
  ws=r.x; print("w*",ws,"rho",-r.fun)
  for N in [50,100,150,200,400]: print(N, rho(a,ws,N))
  eps=1e-4; print("drho/da",(rho(a+eps,ws)-rho(a-eps,ws))/2/eps)
  h=1e-2; print("d2rho/dw2",(rho(a,ws+h)-2*rho(a,ws)+rho(a,ws-h))/h**2)
  v=np.linalg.eigh(M(a,ws,200))[1][:,-1]; print("evec", np.abs(v[[0,1,2,5,10,20,50,100,150]]))
  for w in [0.5,1,2,4,6,7.575,9,12,16,25,40,70,100,150,200,300,500,1000]:
      print(w, rho(a,w,300), rho(a,w,300,-1))
  