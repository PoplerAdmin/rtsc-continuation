import numpy as np
a=4.3315
def hs(w,K=200000):
    k=np.arange(1,K+1); lam=2*a*w/(w**2+(2*np.pi*k)**2)
    H=np.cumsum(1/(2*k-1.0))
    return np.sqrt(4*np.sum(lam**2*H/k))
for w in [10,20,40,60,80,100,150,200,300,500,1000,1e4]:
    print(w, hs(w))
