# float illustration: Bergmann-Rainer functional derivative at the optimal Einstein spectrum.
import numpy as np
def M(a,w,N=400):
    n=np.arange(N); lv=lambda k: 2*a*w/(w**2+(2*np.pi*k)**2)
    K=lv(np.abs(np.subtract.outer(n,n)))+lv(np.add.outer(n,n)+1); np.fill_diagonal(K,lv(2*n+1))
    cum=np.concatenate([[0],np.cumsum(lv(np.arange(1,N)))]); return K-np.diag(2*n+1+2*cum)
a,w=4.33148148714,7.574866
ev,V=np.linalg.eigh(M(a,w)); x=V[:,-1]; print("rho",ev[-1])
# Q_x(a*, w') for w' on a grid:  (affine-in-measure evaluation of x at Einstein a* delta_w')
for wp in [2,4,6,7,7.4,7.5,7.574866,7.65,8,9,12,20,40]:
    print(wp, x@M(a,wp)@x)
