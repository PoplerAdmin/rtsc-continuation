import numpy as np
from fexp import lamv, M, S
def Ghead(a,w,N,L):
    i=np.arange(N)[:,None]; m=np.arange(N,N+L)[None,:]
    C=lamv(m-i,a,w)+lamv(m+i+1,a,w)
    beta=2*np.arange(N,N+L)+1-2*S(a,w)
    return (C/beta)@C.T
a=4.33148; w=7.5748661
for N in [40,60,80,100]:
    B=M(a,w,N); G=Ghead(a,w,N,300)
    print(N, np.linalg.eigvalsh(B)[-1], np.linalg.eigvalsh(B+G)[-1])
