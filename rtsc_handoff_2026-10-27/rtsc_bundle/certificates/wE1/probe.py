import sys
from flint import arb
from cell import certify_cell
a=arb(sys.argv[1]); eta=arb(sys.argv[2]); N=int(sys.argv[3]); K=int(sys.argv[4])
wc=7.57486629486084; h=2.0**-12
w1=round(wc*2**12)/2**12
print(certify_cell(a,eta,w1,w1+h,N=N,K=K))
