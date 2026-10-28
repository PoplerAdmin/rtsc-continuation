import sys, copy, math
sys.path.insert(0,"/home/claude/rtsc"); sys.path.insert(0,".")
import numpy as np, mpmath as mp, indep
from cell_verify import float_margin, push
from rtsc.spectrum import Spectrum
from rtsc.certify import CertifyParams, classify
from rtsc.certify.arbcore import ArbSpectrum
from rtsc.certify.cells import check_cell
from rtsc.certify.rayleigh import check_rayleigh, certified_n_c
sp=Spectrum((300.0, 1500.0),(500.0,700.0),0.2,"fz")
asp=ArbSpectrum(sp.omega_K,sp.weight_K,sp.rel_envelope)
v=classify([sp],CertifyParams(mode="bracket",mu_hi=0.1,omega_c=3000.0,time_budget_s=30))
c=v.certificate["alternatives"][0]["cells"]["cells"][0]
# p=2 variant
c2=copy.deepcopy(c); c2["p"]=2; N=int(c["N"]); c2["C"]=c["C"]/(2*N+1)
for Ne_mul in (1,2,4):
    cc=copy.deepcopy(c2); cc["N_e"]=int(c["N_e"])*Ne_mul; cc["M_v"]=max(int(c["M_v"]),cc["N_e"]+2)
    f,cp=push(asp,cc,"eta",1e3) if check_cell(asp,cc)["ok"] else (0,None)
    print("p=2 N_e x",Ne_mul,"ok" if cp else "fail", check_cell(asp,cc).get("reason"), flush=True)
    if cp: print("   pushed eta margin", float_margin(sp.omega_K,sp.weight_K,0.2,cp,n_rows=8*cc["N_e"]), flush=True)
# Rayleigh T-push
rc=v.certificate["alternatives"][0]["rayleigh"]; T0=float.fromhex(rc["T"])
lo,hi=T0,T0*1.2
def ok(T):
    r=dict(rc); r["T"]=T.hex(); n=certified_n_c(T,rc["omega_c"]); 
    if n is None: return False
    r["N_c"]=n; return check_rayleigh(asp,r)["ok"]
for _ in range(40):
    mid=0.5*(lo+hi)
    if ok(mid): lo=mid
    else: hi=mid
T=lo; x=rc["x"]; Nc=indep.nc_count(T,rc["omega_c"]); e=0.2; p=np.array(sp.weight_K)
mp.mp.dps=40
if len(x)<=200:
    Fs=[indep.mp_atomF(o,T,x) for o in sp.omega_K]
    diag=mp.fsum((2*n+1)*mp.mpf(x[n])**2 for n in range(len(x))); sx=mp.fsum(mp.mpf(t) for t in x[:min(Nc,len(x))])
    bm=mp.fsum(min((1-mp.mpf(e))*mp.mpf(pj)*F,(1+mp.mpf(e))*mp.mpf(pj)*F) for pj,F in zip(p,Fs))-diag-2*mp.mpf(rc["mu_hi"])*sx**2
    print("Rayleigh pushed T",T0,"->",T,"Nx",len(x),"mp boxminQ",mp.nstr(bm,8), "Q/|x|^2", mp.nstr(bm/mp.fsum(mp.mpf(t)**2 for t in x),5))
else:
    print("Nx large",len(x))
