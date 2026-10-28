import sys, copy, math
sys.path.insert(0,"/home/claude/rtsc"); sys.path.insert(0,".")
import numpy as np
from cell_verify import float_margin, push
from rtsc.spectrum import Spectrum
from rtsc.certify import CertifyParams, classify
from rtsc.certify.arbcore import ArbSpectrum
from rtsc.certify.cells import check_cell
rng=np.random.default_rng(1)
sp=Spectrum((300.0, 1500.0),(500.0,700.0),0.2,"fz")
v=classify([sp],CertifyParams(mode="bracket",mu_hi=0.0,time_budget_s=30))
asp=ArbSpectrum(sp.omega_K,sp.weight_K,sp.rel_envelope)
cells=v.certificate["alternatives"][0]["cells"]["cells"]
print(v.status.value, v.T_up_K, len(cells), flush=True)
bad=0; npass=0
for c in cells[:2]:
    # push C up and down
    for key in ("C",):
        f,cp=push(asp,c,key,1e3); m=float_margin(sp.omega_K,sp.weight_K,0.2,cp,n_rows=10*int(c["N_e"]))
        print(" push C x",f," margin",m, flush=True)
        cc=copy.deepcopy(c); 
        lo,hi=1e-3,1.0
        for _ in range(30):
            mid=math.sqrt(lo*hi); cc["C"]=c["C"]*mid
            if check_cell(asp,cc)["ok"]: hi=mid
            else: lo=mid
        cc["C"]=c["C"]*hi; m=float_margin(sp.omega_K,sp.weight_K,0.2,cc,n_rows=10*int(c["N_e"]))
        print(" shrink C x",hi," ok",check_cell(asp,cc)["ok"]," margin",m, flush=True)
    for trial in range(25):
        cc=copy.deepcopy(c)
        cc["u"]=[x*math.exp(rng.normal(0,0.02)) for x in c["u"]]
        if rng.random()<0.5: cc["M_v"]=int(c["N_e"])+2+int(rng.integers(0,50))
        if rng.random()<0.3: cc["p"]=2; cc["C"]=c["C"]*float(rng.uniform(0.2,5))
        cc["eta"]=c["eta"]*float(rng.uniform(0.5,3))
        if check_cell(asp,cc)["ok"]:
            npass+=1
            m=float_margin(sp.omega_K,sp.weight_K,0.2,cc,n_rows=6*int(cc["N_e"]))
            if m>1e-9: bad+=1; print("  BUG? passing cert with float margin",m,cc["p"],cc["M_v"],flush=True)
print("fuzz: passed",npass,"suspicious",bad)
