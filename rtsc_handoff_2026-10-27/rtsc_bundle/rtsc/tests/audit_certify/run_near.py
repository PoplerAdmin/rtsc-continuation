import sys, time, math
sys.path.insert(0, "/home/claude/rtsc"); sys.path.insert(0, ".")
import numpy as np
from run_adversarial import einstein_A_for_tc
from rtsc.spectrum import Spectrum
from rtsc.certify import classify, CertifyParams, recheck
import indep
for Om in (1000.0, 2272.0, 5000.0):
    A_c = einstein_A_for_tc(Om, 300.0, N=400)
    for rel in (-1e-8, -1e-9, 1e-9, 1e-8):
        A = A_c*(1+rel)
        s = Spectrum((Om,), (A,), 0.0)
        for mode in ("decide",):
            v = classify([s], CertifyParams(mode=mode, mu_hi=0.0, time_budget_s=60))
            r300 = indep.rho_N([Om],[A],300.0,400)
            print(Om, rel, mode, v.status.value, v.T_low_K, v.T_up_K, "rho_N(300)=%.3e"%r300, recheck(v,[s])["ok"], flush=True)
