import sys, time, math
sys.path.insert(0, "/mnt/g/ClaudeWorkspace/RTSC_Project/bench_tools")
import analyze_benchmark as ab
from rtsc.spectrum import Spectrum
from pathlib import Path
# simple case: SYNTHETIC Einstein-like (narrow gaussian) spectrum, lambda=1, w=300K
import numpy as np
om = np.linspace(250, 350, 41); g = np.exp(-((om-300)/15)**2); 
lam_target=1.0; w = g/ (2*(g/om).sum()) * lam_target
sp = Spectrum(tuple(om), tuple(w), 0.0, "SYNTHETIC einstein")
for mu in (0.0, 0.10, 0.15):
    t=time.time(); el = ab.eliashberg_tc(sp, mu, 10*om.max())
    ad = ab.allen_dynes(sp.lambda_(), sp.omega_log(), sp.omega2(), mu)
    rf = ab.rtsc_float_tc(sp, mu, 10*om.max()) 
    print("SYNTH mu",mu,"El",el["Tc_K"],"AD",ad["Tc_K"],"rtsc",rf, "t",time.time()-t)
sps = ab.load_qe_dir(Path("/mnt/g/ClaudeWorkspace/RTSC_Project/rtsc_handoff_2026-10-27/rtsc_bundle/rtsc/tests/data_dfpt/qe67_al_smoke"))
s=sps[0]; print(s.label, s.lambda_(), s.omega_log(), s.A(), s.omega_max(), s.provenance.get("smearing_index"), s.rel_envelope)
for mu in (0.0, 0.10, 0.15):
    t=time.time(); el = ab.eliashberg_tc(s, mu, 10*s.omega_max())
    ad = ab.allen_dynes(s.lambda_(), s.omega_log(), s.omega2(), mu)
    print("SMOKE1 mu",mu,"El",el,"AD",ad["Tc_K"],"t",time.time()-t)
