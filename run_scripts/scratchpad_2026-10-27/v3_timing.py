"""A9: cost measurement on Al. One phself pass at K=16^3, Q=16^3 (shifted), sigma=20 meV."""
import json, sys, time
sys.path.insert(0, "/mnt/g/ClaudeWorkspace/RTSC_Project/bench_tools")
import eph_converge as ec
m = json.load(open("/home/user/rtsc_work/epw_Al/material.json"))
K, Q, s = [int(x) for x in sys.argv[1].split(",")], [int(x) for x in sys.argv[2].split(",")], float(sys.argv[3])
t0 = time.time()
d = ec.run_phself(m, K, Q, s, int(sys.argv[4]))
print("run dir", d, "wall", round(time.time() - t0, 1), "s")
sp = ec.load_phself(m, d)
print(f"atoms={len(sp.omega_K)} lambda={sp.lambda_():.5f} omega_log={sp.omega_log():.2f} K A={sp.A():.3f} K")
print({k: v for k, v in sp.provenance.items() if k in ("n_q", "nmodes", "fermi_energy_eV", "dos_ef", "merged_atoms", "counts")})
