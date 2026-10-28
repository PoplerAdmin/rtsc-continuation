"""Wannier variant test for Al (v3.3 candidates), each in its own epw_dir; then G1 on random k."""
import json, os, shutil, subprocess, sys
from pathlib import Path
sys.path.insert(0, "/mnt/g/ClaudeWorkspace/RTSC_Project/bench_tools")
import eph_converge as ec
name = sys.argv[1]
over = json.loads(sys.argv[2])
base = Path("/home/user/rtsc_work/epw_Al")
m = json.load(open(base / "material.json"))
W = base / f"wvar_{name}"
(W / "out").mkdir(parents=True, exist_ok=True)
if not (W / "out" / "pwscf.save").exists():
    shutil.copytree(base / "out" / "pwscf.save", W / "out" / "pwscf.save")
if not (W / "save").exists():
    os.symlink(base / "save", W / "save")
m["epw_dir"] = str(W)
m["wannier_overrides"] = {k: (float(m["E_F_eV"]) + v if k == "scdm_mu" else v) for k, v in over.items()}
d = W / "wannier"; d.mkdir(exist_ok=True)
kw = ec.wannier_keys(m) | {"epwwrite": True, "epwread": False, "epbwrite": False}
(d / "epw.in").write_text(ec.epw_input(m, **kw))
try:
    ec.sh(f"{ec.QE_BIN}/epw.x -nk 4 -in epw.in", d, d / "epw.out", 4)
except RuntimeError as e:
    if not (W / "out" / "pwscf.epmatwp").exists():
        raise
    print("expected post-epwwrite abort (no fine grid given):", str(e)[:80])
wout = (d / "pwscf.wout").read_text()
print([l for l in wout.splitlines() if "Omega Total" in l or "convergence criteria" in l.lower()][-3:])
json.dump(m, open(W / "material.json", "w"))
env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
r = subprocess.run([sys.executable, "-m", "ephproto.bandcheck", "run", str(W / "material.json"),
                    "--wannier", str(d), "--epmatwp", str(W / "out" / "pwscf.epmatwp"),
                    "--label", name, "--np", "4", "--ef-grids", "16", "24", "--kset", "random"],
                   cwd="/mnt/g/ClaudeWorkspace/RTSC_Project/bench_tools", env=env,
                   capture_output=True, text=True)
print(r.stdout[-3000:]); print(r.stderr[-1500:])
