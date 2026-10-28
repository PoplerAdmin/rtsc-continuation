import sys, numpy as np
from ase.io import read
from pathlib import Path
from mattersim.forcefield import MatterSimCalculator
at = read(Path.home()/"rtsc_work/rtsc_bundle/rtsc/examples/discovery/prototypes/H3S_Im-3m.vasp")
at.rattle(0.03, seed=1)
out = {}
for dev in ("cpu","cuda"):
    a = at.copy(); a.calc = MatterSimCalculator(load_path="MatterSim-v1.0.0-5M", device=dev)
    out[dev] = (a.get_potential_energy(), a.get_forces(), a.get_stress())
print("dE", abs(out["cpu"][0]-out["cuda"][0]), "max dF", np.abs(out["cpu"][1]-out["cuda"][1]).max(), "max dS", np.abs(out["cpu"][2]-out["cuda"][2]).max())
