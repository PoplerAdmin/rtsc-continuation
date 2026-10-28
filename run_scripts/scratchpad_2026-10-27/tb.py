import traceback, torch
from ase.build import bulk
from mattersim.forcefield import MatterSimCalculator
at = bulk("Mg","hcp",a=3.209,c=5.211)
at.calc = MatterSimCalculator(load_path="MatterSim-v1.0.0-1M", device="cuda")
try:
    print(at.get_potential_energy(), at.get_stress())
except Exception:
    traceback.print_exc()
