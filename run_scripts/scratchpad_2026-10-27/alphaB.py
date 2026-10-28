"""SIDE CHECK (HEURISTIC, S2): quantify optimism of rtsc's fcc boron fallback using alpha-B12 relaxed with 5M.
alpha-B12 R-3m (hex a=4.908 c=12.567; B1 18h x=0.1187 z=0.8911; B2 18h x=0.1965 z=0.0243) - literature geometry from memory, relaxed anyway."""
import json, re
from pathlib import Path
from pymatgen.core import Structure, Lattice, Composition
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDEntry
from mattersim.forcefield import MatterSimCalculator
from rtsc.mlff.mattersim import relax_at_pressure, RelaxSettings
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
W = Path.home()/"rtsc_work/mattergen_test"
s = Structure.from_spacegroup("R-3m", Lattice.hexagonal(4.908, 12.567), ["B","B"], [[0.1187,-0.1187,0.8911],[0.1965,-0.1965,0.0243]])
s = SpacegroupAnalyzer(s).get_primitive_standard_structure()
print("alpha-B atoms", len(s), "min dist", round(min(s.distance_matrix[s.distance_matrix>0]),3))
calc = MatterSimCalculator(load_path="MatterSim-v1.0.0-5M", device="cuda")
out = {}
for P in (0.0, 150.0):
    ra, h, info = relax_at_pressure(AseAtomsAdaptor.get_atoms(s), calc, P, RelaxSettings())
    sg = SpacegroupAnalyzer(AseAtomsAdaptor.get_structure(ra), 0.1).get_space_group_symbol()
    fcc = json.loads((W/f"hull_cache/B_P{P:g}GPa_MatterSim-v1.0.0-5M.json").read_text())["entries"][0]["enthalpy_eV_per_atom"]
    mg = json.loads((W/f"hull_cache/Mg_P{P:g}GPa_MatterSim-v1.0.0-5M.json").read_text())["entries"][0]["enthalpy_eV_per_atom"]
    mgb2 = json.loads((W/f"hull_cache/B-Mg_P{P:g}GPa_MatterSim-v1.0.0-5M.json").read_text())["entries"][0]["enthalpy_eV_per_atom"]
    print(f"P={P}: alpha-B H={h:.4f} eV/at (sg {sg}, conv {info['converged']}, V/at {info['volume_A3_per_atom']:.3f}); fcc-B fallback H={fcc:.4f}; diff {fcc-h:.3f} eV/at")
    refs = [PDEntry(Composition("Mg"), mg), PDEntry(Composition("MgB2"), 3*mgb2), PDEntry(Composition("B"), h)]
    pd = PhaseDiagram(refs)
    d = json.loads((W/f"prescreen_MatterSim-v1.0.0-5M_{P:g}GPa.json").read_text())
    for r in d["rows"]:
        rel = next(n for n in r["notes"] if n.startswith("relax:"))
        e = float(re.search(r"energy_eV_per_atom=([-\d.e+]+)", rel).group(1)); v = float(re.search(r"volume_A3_per_atom=([-\d.e+]+)", rel).group(1))
        H = e + P*v/160.21766208
        c = Composition(r["formula_in"]); 
        eh = pd.get_e_above_hull(PDEntry(c, H*c.num_atoms), allow_negative=True)
        xb = c.get_atomic_fraction("B")
        if xb > 2/3 + 1e-6:
            print(f"   {r['label']:10} {r['formula_in']:8} xB={xb:.2f} E_hull(rtsc,fccB)={r['e_above_hull_ml_eV']:.3f}  E_hull(alpha-B ref, side check)={eh:.3f}")
