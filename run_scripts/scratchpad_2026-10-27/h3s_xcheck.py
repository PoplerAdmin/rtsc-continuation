"""Cross-check (HEURISTIC, S2) of the H3S ML phonon minimum vs supercell choice. Uses rtsc library functions unchanged."""
import json, numpy as np
from pathlib import Path
from ase.io import read
from mattersim.forcefield import MatterSimCalculator
from rtsc.mlff.mattersim import relax_at_pressure, RelaxSettings
from rtsc.mlff.phonons import ml_phonons, PhononSettings
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
conv = read(Path.home()/"rtsc_work/rtsc_bundle/rtsc/examples/discovery/prototypes/H3S_Im-3m.vasp")
out = []
def special(atoms, calc, sc, prim):
    """frequencies at bcc H and P and N via phonopy with explicit primitive matrix"""
    u = PhonopyAtoms(symbols=atoms.get_chemical_symbols(), cell=np.array(atoms.cell[:]), scaled_positions=atoms.get_scaled_positions())
    ph = Phonopy(u, supercell_matrix=sc, primitive_matrix=prim)
    ph.generate_displacements(distance=0.01)
    from ase import Atoms
    F=[]
    for s in ph.supercells_with_displacements:
        a = Atoms(symbols=s.symbols, cell=s.cell, scaled_positions=s.scaled_positions, pbc=True); a.calc = calc
        f = a.get_forces(); F.append(f - f.mean(0))
    ph.forces = F; ph.produce_force_constants(); ph.symmetrize_force_constants()
    # bcc primitive reduced coords: H=(1/2,-1/2,1/2), N=(0,0,1/2), P=(1/4,1/4,1/4)
    res = {}
    for name,q in [("H",[0.5,-0.5,0.5]),("N",[0,0,0.5]),("P",[0.25,0.25,0.25])]:
        res[name] = float(min(ph.get_frequencies(q))*33.35640951981521)
    return res, len(ph.supercell)
prim = [[-0.5,0.5,0.5],[0.5,-0.5,0.5],[0.5,0.5,-0.5]]
for model in ("MatterSim-v1.0.0-1M","MatterSim-v1.0.0-5M"):
    calc = MatterSimCalculator(load_path=model, device="cuda")
    for P in (150.0, 200.0):
        ra, h, info = relax_at_pressure(conv, calc, P, RelaxSettings())
        row = {"model": model, "P_GPa": P, "a_conv": float(ra.cell.lengths()[0])}
        for sc in ([2,2,2],[3,3,3],[4,4,4]):
            ps = ml_phonons(ra, calc, PhononSettings(supercell=sc, projected_dos=False))
            row[f"mesh_min_{sc[0]}x{sc[0]}x{sc[0]}"] = (round(ps.min_freq_cm1,1), ps.min_freq_qpoint)
        for sc in ([2,2,2],[3,3,3]):
            r, n = special(ra, calc, np.diag(sc), prim)
            row[f"bccHNP_conv{sc[0]}x{sc[0]}x{sc[0]}({n}at)"] = {k: round(v,1) for k,v in r.items()}
        print(json.dumps(row), flush=True); out.append(row)
Path(Path.home()/"rtsc_work/mattergen_test/h3s_supercell_xcheck.json").write_text(json.dumps({"HEURISTIC":"S2 cross-check, ML only","rows":out}, indent=1))
