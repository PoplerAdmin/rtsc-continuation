"""Step 3b: MatterSim prescreen (HEURISTIC, S2) of the real MatterGen Mg-B structures, plus H3S check.
Writes ~/rtsc_work/mattergen_test/prescreen_<model>_<P>GPa.json and prints a table."""
import json, sys, time
from pathlib import Path
import numpy as np
from ase.io import read, write
from ase.build import bulk
from ase import Atoms
from rtsc.generate.mattergen import parse_outputs
from rtsc.discovery.interfaces import GeneratedStructure
from rtsc.mlff.mattersim import MatterSimPrescreener
from rtsc.mlff._atoms import atoms_to_structure
from pymatgen.core import Structure, Lattice
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

W = Path.home()/"rtsc_work/mattergen_test"
B = Path.home()/"rtsc_work/rtsc_bundle/rtsc/examples/discovery/prototypes"
# Reference phases for the ML hull: hcp Mg and the EXPERIMENTAL MgB2 structure (P6/mmm,
# a=3.086, c=3.524 A). Boron uses rtsc's recorded crude fallback (no alpha-B12 here) -> hull optimistic.
ref = W/"refs"; ref.mkdir(exist_ok=True)
write(ref/"Mg_hcp.cif", bulk("Mg", "hcp", a=3.209, c=5.211))
a, c = 3.086, 3.524
mgb2 = Atoms("MgB2", cell=[[a,0,0],[-a/2,a*3**0.5/2,0],[0,0,c]],
             scaled_positions=[[0,0,0],[1/3,2/3,0.5],[2/3,1/3,0.5]], pbc=True)
write(ref/"MgB2_P6mmm_exp.cif", mgb2)

def gs(struct, rid, label):
    return GeneratedStructure(structure=struct, generator="mattergen", request_id=rid,
                              provenance={"label": label, "flags": ["HYPOTHETICAL_STRUCTURE"]})
mg = []
for d, tag in [(W/"MgB_raw", "cli"), *[(p, "s7") for p in (W/"adapter_runs").glob("G-AMB*")]]:
    for i, s in enumerate(parse_outputs(d)[0]):
        mg.append(gs(s, tag, f"{tag}#{i}"))
known = [gs(atoms_to_structure(mgb2), "known", "MgB2-exp(KNOWN)")]
h3s = [gs(atoms_to_structure(read(B/"H3S_Im-3m.vasp")), "known", "H3S-Im-3m(KNOWN)")]

def sg(s):
    st = Structure(Lattice(s["lattice"]), s["species"], s["frac_coords"])
    an = SpacegroupAnalyzer(st, symprec=0.1)
    return st.composition.reduced_formula, f"{an.get_space_group_symbol()}", st
def run(model, P, items, phonons=True):
    cfg = {"model": model, "device": "cuda", "reference_dir": str(ref),
           "cache_dir": str(W/f"hull_cache"), "phonons": phonons}
    ps = MatterSimPrescreener(cfg)
    t0 = time.time(); res = ps.prescreen(items, P); dt = time.time()-t0
    rows = []
    print(f"\n=== {model} at {P} GPa: {len(items)} structures in {dt:.0f} s (HEURISTIC) ===")
    print(f"{'label':18} {'in':>8} {'sg_in':>10} {'relaxed':>10} {'Ehull eV/at':>11} {'min cm-1':>9} {'stable':>6} {'lightHF':>7}  defer?")
    for g, r in zip(items, res):
        f0, s0, _ = sg(g.structure)
        f1, s1, st = sg(r.relaxed_structure) if r.relaxed_structure else ("-", "-", None)
        dfr = "DEFER" if any(n.startswith("ML_SUGGEST_DEFER") for n in r.notes) else ""
        lhf = ps.heuristic_score(r)
        print(f"{g.provenance['label']:18} {f0:>8} {s0:>10} {s1:>10} "
              f"{(r.e_above_hull_ml_eV if r.e_above_hull_ml_eV is not None else float('nan')):11.3f} "
              f"{(r.min_freq_cm1 if r.min_freq_cm1 is not None else float('nan')):9.1f} {str(r.dynamically_stable_ml):>6} "
              f"{(lhf if lhf is not None else float('nan')):7.2f}  {dfr}{'' if r.ok else ' ML-FAIL'}")
        extra = {}
        if st is not None and g.provenance["label"].startswith(("MgB2-exp", "H3S")):
            extra = {"a": st.lattice.a, "c": st.lattice.c, "volume_per_atom": st.volume/len(st)}
            print(f"{'':18} relaxed a={st.lattice.a:.4f} c={st.lattice.c:.4f} V/atom={st.volume/len(st):.3f}")
        rows.append({"label": g.provenance["label"], "formula_in": f0, "sg_in": s0, "sg_relaxed": s1,
                     "e_above_hull_ml_eV": r.e_above_hull_ml_eV, "min_freq_cm1": r.min_freq_cm1,
                     "dynamically_stable_ml": r.dynamically_stable_ml, "light_hf_fraction": lhf,
                     "ok": r.ok, "notes": r.notes, **extra})
    out = W/f"prescreen_{'h3s' if which == 'h3s' else 'mgb'}_{model}_{P:g}GPa.json"
    out.write_text(json.dumps({"HEURISTIC": "ARCH S2: ordering / DEFERRED only", "provenance": ps.provenance(),
                               "wall_s": dt, "rows": rows}, indent=1, default=str))
    return rows

which = sys.argv[1]
if which == "h3s":
    for model in ("MatterSim-v1.0.0-1M", "MatterSim-v1.0.0-5M"):
        for P in (150.0, 200.0):
            run(model, P, h3s)
else:
    for P in (0.0, 150.0):
        run("MatterSim-v1.0.0-5M", P, known + mg)
