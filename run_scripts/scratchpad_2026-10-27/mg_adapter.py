"""Step 3: real MatterGen output through the rtsc adapter, then hand inspection (pymatgen/spglib)."""
import json, sys, zipfile, io
from pathlib import Path
from rtsc.generate.mattergen import parse_outputs, MatterGenGenerator
from rtsc.discovery.interfaces import GenerationRequest
from pymatgen.core import Structure, Lattice
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

raw = Path.home()/"rtsc_work/mattergen_test/MgB_raw"
# (a) parse the CLI's own output: extxyz (preferred) and the CIF zip fallback, and compare
s_xyz, src = parse_outputs(raw)
tmp = Path.home()/"rtsc_work/mattergen_test/cifonly"; tmp.mkdir(exist_ok=True)
(tmp/"generated_crystals_cif.zip").write_bytes((raw/"generated_crystals_cif.zip").read_bytes())
s_cif, src2 = parse_outputs(tmp)
print(f"(a) parsed {len(s_xyz)} from {src}, {len(s_cif)} from {src2}")
def S(d): return Structure(Lattice(d["lattice"]), d["species"], d["frac_coords"])
mism = 0
for i,(a,b) in enumerate(zip(s_xyz, s_cif)):
    A,B = S(a),S(b)
    if A.composition != B.composition or abs(A.volume-B.volume) > 1e-3*A.volume: mism += 1
print(f"    extxyz vs CIF zip: {mism} mismatches (composition/volume)")

# (b) the adapter end to end, calling the real executable as a subprocess (seed 7, new outdir)
gen = MatterGenGenerator({"pretrained_name": "chemical_system_energy_above_hull",
    "executable": [str(Path.home()/"mattergen/.venv/bin/mattergen-generate")],
    "workdir": str(Path.home()/"rtsc_work/mattergen_test/adapter_runs"), "batch_size": 16,
    "timeout_s": 3600})
req = GenerationRequest(chemical_system=("B","Mg"), space_group=None, pressure_GPa=0.0, goal_id="G-AMB",
                        n=16, seed=7, extra={"mattergen_conditioning": {"energy_above_hull": 0.05}})
print("(b) plan:", json.dumps({k: (str(v) if k=="outdir" else v) for k,v in gen.plan(req).items()}, default=str)[:600])
out = gen.generate(req)
print(f"    adapter returned {len(out)} GeneratedStructure; flags of #0: {out[0].provenance.get('flags')}")
print("    provenance keys:", sorted(out[0].provenance))

# (c) hand inspection of both sets
def inspect(label, structs):
    print(f"\n(c) {label}")
    print(f"{'#':>2} {'formula':10} {'sg (0.1)':14} {'sg (0.01)':14} {'nat':>3} {'rho g/cc':>8} {'dmin A':>6} {'dMgB':>5} {'dBB':>5}")
    for i,d in enumerate(structs):
        st = S(d)
        sg1 = SpacegroupAnalyzer(st, symprec=0.1).get_space_group_symbol()
        n1 = SpacegroupAnalyzer(st, symprec=0.1).get_space_group_number()
        sg2 = SpacegroupAnalyzer(st, symprec=0.01).get_space_group_symbol()
        dm = st.distance_matrix + 1e9*__import__("numpy").eye(len(st))
        sp = [x.symbol for x in st.species]
        def pair(a,b):
            v=[dm[i][j] for i in range(len(st)) for j in range(len(st)) if {sp[i],sp[j]}=={a,b} and (a!=b or i<j)]
            return min(v) if v else float('nan')
        print(f"{i:2d} {st.composition.reduced_formula:10} {sg1+' ('+str(n1)+')':14} {sg2:14} {len(st):3d} {st.density:8.3f} {dm.min():6.3f} {pair('Mg','B'):5.2f} {pair('B','B'):5.2f}")
inspect("CLI run (seed default)", s_xyz)
inspect("adapter run (seed 7)", [o.structure for o in out])
