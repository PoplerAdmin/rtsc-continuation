from pathlib import Path
from rtsc.generate.mattergen import parse_outputs
from pymatgen.core import Structure, Lattice
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
W = Path.home()/"rtsc_work/mattergen_test"
for d, tag in [(W/"MgB_raw", "cli"), *[(p, "s7") for p in (W/"adapter_runs").glob("G-AMB*")]]:
    for i, s in enumerate(parse_outputs(d)[0]):
        st = Structure(Lattice(s["lattice"]), s["species"], s["frac_coords"])
        print(f"{tag}#{i}\t{len(st)}\t{st.composition.formula.replace(' ','')}\t{SpacegroupAnalyzer(st,0.1).get_space_group_symbol()}")
