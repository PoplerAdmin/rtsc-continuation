"""ORCH: code paths that use optional packages (skipped when they are not installed)."""
from __future__ import annotations

import pytest

from rtsc.sources import Candidate, dedupe, load_local
from rtsc.stability import HullEntry, e_above_hull

POSCAR = """bcc test
1.0
3.0 0 0
0 3.0 0
0 0 3.0
Nb H
2 2
Direct
0 0 0
0.5 0.5 0.5
0.5 0 0
0 0.5 0.5
"""


def test_pymatgen_phase_diagram_matches_lp():
    pytest.importorskip("pymatgen")
    ents = [HullEntry({"A": 1}, -1.0), HullEntry({"B": 1}, -2.0),
            HullEntry({"A": 1, "B": 1}, -2.5), HullEntry({"A": 1, "B": 2}, -2.4)]
    # use real element symbols for pymatgen
    ren = {"A": "Nb", "B": "H"}
    ents_r = [HullEntry({ren[k]: v for k, v in e.composition.items()}, e.energy_per_atom) for e in ents]
    for comp, e in (({"Nb": 1, "H": 3}, -1.9), ({"Nb": 1, "H": 2}, -2.4), ({"Nb": 2, "H": 1}, -1.5)):
        v1, how1 = e_above_hull(comp, e, ents_r, backend="pymatgen")
        v2, how2 = e_above_hull(comp, e, ents_r, backend="scipy")
        assert "pymatgen" in how1 and "scipy" in how2
        assert v1 == pytest.approx(v2, abs=1e-9)


def test_structure_matcher_dedupe():
    pytest.importorskip("pymatgen")
    lat = [[3.0, 0, 0], [0, 3.0, 0], [0, 0, 3.0]]
    a = Candidate("a", {}, lattice=lat, species=["Nb", "H"], frac_coords=[[0, 0, 0], [.5, .5, .5]],
                  enthalpy_eV_per_atom=-1.0)
    b = Candidate("b", {}, lattice=[[3.03, 0, 0], [0, 3.03, 0], [0, 0, 3.03]], species=["H", "Nb"],
                  frac_coords=[[.1, .1, .1], [.6, .6, .6]], enthalpy_eV_per_atom=-1.2)
    c = Candidate("c", {}, lattice=lat, species=["Nb", "H"], frac_coords=[[0, 0, 0], [.5, 0, 0]])
    for x in (a, b, c):
        x.composition = {"Nb": 1.0, "H": 1.0}
    res = dedupe([a, b, c], method="structure_matcher")
    assert res.duplicate_of == {"a": "b"}
    assert res.method_used["b"] == "structure_matcher"


@pytest.mark.parametrize("backend", ["pymatgen", "ase"])
def test_local_backends(tmp_path, backend):
    pytest.importorskip(backend)
    (tmp_path / "x.vasp").write_text(POSCAR)
    (c,) = load_local(tmp_path, backend=backend)
    assert c.composition == {"Nb": 2.0, "H": 2.0}
    assert c.provenance["reader"] == backend


def test_spglib_spacegroup_fallback():
    pytest.importorskip("spglib")
    from rtsc.sources.dedupe import spacegroup_of
    c = Candidate("x", {}, lattice=[[3, 0, 0], [0, 3, 0], [0, 0, 3]], species=["Nb", "Nb"],
                  frac_coords=[[0, 0, 0], [.5, .5, .5]])
    assert spacegroup_of(c) == "Im-3m"
