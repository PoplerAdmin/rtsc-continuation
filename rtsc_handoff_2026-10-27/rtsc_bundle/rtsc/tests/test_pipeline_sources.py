"""ORCH: sources, readers, dedupe (tiny synthetic inputs; no network, no pymatgen needed)."""
from __future__ import annotations

import json

import pytest

from rtsc.sources import (Candidate, dedupe, load_airss, load_calypso, load_local,
                          load_sources, load_uspex, query_materials_project, substitute)
from rtsc.sources.composition import formula_string, parse_formula
from rtsc.sources.fileformats import read_cif, read_extxyz, read_poscar, read_res

POSCAR = """H3S test cell
1.0
3.0 0.0 0.0
0.0 3.0 0.0
0.0 0.0 3.0
S H
1 3
Direct
0.0 0.0 0.0
0.5 0.0 0.0
0.0 0.5 0.0
0.0 0.0 0.5
"""

CIF = """data_test
_cell_length_a 4.0
_cell_length_b 4.0
_cell_length_c 4.0
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_symmetry_space_group_name_H-M 'I m -3 m'
loop_
_symmetry_equiv_pos_as_xyz
x,y,z
x+1/2,y+1/2,z+1/2
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
S1 S 0 0 0 1.0
"""

RES = """TITL H3S-1234-5-6 200.00 40.0 -123.456 0 0 4 (Im-3m) n - 2
CELL 1.54180 3.0 3.0 3.0 90 90 90
LATT -1
SFAC S H
S 1 0.0 0.0 0.0 1.0
H 2 0.5 0.0 0.0 1.0
H 2 0.0 0.5 0.0 1.0
H 2 0.0 0.0 0.5 1.0
END
"""


def test_formula_parsing():
    assert parse_formula("LaH10") == {"La": 1, "H": 10}
    assert parse_formula("Ca(BH4)2") == {"Ca": 1, "B": 2, "H": 8}
    assert parse_formula("Mg0.5B") == {"Mg": 0.5, "B": 1}
    assert formula_string({"H": 6, "S": 2}) == "SH3"


def test_poscar_and_cartesian():
    raw = read_poscar(POSCAR)
    assert raw.species == ["S", "H", "H", "H"]
    assert raw.frac_coords[1] == [0.5, 0.0, 0.0]
    cart = ("t\n2.0\n1.5 0 0\n0 1.5 0\n0 0 1.5\nS H\n1 1\nCartesian\n0 0 0\n0.75 0 0\n")
    raw2 = read_poscar(cart)                 # scale applies to lattice and cartesian coords
    assert raw2.lattice[0][0] == pytest.approx(3.0)
    assert raw2.frac_coords[1][0] == pytest.approx(0.5)


def test_cif_symmetry_expansion():
    (raw,) = read_cif(CIF)
    assert raw.species == ["S", "S"]           # body-centring doubles the site
    assert raw.info["spacegroup"] == "I m -3 m"


def test_extxyz_and_res():
    text = ('2\nLattice="3 0 0 0 3 0 0 0 3" Properties=species:S:1:pos:R:3 energy=-10.0 '
            'pressure=150\nS 0 0 0\nH 1.5 0 0\n')
    (fr,) = read_extxyz(text)
    assert fr.frac_coords[1][0] == pytest.approx(0.5)
    assert fr.info["energy"] == -10.0
    raw = read_res(RES)
    assert raw.info["pressure_GPa"] == 200.0 and raw.info["enthalpy_eV"] == -123.456
    assert raw.info["spacegroup"] == "Im-3m"


def test_local_dir_with_manifest_and_bad_file(tmp_path):
    (tmp_path / "a.vasp").write_text(POSCAR)
    (tmp_path / "b.cif").write_text(CIF)
    (tmp_path / "broken.cif").write_text("garbage")
    (tmp_path / "a2f.dat").write_text("1 0\n2 1\n3 0\n")
    (tmp_path / "manifest.toml").write_text(
        '[a]\npressure_GPa = 150.0\nproperties = { is_metal = true }\n'
        'spectra_files = [ { path = "a2f.dat", units = "meV" } ]\n')
    cands = load_local(tmp_path, backend="builtin")
    ids = sorted(c.material_id for c in cands)
    assert ids == ["a", "b", "broken"]
    a = next(c for c in cands if c.material_id == "a")
    assert a.pressure_GPa == 150.0 and a.properties["is_metal"] is True
    assert a.spectra_files[0]["path"].endswith("a2f.dat")
    broken = next(c for c in cands if c.material_id == "broken")
    assert "load_error" in broken.properties          # kept, so S4 can account for it
    # JSON round trip
    rt = Candidate.from_dict(json.loads(json.dumps(a.to_dict())))
    assert rt.composition == a.composition and rt.lattice == a.lattice


def test_airss_uspex_calypso(tmp_path):
    d = tmp_path / "airss"
    d.mkdir()
    (d / "H3S-1.res").write_text(RES)
    (c,) = load_airss(d)
    assert c.material_id == "airss:H3S-1"
    assert c.enthalpy_eV_per_atom == pytest.approx(-123.456 / 4)
    assert c.pressure_GPa == 200.0

    u = tmp_path / "uspex"
    u.mkdir()
    block = POSCAR.replace("H3S test cell", "EA7 extra")
    (u / "gatheredPOSCARS").write_text(block + block.replace("EA7", "EA9"))
    (u / "Individuals").write_text(
        "  Gen   ID    Origin   Composition    Enthalpy(eV)   Volume(A^3)  KPOINTS  SYMM\n"
        "    1    7    Random   [ 1  3 ]       -20.0          27.0         [ 4 4 4 ] 229\n"
        "    1    9    Random   [ 1  3 ]       -16.0          27.0         [ 4 4 4 ] 221\n")
    (u / "INPUT.txt").write_text("150 : ExternalPressure\n")
    cs = load_uspex(u)
    assert [x.material_id for x in cs] == ["uspex:EA7", "uspex:EA9"]
    assert cs[0].enthalpy_eV_per_atom == pytest.approx(-5.0) and cs[0].pressure_GPa == 150.0
    assert cs[1].spacegroup == "221"

    k = tmp_path / "calypso"
    k.mkdir()
    (k / "UCell_1_229.vasp").write_text(POSCAR)
    (k / "Analysis_Output.dat").write_text("Index Enthalpy(eV) SpaceGroup\n1 -4.5 229\n")
    (cc,) = load_calypso(k, pressure_GPa=100.0)
    assert cc.enthalpy_eV_per_atom == -4.5 and cc.spacegroup == "229"


def test_seed_list_and_unique_ids(tmp_path):
    p = tmp_path / "seeds.csv"
    p.write_text("material_id,formula,pressure_GPa,is_metal\nx,LaH10,200,true\nx,YH6,150,false\n")
    cands = load_sources([{"type": "seed_list", "path": str(p)}])
    assert [c.material_id for c in cands] == ["x", "x~2"]
    assert cands[0].properties["is_metal"] is True and cands[1].composition == {"Y": 1, "H": 6}


def test_prototype_substitution():
    proto = Candidate("proto", {}, lattice=[[3, 0, 0], [0, 3, 0], [0, 0, 3]],
                      species=["S", "H", "H", "H"],
                      frac_coords=[[0, 0, 0], [.5, 0, 0], [0, .5, 0], [0, 0, .5]])
    out = substitute([proto], {"S": ["Se", "Te"]})
    assert sorted(c.formula for c in out) == ["SeH3", "TeH3"]
    assert all("unrelaxed" in c.provenance["note"] for c in out)


def test_dedupe_fallback_keeps_lowest_enthalpy(monkeypatch):
    import builtins
    real_import = builtins.__import__

    def no_pmg(name, *a, **k):
        if name.startswith(("pymatgen", "spglib")):
            raise ImportError(name)
        return real_import(name, *a, **k)
    monkeypatch.setattr(builtins, "__import__", no_pmg)
    a = Candidate("a", "H3S", spacegroup="Im-3m", pressure_GPa=200, enthalpy_eV_per_atom=-1.0)
    b = Candidate("b", "H6S2", spacegroup="Im-3m", pressure_GPa=200.2, enthalpy_eV_per_atom=-1.1)
    c = Candidate("c", "H3S", spacegroup="R3m", pressure_GPa=200)
    d = Candidate("d", "H3S", spacegroup=None, pressure_GPa=200)
    e = Candidate("e", "H3S", spacegroup="Im-3m", pressure_GPa=100)
    res = dedupe([a, b, c, d, e])
    assert res.duplicate_of == {"a": "b"}
    assert [x.material_id for x in res.kept] == ["b", "c", "d", "e"]


class _FakeSummary:
    def __init__(self, docs):
        self.docs = docs
        self.kw = None

    def search(self, **kw):
        self.kw = kw
        return self.docs


class _FakeRester:
    def __init__(self, docs):
        self.materials = type("M", (), {})()
        self.materials.summary = _FakeSummary(docs)

    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False

    def get_database_version(self):
        return "FAKE-DB"


def test_materials_project_with_fake_client(monkeypatch):
    monkeypatch.setenv("MP_API_KEY", "dummy")
    docs = [{"material_id": "mp-0001", "formula_pretty": "NbH2", "composition": {"Nb": 1, "H": 2},
             "energy_above_hull": 0.02, "is_metal": True, "band_gap": 0.0, "ordering": "NM",
             "total_magnetization": 0.0, "symmetry": {"symbol": "Fm-3m"}}]
    fake = _FakeRester(docs)
    cands = query_materials_project(elements=["H"], energy_above_hull_max=0.05, is_metal=True,
                                    rester_factory=lambda key: fake)
    (c,) = cands
    assert c.material_id == "mp:mp-0001" and c.composition == {"Nb": 1.0, "H": 2.0}
    assert c.properties["energy_above_hull"] == 0.02 and c.spacegroup == "Fm-3m"
    assert c.provenance["database_version"] == "FAKE-DB"
    assert "energy_above_hull" in c.provenance["requested_fields"]
    assert fake.materials.summary.kw["energy_above_hull"] == (0.0, 0.05)
    assert fake.materials.summary.kw["is_metal"] is True


def test_materials_project_requires_key(monkeypatch):
    from rtsc.sources.mp import MPUnavailable
    monkeypatch.delenv("MP_API_KEY", raising=False)
    with pytest.raises(MPUnavailable):
        query_materials_project(elements=["H"])
