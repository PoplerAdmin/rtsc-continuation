"""ML pre-screen adapter tests.  The calculator is ASE's EMT (a cheap empirical potential standing
in for MatterSim so the tests need no torch); numbers are NOT MatterSim results.  The real
MatterSim smoke test is tests/test_mlff_mattersim_real.py (skipped unless mattersim is installed)."""
from __future__ import annotations

from pathlib import Path

import pytest

ase = pytest.importorskip("ase")
from ase.build import bulk  # noqa: E402
from ase.calculators.emt import EMT  # noqa: E402

from rtsc.discovery.interfaces import GeneratedStructure, Prescreener, PrescreenResult  # noqa: E402
from rtsc.mlff._atoms import atoms_to_structure  # noqa: E402
from rtsc.mlff.hull import MLHull, fallback_element, subsystems  # noqa: E402
from rtsc.mlff.mattersim import MatterSimConfig, MatterSimPrescreener, relax_at_pressure  # noqa: E402
from rtsc.mlff.phonons import auto_supercell, light_element_hf_fraction, omega_log_of_dos  # noqa: E402


def gs(atoms, fp="fp-1"):
    return GeneratedStructure(structure=atoms_to_structure(atoms), generator="test", request_id="r",
                              provenance={"fingerprint": fp})


def test_relax_at_pressure_emt():
    a0, h0, i0 = relax_at_pressure(bulk("Al", "fcc", a=4.1), EMT(), 0.0)
    a1, h1, i1 = relax_at_pressure(bulk("Al", "fcc", a=4.1), EMT(), 20.0)
    assert i0["converged"] and i1["converged"]
    assert i1["volume_A3_per_atom"] < i0["volume_A3_per_atom"]
    assert abs(i0["pressure_from_stress_GPa"]) < 0.5 and abs(i1["pressure_from_stress_GPa"] - 20.0) < 0.5
    assert h1 > h0          # enthalpy grows with pressure


def test_hull_cache_and_value(tmp_path):
    refs = tmp_path / "refs"
    refs.mkdir()
    from ase.io import write
    write(refs / "Al.cif", bulk("Al", "fcc", a=4.05))
    write(refs / "Cu.cif", bulk("Cu", "fcc", a=3.6))
    calls = []

    def relax(a, p):
        calls.append(a.get_chemical_formula())
        return relax_at_pressure(a, EMT(), p)
    h = MLHull(relax_fn=relax, model_id="EMT-test", reference_dir=str(refs), cache_dir=str(tmp_path / "c"))
    ra, hh, _ = relax(bulk("Al", "fcc", a=4.05) * (2, 1, 1), 5.0)       # Al itself -> on hull
    ra.symbols[1] = "Cu"
    ra2, hh2, _ = relax(ra, 5.0)
    n0 = len(calls)
    val, how, labels = h.e_above_hull(ra2.get_chemical_symbols(), hh2, 5.0)
    assert val is not None and val >= 0 and len(labels) == 2
    assert len(list((tmp_path / "c").glob("*.json"))) == 3              # Al, Cu, Al-Cu subsystems
    h2 = MLHull(relax_fn=relax, model_id="EMT-test", reference_dir=str(refs), cache_dir=str(tmp_path / "c"))
    n1 = len(calls)
    assert h2.e_above_hull(ra2.get_chemical_symbols(), hh2, 5.0)[0] == val
    assert len(calls) == n1 and n1 - n0 == 2                           # cached: no new relaxations
    assert subsystems(["Cu", "Al"]) == [("Al",), ("Cu",), ("Al", "Cu")]


def test_fallback_elements():
    a, how = fallback_element("H")
    assert how == "fallback:dimer-sc" and len(a) == 2
    a, how = fallback_element("Al")
    assert "bulk" in how


def test_heuristic_helpers():
    import numpy as np
    f = np.linspace(0, 3000, 31)
    p = np.zeros((2, 31))
    p[0, :] = 1.0        # heavy atom everywhere
    p[1, 20:] = 1.0      # light atom only at >= 2000 cm-1
    assert abs(light_element_hf_fraction(f, p, ["La", "H"], ("H",), 500.0) - 11 / 42) < 1e-12
    assert light_element_hf_fraction(f, p, ["La", "La"], ("H",), 500.0) == 0.0
    assert omega_log_of_dos([100.0], [1.0]) == pytest.approx(100 * 1.4387768775039338)
    sc = auto_supercell([[3, 0, 0], [0, 3, 0], [0, 0, 10]], 7.0, 1000, 2)
    assert [sc[i][i] for i in range(3)] == [3, 3, 1]
    sc = auto_supercell([[3, 0, 0], [0, 3, 0], [0, 0, 3]], 7.0, 30, 1)
    assert sc[0][0] * sc[1][1] * sc[2][2] <= 30


def test_prescreener_protocol_and_soundness(tmp_path):
    pytest.importorskip("phonopy")
    cfg = MatterSimConfig.from_dict({"cache_dir": str(tmp_path), "phonon": {"supercell": [2, 2, 2],
                                                                            "mesh": [4, 4, 4]}})
    ps = MatterSimPrescreener(cfg, calculator_factory=EMT, model_id="EMT-test")
    assert isinstance(ps, Prescreener)
    good = gs(bulk("Al", "fcc", a=4.0), "fp-al")
    bad = GeneratedStructure(structure={"lattice": [[0, 0, 0]] * 3, "species": ["Al"],
                                        "frac_coords": [[0, 0, 0]]}, generator="t", request_id="r",
                             provenance={"fingerprint": "fp-bad"})
    r1, r2 = ps.prescreen([good, bad], 0.0)
    assert isinstance(r1, PrescreenResult) and r1.fingerprint == "fp-al" and r1.ok
    assert r1.e_above_hull_ml_eV is not None and r1.e_above_hull_ml_eV < 1e-3
    assert r1.dynamically_stable_ml is True and r1.min_freq_cm1 > -15
    assert r1.relaxed_structure["pressure_GPa"] == 0.0 and r1.phonon_dos["freq_cm1"]
    assert MatterSimPrescreener.heuristic_score(r1) == 0.0          # no light elements in Al
    assert any("HEURISTIC" in n for n in r1.notes)
    assert r2.fingerprint == "fp-bad" and r2.ok is False
    assert not hasattr(r1, "status")


def test_adapters_never_name_proven_or_excluded_statuses():
    root = Path(__file__).resolve().parents[1] / "rtsc"
    for f in list((root / "mlff").glob("*.py")) + list((root / "generate").glob("*.py")):
        text = f.read_text()
        import re
        assert not re.search(r"Status\.[A-Z]", text) and "from ..status" not in text, f
