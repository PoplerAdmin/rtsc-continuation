"""ORCH: existence filters (S3) and validity flags (S6)."""
from __future__ import annotations

import pytest

from rtsc.sources import Candidate
from rtsc.spectrum import Spectrum
from rtsc.stability import FilterSet, HullEntry, e_above_hull
from rtsc.status import Status, Verdict
from rtsc.validity import ValidityConfig, apply_flags, compute_flags


def test_hull_lp_exact():
    ents = [HullEntry({"A": 1}, 0.0), HullEntry({"B": 1}, 0.0), HullEntry({"A": 1, "B": 1}, -1.0)]
    val, how = e_above_hull({"A": 1, "B": 3}, 0.0, ents, backend="scipy")
    # hull at x_B = 3/4 lies halfway between AB (-1.0 eV/atom) and B (0): -0.5 eV/atom
    assert val == pytest.approx(0.5)
    assert "scipy" in how
    val2, _ = e_above_hull({"A": 1, "B": 1}, -1.0, ents, backend="scipy")
    assert val2 == pytest.approx(0.0)
    none, why = e_above_hull({"A": 1, "C": 1}, 0.0, ents)
    assert none is None and "C" in why


def test_hull_filter_from_run_candidates_at_pressure():
    cands = [Candidate("A", {"A": 1}, pressure_GPa=100, enthalpy_eV_per_atom=0.0),
             Candidate("B", {"B": 1}, pressure_GPa=100, enthalpy_eV_per_atom=0.0),
             Candidate("AB", {"A": 1, "B": 1}, pressure_GPa=100, enthalpy_eV_per_atom=-1.0),
             Candidate("AB3", {"A": 1, "B": 3}, pressure_GPa=100, enthalpy_eV_per_atom=0.0),
             Candidate("AB3-lowP", {"A": 1, "B": 3}, pressure_GPa=0, enthalpy_eV_per_atom=0.0)]
    fs = FilterSet.from_config({"hull": {"max_e_above_hull_eV_per_atom": 0.1, "backend": "scipy"}})
    fs.prepare(cands)
    v, res = fs.evaluate_static(cands[3])
    assert v is not None and v.status is Status.EXCLUDED_ASSUMPTION
    assert v.rule_ids == ["X-HULL"]
    r = res[0]
    assert r.value == pytest.approx(0.5) and r.threshold == 0.1 and "<= 0.1" in r.comparison
    assert v.certificate["proved"] is False
    # no references at 0 GPa -> not evaluated, not excluded
    v0, res0 = fs.evaluate_static(cands[4])
    assert v0 is None and res0[0].passed is None


def test_mp_hull_value_not_used_at_other_pressure():
    c = Candidate("mp:x", {"Nb": 1, "H": 2}, source="materials_project", pressure_GPa=0.0,
                  properties={"energy_above_hull": 0.5, "target_pressure_GPa": 150.0})
    fs = FilterSet.from_config({})
    fs.prepare([c])
    v, res = fs.evaluate_static(c)
    hull = res[0]
    assert hull.passed is None and "0.0 GPa" in hull.reason


def test_filters_record_and_can_be_disabled():
    c = Candidate("x", {"Fe": 1, "H": 2}, properties={"band_gap": 1.0, "ordering": "FM",
                                                       "energy_above_hull": 0.0,
                                                       "phonon": {"min_frequency": -2.0, "unit": "THz"}})
    fs = FilterSet.from_config({})
    fs.prepare([c])
    v, res = fs.evaluate_static(c)
    assert v.status is Status.EXCLUDED_ASSUMPTION
    assert set(v.rule_ids) == {"X-METAL", "X-MAGNETIC", "X-DYNAMIC"}
    assert all(r.threshold is not None for r in res)
    off = FilterSet.from_config({k: {"enabled": False} for k in ("hull", "metallic", "magnetic", "dynamic")})
    v2, res2 = off.evaluate_static(c)
    assert v2 is None and res2 == []
    with pytest.raises(ValueError):
        FilterSet.from_config({"hull": {"max_ehull": 1}})


def test_dynamic_filter_ignores_acoustic_gamma():
    fs = FilterSet.from_config({"dynamic": {"imaginary_tol": 10, "imaginary_tol_unit": "cm-1"}})
    c = Candidate("x", "NbH")
    ph = {"frequencies": [[-30.0, -1.0, 0.0, 100.0], [-5.0, 50.0, 60.0, 90.0]],
          "qpoints": [[0, 0, 0], [0.5, 0, 0]], "unit": "cm-1"}
    v, r = fs.evaluate_dynamic(c, ph)
    assert v is None and r.passed and r.value == pytest.approx(-5.0)
    ph["frequencies"][1][0] = -20.0
    v, r = fs.evaluate_dynamic(c, ph)
    assert v.status is Status.EXCLUDED_ASSUMPTION


def _sp(A=1000.0, w=1000.0, **prov):
    return Spectrum((w,), (A,), 0.1, "t", dict(prov))


@pytest.mark.parametrize("status", list(Status))
def test_flags_never_change_status(status):
    c = Candidate("x", {"La": 1, "H": 10}, pressure_GPa=200,
                  lattice=[[3, 0, 0], [0, 3, 0], [0, 0, 9]],
                  properties={"fermi_energy_eV": 0.5, "dos_at_fermi": 3.0,
                              "phonon": {"min_frequency": 0.5, "unit": "THz"}},
                  synthetic=True)
    flags, unchecked = compute_flags(c, [_sp(), _sp(A=1500.0, w=500.0)], ValidityConfig(), mu_hi=0.1)
    codes = {f.code for f in flags}
    assert {"SYNTHETIC_DATA", "NONADIABATIC", "ANHARMONIC_H", "SOFT_MODES", "MAGNETIC_PROXIMITY",
            "MULTIBAND_ANISOTROPY", "MU_STAR_UNCERTAIN", "SPECTRUM_SPREAD", "STRONG_COUPLING"} <= codes
    v = Verdict("x", status, 100.0, 200.0, ["r"], ["R"], ["a"], [], {"c": 1})
    v2 = apply_flags(v, flags)
    assert (v2.status, v2.T_low_K, v2.T_up_K, v2.rule_ids, v2.certificate) == \
        (v.status, v.T_low_K, v.T_up_K, v.rule_ids, v.certificate)
    assert set(v2.flags) == codes


def test_flags_unchecked_without_data():
    flags, unchecked = compute_flags(Candidate("y", "MgB2"), [_sp(A=10.0)], ValidityConfig())
    assert any("NONADIABATIC" in u for u in unchecked)
    assert "NONADIABATIC" not in {f.code for f in flags}
