"""Unit tests for ephproto.bandcheck (gate G1). All fixtures here are SYNTHETIC."""
from fractions import Fraction

import numpy as np
import pytest

from ephproto.bandcheck import (HARTREE_EV, BandCheckError, compare_bands, filkf_text, g1_verdict,
                                g1_windows, parse_band_eig, parse_epw_fermi, parse_pw_xml,
                                pw_bands_input, validation_kpoints, write_kpoint_files, _nexclude)


# ------------------------------------------------------------------------ k points
def test_validation_kpoints_5_on_8():
    keep, drop = validation_kpoints((5, 5, 5), (8, 8, 8))
    half = Fraction(1, 2)
    assert drop == [(half, half, half)]                  # L is the only coarse-grid point
    assert len(keep) == 124
    for p in keep:
        assert not all((c * 8).denominator == 1 for c in p)


def test_write_kpoint_files(tmp_path):
    meta = write_kpoint_files(tmp_path)
    assert meta["N_k"] == 124
    lines = (tmp_path / "kval.txt").read_text().splitlines()
    assert lines[0] == "124 crystal" and len(lines) == 125
    assert [float(x) for x in lines[1].split()] == [0.1, 0.1, 0.1, 1 / 124]
    pw = (tmp_path / "kpoints_pw.txt").read_text().splitlines()
    assert pw[0] == "K_POINTS crystal" and pw[1] == "124" and len(pw) == 126
    assert (tmp_path / "gamma.txt").read_text().split() == ["1", "crystal", "0.0", "0.0", "0.0", "1.0"]


def test_random_kset(tmp_path):
    from ephproto.bandcheck import random_kpoints
    keep, drop = random_kpoints(1000, 20261027)
    assert len(keep) + len(drop) == 1000 and len(keep) >= 990
    assert random_kpoints(1000, 20261027)[0] == keep                     # deterministic
    meta = write_kpoint_files(tmp_path, kset="random", nrand=200, seed=1)
    assert meta["kset"] == "random" and meta["N_k"] >= 195


def test_empty_gate_window_is_inconclusive_not_pass():
    k, E, W = _synthetic()
    E[:, 3] = W[:, 1] = 9.5                         # SYNTHETIC: nothing within E_F +- 80 meV
    cmp = compare_bands(k, E, k, W, EF, 2, WIN)
    g = cmp["windows"][0]
    assert g["n_states"] == 0 and g["inconclusive_no_states"] and not g["pass"]
    cmp["N_k"] = 124
    v = g1_verdict(cmp, EF, EF)
    assert v["gate_inconclusive_no_states"] and not v["G1_pass"]


def test_too_few_kpoints(tmp_path):
    with pytest.raises(BandCheckError):
        write_kpoint_files(tmp_path, n=(4, 4, 4))


def test_filkf_roundtrip_exact():
    keep, _ = validation_kpoints()
    txt = filkf_text(keep)
    vals = [list(map(float, ln.split())) for ln in txt.splitlines()[1:]]
    assert all(v[i] == float(p[i]) for v, p in zip(vals, keep) for i in range(3))


def test_nexclude():
    assert _nexclude("1-4") == 4 and _nexclude(None) == 0
    with pytest.raises(BandCheckError):
        _nexclude("2-5")


# ------------------------------------------------------------------------- parsers
SYNTHETIC_BAND_EIG = (
    "  &plot nbnd=   4, nks=     2 /\n"
    "          0.100000  0.200000 -0.300000\n"
    "     -1.000000000000      9.500000000000     11.100000000000     14.250000000000\n"
    "          0.000000  0.000000  0.000000\n"
    "      2.000000000000      3.000000000000      4.000000000000      5.000000000000\n")


def test_parse_band_eig(tmp_path):
    p = tmp_path / "band.eig"
    p.write_text(SYNTHETIC_BAND_EIG)
    k, E = parse_band_eig(p)
    assert k.shape == (2, 3) and E.shape == (2, 4)
    assert k[0].tolist() == [0.1, 0.2, -0.3] and E[0, 2] == 11.1


def test_parse_band_eig_fixed_width_run_together(tmp_path):
    # SYNTHETIC: f20.12 fields with no separating blank (very large values)
    e = "".join(f"{x:20.12f}" for x in [1.0, 2.0])
    p = tmp_path / "band.eig"
    p.write_text(f"&plot nbnd=   2, nks=     1 /\n          0.0 0.0 0.0\n{e}\n")
    assert parse_band_eig(p)[1].tolist() == [[1.0, 2.0]]


def test_parse_band_eig_bad_count(tmp_path):
    p = tmp_path / "band.eig"
    p.write_text(SYNTHETIC_BAND_EIG.replace("nks=     2", "nks=     3"))
    with pytest.raises(BandCheckError):
        parse_band_eig(p)


SYNTHETIC_XML = """<?xml version="1.0"?>
<qes:espresso xmlns:qes="http://www.quantum-espresso.org/ns/qes/qes-1.0">
 <output><band_structure>
  <ks_energies><k_point weight="1.0">0.1 0.2 -0.3</k_point>
   <npw>100</npw><eigenvalues size="3">0.1 0.2 0.4</eigenvalues>
   <occupations size="3">1 1 0</occupations></ks_energies>
  <ks_energies><k_point weight="1.0">0.0 0.0 0.0</k_point>
   <npw>100</npw><eigenvalues size="3">0.0 0.5 1.0</eigenvalues>
   <occupations size="3">1 1 0</occupations></ks_energies>
 </band_structure></output></qes:espresso>
"""


def test_parse_pw_xml(tmp_path):
    p = tmp_path / "d.xml"
    p.write_text(SYNTHETIC_XML)
    k, E = parse_pw_xml(p)
    assert k.shape == (2, 3) and np.allclose(E[1], [0, 0.5 * HARTREE_EV, HARTREE_EV])


def test_parse_epw_fermi(tmp_path):
    p = tmp_path / "epw.out"
    p.write_text("     Fermi energy coarse grid =  11.109035 eV\n"
                 "     Fermi energy is calculated from the fine k-mesh: Ef =  14.070463 eV\n")
    assert parse_epw_fermi(p) == {"ef_fine_eV": 14.070463, "ef_coarse_eV": 11.109035}
    p.write_text("nothing\n")
    with pytest.raises(BandCheckError):
        parse_epw_fermi(p)


def test_pw_bands_input():
    nscf = ("&control\n calculation = 'nscf'\n outdir = './tmp/'\n tprnfor = .true.\n/\n"
            "&system\n nbnd = 20\n/\nK_POINTS crystal\n2\n0 0 0 1\n0.5 0 0 1\n")
    out = pw_bands_input(nscf, "K_POINTS crystal\n1\n0.1 0.1 0.1 1.0\n")
    assert "calculation = 'bands'" in out and "tprnfor" not in out
    assert out.count("K_POINTS") == 1 and "0.1 0.1 0.1 1.0" in out


# ------------------------------------------------------------------------ matching
EF = 10.0
WIN = [{"name": "gate", "half_width_eV": 0.08, "tol_eV": 0.010, "gating": True},
       {"name": "wide", "half_width_eV": 1.0, "tol_eV": 0.030, "gating": False}]


def _synthetic(nk=3):
    k = np.arange(nk * 3, dtype=float).reshape(nk, 3) / 10
    # DFT: 2 excluded semicore bands, then valence-like bands
    E = np.array([[-80.0, -79.0, 5.0, 9.95, 10.5, 13.0, 14.0]] * nk)
    W = E[:, 2:6].copy()
    return k, E, W


def test_perfect_match():
    k, E, W = _synthetic()
    r = compare_bands(k, E, k, W, EF, 2, WIN)
    g, w = r["windows"]
    assert g["n_states"] == 3 and g["max_abs_dE_meV"] == 0 and g["pass"]
    assert w["n_states"] == 6 and w["pass"]


def test_shift_inside_gate_tolerance():
    k, E, W = _synthetic()
    W[:, 1] += 0.004
    W[1, 2] -= 0.02                  # 10.48: inside wide window, outside gate window
    g, w = compare_bands(k, E, k, W, EF, 2, WIN)["windows"]
    assert g["max_abs_dE_meV"] == pytest.approx(4.0) and g["pass"]
    assert w["max_abs_dE_meV"] == pytest.approx(20.0) and w["pass"]
    assert w["worst"]["ik"] == 2 and w["worst"]["band"] == 3


def test_count_mismatch_spurious_wannier_state():
    k, E, W = _synthetic()
    W[0, 3] = 10.03                  # spurious Wannier state in the gate window (DFT partner at 13)
    g = compare_bands(k, E, k, W, EF, 2, WIN)["windows"][0]
    assert g["n_k_count_mismatch"] == 1 and g["count_mismatches"][0] == {"ik": 1, "n_dft": 1, "n_wann": 2}
    assert g["max_abs_dE_meV"] == pytest.approx(470.0) and not g["pass"]   # sorted: 10.03 pairs with 10.5


def test_missing_wannier_state_beyond_nwann():
    # SYNTHETIC: a 5th non-excluded DFT band in the wide window has no Wannier partner
    k, E, W = _synthetic()
    E[:, 6] = 10.9
    w = compare_bands(k, E, k, W, EF, 2, WIN)["windows"][1]
    assert w["n_k_count_mismatch"] == 3 and not w["pass_count"]
    assert w["max_abs_dE_meV"] == pytest.approx(2100.0)       # sorted pairing: DFT 10.9 vs Wannier 13.0


def test_unsorted_inputs_are_sorted():
    k, E, W = _synthetic()
    W = W[:, ::-1]
    E = np.concatenate([E[:, :2], E[:, 2:][:, ::-1]], axis=1)
    assert compare_bands(k, E, k, W, EF, 2, WIN)["windows"][1]["max_abs_dE_meV"] == 0


def test_k_mismatch_raises():
    k, E, W = _synthetic()
    k2 = k.copy()
    k2[1, 0] += 1e-4
    with pytest.raises(BandCheckError):
        compare_bands(k, E, k2, W, EF, 2, WIN)


def test_too_few_dft_bands():
    k, E, W = _synthetic()
    with pytest.raises(BandCheckError):
        compare_bands(k, E[:, :5], k, W, EF, 2, WIN)


def test_g1_windows_and_verdict():
    gw = g1_windows(0.020)
    assert gw[0]["half_width_eV"] == pytest.approx(0.08) and gw[0]["tol_eV"] == pytest.approx(0.010)
    assert g1_windows(0.100)[0]["tol_eV"] == pytest.approx(0.030)       # capped at 30 meV
    k, E, W = _synthetic()
    cmp = compare_bands(k, E, k, W, EF, 2, WIN)
    v = g1_verdict(cmp, EF + 0.03, EF)
    assert v["dE_ok"] and v["count_ok"] and v["ef_ok"] and not v["N_k_ok"] and not v["G1_pass"]
    cmp["N_k"] = 124
    assert g1_verdict(cmp, EF + 0.03, EF)["G1_pass"]
    assert not g1_verdict(cmp, EF + 0.06, EF)["G1_pass"]
    assert not g1_verdict(cmp, None, EF)["G1_pass"]
