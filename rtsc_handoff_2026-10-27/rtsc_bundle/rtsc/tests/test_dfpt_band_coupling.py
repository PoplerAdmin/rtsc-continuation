"""Tests for rtsc.dfpt.band_coupling (multiband data for rule R-MULTI) and the EPW input writer.

All alpha^2F_ij, lambda_nk and EPW files here are SYNTHETIC (written in the Fortran formats copied
from EPW/src/supercond.f90, qe-6.7.0), except where a test says REAL.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from rtsc.dfpt import (BandCoupling, band_coupling_from_epw_lambda, band_coupling_from_spectra,  # noqa: E402
                       epw_lambda_dos_average, epw_omega_max_meV, nscf_kpoints_card,
                       parse_epw_lambda_aniso, parse_epw_lambda_fs, write_epw_aniso_input)
from rtsc.dfpt import units  # noqa: E402
from rtsc.spectrum import Spectrum  # noqa: E402

REAL_EPW = Path(__file__).parent / "data_dfpt" / "epw53_pb_smoke"


def two_band_synthetic():
    """SYNTHETIC two-band model with detailed balance N_i A_ij = N_j A_ji."""
    N = (0.3, 0.7)
    A = [[400.0, 150.0], [150.0 * 0.3 / 0.7, 80.0]]
    return N, A


def test_band_coupling_basic_and_kappa():
    N, A = two_band_synthetic()
    bc = BandCoupling(("sigma", "pi"), "computed", A, N, label="SYNTHETIC")
    assert bc.n_bands == 2
    assert bc.row_sums_K() == pytest.approx((550.0, 80.0 + 150 * 0.3 / 0.7))
    assert bc.detailed_balance_asymmetry() == pytest.approx(0.0, abs=1e-15)
    K = np.array([[A[0][0], math.sqrt(N[0] / N[1]) * A[0][1]],
                  [math.sqrt(N[1] / N[0]) * A[1][0], A[1][1]]])
    kap = max(np.linalg.eigvalsh(0.5 * (K + K.T)))
    assert bc.kappa_float() == pytest.approx(kap, rel=1e-12)
    # wMB 4.7 ordering: A_avg <= kappa <= max row sum
    A_avg = sum(N[i] * A[i][j] for i in range(2) for j in range(2)) / sum(N)
    assert A_avg <= bc.kappa_float() <= bc.max_row_sum_K()


def test_band_coupling_validation():
    with pytest.raises(ValueError, match=">= 0"):
        BandCoupling(("a", "b"), "computed", [[1.0, -0.1], [0.0, 1.0]])
    with pytest.raises(ValueError, match="2x2"):
        BandCoupling(("a", "b"), "computed", [[1.0]])
    with pytest.raises(ValueError, match="needs A_K"):
        BandCoupling(("a",), "computed")
    with pytest.raises(ValueError, match="needs row_sum_upper_K"):
        BandCoupling(("a",), "row_sum_upper_bound")
    with pytest.raises(ValueError, match="positive"):
        BandCoupling(("a",), "computed", [[1.0]], N=(0.0,))
    with pytest.raises(ValueError, match="exceeds"):
        BandCoupling(("a",), "computed", [[5.0]], row_sum_upper_K=(4.0,))
    with pytest.raises(ValueError, match="kind"):
        BandCoupling(("a",), "guess", [[1.0]])


def test_band_coupling_from_pair_spectra_and_json_roundtrip(tmp_path):
    s00 = Spectrum((100.0, 300.0), (50.0, 100.0), label="SYNTHETIC a2F_00")
    s01 = Spectrum((200.0,), (30.0,), label="SYNTHETIC a2F_01")
    s10 = Spectrum((200.0,), (70.0,), label="SYNTHETIC a2F_10")
    bc = band_coupling_from_spectra(["a", "b"], {(0, 0): s00, (0, 1): s01, (1, 0): s10},
                                    N=(0.7, 0.3), label="SYNTHETIC")
    assert bc.A_K == ((150.0, 30.0), (70.0, 0.0))
    assert bc.detailed_balance_asymmetry() == pytest.approx(0.0, abs=1e-12)
    d = json.loads(bc.to_json())
    p = tmp_path / "bc.json"
    p.write_text(json.dumps({**d, "a2f_pairs": []}))
    bc2 = BandCoupling.from_file(p)
    assert bc2.A_K == bc.A_K and bc2.N == bc.N and bc2.kind == "computed"


def test_band_coupling_toml_units_and_pairs(tmp_path):
    (tmp_path / "p01.dat").write_text("# SYNTHETIC omega[meV] a2F\n10 0.5\n20 0.5\n30 0.5\n")
    toml = tmp_path / "bc.toml"
    toml.write_text("""
# SYNTHETIC band-coupling input
[band_coupling]
labels = ["s", "p"]
unit = "meV"
N = [1.0, 1.0]
rel_envelope = 0.1
[[band_coupling.a2f_pairs]]
i = 0
j = 0
unit = "meV"
omega = [10.0, 20.0, 30.0]
a2f = [1.0, 1.0, 1.0]
[[band_coupling.a2f_pairs]]
i = 0
j = 1
file = "p01.dat"
[[band_coupling.a2f_pairs]]
i = 1
j = 0
file = "p01.dat"
[[band_coupling.a2f_pairs]]
i = 1
j = 1
omega = [10.0, 20.0]
a2f = [0.0, 0.0]
""")
    with pytest.raises(ValueError):          # a pair with zero weight -> empty spectrum
        BandCoupling.from_file(toml, table="band_coupling")
    txt = toml.read_text().replace("a2f = [0.0, 0.0]", "a2f = [0.2, 0.2]")
    toml.write_text(txt)
    bc = BandCoupling.from_file(toml, table="band_coupling")
    # bins rule on a 10 meV grid: A_00 = 3 * 1.0 * 10 meV
    assert bc.A_K[0][0] == pytest.approx(30.0 * units.MEV_TO_K)
    assert bc.A_K[0][1] == pytest.approx(15.0 * units.MEV_TO_K)
    assert bc.A_K[1][1] == pytest.approx(4.0 * units.MEV_TO_K)
    assert bc.rel_envelope == 0.1 and bc.provenance["input_unit"] == "meV"
    assert set(bc.a2f_pairs) == {"0,0", "0,1", "1,0", "1,1"}


def test_band_coupling_from_lambda_fallback():
    d = {"labels": ["a", "b"], "unit": "meV", "lambda": [[1.0, 0.2], [0.3, 0.5]],
         "omega_max": 100.0}
    bc = BandCoupling.from_dict(d)
    assert bc.kind == "entrywise_upper_bound"
    assert bc.A_K[0][0] == pytest.approx(50.0 * units.MEV_TO_K)
    with pytest.raises(ValueError, match="omega_max"):
        BandCoupling.from_dict({"labels": ["a"], "lambda": [[1.0]]})


# ----------------------------------------------------------------------------- EPW synthetic
def write_lambda_fs(path: Path, rows) -> None:
    """SYNTHETIC <prefix>.lambda_FS with EPW's '(3f12.6, i8, f12.6, f24.15)' rows."""
    lines = ["#               k-point                  Band Enk-Ef [eV]            lambda"]
    for kx, ky, kz, b, e, lam in rows:
        lines.append(f"{kx:12.6f}{ky:12.6f}{kz:12.6f}{b:8d}{e:12.6f}{lam:24.15f}")
    path.write_text("\n".join(lines) + "\n")


def test_epw_lambda_fs_row_sum_bound(tmp_path):
    rng = np.random.default_rng(3)
    rows = []
    for i in range(40):
        b = 1 + (i % 2)
        rows.append((rng.uniform(-1, 1), rng.uniform(-1, 1), rng.uniform(-1, 1), b,
                     rng.uniform(-0.3, 0.3), (1.2 if b == 1 else 0.4) + 0.1 * rng.uniform()))
    p = tmp_path / "x.lambda_FS"
    write_lambda_fs(p, rows)
    st = parse_epw_lambda_fs(p)
    assert len(st.lam) == 40 and st.k.shape == (40, 3)
    lam = np.array([r[5] for r in rows])
    bc = band_coupling_from_epw_lambda(st, omega_max=20.0, omega_unit="meV", label="SYNTHETIC")
    assert bc.kind == "row_sum_upper_bound" and bc.n_bands == 1
    assert bc.max_row_sum_K() == pytest.approx(0.5 * 20.0 * units.MEV_TO_K * lam.max(), rel=1e-12)
    assert any("eps_acustic" in c for c in bc.provenance["caveats"])
    # grouping (by EPW window band index, purely illustrative) with max reduction
    g = band_coupling_from_epw_lambda(st, omega_max=20.0, groups=lambda k, b, e: f"b{b}")
    assert g.labels == ("b1", "b2") and g.max_row_sum_K() == pytest.approx(bc.max_row_sum_K())
    # DOS-average (band-isotropic model) is lower and flagged as an assumption
    a = band_coupling_from_epw_lambda(st, omega_max=20.0, groups=lambda k, b, e: f"b{b}",
                                      reduction="dos_average", degaussw_eV=0.1)
    assert a.max_row_sum_K() <= g.max_row_sum_K()
    assert any("M-BAND-ISO" in s for s in a.provenance["assumptions"])
    assert sum(a.N) == pytest.approx(1.0)
    w = np.exp(-(st.e_eV / 0.1) ** 2)
    assert epw_lambda_dos_average(st, 0.1) == pytest.approx(float(np.sum(w * lam) / np.sum(w)), rel=1e-6)
    with pytest.raises(ValueError):
        band_coupling_from_epw_lambda(st, omega_max=-1.0)


def test_epw_lambda_aniso_and_omega_max(tmp_path):
    p = tmp_path / "x.lambda_aniso"
    p.write_text("# Enk-Ef[eV]  lambda_nk  # kpt  # band\n"
                 + "".join(f"{e:12.7f}{l:12.7f}{k:7d}{b:7d}\n"
                           for e, l, k, b in [(-0.1, 0.8, 1, 1), (0.05, 1.3, 2, 1)]))
    st = parse_epw_lambda_aniso(p)
    assert not st.full_grid and list(st.kindex) == [1, 2]
    assert band_coupling_from_epw_lambda(st, omega_max=1.0, omega_unit="K").max_row_sum_K() == pytest.approx(0.65)
    with pytest.raises(ValueError, match="full-grid"):
        band_coupling_from_epw_lambda(st, omega_max=1.0, reduction="dos_average", degaussw_eV=0.1)
    a2f = tmp_path / "x.a2f"
    a2f.write_text(" w[meV] a2f and integrated 2*a2f/w for    1 smearing values\n"
                   "   0.5000000   0.0100000   0.0400000\n   1.0000000   0.0200000   0.0800000\n"
                   " Integrated el-ph coupling\n  #            0.0800000\n")
    assert epw_omega_max_meV(a2f) == 1.0


# ----------------------------------------------------------------------------- EPW input writer
def test_write_epw_aniso_input():
    txt = write_epw_aniso_input(prefix="mgb2", masses=[24.305, 10.811], nk=(6, 6, 6), nq=(3, 3, 3),
                                nkf=(24, 24, 24), nqf=(12, 12, 12), nbndsub=5,
                                proj=["B:pz", "f=0.5,1.0,0.5:s"], dis_froz_max=8.8,
                                temps_K=[10.0, 20.0])
    for frag in ("&inputepw", "laniso = .true.", "eliashberg = .true.", "ephwrite = .true.",
                 "iverbosity = 2", "nkf1 = 24", "nqf3 = 12", "proj(2) = 'f=0.5,1.0,0.5:s'",
                 "amass(2) = 10.811", "temps(2) = 20.0", "nstemp = 2", "dis_froz_max = 8.8"):
        assert frag in txt, frag
    with pytest.raises(ValueError, match="multiple"):
        write_epw_aniso_input(prefix="x", masses=[1.0], nk=(2, 2, 2), nq=(2, 2, 2), nkf=(10, 10, 10),
                              nqf=(4, 4, 4), nbndsub=1, proj=["X:s"])
    with pytest.raises(ValueError, match="namelist"):
        write_epw_aniso_input(prefix="x", masses=[1.0], nk=(2, 2, 2), nq=(2, 2, 2), nkf=(4, 4, 4),
                              nqf=(2, 2, 2), nbndsub=1, proj=["X:s"], extra={"not_a_key": 1})
    with pytest.raises(ValueError, match="required"):
        write_epw_aniso_input(prefix="x", masses=[1.0], nk=(2, 2, 2), nq=(2, 2, 2), nkf=(4, 4, 4),
                              nqf=(2, 2, 2), nbndsub=1, proj=[])


def test_nscf_kpoints_card():
    c = nscf_kpoints_card((2, 2, 3)).splitlines()
    assert c[0] == "K_POINTS crystal" and c[1] == "12" and len(c) == 14
    assert c[3].split()[:3] == ["0.00000000", "0.00000000", "0.33333333"]   # z fastest
    assert sum(float(l.split()[3]) for l in c[2:]) == pytest.approx(1.0, rel=1e-5)


# ----------------------------------------------------------------------------- REAL EPW output
@pytest.mark.skipif(not REAL_EPW.exists(), reason="REAL EPW fixture not present")
def test_real_epw53_pb_smoke():
    """REAL EPW 5.3 output of a deliberately unconverged fcc Pb run (see PROVENANCE.txt)."""
    st = parse_epw_lambda_fs(REAL_EPW / "pb.lambda_FS")
    assert len(st.lam) > 0 and np.all(st.lam >= 0)
    wmax = epw_omega_max_meV(REAL_EPW / "pb.a2f")
    bc = band_coupling_from_epw_lambda(st, omega_max=wmax, omega_unit="meV", label="REAL smoke")
    assert bc.max_row_sum_K() == pytest.approx(0.5 * wmax * units.MEV_TO_K * st.lam.max())
    assert len(st.lam) == 120 and set(st.band.tolist()) == {1, 2}
    # cross-checks against EPW's isotropic numbers (same run)
    from rtsc.dfpt import (parse_epw_a2f, epw_a2f_to_spectra, check_epw_spectrum,
                           parse_epw_summed_lambda)
    tab = parse_epw_a2f(REAL_EPW / "pb.a2f")
    sp = epw_a2f_to_spectra(tab, label_prefix="REAL-smoke:")[0]
    assert all(c.ok for c in check_epw_spectrum(sp))
    # the per-state row-sum bound dominates the isotropic total coupling A (wMB 4.7 ordering)
    assert bc.max_row_sum_K() >= sp.A()
    # DOS-weighted average of lambda_nk vs EPW's summed lambda (normalisation check, ~2 %)
    lam_avg = epw_lambda_dos_average(st, 0.1)
    assert lam_avg == pytest.approx(parse_epw_summed_lambda(REAL_EPW / "pb.a2f"), rel=0.03)
