"""phself -> Spectrum tests on SYNTHETIC fixtures (tests/fixtures/SYNTHETIC_*), plus mutations."""
import math
import shutil

import numpy as np
import pytest

import synthetic_epw as syn
from ephproto.phself import (PhselfError, build_phself_spectrum, parse_epw_stdout,
                             parse_fortran_real, parse_lambda_phself, parse_linewidth_phself)
from ephproto.qlist import sha256_file
from ephproto.tags import EPH_PROTOCOL_TAG

LF, WF, QF = syn.LAMBDA_FILE, syn.LINEWIDTH_FILE, syn.QFILE


def args(d, stdout=True):
    return dict(lambda_path=d / LF, linewidth_path=d / WF, filqf_path=d / QF,
                epw_in_path=d / "epw.in", stdout_path=(d / "epw.out") if stdout else None)


@pytest.fixture
def case(tmp_path):
    """Copy of a committed SYNTHETIC fixture into tmp (so tests can mutate it)."""
    n = [0]

    def _get(name):
        n[0] += 1
        dst = tmp_path / f"SYNTHETIC_{name}_{n[0]}"
        shutil.copytree(syn.fixture_dir(name), dst)
        return dst
    return _get


def edit(path, fn):
    path.write_text(fn(path.read_text()), newline="\n")


# ----------------------------------------------------------------------------- fixtures honest
@pytest.mark.parametrize("name", syn.CASES)
def test_committed_fixtures_match_generator(tmp_path, name):
    d = syn.write_case(tmp_path / name, name)
    for f in (LF, WF, QF, "epw.in", "epw.out"):
        assert (d / f).read_bytes() == (syn.fixture_dir(name) / f).read_bytes(), (name, f)


def test_fixture_byte_format():
    d = syn.fixture_dir("normal")
    lam = (d / LF).read_text().splitlines()
    assert lam[0] == "" and lam[1] == "  #Lambda phonon self-energy" and lam[2] == ""
    assert lam[3].startswith(" #Modes     ") and len(lam[4]) == 9 + 4 + 3 * 18
    lw = (d / WF).read_text().splitlines()
    assert all(len(x) == 9 + 6 + 20 + 22 for x in lw[2:])
    assert "0.15000-119" in (syn.fixture_dir("exp_noE") / LF).read_text()


# ----------------------------------------------------------------------------- Fortran reals
@pytest.mark.parametrize("tok,val", [
    ("0.12345E+01", 1.2345), ("-0.12345E-03", -1.2345e-4), ("0.12345-100", 1.2345e-101),
    ("-0.1234567890+101", -1.23456789e100), ("0.10000D+02", 10.0), ("0.00000E+00", 0.0),
    ("-0.00000E+00", 0.0), ("12.3456", 12.3456), (".5E+00", 0.5), ("0.12345E-399", 0.0)])
def test_parse_fortran_real(tok, val):
    assert parse_fortran_real(tok) == val


@pytest.mark.parametrize("tok", ["***************", "NaN", "Infinity", "12", "0.1E", "1.0e+400", ""])
def test_parse_fortran_real_rejects(tok):
    with pytest.raises(PhselfError):
        parse_fortran_real(tok)


# ----------------------------------------------------------------------------- normal case
def test_normal_spectrum(case):
    d = case("normal")
    sp = build_phself_spectrum(**args(d))
    pts, om, lam = syn.model()
    lam32 = np.vectorize(syn.f32)(lam)
    lam_print = np.vectorize(lambda v: float(syn.fE(v, 15, 5)[:-4] + "e" + syn.fE(v, 15, 5)[-3:]))(lam32)
    om_print = np.vectorize(lambda v: float(syn.fE(v, 20, 8).replace("E", "e")))(om)
    Nq = 8
    # independent rebuild of the atoms from the printed values
    ref = {}
    for r in range(Nq):
        for j in range(3):
            ref[om_print[r, j]] = ref.get(om_print[r, j], 0.0) + lam_print[r, j] * om_print[r, j] / 2 / Nq
    from rtsc.spectrum import unit_factor
    f = unit_factor("meV")
    keys = sorted(ref)
    assert list(sp.omega_K) == [k * f for k in keys]
    assert np.allclose(sp.weight_K, [ref[k] * f for k in keys], rtol=1e-15, atol=0)
    lam_avg = lam_print.sum() / Nq
    assert math.isclose(sp.lambda_(), lam_avg, rel_tol=1e-12)
    pv = sp.provenance
    c = pv["counts"]
    assert c["modes_total"] == 24 and c["atoms_nonzero_before_merge"] == 24
    assert c["atoms_after_merge"] == len(keys) < 24          # degenerate + symmetric q merged
    assert c["atoms_merged_away"] == 24 - len(keys) and c["lambda_zero_modes"] == 0
    assert pv["N_q"] == 8 and pv["nmodes"] == 3 and pv["q_grid"] == [2, 2, 2]
    assert pv["sigma_meV"] == 20.0 and pv["sigma_eV"] == 0.02
    assert pv["stdout_checked"] and pv["stdout_checks"]["need_q"] == [8]
    assert pv["stdout_checks"]["ef_eV"] == 8.123456
    assert pv["assumptions"] == [EPH_PROTOCOL_TAG]
    for role, fn in (("lambda", LF), ("linewidth", WF), ("filqf", QF), ("epw_in", "epw.in"),
                     ("stdout", "epw.out")):
        assert pv["sources"][role]["sha256"] == sha256_file(d / fn)
    assert pv["protocol"] == "EPH_PROTOCOL_DRAFT" and "precision" in pv
    # Spectrum JSON round trip is exact
    from rtsc.spectrum import Spectrum
    sp2 = Spectrum.from_json(sp.to_json())
    assert sp2.omega_K == sp.omega_K and sp2.weight_K == sp.weight_K


def test_partial_zero_lambda_is_counted_not_refused(case):
    d = case("normal")
    L = parse_lambda_phself(d / LF)
    lam = L["lam"].copy()
    lam[2, 0] = 0.0                      # one mode of q=3 has exactly zero coupling
    (d / LF).write_text(syn.lambda_text(lam), newline="\n")
    sp = build_phself_spectrum(**args(d, stdout=False), allow_no_stdout=True)
    c = sp.provenance["counts"]
    assert c["lambda_zero_modes"] == 1 and c["lambda_zero_modes_per_q"] == {3: 1}
    assert c["atoms_nonzero_before_merge"] == 23


# ----------------------------------------------------------------------------- required cases
def test_restart_zeros_detected(case):
    d = case("restart_zeros")
    with pytest.raises(PhselfError, match=r"3 q point\(s\) with lambda = 0 for every mode"):
        build_phself_spectrum(**args(d))
    # even when zero rows are accepted, the stdout shows the restart and the missing blocks
    with pytest.raises(PhselfError, match="restarted run.*without a stdout lambda block"):
        build_phself_spectrum(**args(d), allow_zero_lambda_q=True)
    # without stdout: built only as a diagnostic, counts reported, both directions blocked
    sp = build_phself_spectrum(**args(d, stdout=False), allow_zero_lambda_q=True, allow_no_stdout=True)
    c = sp.provenance["counts"]
    assert c["q_all_lambda_zero"] == [1, 2, 3]
    assert c["lambda_zero_modes_per_q"] == {1: 3, 2: 3, 3: 3} and c["lambda_zero_modes"] == 9
    assert any(a.startswith("A-FILE-INCOMPLETE") for a in sp.provenance["assumptions"])


def test_exponent_without_E(case):
    d = case("exp_noE")
    L = parse_lambda_phself(d / LF)
    assert L["lam"][1, 0] == 1.5e-120
    W = parse_linewidth_phself(d / WF)
    assert W["linewidth_meV"][1 * 3 + 0] == 2.25e-121
    assert W["linewidth_meV"][2 * 3 + 2] == 7.123456789e-105
    sp = build_phself_spectrum(**args(d))
    assert sp.provenance["counts"]["lambda_zero_modes"] == 0
    assert min(sp.weight_K) > 0


def test_negative_lambda_refused(case):
    with pytest.raises(PhselfError, match=r"1 negative lambda_qnu \(first: iq=5 mode=3"):
        build_phself_spectrum(**args(case("negative_lambda")))


def test_nonpositive_omega_refused_with_count(case):
    with pytest.raises(PhselfError, match=r"1 mode\(s\) with omega <= eps.*1 of them imaginary.*iq=5 mode=1"):
        build_phself_spectrum(**args(case("omega_nonpos")))


def test_gamma_present_refused(case):
    d = case("gamma_present")
    with pytest.raises(PhselfError, match="Gamma"):
        build_phself_spectrum(**args(d))
    # same phself files with a valid shifted filqf: the omega = 0 modes are caught and counted
    shutil.copy(syn.fixture_dir("normal") / QF, d / QF)
    with pytest.raises(PhselfError, match=r"3 mode\(s\) with omega <= eps.*0 of them imaginary"):
        build_phself_spectrum(**args(d, stdout=False), allow_no_stdout=True)


def test_wrong_nq_refused(case):
    with pytest.raises(PhselfError, match="lambda file has 8 q rows, filqf has N_q = 27"):
        build_phself_spectrum(**args(case("wrong_nq")))


# ----------------------------------------------------------------------------- more mutations
def test_missing_linewidth_row(case):
    d = case("normal")
    edit(d / WF, lambda t: "\n".join(t.splitlines()[:-1]) + "\n")
    with pytest.raises(PhselfError, match="linewidth file has 23 rows, expected"):
        build_phself_spectrum(**args(d))


def test_truncated_lambda_file(case):
    d = case("normal")
    edit(d / LF, lambda t: t[: t.rstrip().rfind(" ")] + "\n")
    with pytest.raises(PhselfError, match="not a multiple"):
        build_phself_spectrum(**args(d))
    d2 = case("normal")
    edit(d2 / LF, lambda t: "\n".join(t.splitlines()[:-1]) + "\n")
    with pytest.raises(PhselfError, match="7 q rows"):
        build_phself_spectrum(**args(d2))


def test_noncontiguous_q_index(case):
    d = case("normal")
    L = parse_lambda_phself(d / LF)
    (d / LF).write_text(syn.lambda_text(L["lam"], iqs=[1, 2, 3, 4, 5, 6, 7, 9]), newline="\n")
    with pytest.raises(PhselfError, match="q indices are not 1..8 contiguous"):
        build_phself_spectrum(**args(d))
    d2 = case("normal")
    W = parse_linewidth_phself(d2 / WF)
    om = W["omega_meV"].reshape(8, 3)
    (d2 / WF).write_text(syn.linewidth_text(om, om * 0, iqs=[1, 2, 3, 3, 5, 6, 7, 8]), newline="\n")
    with pytest.raises(PhselfError, match="linewidth q indices"):
        build_phself_spectrum(**args(d2))


def test_unsorted_frequencies(case):
    d = case("normal")
    W = parse_linewidth_phself(d / WF)
    om = W["omega_meV"].reshape(8, 3).copy()
    om[6] = om[6][::-1]
    (d / WF).write_text(syn.linewidth_text(om, om * 0), newline="\n")
    with pytest.raises(PhselfError, match="not sorted ascending at iq=7"):
        build_phself_spectrum(**args(d, stdout=False), allow_no_stdout=True)


def test_stdout_lambda_mismatch(case):
    d = case("normal")
    L = parse_lambda_phself(d / LF)
    lam = L["lam"].copy()
    lam[3, 2] *= 1.001                  # 1e-3 relative: far above the 5e-5 print precision
    (d / LF).write_text(syn.lambda_text(lam), newline="\n")
    with pytest.raises(PhselfError, match=r"iq=4 mode=3: file lambda.*iq=4: sum_nu lambda"):
        build_phself_spectrum(**args(d))


def test_stdout_tot_mismatch_alone(case):
    d = case("normal")
    edit(d / "epw.out", lambda t: t.replace("lambda___( tot )=", "lambda___( tot )=9.", 1))
    with pytest.raises(PhselfError, match=r"iq=1: sum_nu lambda"):
        build_phself_spectrum(**args(d))


def test_omega_doubled_like_nimage(case):
    d = case("normal")
    W = parse_linewidth_phself(d / WF)
    om = 2.0 * W["omega_meV"].reshape(8, 3)
    (d / WF).write_text(syn.linewidth_text(om, om * 0), newline="\n")
    with pytest.raises(PhselfError, match="nimage>1 doubling"):
        build_phself_spectrum(**args(d))


def test_stdout_nimage_qwindow_and_ngauss(case):
    d = case("normal")
    edit(d / "epw.out", lambda t: t.replace("     Using q-mesh", "     path-images division:  nimage    =       2\n     Using q-mesh"))
    with pytest.raises(PhselfError, match=r"nimage = \[2\]"):
        build_phself_spectrum(**args(d))
    d = case("normal")
    edit(d / "epw.out", lambda t: t.replace("compute        8 q-points", "compute        7 q-points"))
    with pytest.raises(PhselfError, match=r"qwindow kept \[7\] q points"):
        build_phself_spectrum(**args(d))
    d = case("normal")
    edit(d / "epw.out", lambda t: t.replace("ngauss=   0", "ngauss=   1"))
    with pytest.raises(PhselfError, match="ngauss = 1"):
        build_phself_spectrum(**args(d))


@pytest.mark.parametrize("extra,drop,msg", [
    ({}, ("ngaussw",), "ngaussw must be 0"),
    ({"restart": ".true."}, (), "restart must be .false."),
    ({"mp_mesh_q": ".true."}, (), "mp_mesh_q"),
    ({"a2f": ".true."}, (), "a2f must be .false."),
    ({"nstemp": "2"}, (), "nstemp must be set to 1"),
    ({"efermi_read": ".true."}, (), "efermi_read"),
    ({"selecqread": ".true."}, (), "selecqread"),
    ({"fsthick": "0.1"}, (), "fsthick = 0.1 eV"),
    ({"degaussw": "0.2d0"}, (), r"fsthick = 1.0 eV < max\(1.0 eV, 8.0\*sigma"),
    ({"nqf1": "4"}, (), "nqf1 must be absent"),
    ({}, ("delta_approx",), "delta_approx must be .true."),
    ({"nsmear": "3"}, (), "nsmear must be 1"),
])
def test_epw_pins(case, extra, drop, msg):
    d = case("normal")
    (d / "epw.in").write_text(syn.epw_in_text(extra=extra, drop=drop), newline="\n")
    with pytest.raises(PhselfError, match=msg):
        build_phself_spectrum(**args(d))


def test_filqf_name_must_match_epw_in(case):
    d = case("normal")
    shutil.copy(d / QF, d / "other.txt")
    with pytest.raises(PhselfError, match="filqf"):
        build_phself_spectrum(lambda_path=d / LF, linewidth_path=d / WF, filqf_path=d / "other.txt",
                              epw_in_path=d / "epw.in", stdout_path=d / "epw.out")


def test_stdout_required_by_default(case):
    with pytest.raises(PhselfError, match="stdout is required"):
        build_phself_spectrum(**args(case("normal"), stdout=False))


def test_wrapped_header_and_split_records_parse_identically(case):
    d = case("normal")
    a = parse_lambda_phself(d / LF)
    lines = (d / LF).read_text().splitlines()
    # Intel-style wrapped list-directed header and a data record split after the q index
    hdr = [" #Modes                1           2", "           3"]
    body = lines[4:]
    body[0] = body[0][:13] + "\n" + body[0][13:]
    (d / "alt").write_text("\n".join(lines[:3] + hdr + body) + "\n")
    b = parse_lambda_phself(d / "alt")
    assert b["nmodes"] == 3 and b["iq"] == a["iq"] and np.array_equal(a["lam"], b["lam"])


def test_stdout_parser_fields():
    S = parse_epw_stdout(syn.fixture_dir("normal") / "epw.out")
    assert sorted(S["blocks"]) == list(range(1, 9)) and S["need_q"] == [8]
    b = S["blocks"][1]
    assert b["mode"] == [1, 2, 3] and b["tot"] is not None and b["wt"] == 0.125
    assert S["ngauss"] == 0 and S["gauss_eV"] == 0.02 and not S["restart"]
