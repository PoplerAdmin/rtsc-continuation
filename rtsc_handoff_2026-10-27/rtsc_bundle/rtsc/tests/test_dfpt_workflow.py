"""Workflow tests for rtsc.dfpt: input writers, grids, launcher, pseudo table, structure I/O,
resumability, invalidation, timeouts and error capture.

The workflow tests run FAKE QE executables (a small Python script installed as pw.x, ph.x,
q2r.x, matdyn.x in a temp dir). They write the files the real codes write -- including a
SYNTHETIC a2F.dos in matdyn's Fortran format -- so that the control flow can be tested without
QE. Nothing produced by the fake codes is physical data.
"""
from __future__ import annotations

import json
import stat
import sys
import textwrap
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from rtsc.dfpt import (CellData, EPhSettings, EPhWorkflow, Launcher, PseudoTable,  # noqa: E402
                       WorkflowError, parse_pw_final_structure, select_grids, to_celldata)
from rtsc.dfpt import qe_inputs as qi  # noqa: E402
from rtsc.dfpt.grids import EPhGrids  # noqa: E402

REAL = Path(__file__).parent / "data_dfpt" / "qe67_al_smoke"

FAKE_QE = r'''#!{python}
# FAKE Quantum ESPRESSO executable for tests (writes SYNTHETIC outputs only)
import os, re, sys, time, math
code = os.path.basename(sys.argv[0])
inp = open(sys.argv[sys.argv.index("-in") + 1]).read()
log = os.environ.get("FAKEQE_LOG")
if log:
    open(log, "a").write(code + " " + sys.argv[sys.argv.index("-in") + 1] + "\n")
if os.environ.get("FAKEQE_SLEEP_" + code.replace(".", "_")):
    time.sleep(float(os.environ["FAKEQE_SLEEP_" + code.replace(".", "_")]))
print("     Program %s v.0.0-FAKE starts" % code.upper().replace(".X", ""))
if os.environ.get("FAKEQE_FAIL") == code:
    print(" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    print("     Error in routine fake_routine (1):")
    print("     S matrix not positive definite")
    print(" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    open("CRASH", "w").write("fake crash file\n")
    sys.exit(1)
def nsig():
    m = re.search(r"el_ph_nsigma\s*=\s*(\d+)", inp)
    return int(m.group(1)) if m else 10
if code == "pw.x":
    prefix = re.search(r"prefix\s*=\s*'([^']+)'", inp).group(1)
    os.makedirs("tmp", exist_ok=True)
    if "vc-relax" in inp:
        cell = inp[inp.index("CELL_PARAMETERS"):inp.index("K_POINTS")]
        print("     P=   10.00\n     bfgs converged in   2 scf cycles and   1 bfgs steps")
        print("Begin final coordinates\n" + cell.replace("CELL_PARAMETERS angstrom", "CELL_PARAMETERS (angstrom)")
              .replace("ATOMIC_POSITIONS crystal", "ATOMIC_POSITIONS (crystal)") + "End final coordinates")
        print("     P=   10.00")
    else:
        print("     convergence has been achieved in   5 iterations")
        if "la2F" in inp:
            open("tmp/%s.a2Fsave" % prefix, "w").write("fake\n")
elif code == "ph.x":
    open("dyn0", "w").write("   2   2   2\n   1\n   0.0 0.0 0.0\n")
    open("dyn1", "w").write("fake\n")
    os.makedirs("elph_dir", exist_ok=True)
    for i in range(1, nsig() + 1):
        open("elph_dir/a2Fq2r.%d.1" % (50 + i), "w").write("fake\n")
elif code == "q2r.x":
    open("ifc.fc", "w").write("fake\n")
    for i in range(1, nsig() + 1):
        open("elph_dir/a2Fmatdyn.%d" % (60 + i), "w").write("fake\n")
elif code == "matdyn.x":
    ns = nsig()
    lam_lines = ["", "  Electron-phonon coupling constant, lambda ", ""]
    for isig in range(1, ns + 1):
        delta, nd = 4e-5, 100
        rows, lam, num, den = [], 0.0, 0.0, 0.0
        for n in range(1, nd + 1):
            E = (n - 0.5) * delta
            a = (0.5 + 0.02 * isig) * math.exp(-0.5 * ((E - 0.002) / 0.0003) ** 2)
            lam += 2 * a / E * delta
            num += 2 * a / E * delta * math.log(E)
            rows.append("   %16.6E%16.6E%16.6E" % (E, a, a))
        with open("a2F.dos%d" % isig, "w") as f:
            f.write("\n # Eliashberg function a2F (per both spin)\n #  frequencies in Rydberg  \n"
                    " # DOS normalized to E in Rydberg: a2F_total, a2F(mode) \n # SYNTHETIC (fake matdyn)\n\n")
            f.write("\n".join(rows) + "\n")
            f.write("  lambda =  %r         Delta =    %r\n" % (lam, delta))
        wln = math.exp(num / lam) * 157887.51240203338
        lam_lines.append(" Broadening %8.4f lambda %12.4f dos(Ef)%8.4f omega_ln [K]%12.4f"
                         % (0.005 * isig, lam, 3.0, wln))
    open("lambda", "w").write("\n".join(lam_lines) + "\n")
    open("freq.dos", "w").write(" &plot nbnd=   3, nks=   1 /\n            0.000000  0.000000  0.000000  1.000000\n    0.0000    0.0000    0.0000\n")
print("=" * 78 + "\n   JOB DONE.\n" + "=" * 78)
'''


def install_fake_qe(bindir: Path) -> None:
    bindir.mkdir(parents=True, exist_ok=True)
    src = FAKE_QE.replace("{python}", sys.executable)
    for code in ("pw.x", "ph.x", "q2r.x", "matdyn.x"):
        p = bindir / code
        p.write_text(src)
        p.chmod(p.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


SYNTH_UPF = textwrap.dedent('''\
    <UPF version="2.0.1">
    <PP_INFO> SYNTHETIC header-only pseudopotential for tests </PP_INFO>
    <PP_HEADER generated="test" element="{el}" pseudo_type="NC" functional="PBE"
       wfc_cutoff="{w}" rho_cutoff="{r}" />
    </UPF>
    ''')


def make_pseudos(d: Path, spec: dict[str, tuple[float, float]]) -> PseudoTable:
    d.mkdir(parents=True, exist_ok=True)
    files = {}
    for el, (w, r) in spec.items():
        (d / f"{el}.test.UPF").write_text(SYNTH_UPF.format(el=el, w=w, r=r))
        files[el] = f"{el}.test.UPF"
    return PseudoTable(str(d), files)


def fcc_cell(a: float = 4.05, el: str = "Al", mass: float = 26.98) -> CellData:
    lat = 0.5 * a * np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]], dtype=float)
    return CellData(lat, [el], np.zeros((1, 3)), {el: mass})


def make_wf(tmp_path: Path, **kw) -> EPhWorkflow:
    install_fake_qe(tmp_path / "bin")
    pseudos = make_pseudos(tmp_path / "pp", {"Al": (30.0, 240.0)})
    s = EPhSettings(pressure_GPa=1.0, prefix="t", k_scf=(4, 4, 4), k_fit=(8, 8, 8), q=(2, 2, 2),
                    k_matdyn=(8, 8, 8), el_ph_nsigma=3, **kw)
    lau = Launcher(kind="plain", bin_dir=str(tmp_path / "bin"),
                   env={"FAKEQE_LOG": str(tmp_path / "calls.log")})
    return EPhWorkflow("SYNTHETIC-fake-Al", fcc_cell(), tmp_path / "work", s, pseudos, lau)


def calls(tmp_path: Path) -> list[str]:
    p = tmp_path / "calls.log"
    return [l.split()[0] for l in p.read_text().splitlines()] if p.exists() else []


# ----------------------------------------------------------------------------- writers / grids
def test_pw_inputs_contain_required_settings():
    cell = fcc_cell()
    txt = qi.write_pw_input(calculation="vc-relax", cell=cell, prefix="al", outdir="./tmp/",
                            pseudo_dir="/pp/", pseudo_files={"Al": "Al.UPF"}, kgrid=(4, 4, 4),
                            ecutwfc=30.0, ecutrho=240.0, press_kbar=1000.0)
    for frag in ("calculation = 'vc-relax'", "&cell", "press = 1000.0", "cell_dynamics = 'bfgs'",
                 "ion_dynamics = 'bfgs'", "ibrav = 0", "CELL_PARAMETERS angstrom",
                 "ATOMIC_POSITIONS crystal", "K_POINTS automatic\n 4 4 4 0 0 0",
                 "ecutrho = 240.0", "Al  26.98000000 Al.UPF"):
        assert frag in txt, frag
    fit = qi.write_pw_input(calculation="scf", cell=cell, prefix="al", outdir="./tmp/",
                            pseudo_dir="/pp/", pseudo_files={"Al": "Al.UPF"}, kgrid=(8, 8, 8),
                            ecutwfc=30.0, la2F=True)
    assert "la2F = .true." in fit and "&cell" not in fit and " 8 8 8 0 0 0" in fit


def test_ph_q2r_matdyn_inputs():
    cell = fcc_cell()
    ph = qi.write_ph_input(prefix="al", outdir="./tmp/", cell=cell, qgrid=(2, 3, 4),
                           el_ph_sigma=0.005, el_ph_nsigma=10, tr2_ph=1e-14, recover=True)
    for frag in ("electron_phonon = 'interpolated'", "el_ph_sigma = 0.005", "el_ph_nsigma = 10",
                 "trans = .true.", "ldisp = .true.", "nq1 = 2", "nq2 = 3", "nq3 = 4",
                 "tr2_ph = 1.0d-14", "recover = .true.", "amass(1) = 26.98", "fildvscf = 'dvscf'"):
        assert frag in ph, frag
    q2r = qi.write_q2r_input(el_ph_nsigma=10)
    assert "la2F = .true." in q2r and "el_ph_nsigma = 10" in q2r
    md = qi.write_matdyn_dos_input(cell=cell, kgrid=(16, 16, 16), ndos=300)
    for frag in ("dos = .true.", "la2F = .true.", "nk1 = 16", "ndos = 300", "flfrq = 'freq.dos'"):
        assert frag in md, frag
    md2 = qi.write_matdyn_dos_input(cell=cell, deltaE_cm1=0.5)
    assert "deltaE = 0.5" in md2 and "ndos" not in md2


def test_fortran_value():
    assert qi.fortran_value(1e-14) == "1.0d-14"
    assert qi.fortran_value(True) == ".true."
    assert qi.fortran_value(0.5) == "0.5"
    with pytest.raises(ValueError):
        qi.fortran_value("a'b")


def test_grid_selection_commensurate():
    cell = fcc_cell()
    g = select_grids(cell.reciprocal_lengths_2pi(), k_spacing=0.2, q_spacing=0.6, fit_factor=3)
    for i in range(3):
        assert g.k_fit[i] % g.k_scf[i] == 0 and g.k_scf[i] % g.q[i] == 0
    assert g.k_fit[0] == 3 * g.k_scf[0]
    with pytest.raises(ValueError, match="dense"):
        EPhGrids((4, 4, 4), (6, 6, 6), (4, 4, 4), (8, 8, 8)).validate()
    with pytest.raises(ValueError, match="q grid"):
        select_grids(cell.reciprocal_lengths_2pi(), k_scf=(4, 4, 4), k_fit=(8, 8, 8), q=(3, 3, 3))


def test_launcher_commands():
    pl = Launcher()
    assert pl.command("pw.x", "scf.in") == ["pw.x", "-in", "scf.in"]
    mp = Launcher(kind="mpirun", nprocs=8, npool=4, mpi_args=["--bind-to", "none"], nice=5,
                  bin_dir="/opt/qe/bin")
    assert mp.command("pw.x", "scf.in") == ["nice", "-n", "5", "mpirun", "-np", "8", "--bind-to",
                                             "none", "/opt/qe/bin/pw.x", "-nk", "4", "-in", "scf.in"]
    assert mp.command("matdyn.x", "m.in") == ["nice", "-n", "5", "/opt/qe/bin/matdyn.x", "-in", "m.in"]
    sr = Launcher.from_config({"kind": "srun", "nprocs": 64, "npool": 8, "mpi_args": "--cpu-bind=cores"})
    assert sr.command("ph.x", "ph.in") == ["srun", "-n", "64", "--cpu-bind=cores", "ph.x", "-nk",
                                           "8", "-in", "ph.in"]
    with pytest.raises(ValueError):
        Launcher(kind="pbs")


def test_pseudo_table(tmp_path):
    t = make_pseudos(tmp_path / "pp", {"Mg": (30.0, 240.0), "B": (35.0, 280.0)})
    prov = t.validate(["Mg", "B"])
    assert prov["Mg"]["functional"] == "PBE" and len(prov["B"]["md5"]) == 32
    assert t.cutoffs(["Mg", "B"]) == (35.0, 280.0)
    with pytest.raises(KeyError):
        t.path("H")
    # SSSP-style json
    js = tmp_path / "sssp.json"
    js.write_text(json.dumps({"Mg": {"filename": "Mg.test.UPF", "cutoff_wfc": 40, "cutoff_rho": 320},
                              "B": {"filename": "B.test.UPF", "cutoff_wfc": 45, "cutoff_rho": 360}}))
    t2 = PseudoTable.from_config({"dir": str(tmp_path / "pp"), "sssp_json": str(js)})
    assert t2.cutoffs(["Mg", "B"]) == (45.0, 360.0)
    # wrong element in header is caught
    (tmp_path / "pp" / "bad.UPF").write_text(SYNTH_UPF.format(el="O", w=1, r=2))
    with pytest.raises(ValueError, match="UPF header"):
        PseudoTable(str(tmp_path / "pp"), {"Mg": "bad.UPF"}).validate(["Mg"])


def test_structure_conversion_dict_and_ase():
    d = {"lattice": np.eye(3) * 3.0, "species": ["Nb", "H"], "frac_coords": [[0, 0, 0], [.5, .5, .5]],
         "masses": {"Nb": 92.9, "H": 1.008}}
    c = to_celldata(d)
    assert c.elements == ["Nb", "H"] and c.volume_A3 == pytest.approx(27.0)
    with pytest.raises(ValueError, match="mass"):
        to_celldata({"lattice": np.eye(3), "species": ["X"], "frac_coords": [[0, 0, 0]]})
    ase = pytest.importorskip("ase.build")
    at = ase.bulk("Al", "fcc", a=4.05)
    c2 = to_celldata(at)
    assert c2.species == ["Al"] and c2.masses_amu["Al"] == pytest.approx(26.98, abs=0.01)
    assert c2.volume_A3 == pytest.approx(4.05 ** 3 / 4)


def test_structure_conversion_pymatgen():
    mg = pytest.importorskip("pymatgen.core")
    s = mg.Structure(mg.Lattice.cubic(3.0), ["Nb", "H"], [[0, 0, 0], [.5, .5, .5]])
    c = to_celldata(s)
    assert c.species == ["Nb", "H"] and c.masses_amu["H"] == pytest.approx(1.008, abs=1e-3)


@pytest.mark.skipif(not REAL.exists(), reason="real QE fixture missing")
def test_parse_real_vc_relax_output():
    """REAL QE 6.7 vc-relax output (smoke test) -> relaxed cell."""
    text = (REAL / "relax.out").read_text()
    c = parse_pw_final_structure(text, fcc_cell(mass=26.9815385))
    assert c.lattice_A[0, 1] == pytest.approx(1.949431454)
    assert c.volume_A3 == pytest.approx(99.9886 * 0.529177210903 ** 3, rel=1e-5)


def test_parse_alat_units_block():
    text = textwrap.dedent("""\
         celldm(1)=   7.500000  celldm(2)=   0.000000  celldm(3)=   0.000000
        Begin final coordinates
        CELL_PARAMETERS (alat=  7.50000000)
          -0.500000000   0.000000000   0.500000000
           0.000000000   0.500000000   0.500000000
          -0.500000000   0.500000000   0.000000000
        ATOMIC_POSITIONS (alat)
        Al            0.0000000000        0.0000000000        0.0000000000
        End final coordinates
        """)
    c = parse_pw_final_structure(text, fcc_cell())
    assert c.volume_A3 == pytest.approx((7.5 * 0.529177210903) ** 3 / 4, rel=1e-9)


# ----------------------------------------------------------------------------- workflow
def test_workflow_full_run_and_resume(tmp_path):
    wf = make_wf(tmp_path)
    res = wf.run()
    assert calls(tmp_path) == ["pw.x", "pw.x", "pw.x", "ph.x", "q2r.x", "matdyn.x"]
    assert len(res.spectra) == 3 and res.checks_ok
    assert res.provenance["grids"]["k_fit"] == [8, 8, 8]
    assert res.spectra[0].provenance["pressure_GPa"] == 1.0
    relax_in = (tmp_path / "work" / "relax.in").read_text()
    assert "press = 10.0" in relax_in and "ecutwfc = 30.0" in relax_in
    assert "la2F = .true." in (tmp_path / "work" / "scf_fit.in").read_text()
    assert all(v == "done" for v in wf.status().values())
    # second run: nothing re-executed
    EPhWorkflow("SYNTHETIC-fake-Al", fcc_cell(), tmp_path / "work", wf.settings, wf.pseudos,
                wf.launcher).run()
    assert len(calls(tmp_path)) == 6
    # changing a ph setting reruns ph and everything downstream only
    s2 = EPhSettings(**{**wf.settings.to_dict(), "tr2_ph": 1e-16})
    EPhWorkflow("SYNTHETIC-fake-Al", fcc_cell(), tmp_path / "work", s2, wf.pseudos, wf.launcher).run()
    assert calls(tmp_path)[6:] == ["ph.x", "q2r.x", "matdyn.x"]
    # changing degauss reruns from scf_fit (relax input also changes -> from relax)
    s3 = EPhSettings(**{**s2.to_dict(), "degauss_Ry": 0.03})
    EPhWorkflow("SYNTHETIC-fake-Al", fcc_cell(), tmp_path / "work", s3, wf.pseudos, wf.launcher).run()
    assert calls(tmp_path)[9:] == ["pw.x", "pw.x", "pw.x", "ph.x", "q2r.x", "matdyn.x"]


def test_workflow_adopts_outputs_without_state(tmp_path):
    wf = make_wf(tmp_path)
    wf.run()
    (tmp_path / "work" / "dfpt_state.json").unlink()
    wf2 = EPhWorkflow("SYNTHETIC-fake-Al", fcc_cell(), tmp_path / "work", wf.settings, wf.pseudos,
                      wf.launcher)
    wf2.run()
    assert len(calls(tmp_path)) == 6


def test_workflow_error_capture_and_recover(tmp_path):
    wf = make_wf(tmp_path, do_relax=False)
    wf.launcher.env["FAKEQE_FAIL"] = "ph.x"
    with pytest.raises(WorkflowError) as ei:
        wf.run()
    msg = str(ei.value)
    assert "step ph failed" in msg and "S matrix not positive definite" in msg and "CRASH" in msg
    assert ei.value.record.returncode == 1
    st = json.loads((tmp_path / "work" / "dfpt_state.json").read_text())
    assert st["steps"]["ph"]["status"] == "failed"
    # rerun after the "fix": completed pw steps are skipped; ph restarts with recover if _ph0 exists
    (tmp_path / "work" / "tmp" / "_ph0").mkdir()
    del wf.launcher.env["FAKEQE_FAIL"]
    res = wf.run()
    assert calls(tmp_path) == ["pw.x", "pw.x", "ph.x", "ph.x", "q2r.x", "matdyn.x"]
    assert "recover = .true." in (tmp_path / "work" / "ph.in").read_text()
    assert len(res.spectra) == 3


def test_workflow_timeout(tmp_path):
    wf = make_wf(tmp_path, do_relax=False, timeouts_s={"scf_fit": 1.0})
    wf.launcher.env["FAKEQE_SLEEP_pw_x"] = "30"
    with pytest.raises(WorkflowError) as ei:
        wf.run()
    assert ei.value.record.timed_out and "timeout" in str(ei.value)
    assert ei.value.record.elapsed_s < 25


def test_workflow_missing_binary(tmp_path):
    wf = make_wf(tmp_path, do_relax=False)
    wf.launcher.bin_dir = str(tmp_path / "nonexistent")
    with pytest.raises(WorkflowError, match="failed to start"):
        wf.run()


def test_workflow_use_smearings_subset(tmp_path):
    wf = make_wf(tmp_path, do_relax=False, use_smearings=[2, 3])
    res = wf.run()
    assert [sp.provenance["smearing_index"] for sp in res.spectra] == [2, 3]
    assert any("explicit config" in w for w in res.warnings)
    with pytest.raises(ValueError):
        make_wf(tmp_path / "b", use_smearings=[7])


def test_settings_from_config():
    s = EPhSettings.from_config({"press_kbar": 1500, "q": [3, 3, 3], "el_ph_nsigma": 5})
    assert s.pressure_GPa == 150.0 and s.q == (3, 3, 3)
    with pytest.raises(ValueError, match="unknown"):
        EPhSettings.from_config({"foo": 1})


def test_workflow_from_config(tmp_path):
    from rtsc.dfpt import workflow_from_config
    install_fake_qe(tmp_path / "bin")
    make_pseudos(tmp_path / "pp", {"Al": (30.0, 240.0)})
    cfg = {"press_kbar": 100.0, "k_scf": [4, 4, 4], "k_fit": [8, 8, 8], "q": [2, 2, 2],
           "k_matdyn": [8, 8, 8], "el_ph_nsigma": 2,
           "pseudo": {"dir": str(tmp_path / "pp"), "files": {"Al": "Al.test.UPF"}},
           "launcher": {"kind": "plain", "bin_dir": str(tmp_path / "bin")}}
    wf = workflow_from_config("SYNTHETIC-cfg", fcc_cell(), tmp_path / "w", cfg)
    assert wf.settings.pressure_GPa == 10.0
    res = wf.run()
    assert len(res.spectra) == 2 and res.ok
    json.dumps(res.to_dict())          # fully JSON-serialisable record
