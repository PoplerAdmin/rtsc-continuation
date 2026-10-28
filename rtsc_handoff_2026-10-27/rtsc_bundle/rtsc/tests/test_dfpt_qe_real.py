"""End-to-end REAL Quantum ESPRESSO smoke test (fcc Al, deliberately unconverged).

Runs only when RTSC_RUN_QE=1 and pw.x, ph.x, q2r.x, matdyn.x are on PATH (or in RTSC_QE_BIN).
Environment:
  RTSC_QE_PSEUDO_DIR  directory containing Al.pz-vbc.UPF (default /usr/share/espresso/pseudo)
  RTSC_QE_NP          MPI ranks (default 2; 1 = no mpirun)
  RTSC_QE_SHIM=1      LD_PRELOAD the snprintf_chk shim (needed for Ubuntu 24.04's QE 6.7 package)

The numbers produced are REAL program output from a smoke-test setup (15 Ry, 4x4x4 / 8x8x8 k,
2x2x2 q) and are NOT physical results for aluminium.
"""
from __future__ import annotations

import os
import shutil
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from rtsc.dfpt import (EPhSettings, EPhWorkflow, Launcher, PseudoTable, build_fortify_shim,  # noqa: E402
                       collect_qe_eph, spectrum_lambda)

REAL = Path(__file__).parent / "data_dfpt" / "qe67_al_smoke"
BIN = os.environ.get("RTSC_QE_BIN")
PSEUDO = os.environ.get("RTSC_QE_PSEUDO_DIR", "/usr/share/espresso/pseudo")


def _have_qe() -> bool:
    for c in ("pw.x", "ph.x", "q2r.x", "matdyn.x"):
        if not (Path(BIN, c).exists() if BIN else shutil.which(c)):
            return False
    return (Path(PSEUDO) / "Al.pz-vbc.UPF").exists()


pytestmark = pytest.mark.skipif(os.environ.get("RTSC_RUN_QE") != "1" or not _have_qe(),
                                reason="set RTSC_RUN_QE=1 with QE binaries + Al.pz-vbc.UPF")


def test_al_smoke_end_to_end(tmp_path):
    ase = pytest.importorskip("ase.build")
    atoms = ase.bulk("Al", "fcc", a=4.05)
    env = {}
    if os.environ.get("RTSC_QE_SHIM") == "1":
        env["LD_PRELOAD"] = str(build_fortify_shim(tmp_path / "shim"))
    np_ = int(os.environ.get("RTSC_QE_NP", "2"))
    lau = Launcher(kind="mpirun" if np_ > 1 else "plain", nprocs=np_, npool=np_, nice=10,
                   bin_dir=BIN, env=env,
                   mpi_args=["--allow-run-as-root", "--oversubscribe", "-x", "LD_PRELOAD",
                             "-x", "OMP_NUM_THREADS"] if np_ > 1 else [])
    s = EPhSettings(pressure_GPa=5.0, prefix="al", ecutwfc_Ry=15.0, smearing="mv",
                    degauss_Ry=0.05, k_scf=(4, 4, 4), k_fit=(8, 8, 8), q=(2, 2, 2),
                    k_matdyn=(8, 8, 8), el_ph_sigma_Ry=0.01, el_ph_nsigma=10, tr2_ph=1e-12,
                    ndos=50, conv_thr=1e-9, timeouts_s={k: 1800 for k in
                                                        ("relax", "scf_fit", "scf", "ph", "q2r", "matdyn")})
    wf = EPhWorkflow("Al-fcc-REAL-smoke", atoms, tmp_path / "work", s,
                     PseudoTable(PSEUDO, {"Al": "Al.pz-vbc.UPF"}), lau)
    res = wf.run()
    assert len(res.spectra) == 10 and res.checks_ok, res.flags
    assert res.provenance["code_version"].startswith("MATDYN")
    assert all(min(sp.omega_K) > 0 for sp in res.spectra)
    # regression vs the stored REAL fixture produced by the same setup
    if REAL.exists():
        ref = collect_qe_eph(REAL, "fixture")
        for a, b in zip(res.spectra, ref.spectra):
            assert spectrum_lambda(a) == pytest.approx(spectrum_lambda(b), rel=2e-2)
