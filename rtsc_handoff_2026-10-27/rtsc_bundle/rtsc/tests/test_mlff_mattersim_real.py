"""REAL MatterSim smoke test (MatterSim-v1.0.0-1M on CPU).  Skipped unless ``mattersim`` is
importable; the checkpoint is downloaded by mattersim from raw.githubusercontent.com on first use.
Tolerances are loose sanity bounds, not validation."""
from __future__ import annotations

import pytest

pytest.importorskip("mattersim")
pytest.importorskip("phonopy")
from ase.build import bulk  # noqa: E402

from rtsc.discovery.interfaces import GeneratedStructure  # noqa: E402
from rtsc.mlff._atoms import atoms_to_structure  # noqa: E402
from rtsc.mlff.mattersim import MatterSimConfig, MatterSimPrescreener  # noqa: E402


def test_mattersim_al_relax_and_phonons(tmp_path):
    cfg = MatterSimConfig.from_dict({"device": "cpu", "cache_dir": str(tmp_path),
                                     "phonon": {"supercell": [2, 2, 2], "mesh": [6, 6, 6]}})
    ps = MatterSimPrescreener(cfg)
    g = GeneratedStructure(atoms_to_structure(bulk("Al", "fcc", a=4.1)), "smoke", "r", {"fingerprint": "al"})
    r0, r10 = ps.prescreen([g], 0.0)[0], ps.prescreen([g], 10.0)[0]
    for r in (r0, r10):
        assert r.ok and r.dynamically_stable_ml is True and r.e_above_hull_ml_eV == pytest.approx(0, abs=1e-6)
    vol = lambda r: abs(__import__("numpy").linalg.det(r.relaxed_structure["lattice"]))  # noqa: E731
    a0 = (4 * vol(r0)) ** (1 / 3)
    assert 3.95 < a0 < 4.12 and vol(r10) < vol(r0)
    f = r0.phonon_dos["freq_cm1"]
    top = max(x for x, d in zip(f, r0.phonon_dos["dos"]) if d > 1e-6 * max(r0.phonon_dos["dos"]))
    assert 250 < top < 400          # Al phonons end near 9-10 THz (~300-330 cm^-1)
