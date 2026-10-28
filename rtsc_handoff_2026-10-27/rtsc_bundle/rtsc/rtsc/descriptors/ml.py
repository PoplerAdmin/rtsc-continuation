"""Optional ML-potential hooks (MACE, CHGNet): relax + quick phonon pre-check.  HEURISTIC (S2).

Results are used only to reorder candidates and, if ``defer_if_ml_unstable`` is set, to DEFER
them.  They never exclude or reject: an ML potential's imaginary mode is not evidence that meets
the bar of an EXCLUDED_ASSUMPTION criterion (which needs first-principles phonons).

Requires ``ase`` and the model package (``mace-torch`` or ``chgnet``); phonopy is used for the
phonon check when available, else ASE's finite-difference ``Phonons``.  Install with
``pip install 'rtsc[ml]'``.  Every function returns a plain dict:

    {"ok": bool, "model": str, "dynamically_stable": bool | None, "min_frequency_THz": float | None,
     "energy_eV_per_atom": float | None, "pressure_GPa": float, "error": str | None}
"""
from __future__ import annotations

import math
import tempfile
from typing import Any, Callable

from ..sources.base import Candidate

EV_PER_A3_PER_GPA = 1.0 / 160.21766208


def _atoms(c: Candidate):
    from ase import Atoms
    if c.lattice is None or not c.species or c.frac_coords is None:
        raise ValueError("candidate has no ordered structure")
    return Atoms(symbols=c.species, cell=c.lattice, scaled_positions=c.frac_coords, pbc=True)


def _calculator(model: str, **kw):
    if model == "mace":
        from mace.calculators import mace_mp
        return mace_mp(model=kw.get("model_path", "medium"), default_dtype="float64",
                       device=kw.get("device", "cpu"))
    if model == "chgnet":
        from chgnet.model.dynamics import CHGNetCalculator
        return CHGNetCalculator()
    raise ValueError(f"unknown ML model {model!r} (mace|chgnet)")


def _relax(atoms, pressure_GPa: float, fmax: float, steps: int):
    from ase.optimize import BFGS
    try:
        from ase.filters import FrechetCellFilter as CellFilter
    except ImportError:  # older ASE
        from ase.constraints import ExpCellFilter as CellFilter  # type: ignore
    f = CellFilter(atoms, scalar_pressure=pressure_GPa * EV_PER_A3_PER_GPA)
    BFGS(f, logfile=None).run(fmax=fmax, steps=steps)
    return atoms


def _min_freq_phonopy(atoms, calc, supercell: int, mesh: int) -> float:
    from phonopy import Phonopy
    from phonopy.structure.atoms import PhonopyAtoms
    unit = PhonopyAtoms(symbols=atoms.get_chemical_symbols(), cell=atoms.cell[:],
                        scaled_positions=atoms.get_scaled_positions())
    ph = Phonopy(unit, supercell_matrix=[[supercell, 0, 0], [0, supercell, 0], [0, 0, supercell]])
    ph.generate_displacements(distance=0.01)
    forces = []
    from ase import Atoms
    for sc in ph.supercells_with_displacements:
        a = Atoms(sc.symbols, cell=sc.cell, scaled_positions=sc.scaled_positions, pbc=True)
        a.calc = calc
        forces.append(a.get_forces())
    ph.forces = forces
    ph.produce_force_constants()
    ph.run_mesh([mesh] * 3, with_eigenvectors=False, is_gamma_center=True)
    d = ph.get_mesh_dict()
    freqs, qs = d["frequencies"], d["qpoints"]
    lo = math.inf
    for q, row in zip(qs, freqs):
        vals = sorted(row)
        if math.sqrt(sum(x * x for x in q)) < 1e-8:
            vals = vals[3:]                      # acoustic at Gamma
        if len(vals):
            lo = min(lo, float(vals[0]))
    return lo                                    # THz


def _min_freq_ase(atoms, calc, supercell: int, mesh: int) -> float:
    import numpy as np
    from ase.phonons import Phonons
    with tempfile.TemporaryDirectory() as tmp:
        ph = Phonons(atoms, calc, supercell=(supercell,) * 3, delta=0.01, name=f"{tmp}/ph")
        ph.run()
        ph.read(acoustic=True)
        qs = [(i / mesh, j / mesh, k / mesh) for i in range(mesh) for j in range(mesh)
              for k in range(mesh)]
        e = ph.band_structure(np.array(qs))       # eV, imaginary returned as negative
    thz_per_ev = 241.79892623
    lo = math.inf
    for q, row in zip(qs, e):
        vals = sorted(row)
        if q == (0.0, 0.0, 0.0):
            vals = vals[3:]
        if len(vals):
            lo = min(lo, float(vals[0]) * thz_per_ev)
    return lo


def ml_precheck(c: Candidate, model: str = "mace", *, relax: bool = True, fmax: float = 0.05,
                steps: int = 300, supercell: int = 2, mesh: int = 4,
                imaginary_tol_THz: float = 0.3, **kw) -> dict[str, Any]:
    p = float(c.properties.get("target_pressure_GPa", c.pressure_GPa or 0.0))
    out: dict[str, Any] = {"ok": False, "model": model, "dynamically_stable": None,
                           "min_frequency_THz": None, "energy_eV_per_atom": None,
                           "pressure_GPa": p, "error": None, "heuristic": True}
    try:
        atoms = _atoms(c)
        calc = _calculator(model, **kw)
        atoms.calc = calc
        if relax:
            _relax(atoms, p, fmax, steps)
        out["energy_eV_per_atom"] = float(atoms.get_potential_energy()) / len(atoms)
        try:
            lo = _min_freq_phonopy(atoms, calc, supercell, mesh)
        except ImportError:
            lo = _min_freq_ase(atoms, calc, supercell, mesh)
        out["min_frequency_THz"] = lo
        out["dynamically_stable"] = lo >= -abs(imaginary_tol_THz)
        out["ok"] = True
    except Exception as exc:  # noqa: BLE001 - heuristic stage: failure only loses ordering info
        out["error"] = f"{type(exc).__name__}: {exc}"
    return out


def get_ml_hook(cfg: dict[str, Any] | None) -> Callable[[Candidate], dict[str, Any]] | None:
    """Build the pre-check callable from the ``[ml]`` config table (None if disabled)."""
    cfg = dict(cfg or {})
    if not cfg.get("enabled", False):
        return None
    model = cfg.pop("model", "mace")
    cfg.pop("enabled", None)
    cfg.pop("defer_if_ml_unstable", None)

    def hook(c: Candidate) -> dict[str, Any]:
        return ml_precheck(c, model, **cfg)
    return hook
