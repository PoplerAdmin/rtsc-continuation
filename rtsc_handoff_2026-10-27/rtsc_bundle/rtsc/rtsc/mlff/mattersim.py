"""MatterSim ML pre-screen adapter (``Prescreener``): relax at pressure, ML hull, ML phonons.

HEURISTIC ONLY (ARCH S2, DISCOVERY.md 2.3).  ``PrescreenResult`` values may only reorder candidates
or let the loop mark them DEFERRED with a recorded reason.  This module has no way to set a Status.

Verified against the official sources (github.com/microsoft/mattersim @ 5c84771, 2026-10-21;
PyPI ``mattersim`` 1.2.5, requires Python >= 3.12, MIT licence):

* ASE calculator: ``from mattersim.forcefield import MatterSimCalculator``;
  ``MatterSimCalculator(load_path=..., device="cuda"|"cpu", dtype="float32"|"float64", ...)``;
  implemented properties energy, free_energy, forces, stress (stress in eV/A^3, ASE convention).
* Checkpoints: ``MatterSim-v1.0.0-1M`` (default, faster) and ``MatterSim-v1.0.0-5M`` (more
  accurate), M3GNet architecture.  ``Potential.from_checkpoint`` looks in
  ``~/.local/mattersim/pretrained_models`` and otherwise downloads from
  ``https://raw.githubusercontent.com/microsoft/mattersim/main/pretrained_models/``.  Any other
  ``load_path`` is treated as a file path (e.g. a fine-tuned model).
* Stated range (MODEL_CARD.md): 0-5000 K, pressures up to 1000 GPa; bulk materials only.
* The repo's own relaxer (``mattersim.applications.relax.Relaxer``) wraps ASE optimisers with
  ``ExpCellFilter``/``FrechetCellFilter`` and ``scalar_pressure`` in eV/A^3 (= P[GPa] * ase.units.GPa).
  We call ASE directly with the same filters so the pressure handling is explicit here.

The relaxation minimises the enthalpy H = E + P V at the request's pressure (cell + positions);
energies above hull use H per atom of reference phases relaxed with the same model at the same P.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Callable, Sequence

from ..discovery.interfaces import GeneratedStructure, PrescreenResult
from ._atoms import atoms_to_structure, fallback_fingerprint, structure_to_atoms
from .hull import MLHull, mp_reference_fetcher
from .phonons import PhononSettings, ml_phonons

MATTERSIM_MODELS = ("MatterSim-v1.0.0-1M", "MatterSim-v1.0.0-5M")
MATTERSIM_MAX_PRESSURE_GPA = 1000.0          # MODEL_CARD.md
HEURISTIC_NOTE = "HEURISTIC (ARCH S2): may reorder / DEFER only; never EXCLUDED or proven"


@dataclass
class RelaxSettings:
    fmax: float = 0.02                       # eV/A (forces and filter's virial-derived stress term)
    max_steps: int = 400
    optimizer: str = "FIRE"                  # FIRE | BFGS | LBFGS (ase.optimize)
    cell_filter: str = "FrechetCellFilter"   # FrechetCellFilter | ExpCellFilter (ase.filters)
    constrain_symmetry: bool = False         # ase.constraints.FixSymmetry


@dataclass
class MatterSimConfig:
    model: str = "MatterSim-v1.0.0-1M"       # or a path to a .pth checkpoint
    device: str = "auto"                     # auto | cpu | cuda
    dtype: str = "float32"
    relax: RelaxSettings = field(default_factory=RelaxSettings)
    # hull
    hull: bool = True
    reference_dir: str | None = None
    use_mp: bool = False                     # needs MP_API_KEY and mp-api
    mp_e_hull_max: float | None = 0.05
    cache_dir: str | None = None
    elemental_fallback: bool = True
    # phonons
    phonons: bool = True
    phonon: PhononSettings = field(default_factory=PhononSettings)
    skip_phonons_above_hull_eV: float | None = None   # cost saver; recorded in notes when used
    # thresholds only used to write "ML_SUGGEST_DEFER" notes; the loop owns the decision
    suggest_defer_e_hull_eV: float | None = 0.2
    suggest_defer_unstable: bool = True

    @classmethod
    def from_dict(cls, d: dict[str, Any] | None) -> "MatterSimConfig":
        d = dict(d or {})
        relax = RelaxSettings(**d.pop("relax", {}))
        phon = d.pop("phonon", {})
        if "mesh" in phon:
            phon["mesh"] = tuple(phon["mesh"])
        return cls(relax=relax, phonon=PhononSettings(**phon), **d)


def _device(dev: str) -> str:
    if dev != "auto":
        return dev
    try:
        import torch
        return "cuda" if torch.cuda.is_available() else "cpu"
    except ImportError:
        return "cpu"


def default_calculator_factory(cfg: MatterSimConfig) -> Callable[[], Any]:
    def make():
        from mattersim.forcefield import MatterSimCalculator
        return MatterSimCalculator(load_path=cfg.model, device=_device(cfg.device), dtype=cfg.dtype)
    return make


def relax_at_pressure(atoms, calculator, pressure_GPa: float, settings: RelaxSettings | None = None):
    """Relax cell + positions at hydrostatic ``pressure_GPa``.

    Returns (relaxed atoms, enthalpy eV/atom, info).  ``info`` holds energy per atom, volume per
    atom, the residual pressure from the stress tensor, steps and the convergence flag."""
    import ase.filters as filters
    import ase.optimize as optimize
    import numpy as np
    from ase.units import GPa

    st = settings or RelaxSettings()
    a = atoms.copy()
    a.calc = calculator
    if st.constrain_symmetry:
        from ase.constraints import FixSymmetry
        a.set_constraint(FixSymmetry(a))
    filt_cls = getattr(filters, st.cell_filter)
    p_ev = float(pressure_GPa) * GPa                         # eV/A^3
    flt = filt_cls(a, scalar_pressure=p_ev)
    opt = getattr(optimize, st.optimizer)(flt, logfile=None)
    ret = opt.run(fmax=st.fmax, steps=st.max_steps)
    converged = bool(ret) if isinstance(ret, (bool, np.bool_)) else bool(opt.converged())
    e = float(a.get_potential_energy())
    v = float(a.get_volume())
    n = len(a)
    stress = np.asarray(a.get_stress(voigt=True), dtype=float)   # eV/A^3, ASE sign convention
    p_resid = -float(stress[:3].mean()) / GPa                    # GPa; should approach pressure_GPa
    a.set_constraint()
    info = {"energy_eV_per_atom": e / n, "volume_A3_per_atom": v / n, "pressure_from_stress_GPa": p_resid,
            "steps": int(opt.get_number_of_steps()), "converged": converged,
            "max_force_eV_A": float(np.linalg.norm(a.get_forces(), axis=1).max())}
    return a, (e + p_ev * v) / n, info


class MatterSimPrescreener:
    """``rtsc.discovery.interfaces.Prescreener`` backed by MatterSim (or any ASE calculator).

    ``calculator_factory`` exists for testing / swapping potentials; whatever it returns, the
    reference phases are relaxed with the same calculator, and ``model_id`` keys the hull cache."""

    name = "mattersim"

    def __init__(self, config: MatterSimConfig | dict | None = None, *,
                 calculator_factory: Callable[[], Any] | None = None, model_id: str | None = None,
                 mp_fetch=None):
        self.cfg = config if isinstance(config, MatterSimConfig) else MatterSimConfig.from_dict(config)
        self._factory = calculator_factory or default_calculator_factory(self.cfg)
        self._calc = None
        self.model_id = model_id or self.cfg.model
        if mp_fetch is None and self.cfg.use_mp:
            mp_fetch = mp_reference_fetcher(self.cfg.mp_e_hull_max)
        self.hull = MLHull(relax_fn=self._relax, model_id=self.model_id,
                           reference_dir=self.cfg.reference_dir, mp_fetch=mp_fetch,
                           cache_dir=self.cfg.cache_dir, elemental_fallback=self.cfg.elemental_fallback)

    # ------------------------------------------------------------------ helpers
    @property
    def calculator(self):
        if self._calc is None:
            self._calc = self._factory()
        return self._calc

    def _relax(self, atoms, pressure_GPa: float):
        return relax_at_pressure(atoms, self.calculator, pressure_GPa, self.cfg.relax)

    def provenance(self) -> dict[str, Any]:
        ver = None
        try:
            from importlib.metadata import version
            ver = version("mattersim")
        except Exception:  # noqa: BLE001
            pass
        return {"prescreener": self.name, "model": self.model_id, "mattersim_version": ver,
                "relax": vars(self.cfg.relax), "phonon": {k: (list(v) if isinstance(v, tuple) else v)
                                                         for k, v in vars(self.cfg.phonon).items()}}

    @staticmethod
    def heuristic_score(result: PrescreenResult) -> float | None:
        """Optional cheap ordering score: fraction of the ML phonon DOS carried by light elements
        above ``hf_threshold_cm1`` (HEURISTIC; S2: ordering only).  None when unavailable."""
        if result.phonon_dos is None:
            return None
        return result.phonon_dos.get("light_hf_fraction")

    # ------------------------------------------------------------------ protocol
    def prescreen(self, structures: Sequence[GeneratedStructure], pressure_GPa: float) -> list[PrescreenResult]:
        return [self.prescreen_one(g, pressure_GPa) for g in structures]

    def prescreen_one(self, g: GeneratedStructure, pressure_GPa: float) -> PrescreenResult:
        fp = (g.provenance or {}).get("fingerprint") or fallback_fingerprint(g.structure)
        notes = [HEURISTIC_NOTE, f"model={self.model_id}", f"pressure_GPa={pressure_GPa:g}"]
        if pressure_GPa > MATTERSIM_MAX_PRESSURE_GPA:
            notes.append(f"pressure {pressure_GPa} GPa outside MatterSim stated range (<= 1000 GPa)")
        res = PrescreenResult(fingerprint=fp, notes=notes)
        try:
            atoms = structure_to_atoms(g.structure)
            ra, h, info = self._relax(atoms, pressure_GPa)
        except Exception as exc:  # noqa: BLE001 - ML failure => ok=False, recorded; never a verdict
            res.ok = False
            notes.append(f"relaxation failed: {type(exc).__name__}: {exc}")
            return res
        res.relaxed_structure = atoms_to_structure(ra, pressure_GPa=float(pressure_GPa),
                                                   enthalpy_eV_per_atom=float(h))
        notes.append("relax: " + ", ".join(f"{k}={v:.4g}" if isinstance(v, float) else f"{k}={v}"
                                           for k, v in info.items()))
        if not info["converged"]:
            notes.append(f"relaxation NOT converged in {self.cfg.relax.max_steps} steps")

        if self.cfg.hull:
            try:
                ehull, how, labels = self.hull.e_above_hull(ra.get_chemical_symbols(), h, pressure_GPa)
                res.e_above_hull_ml_eV = ehull
                notes.append(f"ml_hull: {how}; references={labels}")
            except Exception as exc:  # noqa: BLE001
                notes.append(f"ml_hull failed: {type(exc).__name__}: {exc}")
        e = res.e_above_hull_ml_eV
        skip = (self.cfg.skip_phonons_above_hull_eV is not None and e is not None
                and e > self.cfg.skip_phonons_above_hull_eV)
        if self.cfg.phonons and skip:
            notes.append(f"ml_phonons skipped: e_above_hull_ml {e:.3f} > "
                         f"{self.cfg.skip_phonons_above_hull_eV} eV/atom (cost saving)")
        elif self.cfg.phonons:
            try:
                ps = ml_phonons(ra, self.calculator, self.cfg.phonon)
                res.dynamically_stable_ml = ps.dynamically_stable
                res.min_freq_cm1 = ps.min_freq_cm1
                res.phonon_dos = {"freq_cm1": ps.freq_cm1, "dos": ps.dos,
                                  "light_hf_fraction": ps.light_hf_fraction,
                                  "omega_log_dos_K": ps.omega_log_dos_K,
                                  "imaginary_mode_fraction": ps.imaginary_mode_fraction,
                                  "min_freq_qpoint": ps.min_freq_qpoint,
                                  "supercell_matrix": ps.supercell_matrix,
                                  "mesh": list(self.cfg.phonon.mesh),
                                  "imaginary_tol_cm1": self.cfg.phonon.imaginary_tol_cm1}
                notes.append(f"ml_phonons: min {ps.min_freq_cm1:.1f} cm^-1 at q={ps.min_freq_qpoint}, "
                             f"supercell {ps.supercell_matrix}, {ps.n_displacements} displacements")
                notes.extend(ps.notes)
            except Exception as exc:  # noqa: BLE001
                notes.append(f"ml_phonons failed: {type(exc).__name__}: {exc}")
        thr = self.cfg.suggest_defer_e_hull_eV
        if thr is not None and e is not None and e > thr:
            notes.append(f"ML_SUGGEST_DEFER: e_above_hull_ml {e:.3f} eV/atom > {thr} eV/atom")
        if self.cfg.suggest_defer_unstable and res.dynamically_stable_ml is False:
            notes.append(f"ML_SUGGEST_DEFER: ML imaginary modes (min {res.min_freq_cm1:.1f} cm^-1 < "
                         f"-{self.cfg.phonon.imaginary_tol_cm1} cm^-1)")
        if res.min_freq_cm1 is not None and not math.isfinite(res.min_freq_cm1):
            res.min_freq_cm1 = None
        return res


__all__ = ["MatterSimConfig", "RelaxSettings", "MatterSimPrescreener", "relax_at_pressure",
           "default_calculator_factory", "MATTERSIM_MODELS"]
