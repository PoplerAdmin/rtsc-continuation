"""ML phonons (finite displacements via phonopy, forces from an ASE calculator such as MatterSim).

HEURISTIC ONLY (ARCH S2).  The result may reorder candidates or mark them DEFERRED with a recorded
reason; it never yields EXCLUDED_ASSUMPTION (the first-principles dynamic-stability filter in
``rtsc.stability.dynamic`` explicitly refuses ML phonons).

phonopy API used (checked against the installed phonopy, see README "ML adapters"):
``phonopy.Phonopy(unitcell, supercell_matrix=...)``, ``phonopy.structure.atoms.PhonopyAtoms``,
``generate_displacements(distance=)``, ``supercells_with_displacements``, ``forces = ...``,
``produce_force_constants()``, ``symmetrize_force_constants()``,
``run_mesh(mesh, is_gamma_center=True, with_eigenvectors=, is_mesh_symmetry=)`` then ``.mesh``
(qpoints, weights, frequencies), ``run_total_dos()`` then ``.total_dos`` (frequency_points, dos),
``run_projected_dos()`` then ``.projected_dos`` (frequency_points, projected_dos); checked with
phonopy 4.8.6, where the older ``get_*_dict()`` accessors are deprecated.
Frequencies come out in THz for eV/angstrom forces (phonopy's default unit conversion).
"""
from __future__ import annotations

import math
from dataclasses import asdict, dataclass, field
from typing import Any, Sequence

THZ_TO_CM1 = 33.35640951981521          # 1 THz / c  in cm^-1  (c = 299792458 m/s)
CM1_TO_K = 1.4387768775039338           # h c / k_B in K cm

DEFAULT_LIGHT_ELEMENTS = ("H", "Li", "Be", "B", "C", "N", "O", "F")


@dataclass
class PhononSettings:
    supercell: Sequence[int] | Sequence[Sequence[int]] | str = "auto"
    min_supercell_length: float = 7.0       # angstrom, used when supercell == "auto"
    max_supercell_atoms: int = 160          # cap for "auto" (smaller supercell => noisier, recorded)
    displacement: float = 0.01              # angstrom
    mesh: Sequence[int] = (8, 8, 8)
    imaginary_tol_cm1: float = 15.0         # min frequency below -tol => dynamically_stable_ml False
    projected_dos: bool = True              # needed for the light-element heuristic
    light_elements: Sequence[str] = DEFAULT_LIGHT_ELEMENTS
    hf_threshold_cm1: float = 500.0         # "high-frequency" for the heuristic score
    dos_bins: int = 200


@dataclass
class PhononSummary:
    min_freq_cm1: float
    min_freq_qpoint: list[float]
    dynamically_stable: bool
    imaginary_mode_fraction: float          # fraction of mesh modes below -tol (weighted by q weight)
    supercell_matrix: list[list[int]]
    n_supercell_atoms: int
    n_displacements: int
    freq_cm1: list[float]                   # total-DOS frequency points
    dos: list[float]                        # total DOS (states / cm^-1 per unit cell, normalised to 3N)
    light_hf_fraction: float | None         # HEURISTIC score (see light_element_hf_fraction)
    omega_log_dos_K: float | None           # log-average of the PHONON DOS (not alpha^2F), heuristic
    notes: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def auto_supercell(cell, min_length: float, max_atoms: int, n_atoms: int) -> list[list[int]]:
    """Diagonal supercell with every lattice-vector length >= min_length, shrunk to fit max_atoms."""
    import numpy as np
    lengths = np.linalg.norm(np.asarray(cell, dtype=float), axis=1)
    reps = [max(1, int(math.ceil(min_length / L))) for L in lengths]
    while n_atoms * reps[0] * reps[1] * reps[2] > max_atoms and max(reps) > 1:
        k = int(np.argmax([r * L if r > 1 else -1 for r, L in zip(reps, lengths)]))
        reps[k] -= 1
    return [[reps[0], 0, 0], [0, reps[1], 0], [0, 0, reps[2]]]


def _supercell_matrix(settings: PhononSettings, atoms) -> list[list[int]]:
    sc = settings.supercell
    if isinstance(sc, str):
        if sc != "auto":
            raise ValueError(f"supercell must be 'auto', 3 ints or a 3x3 matrix, got {sc!r}")
        return auto_supercell(atoms.cell[:], settings.min_supercell_length,
                              settings.max_supercell_atoms, len(atoms))
    sc = list(sc)
    if len(sc) == 3 and all(isinstance(x, (int, float)) for x in sc):
        return [[int(sc[0]), 0, 0], [0, int(sc[1]), 0], [0, 0, int(sc[2])]]
    return [[int(x) for x in row] for row in sc]


def light_element_hf_fraction(freq_cm1, projected_dos, symbols, light_elements, threshold_cm1) -> float | None:
    """HEURISTIC: fraction of the phonon DOS weight that sits on light-element atoms at
    frequency >= threshold.  Only used for ordering (S2).  ``projected_dos`` has shape
    (n_atoms, n_freq) as returned by phonopy."""
    import numpy as np
    f = np.asarray(freq_cm1, dtype=float)
    p = np.asarray(projected_dos, dtype=float)
    total = float(np.clip(p, 0, None).sum())
    if total <= 0:
        return None
    light = np.array([s in set(light_elements) for s in symbols])
    if not light.any():
        return 0.0
    sel = np.clip(p[light][:, f >= threshold_cm1], 0, None).sum()
    return float(sel / total)


def omega_log_of_dos(freq_cm1, dos) -> float | None:
    """exp(<ln w>) over the positive part of the phonon DOS, in kelvin (HEURISTIC; not alpha^2F)."""
    import numpy as np
    f = np.asarray(freq_cm1, dtype=float)
    g = np.asarray(dos, dtype=float)
    m = (f > 1.0) & (g > 0)
    if not m.any():
        return None
    w = g[m] / f[m]                         # alpha^2F/w-like weighting with alpha^2 = const
    return float(math.exp((w * np.log(f[m])).sum() / w.sum()) * CM1_TO_K)


def ml_phonons(atoms, calculator, settings: PhononSettings | None = None) -> PhononSummary:
    """Finite-displacement harmonic phonons of ``atoms`` (already relaxed) using ``calculator``."""
    import numpy as np
    from phonopy import Phonopy
    from phonopy.structure.atoms import PhonopyAtoms

    st = settings or PhononSettings()
    notes: list[str] = []
    unit = PhonopyAtoms(symbols=atoms.get_chemical_symbols(), cell=np.array(atoms.cell[:]),
                        scaled_positions=atoms.get_scaled_positions())
    scm = _supercell_matrix(st, atoms)
    ph = Phonopy(unit, supercell_matrix=scm)
    ph.generate_displacements(distance=st.displacement)
    scs = ph.supercells_with_displacements
    forces = []
    from ase import Atoms
    for sc in scs:
        a = Atoms(symbols=sc.symbols, cell=sc.cell, scaled_positions=sc.scaled_positions, pbc=True)
        a.calc = calculator
        f = np.array(a.get_forces(), dtype=float)
        forces.append(f - f.mean(axis=0))   # remove net drift (standard practice)
    ph.forces = forces
    ph.produce_force_constants()
    ph.symmetrize_force_constants()
    n_sc = len(ph.supercell)
    lengths = np.linalg.norm(np.array(ph.supercell.cell), axis=1)
    if lengths.min() < st.min_supercell_length:
        why = "capped by max_supercell_atoms" if st.supercell == "auto" else "explicit supercell"
        notes.append(f"supercell shortest length {lengths.min():.2f} A < {st.min_supercell_length} A "
                     f"({why}); instabilities at q not commensurate with it may be missed")

    want_proj = bool(st.projected_dos)
    ph.run_mesh(list(st.mesh), is_gamma_center=True, with_eigenvectors=want_proj,
                is_mesh_symmetry=not want_proj)
    mesh = ph.mesh
    qs = np.asarray(mesh.qpoints, dtype=float)
    wq = np.asarray(mesh.weights, dtype=float)
    fr = np.asarray(mesh.frequencies, dtype=float) * THZ_TO_CM1      # (nq, nbands)

    # ignore the three acoustic modes at Gamma (zero by the acoustic sum rule)
    fmin, qmin, imag_w = math.inf, [0.0, 0.0, 0.0], 0.0
    for q, w, row in zip(qs, wq, fr):
        vals = np.sort(row)
        if np.linalg.norm(q) < 1e-8:
            vals = vals[3:]
        if vals.size and vals[0] < fmin:
            fmin, qmin = float(vals[0]), [float(x) for x in q]
        imag_w += w * float((vals < -abs(st.imaginary_tol_cm1)).sum())
    imag_frac = imag_w / float(wq.sum() * fr.shape[1])

    ph.run_total_dos()
    td = ph.total_dos
    f_cm1 = np.asarray(td.frequency_points, dtype=float) * THZ_TO_CM1
    dos = np.asarray(td.dos, dtype=float) / THZ_TO_CM1                 # per cm^-1

    light = None
    if want_proj:
        ph.run_projected_dos()
        pd = ph.projected_dos
        light = light_element_hf_fraction(np.asarray(pd.frequency_points) * THZ_TO_CM1,
                                          pd.projected_dos, list(ph.primitive.symbols),
                                          st.light_elements, st.hf_threshold_cm1)
    stable = fmin >= -abs(st.imaginary_tol_cm1)
    return PhononSummary(
        min_freq_cm1=float(fmin), min_freq_qpoint=qmin, dynamically_stable=bool(stable),
        imaginary_mode_fraction=float(imag_frac), supercell_matrix=scm, n_supercell_atoms=int(n_sc),
        n_displacements=len(scs), freq_cm1=[float(x) for x in f_cm1], dos=[float(x) for x in dos],
        light_hf_fraction=light, omega_log_dos_K=omega_log_of_dos(f_cm1, dos), notes=notes)
