"""Thermodynamic (convex-hull) existence filter.

Value used, in order of preference:
 1. A convex hull computed from reference entries at the candidate's pressure (enthalpy per atom
    when the pressure is > 0, energy at 0 GPa): pymatgen ``PhaseDiagram`` if installed, otherwise
    an exact linear program (scipy ``linprog``).  Reference entries are the run's own candidates
    that report an enthalpy at the same pressure (``include_run_candidates``) plus user-supplied
    entries (``reference_file``: CSV/JSON with formula, enthalpy_eV_per_atom, pressure_GPa).
    Elemental references for every element are required, otherwise the hull is not evaluated.
 2. A stored ``energy_above_hull`` property (e.g. from Materials Project), but only when it refers
    to the candidate's target pressure (MP values refer to 0 GPa).

Energies from different codes/functionals must not be mixed in one hull; this is the user's
responsibility and recorded as an assumption.
"""
from __future__ import annotations

import csv
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from ..sources.base import Candidate
from ..sources.composition import atomic_fractions, parse_formula
from .base import FilterResult

CRITERION = "X-HULL"


@dataclass
class HullEntry:
    composition: dict[str, float]
    energy_per_atom: float              # enthalpy per atom at pressure_GPa (eV/atom)
    pressure_GPa: float = 0.0
    label: str = ""


def load_reference_entries(path: str | Path) -> list[HullEntry]:
    p = Path(path)
    rows: list[dict[str, Any]]
    if p.suffix == ".csv":
        with open(p, newline="") as fh:
            rows = list(csv.DictReader(fh))
    else:
        data = json.loads(p.read_text())
        rows = data["entries"] if isinstance(data, dict) else data
    out = []
    for r in rows:
        e = r.get("enthalpy_eV_per_atom", r.get("energy_per_atom"))
        out.append(HullEntry(parse_formula(str(r["formula"])), float(e),
                             float(r.get("pressure_GPa") or 0.0), str(r.get("label", r["formula"]))))
    return out


def _lp_hull_energy(target: dict[str, float], entries: list[HullEntry]) -> float | None:
    """min sum_i l_i e_i  s.t. sum_i l_i x_i = x_target, l >= 0 (fractions sum to one)."""
    import numpy as np
    from scipy.optimize import linprog
    els = sorted(set(target))
    xs = atomic_fractions(target)
    usable = [e for e in entries if set(e.composition) <= set(els)]
    if not usable:
        return None
    A = np.array([[atomic_fractions(e.composition).get(el, 0.0) for e in usable] for el in els])
    b = np.array([xs[el] for el in els])
    c = np.array([e.energy_per_atom for e in usable])
    res = linprog(c, A_eq=A, b_eq=b, bounds=[(0, None)] * len(usable), method="highs")
    if not res.success:
        return None
    return float(res.fun)


def _pmg_hull_energy(target: dict[str, float], entries: list[HullEntry]) -> float | None:
    from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram
    from pymatgen.core import Composition
    els = set(target)
    usable = [e for e in entries if set(e.composition) <= els]
    pd_entries = [PDEntry(Composition(e.composition),
                          e.energy_per_atom * sum(e.composition.values())) for e in usable]
    try:
        pd = PhaseDiagram(pd_entries)
    except Exception:  # noqa: BLE001 - e.g. missing terminal (elemental) entries
        return None
    comp = Composition(target)
    try:
        # energy of the hull at this composition (per atom)
        return float(pd.get_hull_energy_per_atom(comp))
    except Exception:  # noqa: BLE001
        return None


def e_above_hull(target: dict[str, float], energy_per_atom: float, entries: list[HullEntry],
                 backend: str = "auto") -> tuple[float | None, str]:
    """Energy above the lower convex hull of ``entries`` (+ the target itself), eV/atom."""
    have = {el for e in entries for el in e.composition if len(e.composition) == 1}
    if not set(target) <= have:
        return None, "missing elemental reference(s): " + ",".join(sorted(set(target) - have))
    allent = list(entries) + [HullEntry(dict(target), energy_per_atom, label="<candidate>")]
    hull = None
    used = ""
    if backend in ("auto", "pymatgen"):
        try:
            hull = _pmg_hull_energy(target, allent)
            used = "pymatgen.PhaseDiagram"
        except ImportError:
            if backend == "pymatgen":
                raise
    if hull is None and backend in ("auto", "scipy"):
        hull = _lp_hull_energy(target, allent)
        used = "scipy.linprog convex hull"
    if hull is None:
        return None, "hull construction failed"
    return max(energy_per_atom - hull, 0.0), used


@dataclass
class HullFilter:
    enabled: bool = True
    max_e_above_hull_eV_per_atom: float = 0.1
    pressure_tol_GPa: float = 1.0
    include_run_candidates: bool = True
    reference_file: str | None = None
    backend: str = "auto"
    use_stored_value: bool = True
    _reference: list[HullEntry] = field(default_factory=list, repr=False)

    def prepare(self, candidates: list[Candidate]) -> None:
        ref: list[HullEntry] = []
        if self.reference_file:
            ref.extend(load_reference_entries(self.reference_file))
        if self.include_run_candidates:
            for c in candidates:
                if c.enthalpy_eV_per_atom is not None and c.composition and not c.synthetic:
                    try:
                        e = HullEntry(dict(c.composition), float(c.enthalpy_eV_per_atom),
                                      float(c.pressure_GPa or 0.0), c.material_id)
                    except (TypeError, ValueError):
                        continue        # malformed: that candidate fails in evaluate() -> ERROR
                    ref.append(e)
        self._reference = ref

    def comparison(self) -> str:
        return f"e_above_hull <= {self.max_e_above_hull_eV_per_atom} eV/atom"

    def evaluate(self, c: Candidate) -> FilterResult:
        thr = self.max_e_above_hull_eV_per_atom
        target_p = float(c.properties.get("target_pressure_GPa", c.pressure_GPa or 0.0))
        res = FilterResult(CRITERION, None, None, thr, self.comparison())
        # 1. computed hull at the candidate's pressure
        if c.enthalpy_eV_per_atom is not None and c.composition and self._reference is not None:
            p = float(c.pressure_GPa or 0.0)
            if abs(p - target_p) <= self.pressure_tol_GPa:
                ents = [e for e in self._reference if abs(e.pressure_GPa - p) <= self.pressure_tol_GPa
                        and e.label != c.material_id]
                val, how = e_above_hull(c.composition, float(c.enthalpy_eV_per_atom), ents,
                                        self.backend)
                if val is not None:
                    res.value = val
                    res.passed = val <= thr
                    res.data_source = f"{how} from {len(ents)} reference entries at {p} GPa"
                    res.details = {"pressure_GPa": p, "n_reference": len(ents),
                                   "enthalpy_eV_per_atom": c.enthalpy_eV_per_atom}
                    return res
                res.reason = how
        # 2. stored value
        if self.use_stored_value and c.properties.get("energy_above_hull") is not None:
            at_p = c.properties.get("energy_above_hull_pressure_GPa")
            if at_p is None:
                at_p = 0.0 if c.source == "materials_project" else (c.pressure_GPa or 0.0)
            if abs(float(at_p) - target_p) <= self.pressure_tol_GPa:
                val = float(c.properties["energy_above_hull"])
                res.value = val
                res.passed = val <= thr
                res.data_source = f"stored energy_above_hull ({c.source}, {at_p} GPa)"
                return res
            res.reason = (f"stored energy_above_hull refers to {at_p} GPa, target {target_p} GPa; "
                          "not used")
        if not res.reason:
            res.reason = "no enthalpy/hull data"
        return res
