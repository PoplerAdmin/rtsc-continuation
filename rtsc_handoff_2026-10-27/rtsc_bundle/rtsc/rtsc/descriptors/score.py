"""Heuristic priority scores.  ORDERING ONLY (contract S2).

Nothing in this package can produce a Status: scores only reorder candidates, and the budget
planner (``budget.py``) can only mark candidates as not computed (DEFERRED by the pipeline).
The components are crude physical proxies for "worth computing first":

  hydrogen      atomic fraction of H                                     (in [0, 1])
  light         atomic fraction of light elements (default H Li Be B C N) (in [0, 1])
  phonon_scale  sqrt(m_H / M_eff), M_eff = harmonic mean atomic mass       (in (0, 1])
                (a harmonic-oscillator estimate: omega ~ 1/sqrt(M) at fixed force constants)
  dos           N(E_F) / (N(E_F) + dos_scale), states/eV/atom if provided  (in [0, 1))
  spectrum      max over attached spectra of A / (A + a0 T*)  (if spectra are already known)
  ml            1 if an ML pre-check found no imaginary modes, 0 if it did, absent otherwise

Missing components are omitted (not zero-filled) and the total is the weighted mean of those
present, so a candidate is not punished for missing optional data.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any

from ..sources.base import Candidate
from ..sources.composition import ATOMIC_MASS, atomic_fractions

DEFAULT_LIGHT = ("H", "Li", "Be", "B", "C", "N")


@dataclass
class PriorityWeights:
    hydrogen: float = 1.0
    light: float = 0.5
    phonon_scale: float = 1.0
    dos: float = 0.5
    spectrum: float = 1.0
    ml: float = 0.5
    light_elements: tuple[str, ...] = DEFAULT_LIGHT
    dos_scale: float = 0.5          # states/eV/atom giving dos component 0.5
    T_target_K: float = 300.0       # only used by the 'spectrum' component
    a0: float = 4.33148145

    @classmethod
    def from_config(cls, cfg: dict[str, Any] | None) -> "PriorityWeights":
        cfg = dict(cfg or {})
        if "light_elements" in cfg:
            cfg["light_elements"] = tuple(cfg["light_elements"])
        known = set(cls.__dataclass_fields__)
        bad = set(cfg) - known
        if bad:
            raise ValueError(f"unknown [priority] keys: {sorted(bad)}")
        return cls(**cfg)


@dataclass
class Score:
    material_id: str
    total: float
    components: dict[str, float] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return {"material_id": self.material_id, "total": self.total,
                "components": dict(self.components), "use": "ordering only (S2)"}


def _spectrum_A(d: dict[str, Any]) -> float:
    return math.fsum(float(w) for w in d.get("weight_K", []))


def components(c: Candidate, w: PriorityWeights, ml: dict[str, Any] | None = None) -> dict[str, float]:
    comp: dict[str, float] = {}
    fr = atomic_fractions(c.composition)
    if fr:
        comp["hydrogen"] = fr.get("H", 0.0)
        comp["light"] = sum(v for e, v in fr.items() if e in w.light_elements)
        masses = [(v, ATOMIC_MASS.get(e)) for e, v in fr.items()]
        if all(m is not None for _, m in masses):
            m_eff = 1.0 / sum(v / m for v, m in masses)       # harmonic mean
            comp["phonon_scale"] = math.sqrt(ATOMIC_MASS["H"] / m_eff)
    dos = c.properties.get("dos_at_fermi")
    if dos is not None:
        x = max(float(dos), 0.0)
        comp["dos"] = x / (x + w.dos_scale)
    if c.spectra:
        a_star = w.a0 * w.T_target_K
        comp["spectrum"] = max(_spectrum_A(s) / (_spectrum_A(s) + a_star) for s in c.spectra)
    if ml is not None and ml.get("ok") and ml.get("dynamically_stable") is not None:
        comp["ml"] = 1.0 if ml["dynamically_stable"] else 0.0
    return comp


def priority(c: Candidate, w: PriorityWeights | None = None, ml: dict[str, Any] | None = None) -> Score:
    w = w or PriorityWeights()
    comp = components(c, w, ml)
    num = den = 0.0
    for k, v in comp.items():
        wk = float(getattr(w, k))
        num += wk * v
        den += wk
    total = num / den if den > 0 else 0.0
    return Score(c.material_id, total, comp)
