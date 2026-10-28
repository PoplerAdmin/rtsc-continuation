"""Dynamic (harmonic phonon) stability filter, from FIRST-PRINCIPLES phonon data only.

Imaginary frequencies are represented as negative numbers (phonopy / QE matdyn convention).
The material is EXCLUDED_ASSUMPTION when the lowest frequency is below ``-imaginary_tol``
(after optionally ignoring the three acoustic branches within ``gamma_q_tol`` of Gamma, where
small negative values are numerical noise).

Accepted data (``candidate.properties['phonon']`` or the ``phonon`` argument):
  {"min_frequency": -0.8, "unit": "THz"}
  {"frequencies": [[...], ...], "qpoints": [[qx,qy,qz], ...], "unit": "cm-1"}
  {"band_yaml": "path/to/band.yaml"}  or {"mesh_yaml": "..."}   (phonopy; needs PyYAML)
ML-potential pre-checks are NOT accepted here: they are heuristic (S2) and live in descriptors.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

from ..sources.base import Candidate
from ..spectrum import unit_factor
from .base import FilterResult

CRITERION = "X-DYNAMIC"


def _load_phonopy_yaml(path: str) -> dict[str, Any]:
    import yaml  # optional
    with open(path) as fh:
        d = yaml.safe_load(fh)
    freqs, qs = [], []
    for ph in d.get("phonon", []):
        qs.append(ph.get("q-position", [0, 0, 0]))
        freqs.append([b["frequency"] for b in ph.get("band", [])])
    return {"frequencies": freqs, "qpoints": qs, "unit": "THz"}   # phonopy default unit


def lowest_frequency(phonon: dict[str, Any], ignore_acoustic_gamma: bool = True,
                     gamma_q_tol: float = 0.05) -> tuple[float, str]:
    """Return (lowest frequency, unit) after the acoustic-at-Gamma exclusion."""
    if "band_yaml" in phonon or "mesh_yaml" in phonon:
        phonon = _load_phonopy_yaml(phonon.get("band_yaml") or phonon["mesh_yaml"])
    unit = phonon.get("unit", "THz")
    if "min_frequency" in phonon and "frequencies" not in phonon:
        return float(phonon["min_frequency"]), unit
    freqs = phonon["frequencies"]
    qs = phonon.get("qpoints")
    lo = math.inf
    for k, row in enumerate(freqs):
        vals = sorted(float(x) for x in row)
        if ignore_acoustic_gamma and qs is not None:
            q = qs[k]
            if math.sqrt(sum(float(x) ** 2 for x in q)) <= gamma_q_tol:
                vals = vals[3:]
        if vals:
            lo = min(lo, vals[0])
    return lo, unit


@dataclass
class DynamicStabilityFilter:
    enabled: bool = True
    imaginary_tol: float = 0.25
    imaginary_tol_unit: str = "THz"
    ignore_acoustic_gamma: bool = True
    gamma_q_tol: float = 0.05

    def comparison(self) -> str:
        return (f"min phonon frequency >= -{self.imaginary_tol} {self.imaginary_tol_unit}"
                + (" (acoustic modes at Gamma ignored)" if self.ignore_acoustic_gamma else ""))

    def evaluate(self, c: Candidate, phonon: dict[str, Any] | None = None) -> FilterResult:
        data = phonon if phonon is not None else c.properties.get("phonon")
        thr = -abs(self.imaginary_tol)
        if not data:
            return FilterResult(CRITERION, None, None, thr, self.comparison(),
                                reason="no first-principles phonon data")
        try:
            lo, unit = lowest_frequency(data, self.ignore_acoustic_gamma, self.gamma_q_tol)
        except ImportError as exc:
            return FilterResult(CRITERION, None, None, thr, self.comparison(),
                                reason=f"cannot read phonon data ({exc})")
        lo_tolunit = lo * unit_factor(unit) / unit_factor(self.imaginary_tol_unit)
        return FilterResult(CRITERION, lo_tolunit >= thr, lo_tolunit, thr, self.comparison(),
                            data_source=str(data.get("source", c.source)),
                            details={"lowest_frequency": lo, "unit": unit})
