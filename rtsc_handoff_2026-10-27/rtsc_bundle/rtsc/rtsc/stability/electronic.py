"""Metallicity and magnetism existence filters (EXCLUDED_ASSUMPTION, S3).

Metallicity: phonon-mediated superconductivity in the Eliashberg model needs a Fermi surface.
  Value: ``band_gap`` (eV) if known, else ``is_metal``.  Fails if band_gap > max_band_gap_eV, or
  (no gap known and) is_metal is False.  DFT (semi)local gaps are usually underestimates, so a
  DFT gap > 0 is a strong (non-proved) indication of an insulator at the computed pressure.
Magnetism: conventional s-wave pairing is suppressed by magnetic order.
  Value: ``ordering`` (MP convention NM/FM/AFM/FiM/Unknown) and |magnetisation| per atom
  (``magnetization_per_atom`` or ``total_magnetization`` / number of sites, in mu_B).
"""
from __future__ import annotations

from dataclasses import dataclass, field

from ..sources.base import Candidate
from .base import FilterResult


@dataclass
class MetallicityFilter:
    enabled: bool = True
    max_band_gap_eV: float = 0.0

    def evaluate(self, c: Candidate) -> FilterResult:
        thr = self.max_band_gap_eV
        comp = f"band_gap <= {thr} eV (or is_metal == True when no gap is known)"
        gap = c.properties.get("band_gap")
        if gap is not None:
            g = float(gap)
            return FilterResult("X-METAL", g <= thr, g, thr, comp,
                                data_source=f"band_gap from {c.source}")
        im = c.properties.get("is_metal")
        if im is not None:
            return FilterResult("X-METAL", bool(im), bool(im), True, comp,
                                data_source=f"is_metal from {c.source}")
        return FilterResult("X-METAL", None, None, thr, comp, reason="no band-gap/is_metal data")


@dataclass
class MagnetismFilter:
    enabled: bool = True
    max_abs_magnetization_muB_per_atom: float = 0.05
    exclude_orderings: tuple[str, ...] = field(default=("FM", "AFM", "FiM"))

    def evaluate(self, c: Candidate) -> FilterResult:
        thr = self.max_abs_magnetization_muB_per_atom
        comp = (f"ordering not in {list(self.exclude_orderings)} and "
                f"|m| <= {thr} mu_B/atom")
        order = c.properties.get("ordering")
        m_pa = c.properties.get("magnetization_per_atom")
        if m_pa is None and c.properties.get("total_magnetization") is not None:
            n = c.properties.get("nsites") or c.n_sites
            if n:
                m_pa = float(c.properties["total_magnetization"]) / float(n)
        if order is None and m_pa is None:
            return FilterResult("X-MAGNETIC", None, None, thr, comp,
                                reason="no magnetic ordering / magnetisation data")
        value = {"ordering": order, "abs_m_per_atom": None if m_pa is None else abs(float(m_pa))}
        bad_order = order is not None and str(order) in self.exclude_orderings
        bad_m = m_pa is not None and abs(float(m_pa)) > thr
        return FilterResult("X-MAGNETIC", not (bad_order or bad_m), value, thr, comp,
                            data_source=f"ordering/magnetisation from {c.source}")
