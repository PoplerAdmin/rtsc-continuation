"""Physical existence filters (ARCH section 6, stage 3; contract S3).

Every filter yields EXCLUDED_ASSUMPTION (never REJECTED_PROVEN) with the exact criterion, value and
threshold, and each can be disabled in the config (``[filters.<name>] enabled = false``).

    hull       X-HULL      e_above_hull <= max_e_above_hull_eV_per_atom
    metallic   X-METAL     band_gap <= max_band_gap_eV  (or is_metal)
    magnetic   X-MAGNETIC  ordering not in exclude_orderings and |m| <= threshold
    dynamic    X-DYNAMIC   min phonon frequency >= -imaginary_tol (first-principles data only)
"""
from __future__ import annotations

from dataclasses import dataclass, field, fields
from typing import Any

from ..sources.base import Candidate
from ..status import Verdict
from .base import FilterResult, exclusion_verdict
from .dynamic import DynamicStabilityFilter
from .electronic import MagnetismFilter, MetallicityFilter
from .hull import HullEntry, HullFilter, e_above_hull, load_reference_entries

__all__ = ["FilterResult", "FilterSet", "HullFilter", "MetallicityFilter", "MagnetismFilter",
           "DynamicStabilityFilter", "HullEntry", "e_above_hull", "load_reference_entries",
           "exclusion_verdict"]


def _build(cls, cfg: dict[str, Any] | None):
    cfg = dict(cfg or {})
    known = {f.name for f in fields(cls) if not f.name.startswith("_")}
    unknown = set(cfg) - known
    if unknown:
        raise ValueError(f"unknown keys for {cls.__name__}: {sorted(unknown)}")
    if "exclude_orderings" in cfg:
        cfg["exclude_orderings"] = tuple(cfg["exclude_orderings"])
    return cls(**cfg)


@dataclass
class FilterSet:
    hull: HullFilter = field(default_factory=HullFilter)
    metallic: MetallicityFilter = field(default_factory=MetallicityFilter)
    magnetic: MagnetismFilter = field(default_factory=MagnetismFilter)
    dynamic: DynamicStabilityFilter = field(default_factory=DynamicStabilityFilter)

    @classmethod
    def from_config(cls, cfg: dict[str, Any] | None) -> "FilterSet":
        cfg = dict(cfg or {})
        unknown = set(cfg) - {"hull", "metallic", "magnetic", "dynamic"}
        if unknown:
            raise ValueError(f"unknown [filters] sections: {sorted(unknown)}")
        return cls(_build(HullFilter, cfg.get("hull")), _build(MetallicityFilter, cfg.get("metallic")),
                   _build(MagnetismFilter, cfg.get("magnetic")),
                   _build(DynamicStabilityFilter, cfg.get("dynamic")))

    def prepare(self, candidates: list[Candidate]) -> None:
        if self.hull.enabled:
            self.hull.prepare(candidates)

    def evaluate_static(self, c: Candidate) -> tuple[Verdict | None, list[FilterResult]]:
        """Stage 3: hull, metallicity, magnetism, and dynamic stability if phonon data is given."""
        results = []
        for f in (self.hull, self.metallic, self.magnetic, self.dynamic):
            if f.enabled:
                results.append(f.evaluate(c))
        return _verdict(c, results), results

    def evaluate_dynamic(self, c: Candidate, phonon: dict[str, Any]) -> tuple[Verdict | None, FilterResult | None]:
        """After DFPT: dynamic stability from the computed phonons (if enabled)."""
        if not self.dynamic.enabled:
            return None, None
        r = self.dynamic.evaluate(c, phonon)
        return _verdict(c, [r]), r


def _verdict(c: Candidate, results: list[FilterResult]) -> Verdict | None:
    if any(r.passed is False for r in results):
        return exclusion_verdict(c.material_id, results)
    return None
