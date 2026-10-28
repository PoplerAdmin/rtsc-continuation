"""Ordering and budget planning (stage 4/5).  ORDERING / DEFERRAL ONLY (contract S2).

``plan_budget`` sorts candidates by heuristic priority (descending, ties by material id) and
returns which ones to compute now and which to defer, each deferral with a reason.  It returns
plain ids and strings: no Status is created here; the pipeline converts deferrals to DEFERRED.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from .score import Score


@dataclass
class BudgetPlan:
    order: list[str]                                   # ids to compute, highest priority first
    deferred: dict[str, str] = field(default_factory=dict)   # id -> reason
    scores: dict[str, Score] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return {"order": list(self.order), "deferred": dict(self.deferred),
                "scores": {k: v.to_dict() for k, v in self.scores.items()}}


def plan_budget(scores: list[Score], *, max_materials: int | None = None,
                min_score: float | None = None,
                ml_results: dict[str, dict[str, Any]] | None = None,
                defer_if_ml_unstable: bool = False) -> BudgetPlan:
    ranked = sorted(scores, key=lambda s: (-s.total, s.material_id))
    plan = BudgetPlan(order=[], scores={s.material_id: s for s in scores})
    for s in ranked:
        ml = (ml_results or {}).get(s.material_id)
        if defer_if_ml_unstable and ml and ml.get("ok") and ml.get("dynamically_stable") is False:
            plan.deferred[s.material_id] = (
                f"heuristic ML pre-check ({ml.get('model', '?')}) suggested imaginary modes "
                f"(min {ml.get('min_frequency_THz')} THz); deferred, NOT excluded (S2)")
            continue
        if min_score is not None and s.total < min_score:
            plan.deferred[s.material_id] = (f"heuristic priority {s.total:.4g} < min_score "
                                            f"{min_score}; not computed (S2)")
            continue
        if max_materials is not None and len(plan.order) >= max_materials:
            plan.deferred[s.material_id] = (f"budget: max_materials={max_materials} reached at "
                                            f"priority {s.total:.4g}; not computed (S2)")
            continue
        plan.order.append(s.material_id)
    return plan
