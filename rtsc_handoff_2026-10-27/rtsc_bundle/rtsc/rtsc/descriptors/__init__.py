"""Heuristic descriptors: priority scores, budget planning, optional ML pre-checks.

CONTRACT S2: everything in this package is heuristic.  It may only REORDER candidates or select
some to be NOT COMPUTED (the pipeline records those as DEFERRED with the reason).  This package
deliberately does not import ``rtsc.status``; a unit test (tests/test_pipeline_s2.py) enforces
that no heuristic input can produce REJECTED_PROVEN or EXCLUDED_ASSUMPTION.
"""
from .budget import BudgetPlan, plan_budget
from .ml import get_ml_hook, ml_precheck
from .score import PriorityWeights, Score, components, priority

__all__ = ["BudgetPlan", "plan_budget", "PriorityWeights", "Score", "components", "priority",
           "get_ml_hook", "ml_precheck"]
