"""Common record for existence filters and the EXCLUDED_ASSUMPTION verdict they produce (S3).

A filter returns a ``FilterResult`` with ``passed``:
  * ``True``  - evaluated, criterion satisfied;
  * ``False`` - evaluated, criterion violated -> the material is EXCLUDED_ASSUMPTION;
  * ``None``  - not evaluated (no applicable data); the material proceeds and the gap is recorded.
Each result records the exact criterion, the value, the threshold and the data source.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any

from ..status import Status, Verdict


@dataclass
class FilterResult:
    criterion: str              # rule id, e.g. "X-HULL"
    passed: bool | None
    value: Any = None
    threshold: Any = None
    comparison: str = ""        # the exact test, e.g. "e_above_hull <= 0.1 eV/atom"
    reason: str = ""
    data_source: str = ""
    details: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    def summary(self) -> str:
        state = {True: "passed", False: "FAILED", None: "not evaluated"}[self.passed]
        v = f"value={self.value!r}" if self.value is not None else "value=n/a"
        return f"{self.criterion}: {self.comparison} [{v}; {state}; source={self.data_source or 'n/a'}]" \
            + (f" -- {self.reason}" if self.reason else "")


def exclusion_verdict(material_id: str, results: list[FilterResult]) -> Verdict:
    """EXCLUDED_ASSUMPTION verdict from the failing filter results (S3: criterion recorded)."""
    failing = [r for r in results if r.passed is False]
    if not failing:
        raise ValueError("exclusion_verdict called without a failing filter")
    return Verdict(
        material_id=material_id,
        status=Status.EXCLUDED_ASSUMPTION,
        T_low_K=None, T_up_K=None,
        reasons=["excluded by physical (non-proved) existence criterion: " + r.summary()
                 for r in failing],
        rule_ids=[r.criterion for r in failing],
        assumptions=[f"{r.criterion}: {r.comparison}" for r in results if r.passed is not None],
        flags=[],
        certificate={"kind": "existence_filters", "proved": False,
                     "filter_results": [r.to_dict() for r in results]},
    )
