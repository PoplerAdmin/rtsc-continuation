"""ORCH: contract S2 -- heuristics (priority scores, budgets, ML pre-checks) can only reorder or
DEFER.  They can never yield REJECTED_PROVEN or EXCLUDED_ASSUMPTION (or CERTIFIED_CANDIDATE)."""
from __future__ import annotations

import ast
import random
from pathlib import Path

import pytest

import rtsc.descriptors as D
import rtsc.pipeline as P
from rtsc.sources import Candidate
from rtsc.sources.synthetic import gaussian_peak_spectrum
from rtsc.status import Status

import test_pipeline_stubs as stubs

HEURISTIC_STATUS = {Status.DEFERRED}


def test_descriptors_cannot_name_a_status():
    """Static: the heuristic package never imports rtsc.status nor mentions a Status member."""
    root = Path(D.__file__).parent
    for f in root.glob("*.py"):
        tree = ast.parse(f.read_text())
        for node in ast.walk(tree):
            if isinstance(node, (ast.Import, ast.ImportFrom)):
                mod = getattr(node, "module", None) or ""
                names = [a.name for a in node.names]
                assert "status" not in mod and not any("status" in n for n in names), f
            if isinstance(node, ast.Name):
                assert node.id not in ("Status", "Verdict"), f
            if isinstance(node, ast.Attribute):
                assert node.attr not in {s.name for s in Status}, f


def _cands(rng: random.Random, n: int) -> list[Candidate]:
    out = []
    for i in range(n):
        A = rng.choice([200.0, 900.0, 1500.0, 3000.0])
        w = rng.choice([500.0, 2000.0, 8000.0])
        comp = rng.choice([{"La": 1, "H": 10}, {"Nb": 1}, {"Mg": 1, "B": 2}, {}])
        out.append(Candidate(f"SYNTHETIC-s2-{i}", comp, source="test",
                             spectra=[gaussian_peak_spectrum(w, w / 20, A, n=6)],
                             properties={"dos_at_fermi": rng.random() * 3,
                                         "is_metal": True, "energy_above_hull": 0.0},
                             synthetic=True))
    return out


def _cfg(tmp_path, budget=None, min_score=None, weights=None, defer_ml=False):
    raw = {"run": {"workdir": str(tmp_path), "resume": False},
           "filters": {k: {"enabled": False} for k in ("hull", "metallic", "magnetic", "dynamic")},
           "budget": {k: v for k, v in (("max_materials", budget), ("min_score", min_score))
                      if v is not None},
           "priority": weights or {},
           "ml": {"enabled": False, "defer_if_ml_unstable": defer_ml}}
    return P.normalize_config(raw)


@pytest.mark.parametrize("seed", range(6))
def test_heuristics_only_defer(tmp_path, monkeypatch, seed):
    """Random scores, budgets, weights and ML outcomes: every non-computed material is DEFERRED,
    every computed one has exactly the status certify.classify gives on its own."""
    stubs.install(monkeypatch, force=True)
    rng = random.Random(seed)
    cands = _cands(rng, 12)
    ml = {c.material_id: {"ok": rng.random() < 0.8, "model": "fake",
                          "dynamically_stable": rng.random() < 0.5, "min_frequency_THz": -1.0}
          for c in cands}
    weights = {k: rng.random() * 3 for k in ("hydrogen", "light", "phonon_scale", "dos", "spectrum", "ml")}
    cfg = _cfg(tmp_path / f"r{seed}", budget=rng.randint(0, 12),
               min_score=rng.choice([None, 0.2, 0.5]), weights=weights, defer_ml=rng.random() < 0.5)
    summ = P.run(cfg, cands, ml_hook=lambda c: ml[c.material_id], write_reports=False)
    recs = {m: P.read_record(cfg["run"]["workdir"], m) for m in summ.inputs}
    for m, v in summ.verdicts.items():
        origin = recs[m]["verdict_origin"]
        if origin == P.CLASSIFY_ORIGIN:
            spectra = [P.Spectrum.from_dict(d) for d in recs[m]["spectra"]]
            assert v.status == stubs.stub_classify(spectra, cfg["certify"]).status
        else:
            assert origin in ("pipeline.budget", "pipeline.ml"), origin
            assert v.status in HEURISTIC_STATUS
        assert v.status is not Status.EXCLUDED_ASSUMPTION


def test_unlimited_budget_statuses_independent_of_heuristics(tmp_path, monkeypatch):
    """With no budget limit, changing every heuristic weight and ML result changes only the
    ORDER, never a status."""
    stubs.install(monkeypatch, force=True)
    cands0 = _cands(random.Random(1), 10)
    results = []
    for k in range(4):
        rng = random.Random(100 + k)
        cands = [Candidate.from_dict(c.to_dict()) for c in cands0]
        weights = {w: rng.random() for w in ("hydrogen", "light", "phonon_scale", "dos", "spectrum", "ml")}
        ml = {c.material_id: {"ok": True, "dynamically_stable": rng.random() < 0.5} for c in cands}
        cfg = _cfg(tmp_path / f"u{k}", weights=weights)
        summ = P.run(cfg, cands, ml_hook=lambda c: ml[c.material_id], write_reports=False)
        results.append({m: v.status for m, v in summ.verdicts.items()})
    assert all(r == results[0] for r in results)


def test_budget_plan_returns_no_status():
    s = [D.Score("a", 0.9), D.Score("b", 0.1), D.Score("c", 0.5)]
    plan = D.plan_budget(s, max_materials=2)
    assert plan.order == ["a", "c"] and set(plan.deferred) == {"b"}
    assert all(isinstance(r, str) for r in plan.deferred.values())


def test_pipeline_guard_rejects_heuristic_exclusion():
    with pytest.raises(AssertionError):
        P._pipeline_verdict("budget", "x", Status.EXCLUDED_ASSUMPTION, ["no"])
    with pytest.raises(AssertionError):
        P._pipeline_verdict("ml", "x", Status.REJECTED_PROVEN, ["no"])
