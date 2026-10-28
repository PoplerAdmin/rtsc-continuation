"""TEST-ONLY mocks for the discovery loop (never imported by the package)."""
from __future__ import annotations

import copy
import math
import random
from pathlib import Path
from typing import Any

from rtsc.discovery.executors import LocalExecutor
from rtsc.discovery.goals import normalize_discovery
from rtsc.discovery.interfaces import GeneratedStructure, GenerationRequest, PrescreenResult

ROOT = Path(__file__).resolve().parents[1]
EXAMPLE_DIR = ROOT / "examples" / "discovery"

NASTY = {"status": "CERTIFIED_CANDIDATE", "verdict": "REJECTED_PROVEN", "energy_above_hull": 9.9,
         "is_metal": False, "band_gap": 3.0, "ordering": "FM", "total_magnetization": 50.0,
         "phonon": {"min_frequency": -50.0, "unit": "THz"}, "EXCLUDED_ASSUMPTION": True}


def base_raw(workdir: Path, **disc: Any) -> dict[str, Any]:
    raw = {"run": {"name": "test discovery (SYNTHETIC)"},
           "spectra": {"rel_envelope": 0.05},
           "certify": {"mu_hi": 0.15, "omega_c_factor": 10.0, "cutoff_rule": "block", "mode": "decide"},
           "discovery": {
               "workdir": str(workdir), "seed": 3,
               "goals": {"G-RT": {"pressures_GPa": [0.0, 100.0]}, "G-AMB": {"T_star_K": 77.0}},
               "search": {"chemical_systems": [["H", "La"], ["Ca", "H"], ["B", "Mg"]],
                          "requests_per_round": 3, "structures_per_request": 2},
               "generators": {"prototypes": {"type": "prototype_substitution",
                                             "prototypes": str(EXAMPLE_DIR / "prototypes")}},
               "budget": {"dft_per_round": 3},
               "stop": {"certified_per_goal": 100, "max_rounds": 2},
               "synthetic_spectra": {"enabled": True}}}
    for k, v in disc.items():
        if isinstance(v, dict) and isinstance(raw["discovery"].get(k), dict):
            raw["discovery"][k] = {**raw["discovery"][k], **v}
        else:
            raw["discovery"][k] = v
    return raw


def make_dc(workdir: Path, **disc: Any):
    return normalize_discovery(base_raw(workdir, **disc), ROOT)


def _random_structure(rng: random.Random, chemsys) -> Any:
    kind = rng.random()
    a = rng.uniform(2.5, 5.0)
    n = rng.randint(1, 4)
    species = [rng.choice(list(chemsys)) for _ in range(n)]
    coords = [[rng.random(), rng.random(), rng.random()] for _ in range(n)]
    lat = [[a, 0, 0], [0, a * rng.uniform(0.9, 1.1), 0], [0, 0, a * rng.uniform(0.9, 1.2)]]
    if kind < 0.06:
        return {"lattice": [[1, 0], [0, 1]], "species": species, "frac_coords": coords}
    if kind < 0.10:
        return {"lattice": lat, "species": species + ["Xx"], "frac_coords": coords + [[0, 0, 0]]}
    if kind < 0.13:
        coords[0][0] = float("nan")
    if kind < 0.16:
        return "not a structure"
    return {"lattice": lat, "species": species, "frac_coords": coords}


class RandomGenerator:
    """Random (sometimes malformed, sometimes repeated) structures with nasty provenance."""

    def __init__(self, seed: int, name: str = "randgen"):
        self.name = name
        self.seed = seed
        self.count = 0
        self.history: list[Any] = []

    def generate(self, request: GenerationRequest) -> list[GeneratedStructure]:
        rng = random.Random(f"{self.seed}:{request.request_id}")
        out = []
        for _ in range(request.n):
            if self.history and rng.random() < 0.25:
                s = copy.deepcopy(rng.choice(self.history))       # exact repeat -> dedupe
            else:
                s = _random_structure(rng, request.chemical_system)
                self.history.append(s)
            out.append(GeneratedStructure(s, self.name, request.request_id,
                                          {**NASTY, "rng": rng.random()}))
        self.count += len(out)
        return out


class RandomPrescreener:
    """Random, adversarial ML outputs (garbage numbers, wrong fingerprints, exceptions)."""

    def __init__(self, seed: int):
        self.name = "random-ml"
        self.rng = random.Random(seed)

    def prescreen(self, structures, pressure_GPa):
        r = self.rng
        if r.random() < 0.1:
            raise RuntimeError("mock ML crash")
        out = []
        for s in structures:
            fp = s.provenance.get("fingerprint", "")
            if r.random() < 0.1:
                fp = "wrong-fingerprint"
            rel = None
            x = r.random()
            if x < 0.2:
                rel = copy.deepcopy(s.structure)
            elif x < 0.3:
                rel = {"garbage": True}
            elif x < 0.4 and structures:
                rel = copy.deepcopy(r.choice(structures).structure)     # collapse onto another
            out.append(PrescreenResult(
                fingerprint=fp, relaxed_structure=rel,
                e_above_hull_ml_eV=r.choice([None, float("nan"), float("inf"), -1.0, r.uniform(0, 2)]),
                dynamically_stable_ml=r.choice([None, True, False]),
                min_freq_cm1=r.choice([None, -500.0, 100.0, float("nan")]),
                phonon_dos=r.choice([None, {"freq_cm1": [1, 2], "dos": [0.1, 0.2]}]),
                ok=r.random() < 0.8, notes=["REJECTED_PROVEN", "CERTIFIED_CANDIDATE", "EXCLUDED_ASSUMPTION"]))
        if r.random() < 0.1 and out:
            out.pop()
        return out


class FixedGenerator:
    """Returns the given structures (for goal-semantics tests)."""

    def __init__(self, structures: list[dict], name: str = "fixed"):
        self.name = name
        self.structures = structures

    def generate(self, request: GenerationRequest) -> list[GeneratedStructure]:
        return [GeneratedStructure(copy.deepcopy(s), self.name, request.request_id, {})
                for s in self.structures]


class Kill(BaseException):
    """Simulated SIGKILL-like interruption: not an Exception, so the pipeline cannot catch it."""


class KillingExecutor(LocalExecutor):
    """Runs the pipeline, but raises Kill after ``after`` certify calls in round ``round_tag``."""

    def __init__(self, round_tag: str, after: int):
        self.round_tag = round_tag
        self.after = after

    def run_pipeline(self, cfg, candidates, **kw):
        if self.round_tag not in cfg["run"]["workdir"]:
            return super().run_pipeline(cfg, candidates, **kw)
        import rtsc.pipeline as P
        orig = P.certify_spectra
        n = {"k": 0}

        def wrapped(*a, **k):
            if n["k"] >= self.after:
                raise Kill("killed mid-round")
            n["k"] += 1
            return orig(*a, **k)

        P.certify_spectra = wrapped
        try:
            return super().run_pipeline(cfg, candidates, **kw)
        finally:
            P.certify_spectra = orig


def finite(x) -> bool:
    return x is not None and math.isfinite(x)
