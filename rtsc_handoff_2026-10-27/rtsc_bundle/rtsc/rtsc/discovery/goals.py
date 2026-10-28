"""Discovery goals and the ``[discovery]`` configuration (DISCOVERY.md sections 1, 2, 5).

A goal is (id, T*, P_max).  A material satisfies a goal ONLY when ``certify.classify`` returned
CERTIFIED_CANDIDATE with that goal's T* for that material at a pressure P <= P_max.  The ambient
goal G-AMB has P_max = 1e-4 GPa (fixed upper limit; a larger value is a configuration error).

The discovery TOML file is an ordinary pipeline config (``[run]`` is ignored, ``[filters]``,
``[spectra]``, ``[dfpt]``, ``[certify]``, ``[validity]``, ``[priority]``, ``[report]`` apply to every
round) plus a ``[discovery]`` table; see README "Discovery loop" for every key.
"""
from __future__ import annotations

import copy
import math
import os
import tomllib
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .. import pipeline
from ..pipeline import ConfigError

AMBIENT_P_MAX_GPa = 1e-4

DEFAULT_GOALS: dict[str, dict[str, Any]] = {
    "G-RT": {"enabled": True, "T_star_K": 300.0, "P_max_GPa": 200.0,
             "pressures_GPa": [0.0, 50.0, 100.0, 150.0, 200.0]},
    "G-AMB": {"enabled": True, "T_star_K": 77.0, "P_max_GPa": AMBIENT_P_MAX_GPa,
              "pressures_GPa": [0.0], "ambient": True},
}

DISCOVERY_DEFAULTS: dict[str, Any] = {
    "workdir": "discovery_run",
    "seed": 0,
    "executor": "local",                     # local | slurm
    "workers": 1,
    "id_prefix": "D",
    "goals": {},
    "search": {"chemical_systems": [], "space_groups": [], "policy": "thompson",
               "requests_per_round": 4, "structures_per_request": 4,
               "exploration_floor": 0.25, "ucb_c": 1.0},
    "generators": {},                        # name -> {type = ..., ...}
    "prescreen": {"type": "null", "defer_if_ml_unstable": True,
                  "max_e_above_hull_ml_eV": 0.2, "use_ml_relaxed": True},
    "budget": {"dft_per_round": 4},
    "stop": {"certified_per_goal": 1, "max_rounds": 10, "max_wall_hours": None,
             "max_dft_jobs": None},
    "dedupe": {"ltol": 0.2, "stol": 0.3, "angle_tol": 5.0, "use_structure_matcher": True,
               "confirm_tol_A": 0.02},
    "slurm": {},
    "synthetic_spectra": {"enabled": False},
}
_TOP_KEYS = set(DISCOVERY_DEFAULTS)
_SUB_KEYS = {
    "search": set(DISCOVERY_DEFAULTS["search"]),
    "prescreen": None,                       # adapter-specific keys allowed
    "budget": {"dft_per_round"},
    "stop": set(DISCOVERY_DEFAULTS["stop"]),
    "dedupe": set(DISCOVERY_DEFAULTS["dedupe"]),
    "synthetic_spectra": {"enabled", "A_range_K", "center_range_K", "width_frac"},
}
_GOAL_KEYS = {"enabled", "T_star_K", "P_max_GPa", "pressures_GPa", "chemical_systems",
              "space_groups", "ambient"}


@dataclass(frozen=True)
class Goal:
    id: str
    T_star_K: float
    P_max_GPa: float
    pressures_GPa: tuple[float, ...]
    chemical_systems: tuple[tuple[str, ...], ...] = ()
    space_groups: tuple[int | None, ...] = (None,)
    ambient: bool = False
    enabled: bool = True

    def in_scope(self, pressure_GPa: float) -> bool:
        return pressure_GPa is not None and float(pressure_GPa) <= self.P_max_GPa

    def to_dict(self) -> dict[str, Any]:
        return {"id": self.id, "T_star_K": self.T_star_K, "P_max_GPa": self.P_max_GPa,
                "pressures_GPa": list(self.pressures_GPa),
                "chemical_systems": [list(c) for c in self.chemical_systems],
                "space_groups": list(self.space_groups), "ambient": self.ambient,
                "enabled": self.enabled}


def _chemsys(x) -> tuple[str, ...]:
    if isinstance(x, str):
        x = [e for e in x.replace("-", " ").split() if e]
    from ..sources.composition import is_element
    out = tuple(sorted(set(str(e) for e in x)))
    if not out or not all(is_element(e) for e in out):
        raise ConfigError(f"bad chemical system {x!r}")
    return out


def parse_goals(gcfg: dict[str, Any], default_chemsys: list) -> list[Goal]:
    merged = copy.deepcopy(DEFAULT_GOALS)
    for gid, g in (gcfg or {}).items():
        if not isinstance(g, dict):
            raise ConfigError(f"[discovery.goals.{gid}] must be a table")
        bad = set(g) - _GOAL_KEYS
        if bad:
            raise ConfigError(f"unknown keys in [discovery.goals.{gid}]: {sorted(bad)}")
        merged[gid] = {**merged.get(gid, {}), **g}
    goals = []
    for gid, g in merged.items():
        if "T_star_K" not in g or "P_max_GPa" not in g:
            raise ConfigError(f"goal {gid}: T_star_K and P_max_GPa are required")
        T = float(g["T_star_K"])
        P = float(g["P_max_GPa"])
        amb = bool(g.get("ambient", gid == "G-AMB"))
        if not (T > 0 and math.isfinite(T)):
            raise ConfigError(f"goal {gid}: T_star_K must be > 0")
        if amb:
            if P > AMBIENT_P_MAX_GPa:
                raise ConfigError(f"goal {gid} is ambient: P_max_GPa must be <= {AMBIENT_P_MAX_GPa}")
            if T > 300.0:
                raise ConfigError(f"goal {gid}: ambient T_star_K is configurable up to 300 K")
        if not (P >= 0 and math.isfinite(P)):
            raise ConfigError(f"goal {gid}: P_max_GPa must be >= 0")
        pressures = tuple(float(p) for p in g.get("pressures_GPa", [0.0]))
        bad_p = [p for p in pressures if p > P or p < 0]
        if bad_p:
            raise ConfigError(f"goal {gid}: pressures {bad_p} outside [0, P_max={P}]")
        cs = g.get("chemical_systems") or default_chemsys
        sgs = tuple(None if (s is None or s == 0) else int(s) for s in (g.get("space_groups") or [None]))
        goals.append(Goal(gid, T, P, pressures, tuple(_chemsys(c) for c in cs), sgs, amb,
                          bool(g.get("enabled", True))))
    if not any(g.enabled for g in goals):
        raise ConfigError("no enabled discovery goal")
    return goals


@dataclass
class DiscoveryConfig:
    d: dict[str, Any]                        # normalized [discovery] table
    pipeline: dict[str, Any]                 # normalized pipeline config (base for every round)
    goals: list[Goal]
    base_dir: str | None = None
    raw: dict[str, Any] = field(default_factory=dict)

    @property
    def workdir(self) -> Path:
        return Path(self.d["workdir"])

    @property
    def enabled_goals(self) -> list[Goal]:
        return [g for g in self.goals if g.enabled]

    def goal(self, gid: str) -> Goal:
        for g in self.goals:
            if g.id == gid:
                return g
        raise KeyError(gid)

    def goal_pipeline_config(self, goal: Goal, workdir: str | os.PathLike) -> dict[str, Any]:
        """The pipeline config of one goal in one round: certify at the goal's T*; the loop has
        already deduplicated and selected, so pipeline dedupe/budget/ML are neutral."""
        c = copy.deepcopy(self.pipeline)
        c["run"]["workdir"] = str(workdir)
        c["run"]["resume"] = True
        c["run"]["workers"] = int(self.d.get("workers", 1))
        c["sources"] = []
        c["dedupe"]["method"] = "none"
        c["budget"] = {}
        c["ml"] = {"enabled": False}
        c["certify"]["T_target_K"] = float(goal.T_star_K)
        pipeline.validate_certify_config(c["certify"])
        return c

    def certify_hash(self, goal: Goal) -> str:
        """Workdir-independent hash of what certify.classify depends on for this goal besides the
        spectra: the [certify] table with T_target_K = goal T*."""
        c = copy.deepcopy(self.pipeline["certify"])
        c["T_target_K"] = float(goal.T_star_K)
        return pipeline.config_hash({"certify": c})

    def dft_hash(self) -> str:
        """Hash of the config parts that determine spectra (everything except [certify])."""
        c = copy.deepcopy(self.pipeline)
        c.pop("certify", None)
        c["run"] = {}
        return pipeline.config_hash(c)

    def to_json(self) -> dict[str, Any]:
        return {"discovery": self.d, "pipeline": self.pipeline,
                "goals": [g.to_dict() for g in self.goals], "base_dir": self.base_dir,
                "raw": self.raw}


def check_dft_pressure_policy(pcfg: dict[str, Any]) -> None:
    """The DFT stage must compute every material AT ITS OWN pressure (the loop sets it per
    material): a global [dfpt.settings] pressure, no relaxation, or a fixed-volume relaxation
    would make the spectra belong to another pressure (audit F3)."""
    st = (pcfg.get("dfpt") or {}).get("settings") or {}
    bad = [k for k in ("pressure_GPa", "press_kbar") if k in st]
    if bad:
        raise ConfigError(f"[dfpt.settings] {bad}: not allowed in a discovery config; each material is "
                          "computed at its own (request) pressure")
    if (pcfg.get("dfpt") or {}).get("enabled") and st.get("do_relax", True) is False:
        raise ConfigError("[dfpt.settings] do_relax = false: not allowed in a discovery config (the "
                          "generated geometry is not at the target pressure)")
    if str(st.get("cell_dofree", "all")) in ("shape", "2Dshape", "x", "y", "z", "xy", "xz", "yz"):
        raise ConfigError(f"[dfpt.settings] cell_dofree = {st['cell_dofree']!r} cannot reach the "
                          "target pressure; use 'all'")


def counts_for(goal: Goal, row: dict[str, Any], *, certify_hash: str, dft_hash: str) -> tuple[bool, str]:
    """THE single predicate (DISCOVERY.md section 3): does this current ledger row make the
    material satisfy ``goal``?  Returns (counts, reason-if-not).  True iff
      status = CERTIFIED_CANDIDATE, stage = pipeline, origin = rtsc.certify.classify,
      T_star_K = goal.T_star_K, material pressure AND DFT pressure in the goal's range (recomputed
      with the CURRENT goal), certify_hash = current certify_hash(goal), dft_hash = current dft_hash.
    """
    if row.get("status") != "CERTIFIED_CANDIDATE":
        return False, "not CERTIFIED_CANDIDATE"
    if row.get("stage") != "pipeline" or row.get("origin") != pipeline.CLASSIFY_ORIGIN:
        return False, "not from certify.classify"
    if row.get("T_star_K") != goal.T_star_K:
        return False, f"certified at T*={row.get('T_star_K')} K, goal T*={goal.T_star_K} K"
    p, pd = row.get("pressure_GPa"), row.get("dft_pressure_GPa")
    if p is None or not goal.in_scope(p):
        return False, f"pressure {p} GPa outside goal range"
    if pd is None or not goal.in_scope(pd) or abs(float(pd) - float(p)) > 1e-9:
        return False, f"DFT pressure {pd} GPa unknown, outside goal range or != material pressure"
    if row.get("certify_hash") != certify_hash:
        return False, "certified under a different [certify] configuration"
    if row.get("dft_hash") != dft_hash:
        return False, "spectra computed under a different DFT/spectra/filter configuration"
    return True, ""


def _merge(a: dict, b: dict) -> dict:
    out = copy.deepcopy(a)
    for k, v in b.items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _merge(out[k], v)
        else:
            out[k] = copy.deepcopy(v)
    return out


def normalize_discovery(raw: dict[str, Any], base_dir: str | os.PathLike | None = None) -> DiscoveryConfig:
    raw = copy.deepcopy(raw)
    draw = raw.pop("discovery", {}) or {}
    bad = set(draw) - _TOP_KEYS
    if bad:
        raise ConfigError(f"unknown keys in [discovery]: {sorted(bad)}")
    d = _merge(DISCOVERY_DEFAULTS, draw)
    for sec, keys in _SUB_KEYS.items():
        if keys is not None:
            bad = set(d[sec]) - keys
            if bad:
                raise ConfigError(f"unknown keys in [discovery.{sec}]: {sorted(bad)}")
    if d["executor"] not in ("local", "slurm"):
        raise ConfigError("discovery.executor must be 'local' or 'slurm'")
    s = d["search"]
    if s["policy"] not in ("thompson", "ucb"):
        raise ConfigError("discovery.search.policy must be 'thompson' or 'ucb'")
    if not (0.0 <= float(s["exploration_floor"]) <= 1.0):
        raise ConfigError("discovery.search.exploration_floor must be in [0, 1]")
    for k in ("requests_per_round", "structures_per_request"):
        if int(s[k]) < 1:
            raise ConfigError(f"discovery.search.{k} must be >= 1")
    if int(d["budget"]["dft_per_round"]) < 0:
        raise ConfigError("discovery.budget.dft_per_round must be >= 0")
    if not d["generators"]:
        raise ConfigError("[discovery.generators.<name>] needs at least one generator")
    raw.pop("sources", None)
    pcfg = pipeline.normalize_config(raw, base_dir)
    check_dft_pressure_policy(pcfg)
    if base_dir is not None:
        b = Path(base_dir)
        if not Path(d["workdir"]).is_absolute():
            d["workdir"] = str((b / d["workdir"]).resolve())
        for g in d["generators"].values():
            for key in ("path", "prototypes"):
                if isinstance(g.get(key), str) and not Path(g[key]).is_absolute():
                    g[key] = str((b / g[key]).resolve())
    goals = parse_goals(d["goals"], s["chemical_systems"])
    for g in goals:
        if g.enabled and not g.chemical_systems:
            raise ConfigError(f"goal {g.id}: no chemical systems (set discovery.search."
                              "chemical_systems or the goal's chemical_systems)")
    return DiscoveryConfig(d, pcfg, goals, None if base_dir is None else str(base_dir),
                           {"discovery": draw, **raw})


def load_discovery_config(path: str | os.PathLike) -> DiscoveryConfig:
    p = Path(path)
    with open(p, "rb") as fh:
        raw = tomllib.load(fh)
    return normalize_discovery(raw, p.parent.resolve())
