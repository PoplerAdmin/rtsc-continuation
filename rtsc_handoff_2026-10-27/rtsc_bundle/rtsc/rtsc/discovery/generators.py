"""Built-in generators / pre-screener usable without ML, and the adapter factory.

* ``PrototypeSubstitutionGenerator``: substitutes the elements of prototype structures (a local
  directory, read with ``rtsc.sources.load_local``) by the request's chemical system, reusing
  ``rtsc.sources.prototype.substitute``.  Hypothetical, unrelaxed (prototype lattice kept).
* ``DirectoryGenerator``: serves structures from a folder (e.g. AIRSS / USPEX / CALYPSO / CIF /
  POSCAR outputs of a crystal-structure-prediction run) whose element set equals the request's.
* ``NullPrescreener``: no ML; every structure passes with no information (ordering by policy only).

External adapters: ``type = "mattergen"`` (``rtsc.generate.mattergen.MatterGenGenerator``; the
other keys of the table form its config dict), prescreen ``type = "mattersim"``
(``rtsc.mlff.mattersim.MatterSimPrescreener``; likewise), or any class via ``type = "python"`` and
``class = "package.module:ClassName"`` (every other key is passed to the constructor as a keyword
argument).
"""
from __future__ import annotations

import importlib
import itertools
from typing import Any

from ..pipeline import ConfigError
from ..sources import load_one
from ..sources.base import Candidate
from ..sources.prototype import substitute
from .interfaces import GeneratedStructure, GenerationRequest, PrescreenResult


def _struct(c: Candidate) -> dict[str, Any]:
    return {"lattice": c.lattice, "species": list(c.species), "frac_coords": c.frac_coords}


def _rotate(items: list, req: GenerationRequest) -> list:
    """Deterministic progress through a finite list: start at arm_pulls * n (cyclic)."""
    if not items:
        return items
    start = (int(req.extra.get("arm_pulls", 0)) * req.n) % len(items)
    return (items[start:] + items[:start])[:req.n]


class PrototypeSubstitutionGenerator:
    def __init__(self, prototypes: str | list[Candidate], *, keep_common: bool = True,
                 name: str = "prototype_substitution", backend: str = "auto",
                 recursive: bool = False):
        self.name = name
        self.keep_common = keep_common
        if isinstance(prototypes, str):
            protos = load_one({"type": "local", "path": prototypes, "backend": backend,
                               "recursive": recursive})
        else:
            protos = list(prototypes)
        self.prototypes = sorted([p for p in protos if p.species and p.lattice is not None
                                  and p.frac_coords is not None and not p.properties.get("load_error")],
                                 key=lambda p: p.material_id)
        if not self.prototypes:
            raise ConfigError(f"{name}: no usable (ordered) prototype structures")

    def generate(self, request: GenerationRequest) -> list[GeneratedStructure]:
        target = tuple(request.chemical_system)
        out: list[GeneratedStructure] = []
        for proto in self.prototypes:
            els = sorted(set(proto.species))
            if len(els) != len(target):
                continue
            sgn = proto.properties.get("spacegroup_number")
            if request.space_group is not None and sgn is not None and int(sgn) != request.space_group:
                continue
            for perm in itertools.permutations(target):
                mapping = dict(zip(els, perm))
                if self.keep_common and any(e in target and mapping[e] != e for e in els):
                    continue
                subs = substitute([proto], {e: [mapping[e]] for e in els}, include_original=True,
                                  id_prefix="")
                for c in subs:
                    out.append(GeneratedStructure(
                        structure=_struct(c), generator=self.name, request_id=request.request_id,
                        provenance={"prototype": proto.material_id, "mapping": mapping,
                                    "prototype_file": proto.provenance.get("file"),
                                    "note": "hypothetical, unrelaxed (prototype lattice kept)"}))
        return _rotate(out, request)


class DirectoryGenerator:
    def __init__(self, path: str, *, format: str = "local", name: str = "directory",
                 recursive: bool = True, pressure_tol_GPa: float = 1.0, **source_kw: Any):
        self.name = name
        self.path = path
        self.tol = float(pressure_tol_GPa)
        entry = {"type": format, "path": path, **source_kw}
        if format == "local":
            entry["recursive"] = recursive
        self.items = sorted([c for c in load_one(entry) if c.species and c.lattice is not None
                             and c.frac_coords is not None and not c.properties.get("load_error")],
                            key=lambda c: c.material_id)

    def generate(self, request: GenerationRequest) -> list[GeneratedStructure]:
        target = set(request.chemical_system)
        out = []
        for c in self.items:
            if set(c.species) != target:
                continue
            if c.pressure_GPa is not None and abs(float(c.pressure_GPa) - request.pressure_GPa) > self.tol:
                continue
            out.append(GeneratedStructure(
                structure=_struct(c), generator=self.name, request_id=request.request_id,
                provenance={"file": c.provenance.get("file"), "source_id": c.material_id,
                            "source_pressure_GPa": c.pressure_GPa,
                            "source_enthalpy_eV_per_atom": c.enthalpy_eV_per_atom}))
        return _rotate(out, request)


class NullPrescreener:
    name = "null"

    def prescreen(self, structures: list[GeneratedStructure], pressure_GPa: float) -> list[PrescreenResult]:
        return [PrescreenResult(fingerprint=str(s.provenance.get("fingerprint", "")), ok=True,
                                notes=["null pre-screener: no ML information"]) for s in structures]


def _import(spec: str):
    if ":" not in spec:
        raise ConfigError(f"class must be 'module:Name', got {spec!r}")
    mod, name = spec.split(":", 1)
    return getattr(importlib.import_module(mod), name)


def make_generator(name: str, gcfg: dict[str, Any]):
    g = dict(gcfg)
    t = g.pop("type", None)
    if t == "prototype_substitution":
        return PrototypeSubstitutionGenerator(g.pop("prototypes", g.pop("path", None)), name=name, **g)
    if t == "directory":
        return DirectoryGenerator(g.pop("path"), name=name, **g)
    if t == "mattergen":       # rtsc.generate.mattergen (other keys -> MatterGenConfig)
        obj = _import("rtsc.generate.mattergen:MatterGenGenerator")(config=g)
        return obj
    if t == "python":
        cls = _import(g.pop("class"))
        obj = cls(**g)
        if getattr(obj, "name", None) in (None, ""):
            obj.name = name
        return obj
    raise ConfigError(f"[discovery.generators.{name}]: unknown type {t!r} "
                      "(prototype_substitution | directory | mattergen | python)")


def make_prescreener(pcfg: dict[str, Any]):
    p = {k: v for k, v in pcfg.items()
         if k not in ("defer_if_ml_unstable", "max_e_above_hull_ml_eV", "use_ml_relaxed")}
    t = p.pop("type", "null")
    if t == "null":
        return NullPrescreener()
    if t == "mattersim":       # rtsc.mlff.mattersim (other keys -> MatterSimConfig)
        return _import("rtsc.mlff.mattersim:MatterSimPrescreener")(config=p)
    if t == "python":
        return _import(p.pop("class"))(**p)
    raise ConfigError(f"[discovery.prescreen]: unknown type {t!r} (null | mattersim | python)")
