"""Candidate sources (ARCH section 6, stages 1-2).

``load_sources(entries)`` dispatches on each entry's ``type``:

=========================  ==============================================================
type                       keys
=========================  ==============================================================
``local``                  path, recursive, pressure_GPa, backend, manifest, id_prefix
``materials_project``      elements, chemsys, exclude_elements, energy_above_hull_max,
                           is_metal, magnetic_ordering, num_elements, material_ids,
                           fields, extra_fields, limit, api_key_env
``airss``                  path, recursive, pressure_GPa, id_prefix
``uspex``                  path, poscars, individuals, enthalpy_unit, pressure_GPa, id_prefix
``calypso``                path, table, enthalpy_unit, pressure_GPa, id_prefix
``seed_list``              path, id_prefix
``prototype_substitution`` path (local dir of prototypes), substitutions, max_per_prototype,
                           include_original, pressure_GPa
``synthetic_demo``         (none) -- clearly labelled SYNTHETIC demo candidates
=========================  ==============================================================
"""
from __future__ import annotations

from typing import Any

from .base import Candidate, safe_id
from .csp import load_airss, load_calypso, load_seed_list, load_uspex
from .dedupe import DedupeResult, dedupe
from .local import load_local
from .mp import fetch_mp_entries, query_materials_project
from .prototype import substitute
from .synthetic import demo_candidates

__all__ = ["Candidate", "safe_id", "load_sources", "load_local", "query_materials_project",
           "fetch_mp_entries", "load_airss", "load_uspex", "load_calypso", "load_seed_list",
           "substitute", "dedupe", "DedupeResult", "demo_candidates", "ensure_unique_ids"]


def _pick(d: dict[str, Any], *keys: str) -> dict[str, Any]:
    return {k: d[k] for k in keys if k in d}


def load_one(entry: dict[str, Any], *, rester_factory=None) -> list[Candidate]:
    t = entry.get("type")
    if t == "local":
        return load_local(entry["path"], **_pick(entry, "recursive", "pressure_GPa", "backend",
                                                   "manifest", "id_prefix"))
    if t == "materials_project":
        kw = _pick(entry, "elements", "chemsys", "exclude_elements", "energy_above_hull_max",
                   "is_metal", "magnetic_ordering", "num_elements", "material_ids", "fields",
                   "extra_fields", "limit", "api_key_env")
        cands = query_materials_project(rester_factory=rester_factory, **kw)
        if entry.get("pressure_GPa"):
            # MP structures are 0 GPa; we only RECORD the requested target pressure (the DFT
            # stage relaxes at it).  MP's energy_above_hull stays a 0 GPa quantity.
            for c in cands:
                c.properties["target_pressure_GPa"] = float(entry["pressure_GPa"])
        return cands
    if t == "airss":
        return load_airss(entry["path"], **_pick(entry, "recursive", "pressure_GPa", "id_prefix"))
    if t == "uspex":
        return load_uspex(entry["path"], **_pick(entry, "poscars", "individuals", "enthalpy_unit",
                                                   "pressure_GPa", "id_prefix"))
    if t == "calypso":
        return load_calypso(entry["path"], **_pick(entry, "table", "enthalpy_unit",
                                                     "pressure_GPa", "id_prefix"))
    if t == "seed_list":
        return load_seed_list(entry["path"], **_pick(entry, "id_prefix"))
    if t == "prototype_substitution":
        protos = load_local(entry["path"], **_pick(entry, "pressure_GPa", "backend"))
        return substitute(protos, entry["substitutions"],
                          **_pick(entry, "max_per_prototype", "include_original"))
    if t == "synthetic_demo":
        return demo_candidates()
    raise ValueError(f"unknown source type {t!r}")


def ensure_unique_ids(cands: list[Candidate]) -> dict[str, str]:
    """Rename repeated material ids (``id``, ``id~2``, ...) so S4 accounting is by unique key.
    Returns {new_id: original_id} for renamed ones (also recorded in provenance)."""
    seen: dict[str, int] = {}
    renamed = {}
    for c in cands:
        n = seen.get(c.material_id, 0) + 1
        seen[c.material_id] = n
        if n > 1:
            new = f"{c.material_id}~{n}"
            while new in seen:
                n += 1
                new = f"{c.material_id}~{n}"
            renamed[new] = c.material_id
            c.provenance["original_material_id"] = c.material_id
            c.material_id = new
            seen[new] = 1
    return renamed


def load_sources(entries: list[dict[str, Any]], *, rester_factory=None) -> list[Candidate]:
    out: list[Candidate] = []
    for i, e in enumerate(entries):
        cands = load_one(e, rester_factory=rester_factory)
        for c in cands:
            c.provenance.setdefault("source_entry", i)
        out.extend(cands)
    ensure_unique_ids(out)
    return out
