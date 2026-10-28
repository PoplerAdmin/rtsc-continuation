"""Materials Project source (optional dependency ``mp-api``; extra ``rtsc[mp]``).

The API key is read from the environment variable ``MP_API_KEY`` (never from the config file, so
configs can be shared).  All MP numbers are DFT (GGA/GGA+U/r2SCAN mixing scheme) values at
zero pressure and zero temperature; they are recorded in ``provenance`` together with the query,
the requested fields, the mp-api version and the database version.

``rester_factory`` lets tests (or users with a custom client) inject an object with the
``MPRester`` interface (context manager exposing ``materials.summary.search`` and optionally
``get_database_version`` / ``get_entries_in_chemsys``).
"""
from __future__ import annotations

import datetime as _dt
import os
from typing import Any, Callable, Iterable

from .base import Candidate, _jsonable

DEFAULT_FIELDS = ("material_id", "formula_pretty", "composition", "structure", "symmetry",
                  "energy_above_hull", "formation_energy_per_atom", "energy_per_atom",
                  "is_metal", "band_gap", "is_magnetic", "ordering", "total_magnetization",
                  "total_magnetization_normalized_formula_units", "theoretical", "nsites",
                  "is_stable", "deprecated")

#: candidate.properties keys copied from MP documents (MP field -> property name)
PROPERTY_MAP = {
    "energy_above_hull": "energy_above_hull", "formation_energy_per_atom": "formation_energy_per_atom",
    "energy_per_atom": "energy_per_atom", "is_metal": "is_metal", "band_gap": "band_gap",
    "is_magnetic": "is_magnetic", "ordering": "ordering",
    "total_magnetization": "total_magnetization", "theoretical": "theoretical",
    "is_stable": "is_stable", "nsites": "nsites",
}


class MPUnavailable(RuntimeError):
    pass


def _default_rester_factory(api_key: str):
    try:
        from mp_api.client import MPRester
    except ImportError as exc:  # pragma: no cover - depends on environment
        raise MPUnavailable("mp-api is not installed; pip install 'rtsc[mp]'") from exc
    return MPRester(api_key)


def _get(doc: Any, key: str, default=None):
    if isinstance(doc, dict):
        return doc.get(key, default)
    return getattr(doc, key, default)


def _api_key(api_key_env: str) -> str:
    key = os.environ.get(api_key_env, "")
    if not key:
        raise MPUnavailable(f"Materials Project API key not found in ${api_key_env}")
    return key


def query_materials_project(*, elements: list[str] | None = None,
                            chemsys: str | list[str] | None = None,
                            exclude_elements: list[str] | None = None,
                            energy_above_hull_max: float | None = 0.1,
                            is_metal: bool | None = None,
                            magnetic_ordering: str | None = None,
                            num_elements: tuple[int, int] | None = None,
                            material_ids: list[str] | None = None,
                            fields: Iterable[str] = DEFAULT_FIELDS,
                            extra_fields: Iterable[str] = (),
                            limit: int | None = None,
                            api_key_env: str = "MP_API_KEY",
                            rester_factory: Callable[[str], Any] | None = None) -> list[Candidate]:
    """Query MP summary documents and return Candidates (pressure 0 GPa).

    ``magnetic_ordering`` is one of MP's ``Ordering`` values ("NM", "FM", "AFM", "FiM", "Unknown").
    NOTE: query-side filters (``is_metal``, hull, ordering) only restrict what is downloaded;
    the recorded EXCLUDED_ASSUMPTION verdicts come from ``rtsc.stability`` on the stored values.
    """
    key = _api_key(api_key_env)
    factory = rester_factory or _default_rester_factory
    fields = list(dict.fromkeys(list(fields) + list(extra_fields)))
    kw: dict[str, Any] = {"fields": fields}
    if elements:
        kw["elements"] = list(elements)
    if chemsys:
        kw["chemsys"] = chemsys
    if exclude_elements:
        kw["exclude_elements"] = list(exclude_elements)
    if energy_above_hull_max is not None:
        kw["energy_above_hull"] = (0.0, float(energy_above_hull_max))
    if is_metal is not None:
        kw["is_metal"] = bool(is_metal)
    if magnetic_ordering:
        try:
            from emmet.core.magnetism import Ordering  # type: ignore
            kw["magnetic_ordering"] = Ordering(magnetic_ordering)
        except Exception:  # noqa: BLE001
            kw["magnetic_ordering"] = magnetic_ordering
    if num_elements:
        kw["num_elements"] = tuple(num_elements)
    if material_ids:
        kw["material_ids"] = list(material_ids)
    query_record = {k: (v if k != "magnetic_ordering" else str(v)) for k, v in kw.items()}
    now = _dt.datetime.now(_dt.timezone.utc).isoformat()
    with factory(key) as mpr:
        docs = list(mpr.materials.summary.search(**kw))
        db_version = None
        if hasattr(mpr, "get_database_version"):
            try:
                db_version = mpr.get_database_version()
            except Exception:  # noqa: BLE001
                db_version = None
    try:
        import mp_api  # type: ignore
        client_version = getattr(mp_api, "__version__", None)
    except ImportError:
        client_version = None
    if limit is not None:
        docs = docs[:limit]
    out = []
    for doc in docs:
        out.append(_doc_to_candidate(doc, fields, query_record, db_version, client_version, now))
    return out


def _doc_to_candidate(doc, fields, query_record, db_version, client_version, now) -> Candidate:
    mid = str(_get(doc, "material_id"))
    comp_obj = _get(doc, "composition")
    structure = _get(doc, "structure")
    comp: dict[str, float] = {}
    if comp_obj is not None:
        if hasattr(comp_obj, "get_el_amt_dict"):
            comp = {str(k): float(v) for k, v in comp_obj.get_el_amt_dict().items()}
        elif isinstance(comp_obj, dict):
            comp = {str(k): float(v) for k, v in comp_obj.items()}
    if not comp and _get(doc, "formula_pretty"):
        comp = str(_get(doc, "formula_pretty"))  # type: ignore[assignment]
    sym = _get(doc, "symmetry")
    sg = None
    if sym is not None:
        sg = _get(sym, "symbol")
    props: dict[str, Any] = {}
    for f, p in PROPERTY_MAP.items():
        v = _get(doc, f)
        if v is not None:
            props[p] = _jsonable(v)
    requested = {f: _jsonable(_get(doc, f)) for f in fields
                 if f not in ("structure", "composition") and _get(doc, f) is not None}
    prov = {"source": "materials_project", "material_id": mid, "query": _jsonable(query_record),
            "requested_fields": list(fields), "values": requested,
            "database_version": db_version, "mp_api_version": client_version,
            "retrieved_utc": now,
            "note": "MP values are DFT, 0 GPa, 0 K; energy_above_hull is relative to MP's hull"}
    c = Candidate(material_id=f"mp:{mid}" if not mid.startswith("mp:") else mid, composition=comp,
                  source="materials_project", spacegroup=sg, pressure_GPa=0.0,
                  enthalpy_eV_per_atom=props.get("energy_per_atom"), properties=props,
                  provenance=prov)
    if structure is not None:
        c.structure = structure
        try:
            if structure.is_ordered:
                c.lattice = structure.lattice.matrix.tolist()
                c.species = [str(s.specie) for s in structure]
                c.frac_coords = structure.frac_coords.tolist()
        except Exception:  # noqa: BLE001
            pass
    return c


def fetch_mp_entries(chemsys: str | list[str], *, api_key_env: str = "MP_API_KEY",
                     rester_factory: Callable[[str], Any] | None = None,
                     compatible_only: bool = True) -> list[Any]:
    """ComputedEntry objects for a chemical system (for hull construction at 0 GPa)."""
    key = _api_key(api_key_env)
    factory = rester_factory or _default_rester_factory
    with factory(key) as mpr:
        return list(mpr.get_entries_in_chemsys(chemsys, compatible_only=compatible_only))
