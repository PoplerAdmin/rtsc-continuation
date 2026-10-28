"""Adapters for crystal-structure-prediction outputs: AIRSS, USPEX, CALYPSO (+ seed lists).

Each adapter returns structures with the enthalpy the CSP code reported, converted to eV/atom,
and records the original value, unit convention and file in ``provenance``.  Unit conventions
differ between codes and versions; the defaults below are documented assumptions and can be
overridden (``enthalpy_unit = "eV/cell" | "eV/atom"``).  Verify them for your version: the
enthalpies feed the convex-hull EXCLUDED_ASSUMPTION filter.

* AIRSS: ``*.res`` files; TITL line ``seed P[GPa] V[A^3] H[eV, whole cell] spin |spin| n_ions
  (sg) n - times``.  Enthalpy per atom = H / n_ions.
* USPEX: ``gatheredPOSCARS`` (titles ``EA<ID> ...``) joined by ID with the ``Individuals`` table
  (column ``Enthalpy``; default unit eV per cell); pressure from ``INPUT.txt``
  ``ExternalPressure`` (GPa) when present, else the ``pressure_GPa`` argument.
* CALYPSO: POSCAR-like files (``*.vasp``, ``POSCAR_*``, ``UCell_*``) in the analysis directory
  written by ``cak.py`` plus an optional table (``Analysis_Output.dat``) with an index column and
  an ``Enthalpy``/``Energy`` column (default unit eV per atom).  Best effort: the format varies.
* Seed lists: CSV/JSON/TOML lists of formulas (no structure) with optional pressure/properties.
"""
from __future__ import annotations

import csv
import json
import re
import tomllib
from pathlib import Path
from typing import Any

from .base import Candidate, candidate_from_structure_arrays
from .fileformats import read_poscar, read_res, split_concatenated_poscars

_UNITS = ("eV/atom", "eV/cell")


def _per_atom(value: float, natoms: int, unit: str) -> float:
    if unit not in _UNITS:
        raise ValueError(f"enthalpy_unit must be one of {_UNITS}, got {unit!r}")
    return value / natoms if unit == "eV/cell" else value


# --------------------------------------------------------------------------- AIRSS
def load_airss(path: str | Path, *, recursive: bool = True, id_prefix: str = "airss:",
               pressure_GPa: float | None = None) -> list[Candidate]:
    root = Path(path)
    files = sorted(root.rglob("*.res") if recursive else root.glob("*.res"))
    out = []
    for f in files:
        mid = id_prefix + f.stem
        try:
            raw = read_res(f.read_text(errors="replace"))
        except Exception as exc:  # noqa: BLE001
            out.append(Candidate(mid, {}, source="airss",
                                 properties={"load_error": f"{type(exc).__name__}: {exc}"},
                                 provenance={"file": str(f.resolve())}))
            continue
        info = raw.info
        n = int(info.get("n_ions") or len(raw.species))
        h = info.get("enthalpy_eV")
        p = info.get("pressure_GPa", pressure_GPa)
        prov = {"source": "airss", "file": str(f.resolve()), "titl": {k: v for k, v in info.items()},
                "enthalpy_convention": "TITL enthalpy is eV per cell; divided by n_ions"}
        out.append(candidate_from_structure_arrays(
            mid, raw.lattice, raw.species, raw.frac_coords, source="airss",
            spacegroup=info.get("spacegroup"), pressure_GPa=p,
            enthalpy_eV_per_atom=None if h is None else float(h) / n, provenance=prov))
    return out


# --------------------------------------------------------------------------- USPEX
def _collapse_brackets(line: str) -> list[str]:
    return re.sub(r"\[[^\]]*\]", lambda m: m.group(0).replace(" ", ""), line).split()


def _read_uspex_individuals(path: Path) -> dict[str, dict[str, Any]]:
    lines = path.read_text(errors="replace").splitlines()
    header = None
    rows: dict[str, dict[str, Any]] = {}
    for ln in lines:
        toks = _collapse_brackets(ln)
        if not toks:
            continue
        if header is None:
            if "ID" in toks and any(t.lower().startswith(("enthalpy", "energy")) for t in toks):
                header = toks
            continue
        if not re.fullmatch(r"-?\d+", toks[0]):
            continue
        rec = dict(zip(header, toks))
        rows[str(rec.get("ID"))] = rec
    if header is None:
        raise ValueError(f"no header with ID/Enthalpy in {path}")
    return rows


def _uspex_pressure(root: Path) -> float | None:
    for cand in (root / "INPUT.txt", root.parent / "INPUT.txt"):
        if cand.is_file():
            txt = cand.read_text(errors="replace")
            m = re.search(r"^\s*([-+0-9.eE]+)\s*:\s*ExternalPressure", txt, re.M)
            if m:
                return float(m.group(1))
    return None


def load_uspex(path: str | Path, *, poscars: str = "gatheredPOSCARS",
               individuals: str = "Individuals", enthalpy_unit: str = "eV/cell",
               id_prefix: str = "uspex:", pressure_GPa: float | None = None) -> list[Candidate]:
    root = Path(path)
    pfile = root / poscars
    if not pfile.is_file():
        hits = sorted(root.rglob(poscars))
        if not hits:
            raise FileNotFoundError(f"{poscars} not found under {root}")
        pfile = hits[0]
    ind_file = pfile.parent / individuals
    table = _read_uspex_individuals(ind_file) if ind_file.is_file() else {}
    p_input = _uspex_pressure(pfile.parent)
    pressure = p_input if p_input is not None else pressure_GPa
    out = []
    for k, block in enumerate(split_concatenated_poscars(pfile.read_text(errors="replace"))):
        title = block.splitlines()[0].strip()
        m = re.match(r"EA(\d+)", title)
        sid = m.group(1) if m else str(k)
        mid = f"{id_prefix}EA{sid}"
        try:
            raw = read_poscar(block)
        except Exception as exc:  # noqa: BLE001
            out.append(Candidate(mid, {}, source="uspex",
                                 properties={"load_error": f"{type(exc).__name__}: {exc}"},
                                 provenance={"file": str(pfile.resolve()), "block": k}))
            continue
        rec = table.get(sid, {})
        h = None
        hkey = next((c for c in rec if c.lower().startswith("enthalpy")), None) or \
            next((c for c in rec if c.lower().startswith("energy")), None)
        if hkey is not None:
            try:
                h = _per_atom(float(rec[hkey]), len(raw.species), enthalpy_unit)
            except ValueError:
                h = None
        sg = None
        for key in ("SYMM", "Symmetry"):
            if key in rec:
                sg = rec[key]
        prov = {"source": "uspex", "file": str(pfile.resolve()), "title": title,
                "individuals_row": rec, "enthalpy_unit_assumed": enthalpy_unit,
                "pressure_from": "INPUT.txt" if p_input is not None else "config"}
        out.append(candidate_from_structure_arrays(
            mid, raw.lattice, raw.species, raw.frac_coords, source="uspex", spacegroup=sg,
            pressure_GPa=pressure, enthalpy_eV_per_atom=h, provenance=prov))
    return out


# --------------------------------------------------------------------------- CALYPSO
_CALYPSO_GLOBS = ("*.vasp", "POSCAR_*", "UCell_*", "CONTCAR_*")


def _read_calypso_table(path: Path) -> dict[str, float]:
    out: dict[str, float] = {}
    header = None
    for ln in path.read_text(errors="replace").splitlines():
        toks = ln.replace(",", " ").split()
        if not toks:
            continue
        low = [t.lower() for t in toks]
        if header is None and any(t.startswith(("enthalpy", "energy")) for t in low):
            header = low
            continue
        if header is None or not re.fullmatch(r"\d+", toks[0]):
            continue
        ie = next(i for i, t in enumerate(header) if t.startswith(("enthalpy", "energy")))
        ii = next((i for i, t in enumerate(header) if t.startswith(("index", "no", "#"))), 0)
        try:
            out[toks[ii]] = float(toks[ie])
        except (IndexError, ValueError):
            continue
    return out


def load_calypso(path: str | Path, *, table: str = "Analysis_Output.dat",
                 enthalpy_unit: str = "eV/atom", id_prefix: str = "calypso:",
                 pressure_GPa: float | None = None) -> list[Candidate]:
    root = Path(path)
    files: set[Path] = set()
    for g in _CALYPSO_GLOBS:
        files.update(root.glob(g))
    tab = _read_calypso_table(root / table) if (root / table).is_file() else {}
    out = []
    for f in sorted(files):
        m = re.search(r"(\d+)", f.stem)
        idx = m.group(1) if m else f.stem
        mid = id_prefix + f.stem
        sgm = re.search(r"_(\d+)$", f.stem)
        try:
            raw = read_poscar(f.read_text(errors="replace"))
        except Exception as exc:  # noqa: BLE001
            out.append(Candidate(mid, {}, source="calypso",
                                 properties={"load_error": f"{type(exc).__name__}: {exc}"},
                                 provenance={"file": str(f.resolve())}))
            continue
        h = tab.get(idx)
        prov = {"source": "calypso", "file": str(f.resolve()), "table": str(root / table)
                if tab else None, "enthalpy_unit_assumed": enthalpy_unit, "table_index": idx}
        out.append(candidate_from_structure_arrays(
            mid, raw.lattice, raw.species, raw.frac_coords, source="calypso",
            spacegroup=sgm.group(1) if (sgm and f.stem.startswith(("UCell", "PCell"))) else None,
            pressure_GPa=pressure_GPa,
            enthalpy_eV_per_atom=None if h is None else _per_atom(h, len(raw.species), enthalpy_unit),
            provenance=prov))
    return out


# --------------------------------------------------------------------------- seed lists
def load_seed_list(path: str | Path, *, id_prefix: str = "") -> list[Candidate]:
    """CSV (columns: material_id, formula, [pressure_GPa], [any property]), JSON or TOML
    (list of Candidate-like dicts under key ``materials`` for TOML)."""
    p = Path(path)
    rows: list[dict[str, Any]]
    if p.suffix == ".csv":
        with open(p, newline="") as fh:
            rows = list(csv.DictReader(fh))
    elif p.suffix == ".toml":
        with open(p, "rb") as fh:
            rows = tomllib.load(fh).get("materials", [])
    else:
        data = json.loads(p.read_text())
        rows = data["materials"] if isinstance(data, dict) else data
    out = []
    for k, r in enumerate(rows):
        r = dict(r)
        mid = id_prefix + str(r.pop("material_id", r.get("formula", f"seed{k}")))
        formula = r.pop("formula", None)
        comp = r.pop("composition", None) or formula or {}
        pressure = r.pop("pressure_GPa", None)
        enth = r.pop("enthalpy_eV_per_atom", None)
        sg = r.pop("spacegroup", None)
        synthetic = str(r.pop("synthetic", "false")).lower() in ("1", "true", "yes")
        props = dict(r.pop("properties", {}) or {})
        spectra_files = r.pop("spectra_files", []) or []
        for key, v in r.items():   # remaining CSV columns become properties
            if v in ("", None):
                continue
            props[key] = _coerce(v)
        for sf in spectra_files:
            sp = Path(sf["path"])
            if not sp.is_absolute():
                sf["path"] = str((p.parent / sp).resolve())
        out.append(Candidate(mid, comp, source="seed_list", spacegroup=sg,
                             pressure_GPa=None if pressure in (None, "") else float(pressure),
                             enthalpy_eV_per_atom=None if enth in (None, "") else float(enth),
                             properties=props, spectra_files=list(spectra_files),
                             synthetic=synthetic,
                             provenance={"source": "seed_list", "file": str(p.resolve()),
                                         "row": k}))
    return out


def _coerce(v: Any) -> Any:
    if not isinstance(v, str):
        return v
    lv = v.strip().lower()
    if lv in ("true", "false"):
        return lv == "true"
    try:
        return int(v)
    except ValueError:
        pass
    try:
        return float(v)
    except ValueError:
        return v
