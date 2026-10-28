"""ML convex hull at pressure: reference phases relaxed with the SAME potential at the SAME pressure.

HEURISTIC ONLY (ARCH S2): the ML energy above hull may reorder candidates or mark them DEFERRED
with a recorded reason; it never produces EXCLUDED_ASSUMPTION (that needs the DFT hull filter in
``rtsc.stability.hull``).

Reference phases for a chemical system ``X-Y-...`` are collected for every non-empty sub-system
(elements, binaries, ...) from

* ``reference_dir``: any file ASE can read (CIF, POSCAR/CONTCAR/*.vasp, extxyz, ...), searched
  recursively; every structure whose element set is a sub-system is used;
* ``mp_fetch``: an optional callable ``subsystem -> list[structure dict]`` (e.g.
  :func:`mp_reference_fetcher`, Materials Project ground states at 0 GPa, re-relaxed here at P);
* an elemental fallback (``ase.build.bulk`` reference states; a crude dimer crystal for
  diatomic elements such as H, N, O) when no elemental reference was supplied.  It is recorded.

Each reference is relaxed with the prescreener's relax function at the target pressure and its
enthalpy per atom H = (E + P V)/N is cached per (sub-system, pressure, model) in ``cache_dir``.

Direction of errors: missing high-pressure polymorphs or crude fallbacks put the hull too HIGH,
so e_above_hull_ml comes out too LOW (optimistic: fewer deferrals).  That is the harmless direction
for a heuristic that may only defer.
"""
from __future__ import annotations

import itertools
import json
import os
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Iterable

from ..stability.hull import HullEntry, e_above_hull
from ._atoms import atoms_to_structure, composition, structure_to_atoms

# relax_fn(atoms, pressure_GPa) -> (relaxed_atoms, enthalpy_eV_per_atom, info dict)
RelaxFn = Callable[[Any, float], tuple[Any, float, dict]]
STRUCT_EXT = (".cif", ".vasp", ".poscar", ".extxyz", ".xyz", ".json", ".res")


def subsystems(elements: Iterable[str]) -> list[tuple[str, ...]]:
    els = sorted(set(elements))
    return [c for k in range(1, len(els) + 1) for c in itertools.combinations(els, k)]


def _p_key(p: float) -> str:
    return f"{float(p):.4g}".replace("+", "")


def fallback_element(symbol: str):
    """Elemental reference when none was supplied (recorded as a fallback)."""
    from ase import Atoms
    from ase.build import bulk
    from ase.data import atomic_numbers, covalent_radii, reference_states
    z = atomic_numbers[symbol]
    ref = reference_states[z] or {}
    if ref.get("symmetry") == "diatom":
        d = float(ref.get("d", 2 * covalent_radii[z]))
        a = max(3.0 * d, 2.5)
        return Atoms(symbol * 2, positions=[[0, 0, 0], [d / 3**0.5] * 3], cell=[a, a, a], pbc=True), \
            "fallback:dimer-sc"
    try:
        return bulk(symbol), f"fallback:ase.build.bulk({ref.get('symmetry')})"
    except Exception:  # noqa: BLE001 - unusual reference state; use fcc at covalent density
        a = 2 * 2**0.5 * covalent_radii[z]
        return bulk(symbol, "fcc", a=a), "fallback:fcc(covalent radius)"


def read_reference_dir(path: str | Path) -> list[tuple[str, Any]]:
    """All periodic structures in ``path`` (recursive) as (label, ase.Atoms)."""
    from ase.io import read
    out = []
    for f in sorted(Path(path).rglob("*")):
        if not f.is_file():
            continue
        name = f.name.lower()
        if not (name.endswith(STRUCT_EXT) or name.startswith(("poscar", "contcar"))):
            continue
        try:
            frames = read(str(f), index=":")
        except Exception:  # noqa: BLE001 - not a structure file
            continue
        for k, a in enumerate(frames if isinstance(frames, list) else [frames]):
            if all(a.pbc):
                out.append((f"{f.name}#{k}", a))
    return out


def mp_reference_fetcher(e_hull_max: float | None = 0.05, api_key_env: str = "MP_API_KEY",
                         rester_factory=None) -> Callable[[tuple[str, ...]], list[tuple[str, dict]]]:
    """Materials Project structures for one exact sub-system (0 GPa data, re-relaxed at P)."""
    from ..sources.mp import query_materials_project

    def fetch(subsys: tuple[str, ...]) -> list[tuple[str, dict]]:
        cands = query_materials_project(chemsys="-".join(subsys), energy_above_hull_max=e_hull_max,
                                        api_key_env=api_key_env, rester_factory=rester_factory)
        return [(f"mp:{c.material_id}", {"lattice": c.lattice, "species": c.species,
                                         "frac_coords": c.frac_coords})
                for c in cands if c.lattice and c.species and c.frac_coords]
    return fetch


@dataclass
class MLHull:
    relax_fn: RelaxFn
    model_id: str                               # part of the cache key (never mix potentials)
    reference_dir: str | None = None
    mp_fetch: Callable[[tuple[str, ...]], list[tuple[str, dict]]] | None = None
    cache_dir: str | None = None
    elemental_fallback: bool = True
    max_reference_atoms: int = 64               # skip larger reference cells (cost)
    _mem: dict[tuple[tuple[str, ...], str], list[dict]] = field(default_factory=dict, repr=False)
    _dir_refs: list[tuple[str, Any]] | None = field(default=None, repr=False)

    # ------------------------------------------------------------------ cache
    def _cache_path(self, subsys: tuple[str, ...], p: float) -> Path | None:
        if not self.cache_dir:
            return None
        mid = re.sub(r"[^A-Za-z0-9_.-]+", "_", self.model_id)
        return Path(self.cache_dir) / f"{'-'.join(subsys)}_P{_p_key(p)}GPa_{mid}.json"

    def _dir_structures(self) -> list[tuple[str, Any]]:
        if self._dir_refs is None:
            self._dir_refs = read_reference_dir(self.reference_dir) if self.reference_dir else []
        return self._dir_refs

    def entries_for_subsystem(self, subsys: tuple[str, ...], pressure_GPa: float) -> list[dict]:
        key = (subsys, _p_key(pressure_GPa))
        if key in self._mem:
            return self._mem[key]
        cp = self._cache_path(subsys, pressure_GPa)
        if cp is not None and cp.exists():
            data = json.loads(cp.read_text())
            if data.get("model_id") == self.model_id:
                self._mem[key] = data["entries"]
                return data["entries"]
        todo: list[tuple[str, Any]] = []
        for label, a in self._dir_structures():
            if tuple(sorted(set(a.get_chemical_symbols()))) == subsys:
                todo.append((f"dir:{label}", a))
        if self.mp_fetch is not None:
            for label, s in self.mp_fetch(subsys):
                todo.append((label, structure_to_atoms(s)))
        if not todo and len(subsys) == 1 and self.elemental_fallback:
            a, how = fallback_element(subsys[0])
            todo.append((how, a))
        entries = []
        for label, a in todo:
            if len(a) > self.max_reference_atoms:
                continue
            try:
                ra, h, info = self.relax_fn(a, pressure_GPa)
            except Exception as exc:  # noqa: BLE001 - a failed reference is skipped, recorded
                entries.append({"label": label, "failed": f"{type(exc).__name__}: {exc}"})
                continue
            entries.append({"label": label, "composition": composition(ra.get_chemical_symbols()),
                            "enthalpy_eV_per_atom": float(h), "converged": bool(info.get("converged")),
                            "structure": atoms_to_structure(ra)})
        self._mem[key] = entries
        if cp is not None:
            cp.parent.mkdir(parents=True, exist_ok=True)
            tmp = cp.with_suffix(".tmp")
            tmp.write_text(json.dumps({"model_id": self.model_id, "pressure_GPa": pressure_GPa,
                                       "subsystem": list(subsys), "entries": entries}))
            os.replace(tmp, cp)
        return entries

    def reference_entries(self, elements: Iterable[str], pressure_GPa: float) -> list[HullEntry]:
        out = []
        for sub in subsystems(elements):
            for e in self.entries_for_subsystem(sub, pressure_GPa):
                if "enthalpy_eV_per_atom" in e:
                    out.append(HullEntry(dict(e["composition"]), float(e["enthalpy_eV_per_atom"]),
                                         float(pressure_GPa), e["label"]))
        return out

    def e_above_hull(self, species: list[str], enthalpy_per_atom: float,
                     pressure_GPa: float) -> tuple[float | None, str, list[str]]:
        """(e_above_hull eV/atom or None, method/reason, labels of reference phases used)."""
        refs = self.reference_entries(species, pressure_GPa)
        val, how = e_above_hull(composition(species), float(enthalpy_per_atom), refs)
        return val, how, [r.label for r in refs]
