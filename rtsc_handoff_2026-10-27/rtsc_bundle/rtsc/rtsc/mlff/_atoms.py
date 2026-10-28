"""Conversions between the discovery structure dict and ASE ``Atoms`` (ASE imported lazily)."""
from __future__ import annotations

import hashlib
import json
from typing import Any


def structure_to_atoms(s: dict[str, Any]):
    """``{"lattice", "species", "frac_coords"}`` -> periodic ``ase.Atoms``."""
    from ase import Atoms
    return Atoms(symbols=list(s["species"]), scaled_positions=[list(map(float, x)) for x in s["frac_coords"]],
                 cell=[list(map(float, r)) for r in s["lattice"]], pbc=True)


def atoms_to_structure(atoms, **extra: Any) -> dict[str, Any]:
    out: dict[str, Any] = {
        "lattice": [[float(x) for x in row] for row in atoms.cell[:]],
        "species": list(atoms.get_chemical_symbols()),
        "frac_coords": [[float(x) % 1.0 for x in row] for row in atoms.get_scaled_positions(wrap=True)],
    }
    out.update(extra)
    return out


def fallback_fingerprint(s: dict[str, Any]) -> str:
    """Only used when the loop did not supply ``provenance['fingerprint']`` (it normally does)."""
    payload = json.dumps({"l": [[round(x, 4) for x in r] for r in s["lattice"]],
                          "s": list(s["species"]),
                          "f": [[round(x, 4) for x in r] for r in s["frac_coords"]]}, sort_keys=True)
    return "raw-" + hashlib.sha256(payload.encode()).hexdigest()[:16]


def composition(species: list[str]) -> dict[str, float]:
    out: dict[str, float] = {}
    for el in species:
        out[el] = out.get(el, 0.0) + 1.0
    return out
