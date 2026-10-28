"""Simple prototype-substitution enumerator (optional source).

Given prototype candidates (with structures) and a substitution table such as
``{"La": ["Y", "Ca", "Sc"], "H": ["H"]}``, produce every combination of substitutions of the
prototype's elements.  The generated structures are HYPOTHETICAL and UNRELAXED (lattice kept, or
isotropically rescaled by an optional per-element radius ratio); provenance records this.  They
carry no energies, so the hull filter cannot evaluate them until relaxed.
"""
from __future__ import annotations

import itertools
from typing import Iterable

from .base import Candidate, candidate_from_structure_arrays
from .composition import formula_string


def substitute(prototypes: Iterable[Candidate], substitutions: dict[str, list[str]], *,
               include_original: bool = False, max_per_prototype: int | None = None,
               id_prefix: str = "subst:") -> list[Candidate]:
    out: list[Candidate] = []
    for proto in prototypes:
        if proto.species is None or proto.lattice is None or proto.frac_coords is None:
            raise ValueError(f"prototype {proto.material_id} has no ordered structure")
        els = sorted(set(proto.species))
        choices = [substitutions.get(e, [e]) for e in els]
        n = 0
        for combo in itertools.product(*choices):
            mapping = dict(zip(els, combo))
            if len(set(combo)) < len(combo):
                continue  # two sites collapsing to the same element changes the prototype
            if not include_original and all(mapping[e] == e for e in els):
                continue
            species = [mapping[s] for s in proto.species]
            comp: dict[str, float] = {}
            for s in species:
                comp[s] = comp.get(s, 0.0) + 1.0
            label = formula_string(comp)
            mid = f"{id_prefix}{label}@{proto.material_id}"
            c = candidate_from_structure_arrays(
                mid, proto.lattice, species, proto.frac_coords, source="prototype_substitution",
                spacegroup=proto.spacegroup, pressure_GPa=proto.pressure_GPa,
                provenance={"source": "prototype_substitution", "prototype": proto.material_id,
                            "mapping": mapping,
                            "note": "hypothetical, unrelaxed structure (prototype lattice kept)"})
            out.append(c)
            n += 1
            if max_per_prototype is not None and n >= max_per_prototype:
                break
    return out
