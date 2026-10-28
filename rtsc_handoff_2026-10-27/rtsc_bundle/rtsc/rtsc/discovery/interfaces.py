"""Protocol interfaces between the discovery loop and generator / ML pre-screen adapters.

These are the ONLY types adapters in ``rtsc.generate`` and ``rtsc.mlff`` need to implement.

Structure dict format (the pipeline's light-weight structure triple, cf. ``rtsc.sources.base.Candidate``)::

    {"lattice": [[ax, ay, az], [bx, by, bz], [cx, cy, cz]],   # angstrom, rows are lattice vectors
     "species": ["La", "H", ...],                              # one element symbol per site
     "frac_coords": [[x, y, z], ...]}                          # fractional coordinates per site

Optional extra keys (ignored by the core, kept in the ledger): ``"spacegroup"`` (str),
``"pressure_GPa"`` (float), ``"enthalpy_eV_per_atom"`` (float).

SOUNDNESS (ARCH S2, DISCOVERY.md section 2.3): everything an adapter returns is HEURISTIC.  The loop
uses a ``PrescreenResult`` only to (a) reorder candidates and (b) mark a material DEFERRED with a
recorded reason.  No adapter output can ever produce REJECTED_PROVEN, CERTIFIED_CANDIDATE or
EXCLUDED_ASSUMPTION; adapters have no API to set a status and must not try.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Protocol, runtime_checkable


@dataclass(frozen=True)
class GenerationRequest:
    """One request issued by the search policy to a generator."""
    chemical_system: tuple[str, ...]      # sorted element symbols, e.g. ("H", "La")
    space_group: int | None               # international number 1..230, or None = unconstrained
    pressure_GPa: float                   # target pressure for the structures
    goal_id: str                          # e.g. "G-RT", "G-AMB"
    n: int                                # number of structures requested (generator may return fewer)
    seed: int                             # RNG seed; same request + seed => same output if possible
    extra: dict[str, Any] = field(default_factory=dict)   # generator-specific conditioning

    @property
    def request_id(self) -> str:
        sg = self.space_group if self.space_group is not None else "any"
        return (f"{self.goal_id}|{'-'.join(self.chemical_system)}|sg{sg}|"
                f"P{self.pressure_GPa:g}|n{self.n}|s{self.seed}")


@dataclass
class GeneratedStructure:
    """One structure produced by a generator.  Always HYPOTHETICAL_STRUCTURE."""
    structure: dict[str, Any]             # {"lattice", "species", "frac_coords"} (see module doc)
    generator: str                        # Generator.name
    request_id: str                       # GenerationRequest.request_id
    provenance: dict[str, Any] = field(default_factory=dict)   # checkpoint, version, cmdline, ...


@dataclass
class PrescreenResult:
    """Heuristic ML pre-screen of one structure (S2: ordering / DEFERRED only).

    ``fingerprint`` must echo the fingerprint the loop passed in via
    ``GeneratedStructure.provenance["fingerprint"]`` (the loop sets it before calling prescreen);
    results are matched back by it, and by position as a fallback.
    """
    fingerprint: str
    relaxed_structure: dict[str, Any] | None = None
    e_above_hull_ml_eV: float | None = None      # eV/atom vs reference phases relaxed with same MLFF
    dynamically_stable_ml: bool | None = None
    min_freq_cm1: float | None = None            # most negative / lowest ML phonon frequency
    phonon_dos: dict[str, Any] | None = None     # e.g. {"freq_cm1": [...], "dos": [...]}
    ok: bool = True                               # False => the ML step itself failed
    notes: list[str] = field(default_factory=list)


@runtime_checkable
class Generator(Protocol):
    name: str

    def generate(self, request: GenerationRequest) -> list[GeneratedStructure]: ...


@runtime_checkable
class Prescreener(Protocol):
    name: str

    def prescreen(self, structures: list[GeneratedStructure],
                  pressure_GPa: float) -> list[PrescreenResult]: ...


__all__ = ["GenerationRequest", "GeneratedStructure", "PrescreenResult", "Generator", "Prescreener"]
