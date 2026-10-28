"""ML force-field pre-screen for the discovery loop (DISCOVERY.md 2.3, 4).  HEURISTIC ONLY (ARCH S2).

``rtsc.mlff.mattersim``: MatterSimPrescreener (relax at pressure, ML hull, ML phonons).
``rtsc.mlff.hull``: reference phases relaxed with the same potential at the same pressure.
``rtsc.mlff.phonons``: phonopy finite displacements + DOS summary + light-element heuristic.
ASE, phonopy, torch and mattersim are imported lazily; importing this package needs none of them.
"""
from .mattersim import MatterSimConfig, MatterSimPrescreener, RelaxSettings
from .phonons import PhononSettings

__all__ = ["MatterSimPrescreener", "MatterSimConfig", "RelaxSettings", "PhononSettings"]
