"""Autonomous discovery loop (DISCOVERY.md): generate -> ML pre-screen -> dedupe -> DFT -> certify.

Soundness (ARCH S1/S2): only ``rtsc.certify.classify`` (reached through the existing pipeline) may
produce REJECTED_PROVEN / CERTIFIED_CANDIDATE.  Generators, ML pre-screens, the search policy and
budgets only order candidates or mark them DEFERRED.
"""
from .interfaces import (GeneratedStructure, GenerationRequest, Generator, Prescreener,
                         PrescreenResult)

__all__ = ["GenerationRequest", "GeneratedStructure", "PrescreenResult", "Generator",
           "Prescreener"]
