"""Provenance tags for protocol-produced spectra, written in the exact form rtsc honours today.

rtsc as-is (rtsc/pipeline.py ``blocking_tags``, applied by ``apply_input_policy`` inside
``certify_spectra``) reads ``Spectrum.provenance["assumptions"]`` and blocks a proved status only for:

  * prefix ``A-FILE-INCOMPLETE``             -> blocks REJECTED and CERTIFIED
  * prefix ``A-FILE-COMPLETE(unverified)``   -> blocks REJECTED
  * regex ``^A-DROP-NONPOS\\(([-+0-9.eE]+)K\\)`` with float(m) > 0 -> blocks REJECTED
  * regex ``^A-NEG-CLIP\\(([-+0-9.eE]+)K\\)``    with float(m) > 0 -> blocks CERTIFIED

Everything else (including any new name such as ``A-DROP-EXCLUDED`` and the protocol tag below) is
only copied into ``Verdict.assumptions`` by ``annotate_model``; it does not change the status.
Note the exact spelling: no space between the number and ``K`` (``A-DROP-NONPOS(5 K)`` does NOT
match and does NOT block), the tag must start the string, ``0K`` / negative masses do not block.
"""
from __future__ import annotations

import math
import re

#: informational, non-blocking; identifies spectra built by this package (protocol still DRAFT)
EPH_PROTOCOL_TAG = ("EPH_PROTOCOL_DRAFT(v3): alpha^2F atoms from EPW phself lambda_qnu, omega_qnu "
                    "(bench_tools.ephproto.phself); benchmark diagnostic, not a pipeline verdict")

# copied verbatim from rtsc/pipeline.py (_TAG_MASS) so tests can detect drift
RTSC_TAG_MASS_RE = re.compile(r"^A-(DROP-NONPOS|NEG-CLIP)\(([-+0-9.eE]+)K\)")


def drop_nonpos_tag(mass_K: float, what: str) -> str:
    """``A-DROP-NONPOS(<m>K): ...`` with m > 0 (blocks REJECTED in rtsc as-is).  Use only when the
    dropped mass is KNOWN (it is printed in reports/report.py as the dropped weight)."""
    m = float(mass_K)
    if not (math.isfinite(m) and m > 0):
        raise ValueError(f"A-DROP-NONPOS needs a finite mass > 0 to block in rtsc; got {mass_K!r} "
                         "(use completeness_unverified_tag for an unknown dropped mass)")
    tag = f"A-DROP-NONPOS({m:.6g}K): {what}"
    mm = RTSC_TAG_MASS_RE.match(tag)
    assert mm and float(mm.group(2)) > 0, tag
    return tag


def completeness_unverified_tag(what: str) -> str:
    """``A-FILE-COMPLETE(unverified): ...`` (blocks REJECTED only; no mass needed).  The existing
    rtsc string for 'weight may be missing, amount unknown'."""
    return f"A-FILE-COMPLETE(unverified): {what}"


def incomplete_tag(what: str) -> str:
    """``A-FILE-INCOMPLETE: ...`` (blocks both directions)."""
    return f"A-FILE-INCOMPLETE: {what}"


# ----------------------------------------------------------------------------- protocol ERROR reasons
# P7 / A3: a gate failure is ERROR ``UNCONVERGED_EPH(<gate>)``; a pin violation is its own reason
# ``EPH_PIN_VIOLATION``.  These are status REASONS (stage pipeline, origin eph_protocol,
# retryable 0), not spectrum assumptions.  Spelling rule (same lesson as ``A-DROP-NONPOS(5K)``):
# no whitespace anywhere, lower-case gate name inside the parentheses.
EPH_PIN_VIOLATION = "EPH_PIN_VIOLATION"
#: gate names used in UNCONVERGED_EPH(<gate>) (spec P1, P3, P4, P5, §B; 'tetra' = gate G4)
UNCONVERGED_EPH_GATES = ("k", "q", "sigma", "envelope", "tetra", "imaginary", "wannier")
EPH_ERROR_TAG_RE = re.compile(r"^(?:UNCONVERGED_EPH\((?:" + "|".join(UNCONVERGED_EPH_GATES)
                              + r")\)|EPH_PIN_VIOLATION)$")


def unconverged_eph_tag(gate: str) -> str:
    """``UNCONVERGED_EPH(<gate>)`` for a known gate name; raises ValueError otherwise."""
    if gate not in UNCONVERGED_EPH_GATES:
        raise ValueError(f"unknown UNCONVERGED_EPH gate {gate!r}; known: {UNCONVERGED_EPH_GATES}")
    tag = f"UNCONVERGED_EPH({gate})"
    assert EPH_ERROR_TAG_RE.match(tag) and not re.search(r"\s", tag), tag
    return tag
