"""Data-path pieces of the alpha^2F convergence protocol (reports/EPH_CONVERGENCE_PROTOCOL.md, v3 decisions).

Modules
-------
qlist   : write/read the EPW ``filqf`` file for a shifted Monkhorst-Pack q grid; minimum-image |q|.
phself  : parse EPW ``lambda.phself.*`` / ``linewidth.phself.*`` (+ stdout, epw.in) and build the
          rtsc ``Spectrum`` from modes (atom construction, P1) with hard assertions.
tags    : provenance tags understood by rtsc as-is (see tags.py for the exact strings), and the
          protocol ERROR reasons UNCONVERGED_EPH(<gate>) / EPH_PIN_VIOLATION.
admit   : P3 Richardson admission (k, sigma, q), the G4 tetrahedron-anchor comparison and the P5
          envelope, on observables of already-parsed spectra (pure Python, no DFT).

Nothing here runs DFT; everything here raises (never warns) on data that could make a spectrum
silently wrong.
"""
from .qlist import (shifted_mp_points, write_filqf, read_filqf, identify_shifted_grid,  # noqa: F401
                    reciprocal_vectors, min_image_qnorm, g_min_norm)
from .tags import (EPH_PROTOCOL_TAG, drop_nonpos_tag, completeness_unverified_tag,  # noqa: F401
                   incomplete_tag, EPH_PIN_VIOLATION, UNCONVERGED_EPH_GATES, EPH_ERROR_TAG_RE,
                   unconverged_eph_tag)

__all__ = ["qlist", "phself", "tags", "admit"]
