"""EPW ``filqf`` files for a shifted uniform Monkhorst-Pack q grid (protocol v3 decision).

EPW 6.1 reads ``filqf`` in ``loadqmesh_serial`` (utilities/bzgrid.f90) as::

    READ(iunqf, *) nqtotf, coordinate_type          ! e.g. "64 crystal"
    DO iq = 1, nqtotf
      READ(iunqf, *) xqf(:, iq), wqf(iq)
    ENDDO

``crystal`` means components along the reciprocal basis b1, b2, b3 (QE ``bg``).  The weights are
stored but NOT used by ``selfen_phon_q`` (lambda_qnu does not depend on wqf); the protocol applies
the uniform weight 1/N_q itself, so the grid must be exactly uniform: that is what
``identify_shifted_grid`` verifies.

Grid: x_i = (i + 1/2)/n_1, y_j = (j + 1/2)/n_2, z_k = (k + 1/2)/n_3, i, j, k = 0..n-1, ordered with
i slowest and k fastest (EPW's own uniform-grid order, iq = i*n2*n3 + j*n3 + k + 1).  A half shift
never contains Gamma for any n (i + 1/2 is never a multiple of n).
"""
from __future__ import annotations

import hashlib
import itertools
import math
from fractions import Fraction
from pathlib import Path
from typing import Sequence

import numpy as np


class QListError(ValueError):
    """filqf content is not the exact shifted uniform grid the protocol needs."""


def _grid(n: Sequence[int]) -> tuple[int, int, int]:
    if len(n) != 3:
        raise ValueError(f"need three grid sizes, got {n!r}")
    out = tuple(int(x) for x in n)
    if any(x < 1 or x != y for x, y in zip(out, n)):
        raise ValueError(f"grid sizes must be positive integers, got {n!r}")
    return out  # type: ignore[return-value]


def shifted_mp_fractions(n: Sequence[int]) -> list[tuple[Fraction, Fraction, Fraction]]:
    """Exact crystal coordinates (Fractions) of the half-shifted n1 x n2 x n3 grid, EPW order."""
    n1, n2, n3 = _grid(n)
    return [(Fraction(2 * i + 1, 2 * n1), Fraction(2 * j + 1, 2 * n2), Fraction(2 * k + 1, 2 * n3))
            for i in range(n1) for j in range(n2) for k in range(n3)]


def shifted_mp_points(n: Sequence[int]) -> np.ndarray:
    """(N, 3) float array of crystal coordinates, EPW order (i slowest)."""
    return np.array([[float(c) for c in p] for p in shifted_mp_fractions(n)], dtype=float)


def _fmt(x: float) -> str:
    # repr is the shortest string that round-trips to the same IEEE double; Fortran list-directed
    # READ into REAL(DP) reproduces that double (correctly rounded decimal conversion).
    r = repr(float(x))
    return r


def filqf_text(n: Sequence[int]) -> str:
    """Text of the filqf file for the half-shifted grid ``n`` (weights 1/N, crystal)."""
    pts = shifted_mp_points(n)
    N = len(pts)
    w = 1.0 / N
    lines = [f"{N} crystal"]
    for p in pts:
        lines.append(f"{_fmt(p[0])} {_fmt(p[1])} {_fmt(p[2])} {_fmt(w)}")
    return "\n".join(lines) + "\n"


def sha256_file(path: str | Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_filqf(path: str | Path, n: Sequence[int]) -> dict:
    """Write the filqf file; return its metadata (grid, N_q, shift, sha256)."""
    n = _grid(n)
    p = Path(path)
    p.write_text(filqf_text(n), encoding="ascii", newline="\n")
    return {"path": str(p), "grid": list(n), "shift": [0.5, 0.5, 0.5], "N_q": n[0] * n[1] * n[2],
            "coordinates": "crystal", "weight": 1.0 / (n[0] * n[1] * n[2]),
            "sha256": sha256_file(p)}


def read_filqf(path: str | Path) -> tuple[np.ndarray, np.ndarray, str]:
    """Read a filqf file as EPW does: (points (N,3), weights (N,), coordinate_type).

    Raises on a count mismatch or trailing non-blank content (EPW itself would silently ignore
    extra lines and crash on missing ones)."""
    lines = Path(path).read_text().splitlines()
    body = [ln for ln in lines if ln.strip()]
    if not body:
        raise QListError(f"{path}: empty filqf")
    head = body[0].split()
    try:
        N = int(head[0])
    except (ValueError, IndexError):
        raise QListError(f"{path}: first line must be 'N crystal|cartesian', got {body[0]!r}")
    ctype = head[1].lower() if len(head) > 1 else "crystal"
    if ctype not in ("crystal", "cartesian"):
        raise QListError(f"{path}: coordinate type must be crystal or cartesian, got {ctype!r}")
    rows = body[1:]
    if len(rows) != N:
        raise QListError(f"{path}: header says {N} points, file has {len(rows)} point lines")
    pts = np.empty((N, 3))
    wts = np.empty(N)
    for i, ln in enumerate(rows):
        t = ln.replace(",", " ").split()
        if len(t) < 4:
            raise QListError(f"{path}: point line {i + 1} has {len(t)} fields, need 4: {ln!r}")
        try:
            vals = [float(x.replace("d", "e").replace("D", "e")) for x in t[:4]]
        except ValueError:
            raise QListError(f"{path}: unparsable point line {i + 1}: {ln!r}")
        pts[i] = vals[:3]
        wts[i] = vals[3]
    if not (np.all(np.isfinite(pts)) and np.all(np.isfinite(wts))):
        raise QListError(f"{path}: non-finite coordinate or weight")
    return pts, wts, ctype


def identify_shifted_grid(points: np.ndarray, weights: np.ndarray, coordinate_type: str = "crystal",
                          *, atol: float = 1e-10) -> tuple[int, int, int]:
    """Return (n1, n2, n3) iff the list is EXACTLY the half-shifted grid in EPW order with uniform
    weights 1/N (to ``atol``); raise QListError otherwise.  Gamma (or any lattice vector) in the list
    is reported explicitly."""
    if coordinate_type.lower() != "crystal":
        raise QListError("only crystal coordinates are accepted (protocol writes crystal)")
    pts = np.asarray(points, dtype=float)
    w = np.asarray(weights, dtype=float)
    N = len(pts)
    if N == 0:
        raise QListError("empty q list")
    red = pts - np.floor(pts + atol)            # into [0, 1)
    gam = np.all(np.minimum(np.abs(red), np.abs(red - 1.0)) <= atol, axis=1)
    if np.any(gam):
        raise QListError(f"q list contains Gamma (or a reciprocal lattice vector) at "
                         f"{int(np.count_nonzero(gam))} point(s), e.g. index {int(np.argmax(gam)) + 1}"
                         "; the protocol requires a half-shifted grid without Gamma")
    ns = []
    for d in range(3):
        vals = np.unique(np.round(red[:, d] / atol).astype(np.int64))
        ns.append(len(vals))
    n1, n2, n3 = ns
    if n1 * n2 * n3 != N:
        raise QListError(f"{N} points do not form a full product grid (distinct coordinates "
                         f"{n1} x {n2} x {n3} = {n1 * n2 * n3})")
    ref = shifted_mp_points((n1, n2, n3))
    if not np.allclose(pts, ref, rtol=0.0, atol=atol):
        bad = int(np.argmax(np.max(np.abs(pts - ref), axis=1) > atol)) + 1
        raise QListError(f"points are not the half-shifted {n1}x{n2}x{n3} grid in EPW order "
                         f"(first mismatch at index {bad}: {pts[bad - 1]} vs {ref[bad - 1]})")
    if not np.allclose(w, 1.0 / N, rtol=0.0, atol=atol / N):
        raise QListError(f"weights are not uniform 1/N = {1.0 / N!r} (min {w.min()!r}, max {w.max()!r})")
    if abs(math.fsum(w) - 1.0) > 1e-9:
        raise QListError(f"weights sum to {math.fsum(w)!r}, not 1")
    return n1, n2, n3


# ----------------------------------------------------------------------------- geometry
def reciprocal_vectors(at: np.ndarray) -> np.ndarray:
    """Rows b_i with a_i . b_j = 2 pi delta_ij, from direct lattice rows ``at`` (any length unit;
    result is in 2 pi / that unit)."""
    a = np.asarray(at, dtype=float).reshape(3, 3)
    return 2.0 * math.pi * np.linalg.inv(a).T


def g_min_norm(bg: np.ndarray, nmax: int = 3) -> float:
    """Length of the shortest non-zero reciprocal lattice vector (rows of ``bg``)."""
    b = np.asarray(bg, dtype=float).reshape(3, 3)
    best = math.inf
    for m in itertools.product(range(-nmax, nmax + 1), repeat=3):
        if m == (0, 0, 0):
            continue
        best = min(best, float(np.linalg.norm(np.asarray(m, dtype=float) @ b)))
    return best


def min_image_qnorm(points_crys: np.ndarray, bg: np.ndarray, nsearch: int = 2) -> np.ndarray:
    """Minimum-image cartesian |q| for crystal coordinates (rows), with reciprocal rows ``bg``:
    min over integer vectors G of |(x - G) . bg|.  The coordinate is first reduced to [-1/2, 1/2)
    and then the 5^3 neighbouring images (``nsearch`` = 2) are searched, which is exact for any
    cell whose reduced basis is not extremely skewed."""
    x = np.atleast_2d(np.asarray(points_crys, dtype=float))
    b = np.asarray(bg, dtype=float).reshape(3, 3)
    red = x - np.floor(x + 0.5)
    shifts = np.array(list(itertools.product(range(-nsearch, nsearch + 1), repeat=3)), dtype=float)
    out = np.empty(len(red))
    for i, r in enumerate(red):
        cand = (r[None, :] - shifts) @ b
        out[i] = float(np.min(np.linalg.norm(cand, axis=1)))
    return out
