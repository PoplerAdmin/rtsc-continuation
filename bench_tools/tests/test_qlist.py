import itertools
import math
from fractions import Fraction

import numpy as np
import pytest

from ephproto.qlist import (QListError, g_min_norm, identify_shifted_grid, min_image_qnorm,
                            read_filqf, reciprocal_vectors, shifted_mp_fractions, shifted_mp_points,
                            write_filqf)


@pytest.mark.parametrize("n", [(1, 1, 1), (2, 2, 2), (3, 3, 3), (4, 4, 4), (2, 3, 5), (6, 6, 6)])
def test_roundtrip_exact_and_gamma_free(tmp_path, n):
    meta = write_filqf(tmp_path / "q.txt", n)
    N = n[0] * n[1] * n[2]
    assert meta["N_q"] == N and len(meta["sha256"]) == 64
    pts, w, ctype = read_filqf(tmp_path / "q.txt")
    assert ctype == "crystal" and pts.shape == (N, 3)
    exact = shifted_mp_fractions(n)
    # coordinates and weights come back as the exact nearest doubles
    assert all(pts[i, k] == float(exact[i][k]) for i in range(N) for k in range(3))
    assert np.all(w == 1.0 / N)
    assert identify_shifted_grid(pts, w, ctype) == tuple(n)
    # no point is a reciprocal lattice vector (no Gamma) -- exact arithmetic
    assert all(any(c.denominator != 1 for c in p) for p in exact)
    # EPW order: iq = i*n2*n3 + j*n3 + k + 1
    n1, n2, n3 = n
    i, j, k = n1 - 1, 0, n3 - 1
    assert exact[i * n2 * n3 + j * n3 + k] == (Fraction(2 * i + 1, 2 * n1), Fraction(1, 2 * n2),
                                               Fraction(2 * k + 1, 2 * n3))


def test_header_line(tmp_path):
    write_filqf(tmp_path / "q.txt", (2, 2, 2))
    lines = (tmp_path / "q.txt").read_text().splitlines()
    assert lines[0] == "8 crystal" and lines[1] == "0.25 0.25 0.25 0.125" and len(lines) == 9


def test_identify_rejects_gamma_grid():
    g = np.array([(i / 2, j / 2, k / 2) for i in range(2) for j in range(2) for k in range(2)])
    with pytest.raises(QListError, match="Gamma"):
        identify_shifted_grid(g, np.full(8, 1 / 8))


def test_identify_rejects_reordered_nonuniform_and_partial():
    pts = shifted_mp_points((2, 2, 2))
    w = np.full(8, 1 / 8)
    with pytest.raises(QListError, match="EPW order"):
        identify_shifted_grid(pts[::-1], w)
    w2 = w.copy()
    w2[0], w2[1] = 0.25, 0.0
    with pytest.raises(QListError, match="uniform"):
        identify_shifted_grid(pts, w2)
    with pytest.raises(QListError, match="product grid"):
        identify_shifted_grid(pts[:7], w[:7])
    with pytest.raises(QListError, match="crystal"):
        identify_shifted_grid(pts, w, "cartesian")


def test_read_count_mismatch(tmp_path):
    write_filqf(tmp_path / "q.txt", (2, 2, 2))
    txt = (tmp_path / "q.txt").read_text().splitlines()
    (tmp_path / "bad.txt").write_text("\n".join(txt[:-1]) + "\n")
    with pytest.raises(QListError, match="header says 8"):
        read_filqf(tmp_path / "bad.txt")


def test_min_image_simple_cubic():
    bg = reciprocal_vectors(np.eye(3))            # 2 pi I
    q = np.array([[0.75, 0.0, 0.0], [0.5, 0.5, 0.5], [0.1, -0.9, 0.0], [1.25, 0, 0]])
    got = min_image_qnorm(q, bg)
    exp = 2 * math.pi * np.array([0.25, math.sqrt(3) / 2, math.sqrt(0.02), 0.25])
    assert np.allclose(got, exp, rtol=1e-14, atol=0)


def test_min_image_fcc_against_bruteforce():
    a = 7.65
    at = a / 2 * np.array([[-1.0, 0.0, 1.0], [0.0, 1.0, 1.0], [-1.0, 1.0, 0.0]])
    bg = reciprocal_vectors(at)
    assert np.allclose(at @ bg.T, 2 * math.pi * np.eye(3))
    pts = shifted_mp_points((4, 4, 4))
    got = min_image_qnorm(pts, bg)
    shifts = np.array(list(itertools.product(range(-4, 5), repeat=3)), dtype=float)
    brute = np.array([np.min(np.linalg.norm((p - shifts) @ bg, axis=1)) for p in pts])
    assert np.allclose(got, brute, rtol=1e-13, atol=0)
    # shortest G of fcc is |(1,1,1)| 2pi/a = sqrt(3) 2pi/a
    assert math.isclose(g_min_norm(bg), math.sqrt(3) * 2 * math.pi / a, rel_tol=1e-13)
