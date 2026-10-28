"""SYNTHETIC EPW 6.1 phself fixtures: byte-level imitations of the EPW writers (no physics).

Formats reproduced from EPW/src/io/io_selfen.f90 (selfen_ph_write) and EPW/src/selfen.f90
(selfen_phon_q stdout FORMATs 100-105 and the per-q header), EPW/src/utilities/bzgrid.f90
(qwindow_wrap 'We only need to compute').  Numbers are invented; every file name carries SYNTHETIC
through its directory (tests/fixtures/SYNTHETIC_<case>/).

Regenerate the committed fixtures with:  python tests/synthetic_epw.py
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from ephproto.qlist import (filqf_text, min_image_qnorm, reciprocal_vectors,  # noqa: E402
                            shifted_mp_points)

RYD2MEV = 13605.693122994
TEMP = "10.000"            # f8.3 of temps(1) = 10 K, adjustl'd
LAMBDA_FILE = f"lambda.phself.{TEMP}K"
LINEWIDTH_FILE = f"linewidth.phself.{TEMP}K"
QFILE = "qlist.txt"
AT_FCC = 7.65 / 2 * np.array([[-1.0, 0.0, 1.0], [0.0, 1.0, 1.0], [-1.0, 1.0, 0.0]])  # QE ibrav=2


# ----------------------------------------------------------------------------- Fortran output
def fE(x: float, w: int, d: int) -> str:
    """gfortran Ew.d: 0.ddddE+ee; |exp| in 100..999 printed without the letter (0.dddd-100)."""
    x = float(x)
    if x == 0.0:
        s = "0." + "0" * d + "E+00"
        return ("-" + s if math.copysign(1, x) < 0 else s).rjust(w)
    sign = "-" if x < 0 else ""
    m, e = f"{abs(x):.{d - 1}e}".split("e")
    digits = m.replace(".", "")
    exp = int(e) + 1
    if abs(exp) <= 99:
        es = f"E{exp:+03d}"
    elif abs(exp) <= 999:
        es = f"{exp:+04d}"
    else:
        raise ValueError("exponent too large for Ew.d")
    s = f"{sign}0.{digits}{es}"
    if len(s) > w:
        return "*" * w
    return s.rjust(w)


def fF(x: float, w: int, d: int) -> str:
    s = f"{float(x):{w}.{d}f}"
    return "*" * w if len(s) > w else s


def f32(x: float) -> float:
    return float(np.float32(x))


# ----------------------------------------------------------------------------- data model
def model(grid=(2, 2, 2)):
    """Invented omega (meV, ascending, modes 1-2 degenerate, symmetric points share omega) and
    lambda (double) for a 1-atom cell (nmodes = 3)."""
    pts = shifted_mp_points(grid)
    bg = reciprocal_vectors(AT_FCC)
    qn = np.round(min_image_qnorm(pts, bg), 6)
    qmax = float(qn.max())
    om, lam = [], []
    for q in qn:
        s = q / qmax
        t = 18.0 + 9.0 * s
        lo = 31.0 + 4.0 * s
        om.append([t, t, lo])
        lam.append([0.11 + 0.03 * s, 0.11 + 0.03 * s, 0.17 + 0.05 * s * s])
    return pts, np.array(om), np.array(lam)


def lambda_text(lam_print: np.ndarray, iqs=None) -> str:
    nq, nm = lam_print.shape
    out = ["", "  #Lambda phonon self-energy", ""]
    out.append(" #Modes     " + "".join(f"{i:12d}" for i in range(1, nm + 1)))
    for r in range(nq):
        iq = r + 1 if iqs is None else iqs[r]
        out.append(f"{iq:9d}    " + "".join("   " + fE(v, 15, 5) for v in lam_print[r]))
    return "\n".join(out) + "\n"


def linewidth_text(om: np.ndarray, gam2: np.ndarray, iqs=None) -> str:
    out = ["# Phonon frequency and phonon lifetime in meV ",
           "# Q-point  Mode   Phonon freq (meV)   Phonon linewidth (meV)"]
    nq, nm = om.shape
    for r in range(nq):
        iq = r + 1 if iqs is None else iqs[r]
        for j in range(nm):
            out.append(f"{iq:9d}{j + 1:6d}{fE(om[r, j], 20, 8)}{fE(gam2[r, j], 22, 10)}")
    return "\n".join(out) + "\n"


def stdout_text(pts, om_stdout, lam, sigma_eV=0.02, need_q=None, blocks=None, restart_from=None,
                nimage=None) -> str:
    nq, nm = lam.shape
    o = ["", "     Program EPW v.6.1 starts on 27Oct2026 at 12: 0: 0 (SYNTHETIC)", ""]
    if nimage:
        o.append(f"     path-images division:  nimage    = {nimage:7d}")
    o.append(f"     Using q-mesh file: {QFILE}")
    o.append(f"     We only need to compute {need_q if need_q is not None else nq:8d} q-points")
    if restart_from:
        o.append(f"     We restart from {restart_from:8d} q-points")
    o.append("")
    o.append("     " + "=" * 67)
    o.append("     Phonon (Imaginary) Self-Energy in the Migdal Approximation")
    o.append("     " + "=" * 67)
    o.append("")
    o.append(f"     Fermi Surface thickness = {fF(1.0, 10, 6)} eV")
    o.append(f"     Golden Rule strictly enforced with T = {fF(10 * 8.617333262e-5, 10, 6)} eV")
    first = True
    for r in range(nq):
        iq = r + 1
        if blocks is not None and iq not in blocks:
            continue
        if first:
            o.append(f"     Gaussian Broadening: {fF(sigma_eV, 10, 6)} eV, ngauss={0:4d}")
            o.append(f"     DOS ={fF(0.391234, 10, 6)} states/spin/eV/Unit Cell at Ef={fF(8.123456, 10, 6)} eV")
            first = False
        x = pts[r]
        o.append("")
        o.append(f"     ismear = {1:5d} iq = {iq:7d} coord.: {fF(x[0], 9, 5)}{fF(x[1], 9, 5)}"
                 f"{fF(x[2], 9, 5)} wt: {fF(1.0 / nq, 9, 5)} Temp: {fF(10.0, 8, 3)}K")
        o.append("     " + "-" * 67)
        for j in range(nm):
            w = om_stdout[r, j]
            o.append(f"     lambda___( {j + 1:3d} )={fF(lam[r, j], 15, 6)}   gamma___={fF(0.5 * lam[r, j] * w, 15, 6)}"
                     f" meV  omega={fF(w, 12, 4)} meV   on_shell_Omega={fF(w, 12, 4)} meV"
                     f"  off_shell_Omega={fF(w, 12, 4)} meV")
            o.append(f"     lambda_tr( {j + 1:3d} )={fF(0.0, 15, 6)}   gamma_tr={fF(0.0, 15, 6)} meV"
                     f"  omega={fF(w, 12, 4)} meV   on_shell_Omega={fF(w, 12, 4)} meV"
                     f"  off_shell_Omega={fF(w, 12, 4)} meV")
        o.append(f"     lambda___( tot )={fF(math.fsum(lam[r]), 15, 6)}")
        o.append(f"     lambda_tr( tot )={fF(0.0, 15, 6)}")
        o.append("     " + "-" * 67)
        o.append("")
        o.append(f"     Number of (k,k+q) pairs on the Fermi surface: {1234:8d} out of {4096:8d}")
    o.append("")
    o.append("     EPW          :      1.00s CPU      1.00s WALL")
    return "\n".join(o) + "\n"


def epw_in_text(sigma_eV=0.02, extra: dict | None = None, drop: tuple = ()) -> str:
    kv = {"prefix": "'synthetic'", "outdir": "'./out/'", "elph": ".true.", "epwread": ".true.",
          "epwwrite": ".false.", "wannierize": ".false.", "nbndsub": "4", "use_ws": ".true.",
          "phonselfen": ".true.", "delta_approx": ".true.", "a2f": ".false.",
          "restart": ".false.", "nstemp": "1", "temps(1)": "10.0", "ngaussw": "0",
          "degaussw": repr(sigma_eV), "fsthick": "1.0d0", "efermi_read": ".false.",
          "selecqread": ".false.", "eps_acoustic": "0.1", "filqf": f"'{QFILE}'",
          "nkf1": "24", "nkf2": "24", "nkf3": "24",
          "nk1": "8", "nk2": "8", "nk3": "8", "nq1": "4", "nq2": "4", "nq3": "4"}
    for k in drop:
        kv.pop(k, None)
    kv.update(extra or {})
    lines = ["--", "&inputepw", "  ! SYNTHETIC fixture input (protocol v3 pins)"]
    lines += [f"  {k} = {v}" for k, v in kv.items()]
    lines += ["/", ""]
    return "\n".join(lines)


# ----------------------------------------------------------------------------- cases
CASES = ("normal", "restart_zeros", "exp_noE", "negative_lambda", "omega_nonpos", "gamma_present",
         "wrong_nq")


def write_case(d: Path, case: str) -> Path:
    d = Path(d)
    d.mkdir(parents=True, exist_ok=True)
    grid = (2, 2, 2)
    pts, om, lam = model(grid)
    nq, nm = lam.shape
    lam_print = np.vectorize(f32)(lam)
    gam2 = 2.0 * 0.5 * lam * om * 1e-3          # invented FWHM (meV); never used
    om_std = om.copy()
    qtext = filqf_text(grid)
    kw: dict = {}
    if case == "restart_zeros":
        lam = lam.copy()
        lam[:3] = 0.0
        lam_print[:3] = 0.0
        gam2[:3] = 0.0
        kw = {"blocks": set(range(4, nq + 1)), "restart_from": 4}
    elif case == "exp_noE":
        lam = lam.copy()
        lam[1, 0] = 1.5e-120          # not representable in single precision: parser test only
        lam_print[1, 0] = 1.5e-120
        gam2[1, 0] = 2.25e-121
        gam2[2, 2] = 7.123456789e-105
    elif case == "negative_lambda":
        lam = lam.copy()
        lam[4, 2] = -1.2e-3
        lam_print[4, 2] = f32(-1.2e-3)
    elif case == "omega_nonpos":
        om = om.copy()
        om[4, 0] = -0.5               # imaginary mode, EPW prints -sqrt|w2| and zero lambda
        om_std = om.copy()
        lam = lam.copy()
        lam[4, 0] = 0.0
        lam_print[4, 0] = 0.0
    elif case == "gamma_present":
        # unshifted 2x2x2 grid: point 1 is Gamma, acoustic omega = 0, lambda = 0
        g = [(i / 2, j / 2, k / 2) for i in range(2) for j in range(2) for k in range(2)]
        qtext = f"{len(g)} crystal\n" + "".join(f"{a!r} {b!r} {c!r} {1 / len(g)!r}\n" for a, b, c in g)
        pts = np.array(g)
        om = om.copy()
        om[0] = 0.0
        om_std = om.copy()
        lam = lam.copy()
        lam[0] = 0.0
        lam_print[0] = 0.0
    elif case == "wrong_nq":
        qtext = filqf_text((3, 3, 3))
    (d / QFILE).write_text(qtext, newline="\n")
    (d / "epw.in").write_text(epw_in_text(), newline="\n")
    (d / LAMBDA_FILE).write_text(lambda_text(lam_print), newline="\n")
    (d / LINEWIDTH_FILE).write_text(linewidth_text(om, gam2), newline="\n")
    (d / "epw.out").write_text(stdout_text(pts, om_std, lam, **kw), newline="\n")
    return d


def fixture_dir(case: str) -> Path:
    return HERE / "fixtures" / f"SYNTHETIC_{case}"


if __name__ == "__main__":
    for c in CASES:
        print(write_case(fixture_dir(c), c))
