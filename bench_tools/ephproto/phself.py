"""EPW 6.1 phonon self-energy files -> rtsc Spectrum by the atom construction (protocol P1, v3).

Writer (EPW/src/io/io_selfen.f90, ``selfen_ph_write``; called once after the q loop)::

    lambda.phself.<T>K     '(/2x,a/)' '#Lambda phonon self-energy'
                           list-directed  '#Modes     ', (imode, imode = 1, nmodes)
                           per q:  '(i9,4x)' iq  then  '(1000(3x,E15.5))' REAL(lambda_all(:, iq, 1, it))
    linewidth.phself.<T>K  two '#' lines, then per (q, mode)
                           '(i9,i6,E20.8,E22.10)' iq, imode, ryd2mev*wf(imode,iq), 2*ryd2mev*gamma

Construction (fine q grid = our half-shifted filqf list, weights exactly 1/N_q)::

    alpha^2F(w) = sum_{q,nu} (1/N_q) lambda_qnu omega_qnu / 2 * delta(w - omega_qnu)

omega from linewidth column 3 (meV), lambda from the lambda file.  The linewidth column 4 is
2*gamma (FWHM) and is never used.  Atoms with identical printed omega are merged exactly.

Every check here RAISES (PhselfError); nothing is silently dropped or clipped.

Writer facts that matter (EPW 6.1 source):
  * lambda is written as ``REAL(lambda_all)`` = default (single) precision, then E15.5: about
    5 significant digits, relative rounding <= 5e-5 per value; magnitudes below ~1.4e-45 print as 0.
    (So an exponent beyond +-99 cannot actually occur in the lambda file; the parser still handles
    the Fortran exponent-without-'E' form, which CAN occur in the linewidth columns, E22.10 of a
    double.)
  * omega is a double printed E20.8: 8 significant digits, relative rounding <= 5e-8.
  * EPW zeroes lambda for omega <= eps_acoustic (g2_tmp = 0) and for imaginary modes; a q point
    outside the qwindow keeps wf = 0; a restarted run (restart=.true.) writes lambda = 0 for every
    q before the restart point.  All of these are turned into hard errors below.
  * stdout per q prints lambda___( i ) (f15.6), omega (f12.4) and lambda___( tot ) (f15.6) from
    the local (pre-mp_sum) arrays: comparing them with the files catches rows that were never
    computed and the nimage>1 omega doubling.
"""
from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Any

import numpy as np

from .qlist import QListError, identify_shifted_grid, read_filqf, sha256_file
from .tags import EPH_PROTOCOL_TAG, incomplete_tag


class PhselfError(ValueError):
    """The EPW phself data cannot be turned into a trustworthy Spectrum."""


# ----------------------------------------------------------------------------- Fortran reals
_FREAL = re.compile(r"^([+-]?(?:\d+\.\d*|\.\d+))(?:[EeDdQq]([+-]?\d+)|([+-]\d+))?$")


def parse_fortran_real(tok: str) -> float:
    """Parse one Fortran E/F/D edit-descriptor output field.

    Accepts ``0.12345E+01``, ``-0.12345D-03``, ``1.5`` and the exponent-without-letter form used
    for |exponent| > 99 (``0.12345-100``, ``0.1234567890+101``).  Rejects asterisk overflow fields,
    NaN/Infinity, integers without a decimal point and anything else (PhselfError)."""
    t = tok.strip()
    m = _FREAL.match(t)
    if not m:
        raise PhselfError(f"not a Fortran real field: {tok!r}")
    mant, e1, e2 = m.groups()
    exp = e1 if e1 is not None else e2
    v = float(mant + ("e" + exp if exp is not None else ""))
    if not math.isfinite(v):
        raise PhselfError(f"non-finite value {tok!r}")
    return v


def _int(tok: str, what: str) -> int:
    if not re.fullmatch(r"[+-]?\d+", tok):
        raise PhselfError(f"{what}: expected an integer, got {tok!r}")
    return int(tok)


# ----------------------------------------------------------------------------- file parsers
def parse_lambda_phself(path: str | Path) -> dict[str, Any]:
    """Return {'nmodes', 'iq' (list), 'lam' (N, nmodes) ndarray, 'raw' (list of token lists)}.

    Token-stream parsing: robust to the list-directed '#Modes' header being wrapped (Intel
    Fortran wraps list-directed output at 80 columns) and to a data record being split."""
    p = Path(path)
    lines = p.read_text().splitlines()
    i = 0
    while i < len(lines) and not lines[i].strip():
        i += 1
    if i >= len(lines) or lines[i].strip() != "#Lambda phonon self-energy":
        raise PhselfError(f"{p}: missing '#Lambda phonon self-energy' header")
    i += 1
    while i < len(lines) and not lines[i].strip():
        i += 1
    if i >= len(lines) or not lines[i].strip().startswith("#Modes"):
        raise PhselfError(f"{p}: missing '#Modes' header line")
    modes = [_int(t, f"{p} #Modes header") for t in lines[i].strip()[len("#Modes"):].split()]
    i += 1
    # continuation lines of a wrapped list-directed header: all integers continuing 1..nmodes
    while i < len(lines):
        toks = lines[i].split()
        if toks and all(re.fullmatch(r"\d+", t) for t in toks) and \
                [int(t) for t in toks] == list(range(len(modes) + 1, len(modes) + 1 + len(toks))):
            modes += [int(t) for t in toks]
            i += 1
        else:
            break
    nmodes = len(modes)
    if nmodes == 0 or modes != list(range(1, nmodes + 1)):
        raise PhselfError(f"{p}: '#Modes' header must list 1..nmodes, got {modes}")
    toks = " ".join(lines[i:]).split()
    rec = 1 + nmodes
    if len(toks) % rec:
        raise PhselfError(f"{p}: {len(toks)} data tokens is not a multiple of 1+nmodes = {rec} "
                          "(truncated or malformed file)")
    nq = len(toks) // rec
    iq, lam, raw = [], np.empty((nq, nmodes)), []
    for r in range(nq):
        chunk = toks[r * rec:(r + 1) * rec]
        iq.append(_int(chunk[0], f"{p} record {r + 1} q index"))
        for j, t in enumerate(chunk[1:]):
            lam[r, j] = parse_fortran_real(t)
        raw.append(chunk[1:])
    return {"path": str(p), "nmodes": nmodes, "iq": iq, "lam": lam, "raw": raw}


def parse_linewidth_phself(path: str | Path) -> dict[str, Any]:
    """Return {'iq', 'imode' (int arrays), 'omega_meV', 'linewidth_meV' (float arrays),
    'omega_raw' (printed strings)}."""
    p = Path(path)
    lines = p.read_text().splitlines()
    hdr = [ln for ln in lines[:2]]
    if len(hdr) < 2 or not hdr[0].startswith("# Phonon frequency") or not hdr[1].startswith("# Q-point"):
        raise PhselfError(f"{p}: missing the two '#' header lines of linewidth.phself")
    iq, im, om, lw, oraw = [], [], [], [], []
    for n, ln in enumerate(lines[2:], start=3):
        if not ln.strip():
            continue
        t = ln.split()
        if len(t) != 4:
            raise PhselfError(f"{p}:{n}: expected 4 fields (iq imode omega linewidth), got {ln!r}")
        iq.append(_int(t[0], f"{p}:{n} iq"))
        im.append(_int(t[1], f"{p}:{n} imode"))
        om.append(parse_fortran_real(t[2]))
        lw.append(parse_fortran_real(t[3]))
        oraw.append(t[2])
    return {"path": str(p), "iq": np.array(iq, dtype=int), "imode": np.array(im, dtype=int),
            "omega_meV": np.array(om), "linewidth_meV": np.array(lw), "omega_raw": oraw}


_RE_BLOCK = re.compile(r"ismear\s*=\s*(\d+)\s+iq\s*=\s*(\d+)\s+coord\.:\s*(\S+)\s+(\S+)\s+(\S+)\s+"
                       r"wt:\s*(\S+)\s+Temp:\s*(\S+?)\s*K")
_RE_MODE = re.compile(r"lambda___\(\s*(\S+)\s*\)=\s*(\S+)\s+gamma___=\s*(\S+)\s*meV\s+omega=\s*(\S+)\s*meV")
_RE_TOT = re.compile(r"lambda_*\(\s*(?:tot|sum)\s*\)\s*=\s*(\S+)")
_RE_NEED = re.compile(r"We only need to compute\s+(\d+)\s+q-points")
_RE_DOS = re.compile(r"DOS\s*=\s*(\S+)\s+states/spin/eV/Unit Cell at Ef=\s*(\S+)\s*eV")
_RE_GAUSS = re.compile(r"Gaussian Broadening:\s*(\S+)\s*eV,\s*ngauss=\s*(\S+)")
_RE_NIMAGE = re.compile(r"nimage\s*=\s*(\d+)", re.I)
_RE_RESTART = re.compile(r"We restart from", re.I)


def _f(tok: str, what: str) -> float:
    try:
        return parse_fortran_real(tok)
    except PhselfError as exc:
        raise PhselfError(f"EPW stdout {what}: {exc}") from None


def parse_epw_stdout(path: str | Path) -> dict[str, Any]:
    """Parse the phonon self-energy part of an epw.x stdout.

    Returns {'blocks': {iq: {'ismear', 'coord', 'wt', 'temp_K', 'lam' [per mode], 'omega_meV'
    [per mode], 'tot'}}, 'need_q', 'dos', 'ef_eV', 'gauss_eV', 'ngauss', 'nimage', 'restart'}."""
    p = Path(path)
    blocks: dict[int, dict[str, Any]] = {}
    cur: dict[str, Any] | None = None
    out: dict[str, Any] = {"path": str(p), "blocks": blocks, "need_q": [], "dos": None,
                           "ef_eV": None, "gauss_eV": None, "ngauss": None, "nimage": [],
                           "restart": False, "duplicate_blocks": []}
    for n, ln in enumerate(p.read_text(errors="replace").splitlines(), start=1):
        if _RE_RESTART.search(ln):
            out["restart"] = True
        m = _RE_NIMAGE.search(ln)
        if m:
            out["nimage"].append(int(m.group(1)))
        m = _RE_NEED.search(ln)
        if m:
            out["need_q"].append(int(m.group(1)))
            continue
        m = _RE_GAUSS.search(ln)
        if m and out["gauss_eV"] is None:
            out["gauss_eV"] = _f(m.group(1), f"line {n} Gaussian Broadening")
            out["ngauss"] = _int(m.group(2), f"{p}:{n} ngauss")
            continue
        m = _RE_DOS.search(ln)
        if m and out["dos"] is None:
            out["dos"] = _f(m.group(1), f"line {n} DOS")
            out["ef_eV"] = _f(m.group(2), f"line {n} Ef")
            continue
        m = _RE_BLOCK.search(ln)
        if m:
            iq = int(m.group(2))
            cur = {"ismear": int(m.group(1)), "line": n,
                   "coord": [_f(m.group(k), f"line {n} coord") for k in (3, 4, 5)],
                   "wt": _f(m.group(6), f"line {n} wt"), "temp_K": _f(m.group(7), f"line {n} Temp"),
                   "mode": [], "lam": [], "omega_meV": [], "tot": None}
            if iq in blocks:
                out["duplicate_blocks"].append(iq)
            blocks[iq] = cur
            continue
        if cur is None:
            continue
        m = _RE_MODE.search(ln)
        if m:
            cur["mode"].append(_int(m.group(1), f"{p}:{n} mode index"))
            cur["lam"].append(_f(m.group(2), f"line {n} lambda___"))
            cur["omega_meV"].append(_f(m.group(4), f"line {n} omega"))
            continue
        m = _RE_TOT.search(ln)
        if m:
            cur["tot"] = _f(m.group(1), f"line {n} lambda___( tot )")
            cur = None
    return out


# ----------------------------------------------------------------------------- epw.in
_NML_PAIR = re.compile(r"([A-Za-z_][A-Za-z0-9_]*(?:\(\s*\d+\s*\))?)\s*=\s*"
                       r"('(?:[^']|'')*'|\"[^\"]*\"|[^,\s!/]+)")


def _strip_comment(line: str) -> str:
    q = None
    for i, ch in enumerate(line):
        if q:
            if ch == q:
                q = None
        elif ch in "'\"":
            q = ch
        elif ch == "!":
            return line[:i]
    return line


def parse_epw_namelist(path: str | Path) -> dict[str, str]:
    """Raw values of the &inputepw namelist, keys lower-cased (``temps(1)`` kept as is)."""
    text = Path(path).read_text()
    lines = [_strip_comment(ln) for ln in text.splitlines()]
    body, inside = [], False
    for ln in lines:
        s = ln.strip()
        if not inside:
            if s.lower().startswith("&inputepw"):
                inside = True
                body.append(s[len("&inputepw"):])
            continue
        if s.startswith("/"):
            break
        body.append(ln)
    if not inside:
        raise PhselfError(f"{path}: no &inputepw namelist")
    out: dict[str, str] = {}
    for m in _NML_PAIR.finditer("\n".join(body)):
        k = re.sub(r"\s+", "", m.group(1).lower())
        out[k] = m.group(2)
    return out


def _nml_bool(v: str | None, default: bool) -> bool:
    if v is None:
        return default
    s = v.strip().lower().strip(".")
    if s in ("true", "t"):
        return True
    if s in ("false", "f"):
        return False
    raise PhselfError(f"not a Fortran logical: {v!r}")


def _nml_real(v: str | None, default: float | None) -> float | None:
    if v is None:
        return default
    s = v.strip()
    if re.fullmatch(r"[+-]?\d+", s):
        return float(int(s))
    return parse_fortran_real(s)


def _nml_int(v: str | None, default: int | None) -> int | None:
    if v is None:
        return default
    return _int(v.strip(), "namelist integer")


def _nml_str(v: str | None) -> str | None:
    if v is None:
        return None
    s = v.strip()
    if len(s) >= 2 and s[0] == s[-1] and s[0] in "'\"":
        s = s[1:-1]
    return s.strip()


def check_epw_pins(nml: dict[str, str], *, fsthick_min_eV: float = 1.0,
                   fsthick_sigma_factor: float = 8.0) -> dict[str, Any]:
    """Enforce the v3 EPW settings; raise PhselfError listing EVERY violation.

    Defaults are EPW 6.1 readin defaults (ngaussw=1, eps_acoustic=0.1 cm-1, fsthick=1e10 eV)."""
    bad: list[str] = []

    def need(cond: bool, msg: str) -> None:
        if not cond:
            bad.append(msg)

    g = nml.get
    need(_nml_bool(g("phonselfen"), False) is True, "phonselfen must be .true.")
    need(_nml_bool(g("delta_approx"), False) is True, "delta_approx must be .true.")
    need(_nml_bool(g("a2f"), False) is False, "a2f must be .false. (EPW's binned a2F is not used)")
    need(_nml_bool(g("restart"), False) is False, "restart must be .false. (restart writes lambda=0 "
         "for q before the restart point)")
    need(_nml_int(g("nstemp"), None) == 1, "nstemp must be set to 1")
    need(_nml_int(g("ngaussw"), 1) == 0, "ngaussw must be 0 (Gaussian; EPW default is 1 = M-P, "
         "which gives negative weights)")
    need(_nml_bool(g("efermi_read"), False) is False, "efermi_read must be .false.")
    need(_nml_bool(g("selecqread"), False) is False, "selecqread must be .false.")
    need(_nml_bool(g("mp_mesh_q"), False) is False, "mp_mesh_q must be .false.")
    need(_nml_bool(g("mp_mesh_k"), False) is False, "mp_mesh_k must be .false.")
    need(_nml_bool(g("rand_q"), False) is False, "rand_q must be .false.")
    need(_nml_bool(g("rand_k"), False) is False, "rand_k must be .false.")
    need(_nml_bool(g("shortrange"), False) is False, "shortrange must be .false.")
    need(_nml_int(g("nsmear"), 1) == 1, "nsmear must be 1 (the file holds smearing index 1 only)")
    filqf = _nml_str(g("filqf"))
    need(bool(filqf), "filqf must name the protocol q list")
    for k in ("nqf1", "nqf2", "nqf3"):
        need(_nml_int(g(k), 0) == 0, f"{k} must be absent/0 when filqf is used")
    need(not _nml_str(g("filkf")), "filkf must not be used (fine k = uniform nkf grid)")
    nkf = [_nml_int(g(k), 0) for k in ("nkf1", "nkf2", "nkf3")]
    need(all(x and x > 0 for x in nkf), "nkf1..3 must be > 0")
    need(_nml_int(g("etf_mem"), 1) != 3, "etf_mem=3 (lfast_kmesh) is not allowed")
    sigma = _nml_real(g("degaussw"), 0.025)
    need(sigma is not None and sigma > 0, "degaussw must be > 0 (0 = adaptive smearing)")
    fst = _nml_real(g("fsthick"), 1e10)
    if sigma and sigma > 0:
        need(fst >= max(fsthick_min_eV, fsthick_sigma_factor * sigma) * (1 - 1e-12),
             f"fsthick = {fst} eV < max({fsthick_min_eV} eV, {fsthick_sigma_factor}*sigma = "
             f"{fsthick_sigma_factor * sigma} eV)")
    for k in ("lscreen", "plselfen", "specfun_pl", "system_2d"):
        if k == "system_2d":
            v = _nml_str(g(k))
            need(v in (None, "no"), f"{k} must be 'no'")
        else:
            need(_nml_bool(g(k), False) is False, f"{k} must be .false.")
    eps = _nml_real(g("eps_acoustic"), 0.1)
    temps = _nml_real(g("temps(1)"), None)
    if temps is None and "temps" in nml:
        temps = _nml_real(g("temps"), None)
    if bad:
        raise PhselfError("epw.in violates the protocol pins: " + "; ".join(bad))
    return {"sigma_eV": sigma, "sigma_meV": sigma * 1e3, "fsthick_eV": fst,
            "eps_acoustic_cm1": eps, "temps_K": temps, "nkf": nkf, "filqf": filqf,
            "ngaussw": 0, "nstemp": 1, "phonselfen": True, "delta_approx": True, "a2f": False,
            "restart": False, "efermi_read": False, "selecqread": False}


# ----------------------------------------------------------------------------- builder
LAMBDA_REL_ROUND = 5e-5      # E15.5 mantissa 0.ddddd: half unit / 0.1
OMEGA_REL_ROUND = 5e-8       # E20.8
STDOUT_ABS_F15_6 = 5e-7      # f15.6 half unit
STDOUT_ABS_F12_4 = 5e-5      # f12.4 half unit
STDOUT_ABS_F9_5 = 5e-6       # f9.5 half unit
CM1_TO_MEV = 0.12398419843320026   # h c / e * 100 cm/m * 1000 meV/eV  (CODATA exact constants)


def _lam_tol(sum_abs: float) -> float:
    return 1.001 * (LAMBDA_REL_ROUND * sum_abs + STDOUT_ABS_F15_6) + 1e-12


def build_phself_spectrum(lambda_path: str | Path, linewidth_path: str | Path,
                          filqf_path: str | Path, epw_in_path: str | Path,
                          stdout_path: str | Path | None = None, *, label: str | None = None,
                          protocol_eps_acoustic_cm1: float = 0.1, allow_zero_lambda_q: bool = False,
                          allow_no_stdout: bool = False, rel_envelope: float = 0.0,
                          extra_provenance: dict | None = None):
    """Parse, check and build the rtsc Spectrum (kelvin) of one EPW phself run.

    Raises PhselfError on any failed assertion.  ``allow_zero_lambda_q`` accepts q points whose
    every lambda is exactly 0 (default: refuse; the restart-zeros signature).  Without stdout the
    restart/qwindow/nimage/omega checks cannot be run: refused unless ``allow_no_stdout``, in which
    case the spectrum carries ``A-FILE-INCOMPLETE`` (rtsc blocks both proved directions)."""
    from rtsc.spectrum import Spectrum, unit_factor

    paths = {"lambda": Path(lambda_path), "linewidth": Path(linewidth_path),
             "filqf": Path(filqf_path), "epw_in": Path(epw_in_path)}
    if stdout_path is not None:
        paths["stdout"] = Path(stdout_path)
    elif not allow_no_stdout:
        raise PhselfError("EPW stdout is required (restart, qwindow, nimage and per-q lambda "
                          "checks); pass allow_no_stdout=True for a diagnostic-only spectrum")
    for k in ("lambda", "linewidth"):
        if paths[k].name.endswith(".down"):
            raise PhselfError(f"{paths[k]}: spin-polarised (lsda) output is not supported")
    tl = paths["lambda"].name.removeprefix("lambda.phself.")
    tw = paths["linewidth"].name.removeprefix("linewidth.phself.")
    if tl != tw:
        raise PhselfError(f"lambda/linewidth files are from different temperatures: {tl!r} vs {tw!r}")

    # --- epw.in pins and sigma
    nml = parse_epw_namelist(paths["epw_in"])
    settings = check_epw_pins(nml)
    if Path(settings["filqf"]).name != paths["filqf"].name:
        raise PhselfError(f"epw.in filqf = {settings['filqf']!r} but q list given is "
                          f"{paths['filqf'].name!r}")

    # --- q list
    try:
        pts, wts, ctype = read_filqf(paths["filqf"])
        grid = identify_shifted_grid(pts, wts, ctype)
    except QListError as exc:
        raise PhselfError(f"filqf: {exc}") from None
    Nq = len(pts)

    # --- files
    L = parse_lambda_phself(paths["lambda"])
    W = parse_linewidth_phself(paths["linewidth"])
    nm = L["nmodes"]
    lam = L["lam"]
    if lam.shape[0] != Nq:
        raise PhselfError(f"lambda file has {lam.shape[0]} q rows, filqf has N_q = {Nq}")
    if L["iq"] != list(range(1, Nq + 1)):
        raise PhselfError(f"lambda file q indices are not 1..{Nq} contiguous")
    if len(W["iq"]) != Nq * nm:
        raise PhselfError(f"linewidth file has {len(W['iq'])} rows, expected N_q*nmodes = "
                          f"{Nq}*{nm} = {Nq * nm}")
    exp_iq = np.repeat(np.arange(1, Nq + 1), nm)
    exp_im = np.tile(np.arange(1, nm + 1), Nq)
    if not np.array_equal(W["iq"], exp_iq):
        raise PhselfError(f"linewidth q indices are not 1..{Nq} contiguous (each repeated nmodes)")
    if not np.array_equal(W["imode"], exp_im):
        raise PhselfError(f"linewidth mode indices are not 1..{nm} within each q")
    om = W["omega_meV"].reshape(Nq, nm)

    # lambda: finite and >= 0
    if not np.all(np.isfinite(lam)):
        raise PhselfError("non-finite lambda")
    nneg = int(np.count_nonzero(lam < 0))
    if nneg:
        r, c = np.argwhere(lam < 0)[0]
        raise PhselfError(f"{nneg} negative lambda_qnu (first: iq={r + 1} mode={c + 1} "
                          f"lambda={lam[r, c]!r}); Gaussian double-delta weights are >= 0")
    # omega: > eps (no Gamma in a shifted grid; EPW zeroes lambda at omega <= eps_acoustic)
    eps_cm1 = max(float(settings["eps_acoustic_cm1"]), float(protocol_eps_acoustic_cm1))
    eps_meV = eps_cm1 * CM1_TO_MEV
    low = om <= eps_meV
    if np.any(low):
        nl = int(np.count_nonzero(low))
        nimag = int(np.count_nonzero(om < 0))
        r, c = np.argwhere(low)[0]
        raise PhselfError(f"{nl} mode(s) with omega <= eps = {eps_cm1} cm-1 ({eps_meV:.6g} meV), "
                          f"{nimag} of them imaginary (first: iq={r + 1} mode={c + 1} omega="
                          f"{om[r, c]!r} meV); no exclusions are allowed on a shifted grid "
                          "(EPW zeroes their lambda; a q outside the qwindow has omega = 0)")
    if np.any(np.diff(om, axis=1) < 0):
        r = int(np.argmax(np.any(np.diff(om, axis=1) < 0, axis=1)))
        raise PhselfError(f"frequencies not sorted ascending at iq={r + 1}: {om[r].tolist()}")

    zero_per_q = (lam == 0).sum(axis=1)
    q_all_zero = [int(i) + 1 for i in np.nonzero(zero_per_q == nm)[0]]
    if q_all_zero and not allow_zero_lambda_q:
        raise PhselfError(f"{len(q_all_zero)} q point(s) with lambda = 0 for every mode (e.g. iq="
                          f"{q_all_zero[:10]}): restart-zeros signature or no Fermi-surface pairs; "
                          "refusing (allow_zero_lambda_q=True to accept)")

    # --- stdout cross-checks
    sinfo: dict[str, Any] | None = None
    assumptions = [EPH_PROTOCOL_TAG]
    if "stdout" in paths:
        S = parse_epw_stdout(paths["stdout"])
        prob: list[str] = []
        if S["restart"]:
            prob.append("stdout shows a restarted run ('We restart from')")
        if any(x != 1 for x in S["nimage"]):
            prob.append(f"nimage = {S['nimage']} (must be 1)")
        if not S["need_q"]:
            prob.append("no 'We only need to compute N q-points' line")
        elif any(x != Nq for x in S["need_q"]):
            prob.append(f"qwindow kept {S['need_q']} q points, N_q = {Nq} (rows for dropped q are "
                        "not computed)")
        if S["duplicate_blocks"]:
            prob.append(f"duplicate stdout blocks for iq {sorted(set(S['duplicate_blocks']))[:10]}")
        if S["gauss_eV"] is None:
            prob.append("no 'Gaussian Broadening' line")
        else:
            if abs(S["gauss_eV"] - settings["sigma_eV"]) > 5e-7 * 1.001:
                prob.append(f"stdout broadening {S['gauss_eV']} eV != degaussw {settings['sigma_eV']} eV")
            if S["ngauss"] != 0:
                prob.append(f"stdout ngauss = {S['ngauss']} (must be 0)")
        if S["ef_eV"] is None:
            prob.append("no 'DOS = ... at Ef=' line")
        missing = [i for i in range(1, Nq + 1) if i not in S["blocks"]]
        extra = sorted(set(S["blocks"]) - set(range(1, Nq + 1)))
        if missing:
            prob.append(f"{len(missing)} q point(s) without a stdout lambda block (e.g. {missing[:10]})")
        if extra:
            prob.append(f"stdout blocks for q indices outside 1..{Nq}: {extra[:10]}")
        max_dl = max_dtot = max_dom = max_dx = 0.0
        for iq in range(1, Nq + 1):
            b = S["blocks"].get(iq)
            if b is None:
                continue
            if b["ismear"] != 1:
                prob.append(f"iq={iq}: ismear = {b['ismear']}")
            if b["mode"] != list(range(1, nm + 1)) or b["tot"] is None:
                prob.append(f"iq={iq}: incomplete stdout block (modes {b['mode'][:5]}..., tot "
                            f"{b['tot']})")
                continue
            row = lam[iq - 1]
            for j in range(nm):
                d = abs(row[j] - b["lam"][j])
                max_dl = max(max_dl, d)
                if d > _lam_tol(abs(row[j])):
                    prob.append(f"iq={iq} mode={j + 1}: file lambda {row[j]!r} vs stdout {b['lam'][j]!r}")
                dw = abs(om[iq - 1, j] - b["omega_meV"][j])
                max_dom = max(max_dom, dw)
                if dw > 1.001 * (STDOUT_ABS_F12_4 + OMEGA_REL_ROUND * abs(om[iq - 1, j])) + 1e-12:
                    prob.append(f"iq={iq} mode={j + 1}: file omega {om[iq - 1, j]!r} meV vs stdout "
                                f"{b['omega_meV'][j]!r} meV (nimage>1 doubling?)")
            s = math.fsum(row)
            dt = abs(s - b["tot"])
            max_dtot = max(max_dtot, dt)
            if dt > _lam_tol(math.fsum(abs(x) for x in row)):
                prob.append(f"iq={iq}: sum_nu lambda (file) = {s!r} vs stdout lambda( tot ) {b['tot']!r}")
            dx = max(abs(((b["coord"][k] - pts[iq - 1, k]) + 0.5) % 1.0 - 0.5) for k in range(3))
            max_dx = max(max_dx, dx)
            if dx > 1.001 * STDOUT_ABS_F9_5:
                prob.append(f"iq={iq}: stdout coord {b['coord']} != filqf {pts[iq - 1].tolist()}")
            if abs(b["wt"] - wts[iq - 1]) > 1.001 * STDOUT_ABS_F9_5:
                prob.append(f"iq={iq}: stdout wt {b['wt']} != filqf {wts[iq - 1]}")
        if prob:
            raise PhselfError("EPW stdout consistency failed: " + "; ".join(prob[:20])
                              + (f"; ... ({len(prob)} problems)" if len(prob) > 20 else ""))
        sinfo = {"ef_eV": S["ef_eV"], "dos_states_per_spin_eV": S["dos"], "gauss_eV": S["gauss_eV"],
                 "need_q": S["need_q"], "max_abs_dlambda_mode": max_dl, "max_abs_dlambda_tot": max_dtot,
                 "max_abs_domega_meV": max_dom, "max_abs_dcoord": max_dx}
    else:
        assumptions.append(incomplete_tag(
            "EPW stdout not supplied: restart, qwindow, nimage and per-q lambda checks of the "
            "phself files were not run (diagnostic only)"))

    # --- atoms: weight = lambda * omega / 2 / N_q  (meV), merged on identical printed omega
    groups: dict[float, list[float]] = {}
    order: list[float] = []
    n_raw = 0
    for r in range(Nq):
        for j in range(nm):
            if lam[r, j] == 0.0:
                continue
            n_raw += 1
            w = om[r, j]
            if w not in groups:
                groups[w] = []
                order.append(w)
            groups[w].append(lam[r, j] * w / 2.0 / Nq)
    if not groups:
        raise PhselfError("every lambda_qnu is zero: no spectrum")
    order.sort()
    om_meV = order
    wt_meV = [math.fsum(groups[w]) for w in order]

    lam_qavg = math.fsum(float(x) for x in lam.ravel()) / Nq
    f = unit_factor("meV")
    A_meV = math.fsum(wt_meV)
    prec = {
        "lambda_print": "E15.5 of REAL(lambda_all) (single precision): <= 5e-5 relative per value; "
                        "|lambda| < ~1.4e-45 prints as 0",
        "omega_print": "E20.8 of a double (meV): <= 5e-8 relative per value (position uncertainty, "
                       "outside the weight envelope)",
        "stdout_lambda_print": "f15.6 absolute (5e-7)",
        "atom_weight_rel_rounding_max": LAMBDA_REL_ROUND + OMEGA_REL_ROUND,
        "merge_rule": "atoms with identical printed omega merged exactly (fsum of weights)",
    }
    prov: dict[str, Any] = {
        "kind": "eph_protocol_phself_atoms",
        "protocol": "EPH_PROTOCOL_DRAFT",
        "protocol_version": "v3-draft",
        "construction": "alpha2F = sum_{q,nu} (1/N_q) lambda_qnu omega_qnu/2 delta(omega - omega_qnu); "
                        "omega = linewidth.phself col 3 (meV), lambda = lambda.phself",
        "sources": {k: {"path": str(v), "sha256": sha256_file(v)} for k, v in paths.items()},
        "N_q": Nq, "nmodes": nm, "q_grid": list(grid), "q_shift": [0.5, 0.5, 0.5],
        "q_weight": 1.0 / Nq,
        "sigma_eV": settings["sigma_eV"], "sigma_meV": settings["sigma_meV"],
        "fsthick_eV": settings["fsthick_eV"], "eps_acoustic_cm1_epw": settings["eps_acoustic_cm1"],
        "eps_acoustic_cm1_applied": eps_cm1, "temperature_label": tl,
        "epw_settings": settings, "precision": prec,
        "counts": {
            "modes_total": Nq * nm,
            "excluded_modes": 0,
            "lambda_zero_modes": int((lam == 0).sum()),
            "lambda_zero_modes_per_q": {int(i) + 1: int(c) for i, c in enumerate(zero_per_q) if c},
            "q_all_lambda_zero": q_all_zero,
            "atoms_nonzero_before_merge": n_raw,
            "atoms_after_merge": len(om_meV),
            "atoms_merged_away": n_raw - len(om_meV),
        },
        "lambda_q_average": lam_qavg,
        "A_meV": A_meV,
        "stdout_checks": sinfo,
        "stdout_checked": sinfo is not None,
        "assumptions": assumptions,
    }
    prov.update(extra_provenance or {})
    lab = label if label is not None else (
        f"phself:sigma={settings['sigma_meV']:g}meV:Nq={Nq}({grid[0]}x{grid[1]}x{grid[2]})")
    sp = Spectrum.from_atoms(om_meV, wt_meV, unit="meV", rel_envelope=rel_envelope, label=lab,
                             provenance=prov, drop_zero=False)
    lam_spec = sp.lambda_()
    if not (lam_qavg > 0 and abs(lam_spec - lam_qavg) <= 1e-6 * lam_qavg):
        raise PhselfError(f"lambda of the spectrum {lam_spec!r} != q-average of sum_nu lambda_qnu "
                          f"{lam_qavg!r} (rel tol 1e-6)")
    sp.provenance["lambda_spectrum"] = lam_spec
    sp.provenance["A_K"] = sp.A()
    sp.provenance["omega_log_K"] = sp.omega_log()
    assert abs(sp.A() - A_meV * f) <= 1e-12 * max(1.0, sp.A())
    return sp
