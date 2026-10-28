"""Gate G1 (Wannier band quality) of the alpha^2F convergence protocol (v3 P4 + v3.1 A7).

G1 passes when, on >= 100 validation k points that are NOT on the coarse (nk1 x nk2 x nk3) grid,

  (a) max |E_Wannier - E_DFT| <= min(30 meV, sigma_min,adm / 2) over the states in
      W = [E_F - 4 sigma_min,adm, E_F + 4 sigma_min,adm]            (P4)
  (b) every validation k has the same number of DFT and Wannier states in W   (A7)
  (c) EPW's fine-grid Fermi energy lies within 50 meV of the DFT E_F           (A7)

Validation k rule (pre-registered): the half-shifted 5x5x5 Monkhorst-Pack grid in crystal
coordinates, x_i = (i + 1/2)/5, minus any point that lies on the coarse grid (exact rational test).
For an 8^3 coarse grid exactly one point is removed, (1/2, 1/2, 1/2) = L, leaving 124.

Matching rule. The Wannier bands are compared with the DFT bands above the excluded (semicore)
bands.  Both lists are sorted per k and paired by index (Wannier j <-> j-th non-excluded DFT
band).  A pair enters the statistics if EITHER member lies in W, so a DFT state with no Wannier
partner and a spurious Wannier state both show up as a large |dE|.  The state counts in W are
compared separately (A7).  A nearest-neighbour distance (each in-window state to the closest state
of the other set, any energy) is reported as a secondary, ordering-independent diagnostic.

Files
  band.eig (EPW plot_band): '&plot nbnd=NW, nks=NK /', then per k a line (10x,3f10.6) with k in
      cartesian 2pi/alat and a line (1000f20.12) with absolute energies in eV.
  pw.x data-file-schema.xml: ks_energies/k_point (cartesian 2pi/alat), eigenvalues in Hartree.

CLI (WSL, venv ~/rtsc_work/.venv):
  python -m ephproto.bandcheck kpoints OUTDIR [--n 5] [--coarse 8 8 8]
  python -m ephproto.bandcheck compare --pw-xml X --band-eig B --ef 11.109 --nexclude 4
         [--epw-out EPW_EF_OUT] [--json OUT.json]
  python -m ephproto.bandcheck run material.json --wannier DIR --epmatwp FILE --label NAME [--np 4]
"""
from __future__ import annotations

import argparse
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time
import xml.etree.ElementTree as ET
from fractions import Fraction
from pathlib import Path
from typing import Sequence

import numpy as np

HARTREE_EV = 27.211386245988       # QE 7.x AUTOEV (= 2 * EPW ryd2ev)
K_TOL = 2.0e-6                     # band.eig writes k with f10.6
EF_TOL_EV = 0.050                  # A7
SIGMA_MIN_ADM_EV = 0.020           # evaluation value for now (protocol: sigma_min,adm <= 20 meV)


class BandCheckError(ValueError):
    pass


# ------------------------------------------------------------------------------ windows
def g1_windows(sigma_min_adm_eV: float = SIGMA_MIN_ADM_EV) -> list[dict]:
    """The gate window (P4) plus the wider diagnostic window (spec section 6)."""
    return [
        {"name": "gate", "half_width_eV": 4.0 * sigma_min_adm_eV,
         "tol_eV": min(0.030, sigma_min_adm_eV / 2.0), "gating": True},
        {"name": "wide", "half_width_eV": 1.0, "tol_eV": 0.030, "gating": False},
    ]


# ----------------------------------------------------------------------------- k points
def on_grid(p: Sequence[Fraction], coarse: Sequence[int]) -> bool:
    return all((Fraction(c) * n).denominator == 1 for c, n in zip(p, coarse))


def validation_kpoints(n: Sequence[int] = (5, 5, 5), coarse: Sequence[int] = (8, 8, 8)):
    """(kept Fractions, removed Fractions): half-shifted n-grid minus coarse-grid points."""
    pts = [(Fraction(2 * i + 1, 2 * n[0]), Fraction(2 * j + 1, 2 * n[1]), Fraction(2 * k + 1, 2 * n[2]))
           for i in range(n[0]) for j in range(n[1]) for k in range(n[2])]
    keep = [p for p in pts if not on_grid(p, coarse)]
    drop = [p for p in pts if on_grid(p, coarse)]
    return keep, drop


def filkf_text(pts) -> str:
    w = 1.0 / len(pts)
    lines = [f"{len(pts)} crystal"]
    lines += [" ".join(repr(float(c)) for c in p) + f" {w!r}" for p in pts]
    return "\n".join(lines) + "\n"


def pw_kpoints_text(pts) -> str:
    lines = ["K_POINTS crystal", str(len(pts))]
    lines += [" ".join(f"{float(c):.12f}" for c in p) + " 1.0" for p in pts]
    return "\n".join(lines) + "\n"


def gamma_filqf_text() -> str:
    return "1 crystal\n0.0 0.0 0.0 1.0\n"


def random_kpoints(nrand: int, seed: int, coarse: Sequence[int] = (8, 8, 8)):
    """SUPPLEMENTARY (post hoc, not the pre-registered rule): nrand uniform crystal points in
    [0,1)^3 from numpy default_rng(seed), rounded to 1e-6 (filkf/band.eig precision).  Points within
    1e-4 of a coarse-grid point (in any component combination) are dropped."""
    rng = np.random.default_rng(seed)
    pts = np.round(rng.random((nrand, 3)), 6)
    c = np.asarray(coarse, float)
    near = np.all(np.abs(pts * c - np.round(pts * c)) < 1e-4 * c, axis=1)
    keep = [tuple(Fraction(repr(float(x))) for x in p) for p in pts[~near]]
    drop = [tuple(Fraction(repr(float(x))) for x in p) for p in pts[near]]
    return keep, drop


def write_kpoint_files(outdir: str | Path, n=(5, 5, 5), coarse=(8, 8, 8),
                       kset: str = "mp5", nrand: int = 1000, seed: int = 20261027) -> dict:
    d = Path(outdir)
    d.mkdir(parents=True, exist_ok=True)
    if kset == "mp5":
        keep, drop = validation_kpoints(n, coarse)
    elif kset == "random":
        keep, drop = random_kpoints(nrand, seed, coarse)
    else:
        raise BandCheckError(f"unknown kset {kset!r}")
    if len(keep) < 100:
        raise BandCheckError(f"only {len(keep)} validation k points (< 100)")
    (d / "kval.txt").write_text(filkf_text(keep))
    (d / "kpoints_pw.txt").write_text(pw_kpoints_text(keep))
    (d / "gamma.txt").write_text(gamma_filqf_text())
    return {"kset": kset, "n": list(n) if kset == "mp5" else None,
            "nrand": nrand if kset == "random" else None, "seed": seed if kset == "random" else None,
            "coarse": list(coarse), "N_k": len(keep),
            "removed_on_coarse_grid": [[str(c) for c in p] for p in drop]}


# ------------------------------------------------------------------------------ parsers
_RE_PLOT = re.compile(r"&plot\s+nbnd\s*=\s*(\d+)\s*,\s*nks\s*=\s*(\d+)", re.I)


def _fixed(line: str, width: int) -> list[float]:
    s = line.rstrip("\n")
    if not s.strip():
        return []
    if len(s) % width:                       # tolerate stripped leading blanks
        s = s.rjust(len(s) + width - len(s) % width)
    out = []
    for i in range(0, len(s), width):
        tok = s[i:i + width].strip()
        if tok:
            if "*" in tok:
                raise BandCheckError(f"Fortran overflow field {tok!r}")
            out.append(float(tok.replace("D", "E").replace("d", "e")))
    return out


def parse_band_eig(path: str | Path) -> tuple[np.ndarray, np.ndarray]:
    """band.eig -> (k_cart (NK,3) in 2pi/alat, E (NK, NW) in eV, absolute)."""
    lines = [ln for ln in Path(path).read_text().splitlines() if ln.strip()]
    if not lines:
        raise BandCheckError(f"{path}: empty")
    m = _RE_PLOT.search(lines[0])
    if not m:
        raise BandCheckError(f"{path}: bad header {lines[0]!r}")
    nw, nk = int(m.group(1)), int(m.group(2))
    body = lines[1:]
    if len(body) != 2 * nk:
        raise BandCheckError(f"{path}: expected {2 * nk} data lines, got {len(body)}")
    ks, es = [], []
    for ik in range(nk):
        kl, el = body[2 * ik], body[2 * ik + 1]
        k = [float(t) for t in kl.split()]
        if len(k) != 3:
            raise BandCheckError(f"{path}: k line {ik + 1} has {len(k)} values")
        e = _fixed(el, 20) if len(el.strip().split()) != nw else [float(t) for t in el.split()]
        if len(e) != nw:
            raise BandCheckError(f"{path}: k {ik + 1}: {len(e)} energies, header says {nw}")
        ks.append(k)
        es.append(e)
    return np.array(ks, float), np.array(es, float)


def parse_pw_xml(path: str | Path) -> tuple[np.ndarray, np.ndarray]:
    """pw.x data-file-schema.xml -> (k_cart (NK,3) in 2pi/alat, E (NK, nbnd) in eV)."""
    root = ET.parse(path).getroot()
    ks, es = [], []
    for kse in root.iter("ks_energies"):
        kp = kse.find("k_point")
        ev = kse.find("eigenvalues")
        if kp is None or ev is None:
            raise BandCheckError(f"{path}: ks_energies without k_point/eigenvalues")
        ks.append([float(t) for t in kp.text.split()])
        es.append([float(t) * HARTREE_EV for t in ev.text.split()])
    if not ks:
        raise BandCheckError(f"{path}: no ks_energies")
    if len({len(e) for e in es}) != 1:
        raise BandCheckError(f"{path}: ragged eigenvalue lists")
    return np.array(ks, float), np.array(es, float)


_RE_EF_FINE = re.compile(r"Fermi energy is calculated from the fine k-mesh:\s*Ef\s*=\s*([-+0-9.EeDd]+)\s*eV")
_RE_EF_COARSE = re.compile(r"Fermi energy coarse grid\s*=\s*([-+0-9.EeDd]+)\s*eV")


def parse_epw_fermi(path: str | Path) -> dict:
    txt = Path(path).read_text(errors="replace")
    fine = _RE_EF_FINE.findall(txt)
    coarse = _RE_EF_COARSE.findall(txt)
    if not fine:
        raise BandCheckError(f"{path}: no 'Fermi energy is calculated from the fine k-mesh' line")
    if len(set(fine)) != 1:
        raise BandCheckError(f"{path}: inconsistent fine-mesh E_F lines {fine}")
    return {"ef_fine_eV": float(fine[0].replace("D", "E")),
            "ef_coarse_eV": float(coarse[0].replace("D", "E")) if coarse else None}


# ----------------------------------------------------------------------------- matching
def compare_bands(k_dft, E_dft, k_wan, E_wan, ef: float, nexclude: int, windows: list[dict],
                  k_tol: float = K_TOL) -> dict:
    k_dft, E_dft = np.asarray(k_dft, float), np.asarray(E_dft, float)
    k_wan, E_wan = np.asarray(k_wan, float), np.asarray(E_wan, float)
    if k_dft.shape != k_wan.shape:
        raise BandCheckError(f"k lists differ in shape: DFT {k_dft.shape}, Wannier {k_wan.shape}")
    dk = float(np.max(np.abs(k_dft - k_wan))) if len(k_dft) else 0.0
    if dk > k_tol:
        bad = int(np.argmax(np.max(np.abs(k_dft - k_wan), axis=1)))
        raise BandCheckError(f"k lists differ (max |dk| = {dk:.2e} at k #{bad + 1}: "
                             f"{k_dft[bad]} vs {k_wan[bad]})")
    D = np.sort(E_dft[:, nexclude:], axis=1)
    Wn = np.sort(E_wan, axis=1)
    nw = Wn.shape[1]
    if D.shape[1] < nw:
        raise BandCheckError(f"only {D.shape[1]} DFT bands after excluding {nexclude}; need >= {nw}")
    res = {"N_k": int(len(k_dft)), "max_dk": dk, "nexclude": nexclude, "nwann": nw,
           "nbnd_dft": int(E_dft.shape[1]), "ef_eV": ef, "windows": []}
    for w in windows:
        lo, hi = ef - w["half_width_eV"], ef + w["half_width_eV"]
        dEs, nn, mism, worst = [], [], [], None
        for ik in range(len(D)):
            d_in = (D[ik] >= lo) & (D[ik] <= hi)
            w_in = (Wn[ik] >= lo) & (Wn[ik] <= hi)
            nd, nwi = int(d_in.sum()), int(w_in.sum())
            if nd != nwi:
                mism.append({"ik": ik + 1, "n_dft": nd, "n_wann": nwi})
            for j in range(nw):                       # index pairing, pair counted if either in W
                if d_in[j] or w_in[j]:
                    de = float(Wn[ik, j] - D[ik, j])
                    dEs.append(de)
                    if worst is None or abs(de) > abs(worst["dE_eV"]):
                        worst = {"ik": ik + 1, "band": j + 1, "E_dft": float(D[ik, j]),
                                 "E_wann": float(Wn[ik, j]), "dE_eV": de}
            for j in range(D.shape[1]):               # DFT states in W with no pair partner
                if j >= nw and d_in[j]:
                    de = float(np.min(np.abs(Wn[ik] - D[ik, j])))
                    dEs.append(de)
                    if worst is None or abs(de) > abs(worst["dE_eV"]):
                        worst = {"ik": ik + 1, "band": j + 1, "E_dft": float(D[ik, j]),
                                 "E_wann": None, "dE_eV": de}
            for e in D[ik][d_in]:
                nn.append(float(np.min(np.abs(Wn[ik] - e))))
            for e in Wn[ik][w_in]:
                nn.append(float(np.min(np.abs(D[ik] - e))))
        a = np.abs(np.array(dEs)) if dEs else np.zeros(0)
        mx = float(a.max()) if a.size else 0.0
        out = {**w, "E_lo": lo, "E_hi": hi, "n_states": int(a.size),
               "max_abs_dE_meV": 1000 * mx,
               "rms_dE_meV": 1000 * float(np.sqrt(np.mean(a ** 2))) if a.size else 0.0,
               "max_nn_dE_meV": 1000 * (max(nn) if nn else 0.0),
               "n_k_count_mismatch": len(mism), "count_mismatches": mism[:20], "worst": worst,
               "n_k_with_states": int(sum(1 for ik in range(len(D)) if
                                      np.any((D[ik] >= lo) & (D[ik] <= hi)) or
                                      np.any((Wn[ik] >= lo) & (Wn[ik] <= hi)))),
               "inconclusive_no_states": bool(a.size == 0),
               "pass_dE": bool(a.size > 0 and mx <= w["tol_eV"]),
               "pass_count": len(mism) == 0}
        out["pass"] = out["pass_dE"] and out["pass_count"]
        res["windows"].append(out)
    return res


def g1_verdict(cmp: dict, ef_epw: float | None, ef_dft: float) -> dict:
    gate = next(w for w in cmp["windows"] if w["gating"])
    v = {"N_k": cmp["N_k"], "N_k_ok": cmp["N_k"] >= 100,
         "dE_ok": gate["pass_dE"], "count_ok": gate["pass_count"],
         # no state in the gate window at any validation k: the dE test is vacuous, not passed
         "gate_inconclusive_no_states": gate["inconclusive_no_states"],
         "ef_epw_eV": ef_epw, "ef_dft_eV": ef_dft,
         "dEF_meV": None if ef_epw is None else 1000 * (ef_epw - ef_dft)}
    v["ef_ok"] = ef_epw is not None and abs(ef_epw - ef_dft) <= EF_TOL_EV
    v["G1_pass"] = all(v[k] for k in ("N_k_ok", "dE_ok", "count_ok", "ef_ok"))
    return v


# ---------------------------------------------------------------------------- run stage
QE_BIN = Path.home() / "micromamba/envs/qe/bin"


def _sh(cmd: str, cwd: Path, log: Path, np_: int) -> float:
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    if np_ > 1:
        cmd = (f"{QE_BIN}/mpirun -np {np_} --bind-to none -x OMP_NUM_THREADS -x OPENBLAS_NUM_THREADS "
               + cmd)
    t0 = time.time()
    with open(log, "w") as fh:
        rc = subprocess.call(cmd, shell=True, cwd=cwd, stdout=fh, stderr=subprocess.STDOUT, env=env)
    if rc != 0:
        raise RuntimeError(f"rc={rc}: {cmd} (see {log})")
    return time.time() - t0


def pw_bands_input(nscf_text: str, kpts_block: str) -> str:
    head = re.split(r"^\s*K_POINTS", nscf_text, flags=re.M | re.I)[0]
    head = re.sub(r"calculation\s*=\s*'[^']*'", "calculation = 'bands'", head)
    head = re.sub(r"outdir\s*=\s*'[^']*'", "outdir = './tmp/'", head)
    head = re.sub(r"^\s*(tprnfor|tstress)\s*=.*\n", "", head, flags=re.M)
    if "calculation = 'bands'" not in head:
        raise BandCheckError("nscf.in has no calculation key to replace")
    return head.rstrip() + "\n" + kpts_block


def run_pw_bands(m: dict, root: Path, np_: int) -> Path:
    """pw.x 'bands' on the validation k (shared by every Wannierisation)."""
    d = root / "pw_bands"
    xml = d / "tmp" / f"{m['prefix']}.save" / "data-file-schema.xml"
    if (d / "DONE").exists():
        return xml
    if d.exists():
        shutil.rmtree(d)
    sv = d / "tmp" / f"{m['prefix']}.save"
    sv.mkdir(parents=True)
    src = Path(m["epw_dir"]) / "tmp" / f"{m['prefix']}.save"
    for f in src.iterdir():                     # copy of the scf density only, never the original
        if f.name.startswith("charge-density") or f.name == "data-file-schema.xml" \
                or f.suffix.lower() == ".upf":
            shutil.copy(f, sv / f.name)
    nscf = (Path(m["epw_dir"]) / "nscf.in").read_text()
    (d / "bands.in").write_text(pw_bands_input(nscf, (root / "kpoints_pw.txt").read_text()))
    t = _sh(f"{QE_BIN}/pw.x -nk {np_} -in bands.in", d, d / "bands.out", np_)
    (d / "DONE").write_text(json.dumps({"wall_s": round(t, 1), "np": np_}))
    return xml


EPWREAD_FILES = ("crystal.fmt", "epwdata.fmt", "wigner.fmt", "vmedata.fmt", "dmedata.fmt",
                 "{prefix}.ukk", "{prefix}.kgmap", "{prefix}_0.kmap")


def _epw_dir(m: dict, d: Path, wannier: Path, epmatwp: Path) -> Path:
    """Fresh run dir with the epwread files and a private outdir (symlinks, never writes out/)."""
    if d.exists():
        shutil.rmtree(d)
    d.mkdir(parents=True)
    for f in EPWREAD_FILES:
        f = f.format(prefix=m["prefix"])
        if not (wannier / f).exists():
            raise BandCheckError(f"{wannier / f} missing")
        shutil.copy(wannier / f, d / f)
    od = d / "out"
    od.mkdir()
    os.symlink(Path(m["epw_dir"]) / "out" / f"{m['prefix']}.save", od / f"{m['prefix']}.save")
    os.symlink(epmatwp.resolve(), od / f"{m['prefix']}.epmatwp")
    return od


def run_epw_bands(m: dict, root: Path, label: str, wannier: Path, epmatwp: Path, np_: int) -> Path:
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from eph_converge import epw_input
    d = root / label / "epw_bands"
    if (d / "DONE").exists():
        return d / "band.eig"
    od = _epw_dir(m, d, wannier, epmatwp)
    shutil.copy(root / "kval.txt", d / "kval.txt")
    shutil.copy(root / "gamma.txt", d / "gamma.txt")
    kw = {"outdir": str(od) + "/", "wannierize": False, "nbndsub": int(m["nbndsub"]),
          "epwwrite": False, "epwread": True, "elph": True, "band_plot": True,
          "filkf": "kval.txt", "filqf": "gamma.txt", "efermi_read": True,
          "fermi_energy": float(m["E_F_eV"])}
    if m.get("exclude_bands"):
        kw["bands_skipped"] = f"exclude_bands = {m['exclude_bands']}"
    (d / "epw.in").write_text(epw_input(m, **kw))
    t = _sh(f"{QE_BIN}/epw.x -nk {np_} -in epw.in", d, d / "epw.out", np_)
    (d / "DONE").write_text(json.dumps({"wall_s": round(t, 1), "np": np_}))
    return d / "band.eig"


def run_epw_fermi(m: dict, root: Path, label: str, wannier: Path, epmatwp: Path, np_: int,
                  K=(16, 16, 16), sigma_meV: float = 20.0) -> Path:
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from eph_converge import epw_input, phself_keys
    from ephproto.qlist import write_filqf
    from ephproto.phself import check_epw_pins, parse_epw_namelist
    d = Path(m["epw_dir"]) / "g1" / label / f"epw_ef_K{K[0]}"      # independent of the k set
    if (d / "DONE").exists():
        return d / "epw.out"
    od = _epw_dir(m, d, wannier, epmatwp)
    write_filqf(d / "qlist.txt", [1, 1, 1])
    kw = phself_keys(m, sigma_meV, list(K))
    kw["outdir"] = str(od) + "/"
    (d / "epw.in").write_text(epw_input(m, **kw))
    check_epw_pins(parse_epw_namelist(d / "epw.in"))
    t = _sh(f"{QE_BIN}/epw.x -nk {np_} -in epw.in", d, d / "epw.out", np_)
    (d / "DONE").write_text(json.dumps({"wall_s": round(t, 1), "np": np_}))
    return d / "epw.out"


def stage_run(m: dict, wannier: Path, epmatwp: Path, label: str, np_: int,
              ef_grids: Sequence[int] = (16,), kset: str = "mp5", nrand: int = 1000,
              seed: int = 20261027) -> dict:
    root = Path(m["epw_dir"]) / ("g1" if kset == "mp5" else f"g1_{kset}{nrand}_s{seed}")
    kmeta = write_kpoint_files(root, coarse=m["nk"], kset=kset, nrand=nrand, seed=seed)
    xml = run_pw_bands(m, root, np_)
    beig = run_epw_bands(m, root, label, wannier, epmatwp, np_)
    kd, Ed = parse_pw_xml(xml)
    kw_, Ew = parse_band_eig(beig)
    nexcl = _nexclude(m.get("exclude_bands"))
    cmp = compare_bands(kd, Ed, kw_, Ew, float(m["E_F_eV"]), nexcl, g1_windows())
    efs = {}
    for n in ef_grids:
        efs[n] = parse_epw_fermi(run_epw_fermi(m, root, label, wannier, epmatwp, np_, K=(n, n, n)))
    ef_gate = efs[ef_grids[0]]["ef_fine_eV"]
    res = {"label": label, "wannier": str(wannier), "epmatwp": str(epmatwp), "kpoints": kmeta,
           "compare": cmp, "ef_runs": {f"{n}^3": v for n, v in efs.items()},
           "verdict": g1_verdict(cmp, ef_gate, float(m["E_F_eV"]))}
    (root / label / "g1_result.json").write_text(json.dumps(res, indent=1))
    return res


def _nexclude(s: str | None) -> int:
    """'1-4' -> 4. Only a leading contiguous block 1..n is supported."""
    if not s:
        return 0
    m = re.fullmatch(r"\s*1\s*-\s*(\d+)\s*", s)
    if not m:
        raise BandCheckError(f"exclude_bands {s!r}: only '1-n' supported")
    return int(m.group(1))


# ---------------------------------------------------------------------------------- CLI
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="bandcheck", description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("kpoints")
    a.add_argument("outdir")
    a.add_argument("--n", type=int, nargs=3, default=[5, 5, 5])
    a.add_argument("--coarse", type=int, nargs=3, default=[8, 8, 8])
    c = sub.add_parser("compare")
    c.add_argument("--pw-xml", required=True)
    c.add_argument("--band-eig", required=True)
    c.add_argument("--ef", type=float, required=True)
    c.add_argument("--nexclude", type=int, default=0)
    c.add_argument("--sigma-min-adm", type=float, default=SIGMA_MIN_ADM_EV, help="eV")
    c.add_argument("--epw-out", help="epw.x stdout with the fine-mesh E_F line")
    c.add_argument("--json")
    r = sub.add_parser("run")
    r.add_argument("material")
    r.add_argument("--wannier", required=True)
    r.add_argument("--epmatwp", required=True)
    r.add_argument("--label", required=True)
    r.add_argument("--np", type=int, default=4)
    r.add_argument("--ef-grids", type=int, nargs="+", default=[16])
    r.add_argument("--kset", choices=["mp5", "random"], default="mp5")
    r.add_argument("--nrand", type=int, default=1000)
    r.add_argument("--seed", type=int, default=20261027)
    args = ap.parse_args(argv)

    if args.cmd == "kpoints":
        print(json.dumps(write_kpoint_files(args.outdir, args.n, args.coarse), indent=1))
        return 0
    if args.cmd == "compare":
        kd, Ed = parse_pw_xml(args.pw_xml)
        kw_, Ew = parse_band_eig(args.band_eig)
        cmp = compare_bands(kd, Ed, kw_, Ew, args.ef, args.nexclude, g1_windows(args.sigma_min_adm))
        ef = parse_epw_fermi(args.epw_out)["ef_fine_eV"] if args.epw_out else None
        res = {"compare": cmp, "verdict": g1_verdict(cmp, ef, args.ef)}
    else:
        m = json.loads(Path(args.material).read_text())
        if args.np > 4:
            raise SystemExit("at most 4 ranks")
        res = stage_run(m, Path(args.wannier), Path(args.epmatwp), args.label, args.np,
                        tuple(args.ef_grids), args.kset, args.nrand, args.seed)
    txt = json.dumps(res, indent=1)
    if getattr(args, "json", None):
        Path(args.json).write_text(txt)
    for w in res["compare"]["windows"]:
        print(f"{w['name']:5s} +-{1000 * w['half_width_eV']:.0f} meV tol {1000 * w['tol_eV']:.0f} meV: "
              f"max|dE| {w['max_abs_dE_meV']:.2f} meV, RMS {w['rms_dE_meV']:.2f} meV, "
              f"nn-max {w['max_nn_dE_meV']:.2f} meV, states {w['n_states']}, "
              f"count-mismatch k {w['n_k_count_mismatch']}, pass={w['pass']}")
    print(json.dumps(res["verdict"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
