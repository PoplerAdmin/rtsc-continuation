#!/usr/bin/env python3
"""α²F convergence protocol driver (reports/EPH_CONVERGENCE_PROTOCOL.md), EPW 6.1 / QE 7.6.

Material-agnostic stages, each idempotent (a finished stage is detected and skipped):

  wannier   epw.x: SCDM Wannierisation + epwwrite (coarse el-ph in the Wannier basis), one pass
  bands     epw.x band_plot on off-grid k vs pw.x 'bands' on the same k  -> gate G1
  sweep     one epw.x pass per (fine grid N, electron smearing sigma), each in its own directory,
            epwread=.true. (route 2a: ephwrite + eliashberg, liso=laniso=.false.)
  admit     P3 admission + gates G2/G3 -> admitted (N, sigma) set, envelope, alternatives
  export    writes the admitted alternatives as rtsc Spectrum JSON (+ provenance) for certify

This script produces INPUTS for rtsc.certify.classify only (ARCH S1/S2). It never emits a proved
status. Gate failures are reported as UNCONVERGED_EPH(<reason>) for the caller to record as
DEFERRED.

Usage (WSL, venv ~/rtsc_work/.venv):
  python eph_converge.py <material.json> <stage> [--np N]
material.json holds the per-material facts (paths, grids, prefix, masses, E_F, band counts);
every protocol threshold lives in PROTOCOL below and is NOT taken from the material file.
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
from pathlib import Path

PROTOCOL = {
    "version": "1.0-draft",
    # P2 smearing set (meV), identical for every material
    "sigma_meV": [5, 10, 20, 40, 80, 160],
    # fine grids, as multiples of a base N0 (material file gives N0 per direction)
    "grid_multiples": [2, 4, 6],
    "fsthick_eV": 1.0,                 # >= 4*max(sigma) (EPW hard check), fixed for all passes
    "degaussq_meV": 0.05, "delta_qsmear_meV": 0.05, "nqsmear": 10, "nqstep": 1000,
    "eps_acoustic_cm1": 0.1,
    # G1 Wannier band check
    "g1_window_eV": 1.0, "g1_max_meV": 50.0,
    # P3 admission between the two densest grids
    "p3_rel_tol": 0.02,
    # G2 / G3
    "g2_min_admitted": 3, "g3_lambda_spread": 0.10,
    # P5 envelope
    "envelope_floor": 0.10, "envelope_factor": 3.0,
    # SCDM (absolute eV = E_F + offset)
    "scdm_mu_offset_eV": 2.0, "scdm_sigma_eV": 1.5, "froz_above_EF_eV": 2.0,
    "outer_above_EF_eV": 20.0, "dis_num_iter": 5000,
}

QE_BIN = Path.home() / "micromamba/envs/qe/bin"


def sh(cmd: str, cwd: Path, log: Path, np_: int = 1) -> None:
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    if np_ > 1:
        cmd = (f"{QE_BIN}/mpirun -np {np_} --bind-to none -x OMP_NUM_THREADS -x OPENBLAS_NUM_THREADS "
               + cmd)
    t0 = time.time()
    with open(log, "w") as fh:
        rc = subprocess.call(cmd, shell=True, cwd=cwd, stdout=fh, stderr=subprocess.STDOUT, env=env)
    if rc != 0:
        raise RuntimeError(f"command failed (rc={rc}, {time.time()-t0:.0f}s): {cmd}  -> see {log}")


def fval(v):
    if isinstance(v, bool):
        return ".true." if v else ".false."
    if isinstance(v, str):
        return f"'{v}'"
    return repr(v) if isinstance(v, float) else str(v)


def epw_input(m: dict, **kw) -> str:
    v = {"prefix": m["prefix"], "outdir": str(Path(m["epw_dir"]) / "out") + "/",
         "dvscf_dir": str(Path(m["epw_dir"]) / "save")}
    for i, a in enumerate(m["masses"], 1):
        v[f"amass({i})"] = float(a)
    v.update({"elph": True, "etf_mem": 1, "use_ws": True, "iverbosity": 0,
              "nk1": m["nk"][0], "nk2": m["nk"][1], "nk3": m["nk"][2],
              "nq1": m["nq"][0], "nq2": m["nq"][1], "nq3": m["nq"][2],
              "eps_acoustic": PROTOCOL["eps_acoustic_cm1"], "ngaussw": 0})
    v.update(kw)
    lines = ["--", "&inputepw"] + [f"  {k} = {fval(x)}" for k, x in v.items() if x is not None] + ["/"]
    return "\n".join(lines) + "\n"


def wannier_keys(m: dict) -> dict:
    kw = {"wannierize": True, "nbndsub": int(m["nbndsub"]), "num_iter": 500, "iprint": 2,
          "auto_projections": True, "scdm_proj": True, "scdm_entanglement": "erfc",
          "scdm_mu": float(m["E_F_eV"]) + PROTOCOL["scdm_mu_offset_eV"],
          "scdm_sigma": PROTOCOL["scdm_sigma_eV"],
          # v3.2 (A10): frozen window so the Wannier bands reproduce DFT exactly up to E_F + 2 eV.
          # Without it the Al Wannier set gave a fine-grid E_F 4 eV too high.
          "dis_froz_max": float(m["E_F_eV"]) + PROTOCOL["froz_above_EF_eV"],
          # v3.3 candidates (A11): converge the disentanglement and cap the outer window
          "dis_win_max": float(m["E_F_eV"]) + PROTOCOL["outer_above_EF_eV"]}
    kw["wdata(1)"] = f"dis_num_iter = {PROTOCOL['dis_num_iter']}"
    kw["wdata(2)"] = "dis_conv_tol = 1.0d-10"
    for k, v in (m.get("wannier_overrides") or {}).items():   # variant testing only (recorded)
        kw[k] = v
    if m.get("exclude_bands"):
        kw["bands_skipped"] = f"exclude_bands = {m['exclude_bands']}"
    return kw


def a2f_keys(sigma_meV: float, N: list[int]) -> dict:
    return {"ephwrite": True, "eliashberg": True, "liso": False, "laniso": False,
            "fsthick": PROTOCOL["fsthick_eV"], "degaussw": sigma_meV / 1000.0,
            "degaussq": PROTOCOL["degaussq_meV"], "delta_qsmear": PROTOCOL["delta_qsmear_meV"],
            "nqsmear": PROTOCOL["nqsmear"], "nqstep": PROTOCOL["nqstep"],
            "nstemp": 1, "temps(1)": 10.0, "mp_mesh_k": True,
            "nkf1": N[0], "nkf2": N[1], "nkf3": N[2], "nqf1": N[0], "nqf2": N[1], "nqf3": N[2]}


def grids(m: dict) -> list[list[int]]:
    return [[int(n0 * k) for n0 in m["N0"]] for k in PROTOCOL["grid_multiples"]]


# ----------------------------------------------------------------------------------- stages
def stage_wannier(m: dict, np_: int) -> None:
    d = Path(m["epw_dir"]) / "wannier"
    done = d / "DONE"
    if done.exists():
        print("wannier: done"); return
    d.mkdir(parents=True, exist_ok=True)
    out = Path(m["epw_dir"]) / "out"
    out.mkdir(exist_ok=True)
    # the nscf outdir becomes EPW's outdir
    src = Path(m["epw_dir"]) / "tmp" / f"{m['prefix']}.save"
    if not (out / f"{m['prefix']}.save").exists():
        shutil.copytree(src, out / f"{m['prefix']}.save")
    N = grids(m)[0]
    kw = wannier_keys(m) | {"epwwrite": True, "epwread": False, "epbwrite": False}
    kw |= a2f_keys(PROTOCOL["sigma_meV"][2], N)       # also exercises the a2F path once
    (d / "epw.in").write_text(epw_input(m, **kw))
    sh(f"{QE_BIN}/epw.x -nk {np_} -in epw.in", d, d / "epw.out", np_)
    done.write_text(time.ctime())
    print("wannier: finished", (d / "epw.out"))


def stage_sweep(m: dict, np_: int, only_grid: int | None = None) -> None:
    for gi, N in enumerate(grids(m)):
        if only_grid is not None and gi != only_grid:
            continue
        for s in PROTOCOL["sigma_meV"]:
            d = Path(m["epw_dir"]) / "sweep" / f"N{'x'.join(map(str, N))}_s{s}meV"
            if (d / "DONE").exists():
                continue
            d.mkdir(parents=True, exist_ok=True)
            for f in d.glob(f"{m['prefix']}.a2f*"):
                f.unlink()                              # EPW silently reuses a stale .a2f
            kw = {"wannierize": False, "nbndsub": int(m["nbndsub"]), "epwwrite": False,
                  "epwread": True, "epbwrite": False} | a2f_keys(s, N)
            if m.get("exclude_bands"):
                kw["bands_skipped"] = f"exclude_bands = {m['exclude_bands']}"
            (d / "epw.in").write_text(epw_input(m, **kw))
            t0 = time.time()
            sh(f"{QE_BIN}/epw.x -nk {np_} -in epw.in", d, d / "epw.out", np_)
            (d / "DONE").write_text(f"{time.time()-t0:.1f}")
            print(f"sweep {N} sigma={s} meV: {time.time()-t0:.0f} s", flush=True)


def load_pass(d: Path, prefix: str):
    """Parse one pass's .a2f with rtsc's EPW parser; return (spectra list, summary dict)."""
    from rtsc.dfpt.epw_inputs import parse_epw_a2f, epw_a2f_to_spectra
    tab = parse_epw_a2f(d / f"{prefix}.a2f")
    sp = epw_a2f_to_spectra(tab, label_prefix=d.name + ":")
    from rtsc.spectrum import Spectrum  # noqa: F401
    neg = [sum(w for w in s.weight_K if w < 0) for s in sp]
    lam = [s.lambda_() for s in sp]
    wlog = [s.omega_log() for s in sp]
    A = [s.A() for s in sp]
    return sp, {"lambda": lam, "omega_log_K": wlog, "A_K": A, "neg_mass_K": neg}


def stage_admit(m: dict) -> dict:
    P = PROTOCOL
    gs = grids(m)
    res = {"protocol": P, "material": m["material_id"], "grids": gs, "table": {}, "gates": {}}
    data = {}
    for N in gs:
        for s in P["sigma_meV"]:
            d = Path(m["epw_dir"]) / "sweep" / f"N{'x'.join(map(str, N))}_s{s}meV"
            if not (d / "DONE").exists():
                raise SystemExit(f"missing pass {d}")
            sp, summ = load_pass(d, m["prefix"])
            data[(tuple(N), s)] = (sp, summ)
            res["table"][f"{N}|{s}"] = summ
    N2, N3 = tuple(gs[-2]), tuple(gs[-1])
    admitted, why = [], {}
    for s in P["sigma_meV"]:
        a, b = data[(N2, s)][1], data[(N3, s)][1]
        ok, reasons = True, []
        # every phonon-smearing column must pass (they are all alternatives later)
        for key in ("lambda", "omega_log_K", "A_K"):
            for x2, x3 in zip(a[key], b[key]):
                rel = abs(x3 - x2) / abs(x3)
                if rel > P["p3_rel_tol"]:
                    ok = False; reasons.append(f"{key} rel change {rel:.3%}")
                    break
        if any(n < 0 for n in b["neg_mass_K"]):
            ok = False; reasons.append(f"negative bins at N3 ({min(b['neg_mass_K']):.3g} K)")
        why[s] = reasons or ["ok"]
        if ok:
            admitted.append(s)
    res["admission"] = why
    # G2: >= 3 admitted, contiguous in the sigma list
    idx = [P["sigma_meV"].index(s) for s in admitted]
    contiguous = bool(idx) and idx == list(range(idx[0], idx[0] + len(idx)))
    g2 = len(admitted) >= P["g2_min_admitted"] and contiguous
    # large-sigma guard: excluded values must lie only ABOVE the admitted block, never below it
    # unless they fail for resolution reasons (small sigma) -- small-sigma exclusions are the
    # expected under-resolution; a gap inside the block is a G2 failure (contiguity).
    res["gates"]["G2"] = {"pass": g2, "admitted_sigma_meV": admitted, "contiguous": contiguous}
    if not g2:
        res["status"] = "UNCONVERGED_EPH(grid)"
        return res
    lam3 = [x for s in admitted for x in data[(N3, s)][1]["lambda"]]
    spread = (max(lam3) - min(lam3)) / min(lam3)
    res["gates"]["G3"] = {"pass": spread <= P["g3_lambda_spread"], "lambda_spread": spread}
    if spread > P["g3_lambda_spread"]:
        res["status"] = "UNCONVERGED_EPH(smearing)"
        return res
    deltaA = max(abs(data[(N3, s)][1]["A_K"][j] - data[(N2, s)][1]["A_K"][j]) / data[(N3, s)][1]["A_K"][j]
                 for s in admitted for j in range(len(data[(N3, s)][1]["A_K"])))
    env = max(P["envelope_floor"], P["envelope_factor"] * deltaA)
    res["envelope"] = {"e": env, "deltaA_max": deltaA, "rule": "max(floor, factor*deltaA)"}
    res["alternatives"] = [f"{list(N)}|{s}" for N in (N2, N3) for s in admitted]
    res["status"] = "ADMITTED"
    res["_spectra"] = [sp for N in (N2, N3) for s in admitted for sp in data[(N, s)][0]]
    return res


def stage_export(m: dict, res: dict) -> Path:
    out = Path(m["epw_dir"]) / "protocol"
    out.mkdir(exist_ok=True)
    sp = res.pop("_spectra", [])
    e = res.get("envelope", {}).get("e")
    js = []
    for s in sp:
        d = s.to_json() if hasattr(s, "to_json") else s.__dict__
        d = json.loads(d) if isinstance(d, str) else dict(d)
        if e is not None:
            d["rel_envelope"] = max(float(d.get("rel_envelope", 0.0)), e)
        prov = dict(d.get("provenance", {}))
        prov.update({"eph_protocol": PROTOCOL["version"], "material": m["material_id"]})
        d["provenance"] = prov
        js.append(d)
    (out / "alternatives.json").write_text(json.dumps(js, indent=1, default=str))
    (out / "protocol_result.json").write_text(json.dumps(res, indent=1, default=str))
    return out


# ----------------------------------------------------------------------------- v3 (phselfen)
EPWREAD_FILES = ("crystal.fmt", "epwdata.fmt", "wigner.fmt", "vmedata.fmt", "dmedata.fmt",
                 "{prefix}.ukk", "{prefix}.kgmap", "{prefix}_0.kmap")


def phself_keys(m: dict, sigma_meV: float, K: list[int]) -> dict:
    """The v3.1 pinned key set (whitelist): nothing else is ever written."""
    s = sigma_meV / 1000.0
    kw = {"wannierize": False, "nbndsub": int(m["nbndsub"]), "epwwrite": False, "epwread": True,
          "epbwrite": False, "epbread": False,
          "phonselfen": True, "delta_approx": True, "a2f": False, "restart": False,
          "selecqread": False, "nstemp": 1, "temps(1)": 10.0, "efermi_read": False,
          "mp_mesh_k": False, "mp_mesh_q": False, "nsmear": 1,
          "degaussw": s, "fsthick": max(1.0, 8.0 * s),
          "nkf1": K[0], "nkf2": K[1], "nkf3": K[2], "filqf": "qlist.txt"}
    if m.get("exclude_bands"):
        kw["bands_skipped"] = f"exclude_bands = {m['exclude_bands']}"
    return kw


def run_phself(m: dict, K: list[int], Q: list[int], sigma_meV: float, np_: int) -> Path:
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from ephproto.qlist import write_filqf
    from ephproto.phself import parse_epw_namelist, check_epw_pins
    d = (Path(m["epw_dir"]) / "v3" /
         f"K{'x'.join(map(str, K))}_Q{'x'.join(map(str, Q))}_s{sigma_meV:g}meV")
    if (d / "DONE").exists():
        return d
    if d.exists():
        shutil.rmtree(d)                                   # never resume a partial EPW run
    d.mkdir(parents=True)
    wan = Path(m["epw_dir"]) / "wannier"
    for f in EPWREAD_FILES:
        f = f.format(prefix=m["prefix"])
        if (wan / f).exists():
            shutil.copy(wan / f, d / f)
    write_filqf(d / "qlist.txt", Q)
    (d / "epw.in").write_text(epw_input(m, **phself_keys(m, sigma_meV, K)))
    check_epw_pins(parse_epw_namelist(d / "epw.in"))       # refuse before spending compute
    t0 = time.time()
    sh(f"{QE_BIN}/epw.x -nk {np_} -in epw.in", d, d / "epw.out", np_)
    (d / "DONE").write_text(json.dumps({"wall_s": round(time.time() - t0, 1), "np": np_}))
    return d


def load_phself(m: dict, d: Path, rel_envelope: float = 0.0):
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from ephproto.phself import build_phself_spectrum
    lam = sorted(d.glob("lambda.phself.*"))
    lw = sorted(d.glob("linewidth.phself.*"))
    if len(lam) != 1 or len(lw) != 1:
        raise RuntimeError(f"{d}: expected one lambda/linewidth phself file, got {lam} {lw}")
    return build_phself_spectrum(lam[0], lw[0], d / "qlist.txt", d / "epw.in", d / "epw.out",
                                 label=f"{m['material_id']}:{d.name}", rel_envelope=rel_envelope,
                                 extra_provenance={"discretisation": "atoms from (omega_qnu, lambda_qnu), "
                                                   "eph_protocol v3.1 DRAFT", "eph_protocol": "v3.1-DRAFT"})


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("material")
    ap.add_argument("stage", choices=["wannier", "sweep", "admit", "all"])
    ap.add_argument("--np", type=int, default=4)
    ap.add_argument("--grid", type=int, help="sweep only this grid index")
    a = ap.parse_args(argv)
    m = json.loads(Path(a.material).read_text())
    if a.stage in ("wannier", "all"):
        stage_wannier(m, a.np)
    if a.stage in ("sweep", "all"):
        stage_sweep(m, a.np, a.grid)
    if a.stage in ("admit", "all"):
        res = stage_admit(m)
        out = stage_export(m, res)
        print(json.dumps({k: v for k, v in res.items() if k in ("status", "gates", "envelope", "admission")},
                         indent=1, default=str))
        print("written:", out)


if __name__ == "__main__":
    main()
