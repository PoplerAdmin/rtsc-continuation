#!/usr/bin/env python3
"""extract_perf.py -- read-only performance-metrics extractor for the rtsc benchmark campaign.

HOW TO RERUN (idempotent; it only reads run files and rewrites the one JSON output):

  From Windows (Git Bash):
    MSYS_NO_PATHCONV=1 wsl.exe -d Ubuntu --exec bash -lc \
      '~/rtsc_work/.venv/bin/python /mnt/g/ClaudeWorkspace/RTSC_Project/bench_tools/perf/extract_perf.py' \
      2>&1 | tr -d '\\000'

  Inside WSL:
    ~/rtsc_work/.venv/bin/python /mnt/g/ClaudeWorkspace/RTSC_Project/bench_tools/perf/extract_perf.py \
        [--work ~/rtsc_work] [--repo /mnt/g/ClaudeWorkspace/RTSC_Project] [--out <json>] [--no-proc]

  Output: <repo>/bench_tools/perf/perf_metrics.json (overwritten each run).
  Rerun whenever Nb / H3S / MgB2 (or the validation chain, tetra, EPW runs) progress; missing or
  partial runs are reported as such (status "missing" / "running" / "incomplete"), never as errors.

What it reads (never writes, never starts/stops processes, never runs QE):
  ~/rtsc_work/bench/<M>/{config.toml, run*.log, run/<id>/{record.json,verdict.json,dfpt/*}}
  ~/rtsc_work/bench/queue.log
  ~/rtsc_work/validate/<M>/<tag>/{config.toml, run.log, run/<id>/record.json}, validate/chain.log
  ~/rtsc_work/tetra_Al/**.out*, ~/rtsc_work/epw_Al/**.out, *.log, timing.txt
  <repo>/reports/data/pruned_manifests/*.tsv  (bytes of scratch files that were deleted)
  /proc/loadavg and /proc/<pid>/{comm,status} (one light snapshot of what is running now; --no-proc skips)

Conventions:
  * QE "CPU" in timing tables is the CPU time of MPI rank 0 only; WALL is wall-clock.
    cpu_wall_ratio = CPU/WALL of rank 0. ~1.0 = rank 0 was computing (or MPI busy-polling) the whole
    time; well below 1 = rank 0 waited (disk I/O, descheduled because the machine was oversubscribed).
  * `time` user/sys in run*.log are summed over ALL processes of `rtsc run` (all MPI ranks).
  * Per-q ph.x wall: from the cumulative "PHONON : ... WALL" line printed after each q (exact to 0.01 s
    below 1 h; QE prints only h+m above 1 h, so resolution is 60 s there) and, independently, from the
    mtimes of the dyn<i> files (finish time of q i). Both are given.
  * All paths are written with ~/ ; /home/<user> is never stored.
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import re
import sys
import time
from pathlib import Path

MATERIALS = ["Al", "Pb", "Nb", "H3S", "MgB2"]
STEPS = ["relax", "scf_fit", "scf", "ph", "q2r", "matdyn"]
SCHEMA = "rtsc.perf_metrics/1"
RUNNING_GRACE_S = 900          # a QE .out without JOB DONE modified within this window counts as running

HOME = os.path.expanduser("~")
_HOME_RE = re.compile(r"/home/[^/\s\"']+")


# ----------------------------------------------------------------------------- small helpers
def sanitize(obj):
    """Replace /home/<user> (and $HOME) by ~ everywhere in the output."""
    if isinstance(obj, str):
        s = obj.replace(HOME, "~") if HOME not in ("", "/") else obj
        return _HOME_RE.sub("~", s)
    if isinstance(obj, dict):
        return {sanitize(k): sanitize(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [sanitize(v) for v in obj]
    return obj


def tpath(p) -> str:
    return sanitize(str(p))


def iso(ts):
    if ts is None:
        return None
    return dt.datetime.fromtimestamp(ts).astimezone().isoformat(timespec="seconds")


def read_text(p: Path, limit: int | None = None) -> str | None:
    try:
        with open(p, "r", errors="replace") as f:
            return f.read(limit) if limit else f.read()
    except (OSError, ValueError):
        return None


def load_json(p: Path):
    try:
        with open(p) as f:
            return json.load(f)
    except (OSError, ValueError):
        return None


def mtime(p: Path):
    try:
        return p.stat().st_mtime
    except OSError:
        return None


def rnd(x, n=3):
    return None if x is None else round(float(x), n)


def ratio(a, b, n=3):
    try:
        if a is None or b in (None, 0):
            return None
        return round(float(a) / float(b), n)
    except (TypeError, ValueError, ZeroDivisionError):
        return None


_MON = {m: i for i, m in enumerate(
    ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"], 1)}


def qe_time_to_s(s: str):
    """'53m23.30s' / '1h 3m' / '2m 8.32s' / '0.12s' / '1d 2h' -> (seconds, resolution_s)."""
    s = s.strip()
    m = re.fullmatch(r"(?:(\d+)d)?\s*(?:(\d+)h)?\s*(?:(\d+)m)?\s*(?:([\d.]+)s)?", s)
    if not m or not any(m.groups()):
        return None, None
    d, h, mi, se = m.groups()
    tot = (int(d or 0) * 86400 + int(h or 0) * 3600 + int(mi or 0) * 60 + float(se or 0.0))
    res = 0.01 if se is not None else (60.0 if mi is not None or h is not None else 3600.0)
    return tot, res


_TIMING_RE = re.compile(
    r"^\s+(\S.*?)\s+:\s+([0-9dhms. ]+?)\s*CPU\s+([0-9dhms. ]+?)\s*WALL(?:\s*\(\s*(\d+)\s*calls\))?\s*$")
_START_RE = re.compile(r"Program\s+(\S+)\s+v\.(\S+)\s+starts on\s+(\d+)(\w{3})\s*(\d{4})\s+at\s+([\d ]+):([\d ]+):([\d ]+)")
_END_RE = re.compile(r"terminated on:\s+([\d ]+):([\d ]+):([\d ]+)\s+(\d+)(\w{3})\s*(\d{4})")


def _local_ts(y, mon, d, hh, mm, ss):
    try:
        return dt.datetime(int(y), _MON[mon], int(d), int(hh), int(mm), int(ss)).timestamp()
    except (KeyError, ValueError):
        return None


# ----------------------------------------------------------------------------- QE output parser
def parse_qe_out(path: Path, now: float) -> dict | None:
    txt = read_text(path)
    if txt is None:
        return None
    info: dict = {"path": tpath(path), "bytes": len(txt.encode(errors="replace")),
                  "mtime": iso(mtime(path))}
    m = _START_RE.search(txt)
    if m:
        info["program"], info["version"] = m.group(1), m.group(2)
        info["start_ts"] = _local_ts(m.group(5), m.group(4), m.group(3), m.group(6), m.group(7), m.group(8))
    else:
        info["program"], info["start_ts"] = None, None
    ends = _END_RE.findall(txt)
    info["end_ts"] = (_local_ts(ends[-1][5], ends[-1][4], ends[-1][3], *ends[-1][:3]) if ends else None)
    info["job_done"] = "JOB DONE" in txt
    mt = mtime(path)
    if info["job_done"]:
        info["status"] = "done"
    elif mt is not None and now - mt < RUNNING_GRACE_S:
        info["status"] = "running"
    else:
        info["status"] = "incomplete"
    info["end_or_last_write_ts"] = info["end_ts"] or mt
    for key, rx, typ in [
        ("mpi_procs", r"Number of MPI processes:\s+(\d+)", int),
        ("threads_per_mpi", r"Threads/MPI process:\s+(\d+)", int),
        ("npool", r"K-points division:\s+npool\s+=\s+(\d+)", int),
        ("n_ks_states", r"number of Kohn-Sham states=\s+(\d+)", int),
        ("n_k_points", r"number of k points=\s+(\d+)", int),
        ("ecutwfc_Ry", r"kinetic-energy cut-?off\s+=\s+([\d.]+)", float),
        ("ecutrho_Ry", r"charge density cut-?off\s+=\s+([\d.]+)", float),
        ("n_electrons", r"number of electrons\s+=\s+([\d.]+)", float),
        ("mem_available_MiB", r"(\d+) MiB available memory", int),
    ]:
        mm = re.search(rx, txt)
        info[key] = typ(mm.group(1)) if mm else None
    if info.get("threads_per_mpi") is None and info.get("mpi_procs") is not None:
        info["threads_per_mpi"] = 1
    mm = re.search(r"Dense\s+grid:\s+(\d+) G-vectors\s+FFT dimensions: \(\s*(\d+),\s*(\d+),\s*(\d+)\)", txt)
    if mm:
        info["dense_G_vectors"] = int(mm.group(1))
        info["fft_dims"] = [int(mm.group(i)) for i in (2, 3, 4)]
    mm = re.search(r"Estimated max dynamical RAM per process >\s+([\d.]+)\s*(\w+)", txt)
    if mm:
        info["est_ram_per_process"] = f"{mm.group(1)} {mm.group(2)}"
    scf_it = [int(x) for x in re.findall(r"convergence has been achieved in\s+(\d+) iterations", txt)]
    if scf_it:
        info["scf_iterations"] = scf_it
        info["scf_iterations_total"] = sum(scf_it)
    mm = re.search(r"bfgs converged in\s+(\d+) scf cycles and\s+(\d+) bfgs steps", txt)
    if mm:
        info["bfgs"] = {"scf_cycles": int(mm.group(1)), "bfgs_steps": int(mm.group(2))}
    # timing lines
    routines: dict = {}
    prog_lines: list = []
    prog = (info.get("program") or "").upper()
    for line in txt.splitlines():
        t = _TIMING_RE.match(line)
        if not t:
            continue
        name, cpu_s, wall_s, calls = t.group(1).strip(), t.group(2), t.group(3), t.group(4)
        cpu, cres = qe_time_to_s(cpu_s)
        wall, wres = qe_time_to_s(wall_s)
        if cpu is None or wall is None:
            continue
        if calls is None and (name.upper() == prog or (not prog and name.isupper())):
            prog_lines.append({"cpu_s": cpu, "wall_s": wall, "resolution_s": max(cres, wres)})
        elif calls is not None:
            routines[name] = {"cpu_s": cpu, "wall_s": wall, "calls": int(calls)}
    info["program_cumulative"] = prog_lines
    if prog_lines:
        last = prog_lines[-1]
        info["program_cpu_s"], info["program_wall_s"] = last["cpu_s"], last["wall_s"]
        info["program_time_resolution_s"] = last["resolution_s"]
        info["cpu_wall_ratio_rank0"] = ratio(last["cpu_s"], last["wall_s"])
    else:
        info["program_cpu_s"] = info["program_wall_s"] = None
        info["cpu_wall_ratio_rank0"] = None
    info["routines"] = routines
    info["_text"] = txt          # removed before output
    return info


def top_routines(routines: dict, n=15, exclude=()):
    items = [(k, v) for k, v in routines.items() if k not in exclude]
    items.sort(key=lambda kv: -kv[1]["wall_s"])
    out = []
    for k, v in items[:n]:
        out.append({"routine": k, "wall_s": rnd(v["wall_s"], 2), "cpu_s": rnd(v["cpu_s"], 2),
                    "calls": v["calls"], "cpu_wall_ratio": ratio(v["cpu_s"], v["wall_s"])})
    return out


# ----------------------------------------------------------------------------- ph.x specifics
def ph_breakdown(qe: dict) -> dict | None:
    """Partition of ph.x wall time into its top-level parts (rank 0 timers)."""
    r = qe.get("routines") or {}
    total = qe.get("program_wall_s")
    if not r or not total:
        return None
    w = lambda k: r.get(k, {}).get("wall_s")
    parts = {
        "nscf_per_q (electrons)": w("electrons"),
        "linear_response (phqscf)": w("phqscf"),
        "elphsum (e-ph interpolation/sum)": w("elphsum"),
        "phq_init": w("phq_init"),
        "phq_setup": w("phq_setup"),
        "dynmat0": w("dynmat0"),
        "init_run (nscf init)": w("init_run"),
    }
    parts = {k: v for k, v in parts.items() if v is not None}
    known = sum(parts.values())
    parts["other (total - listed)"] = total - known
    out = {"program_wall_s": total,
           "parts": {k: {"wall_s": rnd(v, 2), "frac_of_ph_wall": ratio(v, total, 4)} for k, v in parts.items()}}
    # nested pieces of interest (fractions of ph wall; they overlap with the parts above)
    nested = {}
    for k in ["solve_linter", "cgsolve", "ch_psi", "h_psi", "vloc_psi", "fftw", "fft", "ffts",
              "incdrhoscf", "apply_dpot_b", "dvqpsi_us", "elphel", "c_bands", "cegterg",
              "davcio", "write_rec", "init_us_2", "drhodv", "psymdvscf", "mix_pot", "dv_of_drho"]:
        if k in r:
            nested[k] = {"wall_s": rnd(r[k]["wall_s"], 2), "cpu_s": rnd(r[k]["cpu_s"], 2),
                         "calls": r[k]["calls"], "frac_of_ph_wall": ratio(r[k]["wall_s"], total, 4),
                         "cpu_wall_ratio": ratio(r[k]["cpu_s"], r[k]["wall_s"])}
    out["nested_routines"] = nested
    io = [k for k in ("davcio", "write_rec") if k in r]
    io_wall = sum(r[k]["wall_s"] for k in io)
    io_cpu = sum(r[k]["cpu_s"] for k in io)
    out["io_timed"] = {"routines": io, "wall_s": rnd(io_wall, 2), "cpu_s": rnd(io_cpu, 2),
                       "frac_of_ph_wall": ratio(io_wall, total, 4),
                       "wall_minus_cpu_s": rnd(io_wall - io_cpu, 2),
                       "note": "only I/O that QE times explicitly on rank 0; page-cache writes are mostly hidden"}
    if "elphel" in r and "phqscf" in r:
        s_all = known + r["elphel"]["wall_s"]
        out["elphel_nesting_check"] = {
            "sum_parts_plus_elphel_s": rnd(s_all, 2), "ph_wall_s": rnd(total, 2),
            "inference": ("elphel is nested inside phqscf (adding it would exceed the ph total)"
                          if s_all > total * 1.002 else "elphel may be outside phqscf")}
    return out


def ph_q_list(txt: str):
    """q list printed in the ph.x header (fresh run: '( 29 q-points):'; recover run with a q range:
    '12 /  29 q-points for this run, from 18 to 29:')."""
    head = txt[:200000]
    g = re.search(r"Dynamical matrices for \(\s*(\d+),\s*(\d+),\s*(\d+)\)", head)
    grid = [int(g.group(i)) for i in (1, 2, 3)] if g else None
    n = re.search(r"\(\s*(\d+) q-points\):", head) or re.search(r"/\s*(\d+) q-points for this run", head)
    nq = int(n.group(1)) if n else None
    m = re.search(r"\n\s*N\s+xq\(1\)\s+xq\(2\)\s+xq\(3\)\s*\n((?:\s+\d+\s+\S+\s+\S+\s+\S+[ \t]*\n)+)", head)
    qs = []
    if m:
        for line in m.group(1).strip().splitlines():
            p = line.split()
            try:
                qs.append((int(p[0]), [float(x) for x in p[1:4]]))
            except ValueError:
                break
    return grid, nq, qs


def ph_qblocks(qe: dict, marks: dict | None = None) -> dict:
    """Per-q structure of a ph.x output.

    Timing model (verified on the Al ph.out): ph.x prints the cumulative "PHONON : cpu wall" line inside
    each q block right AFTER that q's nscf and BEFORE its linear response. `marks` maps q index -> mtime
    of a file written when that q is completely finished (elph_dir/elph.inp_lambda.<i> for the e-ph
    runs, <fildyn><i> otherwise). With both:
        wall_s       = mark(q_i) - mark(q_{i-1})     (first q of this output: - program start)
        nscf_wall_s  = (start + PHONON_i) - mark(q_{i-1})
        dfpt_elph_wall_s = mark(q_i) - (start + PHONON_i)
    The rank-0 CPU/WALL ratio is given for the window between consecutive PHONON lines.
    """
    marks = marks or {}
    start_ts = qe.get("start_ts")
    txt = qe.get("_text") or ""
    grid, nq, qlist = ph_q_list(txt)
    if nq is None and qlist:
        nq = len(qlist)
    rng = re.search(r"q-points for this run, from\s+(\d+) to\s+(\d+)", txt[:200000])
    pos = [m.start() for m in re.finditer(r"Calculation of q =", txt)]
    # cumulative PHONON lines with their text positions
    cum = []
    for m in re.finditer(r"^\s+PHONON\s+:\s+([0-9dhms. ]+?)\s*CPU\s+([0-9dhms. ]+?)\s*WALL\s*$", txt, re.M):
        c, cr = qe_time_to_s(m.group(1))
        w, wr = qe_time_to_s(m.group(2))
        cum.append((m.start(), c, w, max(cr or 0, wr or 0)))
    blocks = []
    bounds = pos + [len(txt)]
    own_lines = []
    prev_end = start_ts
    for bi, start in enumerate(pos):
        end = bounds[bi + 1]
        b = txt[start:end]
        mq = re.match(r"Calculation of q =\s+(\S+)\s+(\S+)\s+(\S+)", b)
        xq = [float(mq.group(i)) for i in (1, 2, 3)] if mq else None
        idx = None
        if xq and qlist:
            best = min(qlist, key=lambda t: sum((a - c) ** 2 for a, c in zip(t[1], xq)))
            if sum((a - c) ** 2 for a, c in zip(best[1], xq)) < 1e-10:
                idx = best[0]
        reps = re.findall(r"Representation\s+(\d+)\s+(\d+) modes?\s+-\s+(To be done|Done|Not done in this run)", b)
        mk = re.search(r"number of k points=\s+(\d+)", b)
        md = re.search(r"avg # of iterations =\s+([\d.]+)", b)
        iters = re.findall(r"iter #\s+(\d+) total cpu time :\s+([\d.]+) secs\s+av\.it\.:\s+([\d.]+)", b)
        conv = len(re.findall(r"Convergence has been achieved", b))
        star = re.search(r"Number of q in the star =\s+(\d+)", b)
        bcum = [c for c in cum if start <= c[0] < end]
        blk = {"q_index": idx, "xq": xq,
               "n_irreps": len(reps) or None,
               "n_modes": sum(int(r[1]) for r in reps) if reps else None,
               "irreps_already_done_at_start": sum(1 for r in reps if r[2] == "Done"),
               "nscf_k_points": int(mk.group(1)) if mk else None,
               "nscf_davidson_avg_iter": float(md.group(1)) if md else None,
               "dfpt_scf_iterations": len(iters),
               "irreps_converged": conv,
               "avg_linear_solver_iter": rnd(sum(float(i[2]) for i in iters) / len(iters), 2) if iters else None,
               "star_size": int(star.group(1)) if star else None,
               "has_elph": "Electron-phonon interaction" in b}
        own = bcum[0] if bcum else None
        own_lines.append((own, bcum[-1] if len(bcum) > 1 else None))
        t_end = marks.get(idx)
        if t_end is not None and start_ts is not None and t_end < start_ts - 1:
            t_end = None        # finished in an earlier (pre-restart) run
        blk["complete"] = t_end is not None
        if t_end is not None and prev_end is not None:
            blk["wall_s"] = rnd(t_end - prev_end, 1)
            if own is not None and start_ts is not None:
                t_dfpt = start_ts + own[2]
                blk["nscf_wall_s"] = rnd(t_dfpt - prev_end, 1)
                blk["dfpt_elph_wall_s"] = rnd(t_end - t_dfpt, 1)
                blk["time_resolution_s"] = max(1.0, own[3])
        if blk.get("wall_s") and blk["nscf_k_points"] and blk["dfpt_scf_iterations"]:
            blk["wall_ms_per_kpoint_per_dfpt_iter"] = rnd(
                1000.0 * blk["wall_s"] / (blk["nscf_k_points"] * blk["dfpt_scf_iterations"]), 3)
        if blk["irreps_already_done_at_start"]:
            blk["partial_resumed_q"] = True    # some irreps were computed before the restart
        prev_end = t_end if t_end is not None else None
        blocks.append(blk)
    # CPU/WALL ratio windows between consecutive PHONON lines
    for i, blk in enumerate(blocks):
        own, final = own_lines[i]
        nxt = own_lines[i + 1][0] if i + 1 < len(blocks) else final
        if own is not None and nxt is not None and nxt[2] > own[2]:
            blk["window_cpu_s_rank0"] = rnd(nxt[1] - own[1], 2)
            blk["window_wall_s"] = rnd(nxt[2] - own[2], 2)
            blk["window_cpu_wall_ratio_rank0"] = ratio(nxt[1] - own[1], nxt[2] - own[2])
    restart = bool(re.search(r"Restart in Phonon calculation", txt))
    return {"q_grid": grid, "n_irreducible_q": nq, "q_blocks_in_this_output": len(blocks),
            "restarted_output (recover)": restart,
            "q_range_this_run": [int(rng.group(1)), int(rng.group(2))] if rng else None,
            "per_q_definitions": {
                "wall_s": "finish(q_i) - finish(q_{i-1}) from completion-file mtimes (first q: - ph.x start)",
                "nscf_wall_s": "nscf part of q_i (from the PHONON line printed after the nscf)",
                "dfpt_elph_wall_s": "linear response + e-ph part of q_i",
                "window_*": "between this q's PHONON line and the next one (= DFPT+e-ph of q_i + nscf of q_{i+1})"},
            "per_q": blocks}


def _fit_origin(xs, ys):
    sxx = sum(x * x for x in xs)
    if not xs or sxx <= 0:
        return None
    a = sum(x * y for x, y in zip(xs, ys)) / sxx
    ybar = sum(ys) / len(ys)
    ss_tot = sum((y - ybar) ** 2 for y in ys) or 1.0
    ss_res = sum((y - a * x) ** 2 for x, y in zip(xs, ys))
    return a, 1 - ss_res / ss_tot


def per_q_summary(per_q: list) -> dict | None:
    pq = [b for b in per_q if b.get("wall_s")]
    if not pq:
        return None
    ws = [b["wall_s"] for b in pq]
    out = {"n_q_timed": len(pq), "wall_s_sum": rnd(sum(ws), 1), "wall_s_min": min(ws), "wall_s_max": max(ws),
           "wall_s_mean": rnd(sum(ws) / len(ws), 1),
           "nscf_wall_s_sum": rnd(sum(b.get("nscf_wall_s") or 0 for b in pq), 1),
           "dfpt_elph_wall_s_sum": rnd(sum(b.get("dfpt_elph_wall_s") or 0 for b in pq), 1),
           "dfpt_iterations_total": sum(b["dfpt_scf_iterations"] for b in pq),
           "irreps_total": sum(b["n_irreps"] or 0 for b in pq),
           "nscf_k_points_range": [min(b["nscf_k_points"] or 0 for b in pq), max(b["nscf_k_points"] or 0 for b in pq)]}
    fits = {}
    for name, f in [("nscf_k_points", lambda b: (b["nscf_k_points"] or 0)),
                    ("nscf_k_points*dfpt_iterations", lambda b: (b["nscf_k_points"] or 0) * b["dfpt_scf_iterations"])]:
        r = _fit_origin([f(b) for b in pq], ws)
        if r:
            fits[f"wall_s ~ a*{name}"] = {"a_ms": rnd(1000 * r[0], 4), "R2": rnd(r[1], 3)}
    dq = [b for b in pq if b.get("dfpt_elph_wall_s") and b["dfpt_scf_iterations"]]
    r = _fit_origin([(b["nscf_k_points"] or 0) * b["dfpt_scf_iterations"] for b in dq],
                    [b["dfpt_elph_wall_s"] for b in dq])
    if r:
        fits["dfpt_elph_wall_s ~ a*nscf_k_points*dfpt_iterations"] = {"a_ms": rnd(1000 * r[0], 4), "R2": rnd(r[1], 3)}
    nq = [b for b in pq if b.get("nscf_wall_s") and b["nscf_k_points"]]
    r = _fit_origin([b["nscf_k_points"] for b in nq], [b["nscf_wall_s"] for b in nq])
    if r:
        fits["nscf_wall_s ~ a*nscf_k_points"] = {"a_ms": rnd(1000 * r[0], 4), "R2": rnd(r[1], 3)}
    out["scaling_fits (through origin)"] = fits
    low = [{"q_index": b["q_index"], "window_cpu_wall_ratio_rank0": b.get("window_cpu_wall_ratio_rank0"),
            "window_wall_s": b.get("window_wall_s")} for b in per_q
           if b.get("window_cpu_wall_ratio_rank0") is not None and b["window_cpu_wall_ratio_rank0"] < 0.7]
    out["low_cpu_wall_windows (<0.7; stalls)"] = low
    return out


def completion_marks(d: Path, patterns: list) -> dict:
    """q index -> mtime of the first existing completion file per pattern ('{i}' placeholder)."""
    out = {}
    try:
        names = os.listdir(d)
    except OSError:
        return out
    for pat in patterns:
        rx = re.compile("^" + re.escape(pat).replace(r"\{i\}", r"(\d+)") + "$")
        for n in names:
            m = rx.match(n)
            if m and int(m.group(1)) not in out:
                try:
                    st = (d / n).stat()
                except OSError:
                    continue
                if st.st_size > 0:
                    out[int(m.group(1))] = st.st_mtime
        if out:
            break
    return out


def marks_timeline(marks: dict, restarts_ts: list) -> list:
    out = []
    keys = sorted(marks)
    for a, b in zip([None] + keys, keys):
        d = {"q_index": b, "finished": iso(marks[b])}
        if a is not None:
            d["since_prev_s"] = rnd(marks[b] - marks[a], 1)
            if any(marks[a] < r <= marks[b] for r in restarts_ts):
                d["spans_restart"] = True
        out.append(d)
    return out


def dyn_timeline(dfpt: Path, prefix: str = "dyn"):
    out = []
    try:
        names = os.listdir(dfpt)
    except OSError:
        return out
    rx = re.compile(rf"^{re.escape(prefix)}(\d+)$")
    for n in names:
        m = rx.match(n)
        if m:
            p = dfpt / n
            try:
                st = p.stat()
            except OSError:
                continue
            out.append({"i": int(m.group(1)), "mtime_ts": st.st_mtime, "bytes": st.st_size})
    out.sort(key=lambda d: d["i"])
    for a, b in zip(out, out[1:]):
        if b["bytes"] > 0 and a["mtime_ts"] is not None:
            b["wall_since_prev_s"] = rnd(b["mtime_ts"] - a["mtime_ts"], 1)
    for d in out:
        d["mtime"] = iso(d["mtime_ts"])
        d["finished"] = d["bytes"] > 0
    return out


# ----------------------------------------------------------------------------- run logs, configs, queue
def parse_time_block(txt: str) -> dict:
    out = {}
    for k in ("real", "user", "sys"):
        m = re.search(rf"^{k}\s+(\d+)m([\d.]+)s", txt, re.M)
        if m:
            out[f"{k}_s"] = round(int(m.group(1)) * 60 + float(m.group(2)), 3)
    m = re.search(r"^wall_s\s+([\d.]+)", txt, re.M)
    if m:
        out["wall_s"] = float(m.group(1))
    counts = dict((k, int(v)) for k, v in re.findall(
        r"^\s+(REJECTED_PROVEN|CERTIFIED_CANDIDATE|UNDECIDED|EXCLUDED_ASSUMPTION|DEFERRED|ERROR|total)\s+(\d+)", txt, re.M))
    if counts:
        out["verdict_counts"] = counts
    return out


def load_toml(p: Path):
    try:
        import tomllib
        with open(p, "rb") as f:
            return tomllib.load(f)
    except Exception:
        return None


def config_summary(cfg: dict | None) -> dict | None:
    if not cfg:
        return None
    d = cfg.get("dfpt", {}) or {}
    la, se = d.get("launcher", {}) or {}, d.get("settings", {}) or {}
    ce = cfg.get("certify", {}) or {}
    env = la.get("env", {}) or {}
    return {
        "launcher": {"kind": la.get("kind"), "nprocs": la.get("nprocs"), "npool": la.get("npool"),
                     "OMP_NUM_THREADS": env.get("OMP_NUM_THREADS"), "mpi_args": la.get("mpi_args")},
        "settings": {k: se.get(k) for k in ("ecutwfc_Ry", "ecutrho_Ry", "smearing", "degauss_Ry", "k_scf",
                                            "k_fit", "q", "k_matdyn", "el_ph_sigma_Ry", "el_ph_nsigma",
                                            "tr2_ph", "do_relax", "fildvscf", "smearing_admission",
                                            "timeouts_s", "extra_ph", "extra_pw_control") if k in se},
        "pseudos": (d.get("pseudos", {}) or {}).get("files"),
        "certify": {k: ce.get(k) for k in ("T_target_K", "mu_hi", "mode", "cutoff_rule", "omega_c_factor",
                                           "time_budget_s", "rel_tol") if k in ce},
        "run": {k: (cfg.get("run", {}) or {}).get(k) for k in ("workers", "resume")},
    }


def parse_queue_log(p: Path) -> list:
    txt = read_text(p) or ""
    ev = []
    for line in txt.splitlines():
        m = re.match(r"^(\w{3} \w{3}\s+\d+ \d\d:\d\d:\d\d \w+ \d{4})\s+(START|END|RESUME)\s+(.*)$", line)
        if not m:
            continue
        stamp = m.group(1)
        try:
            parts = stamp.split()
            ts = dt.datetime.strptime(" ".join(parts[:4] + parts[5:]), "%a %b %d %H:%M:%S %Y").timestamp()
        except ValueError:
            ts = None
        rest = m.group(3).strip()
        mat = re.match(r"(?:after [^)]*\)\s*)?(\w+)", rest)
        mm = re.match(r"^(\w+)", rest)
        ev.append({"ts": ts, "time": iso(ts), "event": m.group(2), "text": rest,
                   "material": (mm.group(1).rstrip(":") if mm and mm.group(1) in MATERIALS else
                                next((x for x in MATERIALS if re.search(rf"\b{x}\b", rest)), None))})
    return ev


def parse_chain_log(p: Path) -> list:
    txt = read_text(p) or ""
    ev = []
    for line in txt.splitlines():
        m = re.match(r"^(\d{4}-\d\d-\d\d \d\d:\d\d:\d\d)\s+(START|END)\s+(\S+)\s+(.*)$", line)
        if m:
            ts = dt.datetime.strptime(m.group(1), "%Y-%m-%d %H:%M:%S").timestamp()
            ev.append({"ts": ts, "time": iso(ts), "event": m.group(2), "material": m.group(3), "text": m.group(4)})
    return ev


# ----------------------------------------------------------------------------- disk accounting
def categorize(path: str) -> str:
    n = path.rsplit("/", 1)[-1]
    if re.search(r"\.wfc\d*$", n) or ".wfc" in n:
        return "wavefunctions (wfc)"
    if "dvscf" in n:
        return "dvscf"
    if re.search(r"\.(bar|dwf|mixd|prd|ebar|dbecsum|recover)\d*$", n) or ".bar" in n or ".dwf" in n or ".mixd" in n:
        return "ph scratch (bar/dwf/mixd/prd)"
    if "/elph_dir/" in path or "a2Fsave" in path or "a2Fq2r" in n:
        return "e-ph (elph_dir/a2F)"
    if "charge-density" in n or n.endswith(".xml") or "/pwscf.save/" in path or ".save/" in path:
        return "save (xml/charge/pseudo)"
    if "phsave" in path:
        return "phsave"
    if "epmat" in n or n.endswith(".fmt") or n in ("pwscf.mmn", "pwscf.ukk", "pwscf.chk"):
        return "epw/wannier"
    return "other"


def walk_sizes(root: Path, max_files=200000) -> dict:
    tot, by_cat, per_q, n = 0, {}, {}, 0
    if not root.exists():
        return {"exists": False}
    for dp, dn, fn in os.walk(root, followlinks=False):
        for f in fn:
            p = os.path.join(dp, f)
            try:
                st = os.lstat(p)
            except OSError:
                continue
            if not os.path.isfile(p) or os.path.islink(p):
                continue
            n += 1
            tot += st.st_size
            c = categorize(p)
            by_cat[c] = by_cat.get(c, 0) + st.st_size
            m = re.search(r"\.q_(\d+)/", p)
            if m:
                per_q[int(m.group(1))] = per_q.get(int(m.group(1)), 0) + st.st_size
            if n >= max_files:
                break
    return {"exists": True, "files": n, "bytes": tot, "GB": round(tot / 1e9, 3),
            "by_category_bytes": dict(sorted(by_cat.items(), key=lambda kv: -kv[1])),
            "per_q_dir_bytes": dict(sorted(per_q.items()))}


def parse_manifests(mdir: Path) -> list:
    out = []
    if not mdir.is_dir():
        return out
    for p in sorted(mdir.glob("*.tsv")):
        txt = read_text(p) or ""
        lines = txt.splitlines()
        if not lines:
            continue
        hdr = lines[0].split("\t")
        try:
            ip, ib = hdr.index("path"), hdr.index("bytes")
        except ValueError:
            continue
        kind = "kept_listing" if "status" in hdr else "deleted"
        ir = hdr.index("reason") if "reason" in hdr else (hdr.index("status") if "status" in hdr else None)
        rows = []
        for line in lines[1:]:
            c = line.split("\t")
            if len(c) <= max(ip, ib):
                continue
            try:
                b = int(c[ib])
            except ValueError:
                continue
            rows.append((c[ip], b, c[ir] if ir is not None and ir < len(c) else ""))
        tot = sum(r[1] for r in rows)
        by_cat, per_q, reasons = {}, {}, {}
        for path, b, why in rows:
            cat = categorize(path)
            by_cat[cat] = by_cat.get(cat, 0) + b
            m = re.search(r"\.q_(\d+)/", path)
            if m:
                per_q[int(m.group(1))] = per_q.get(int(m.group(1)), 0) + b
            reasons[why] = reasons.get(why, 0) + b
        prefix = os.path.commonpath([r[0] for r in rows]) if rows else None
        out.append({"file": tpath(p), "kind": kind, "rows": len(rows), "bytes": tot, "GB": round(tot / 1e9, 3),
                    "common_path": tpath(prefix) if prefix else None,
                    "by_category_bytes": dict(sorted(by_cat.items(), key=lambda kv: -kv[1])),
                    "per_q_dir_bytes": dict(sorted(per_q.items())),
                    "bytes_by_reason": reasons,
                    "_rows": [(r[0], r[1]) for r in rows]})
    return out


def deleted_bytes_under(manifests: list, prefix_tilde: str) -> dict:
    """Bytes of files listed as deleted (pruned) under a ~/ prefix that no longer exist on disk."""
    tot, cats, per_q, files = 0, {}, {}, set()
    for m in manifests:
        if m["kind"] != "deleted":
            continue
        for p, b in m["_rows"]:
            if not p.startswith(prefix_tilde):
                continue
            if os.path.exists(p.replace("~", HOME, 1)):
                continue
            files.add(m["file"])
            tot += b
            cat = categorize(p)
            cats[cat] = cats.get(cat, 0) + b
            mq = re.search(r"\.q_(\d+)/", p)
            if mq:
                per_q[int(mq.group(1))] = per_q.get(int(mq.group(1)), 0) + b
    return {"manifests": sorted(files), "bytes": tot, "GB": round(tot / 1e9, 3),
            "by_category_bytes": dict(sorted(cats.items(), key=lambda kv: -kv[1])),
            "per_q_dir_bytes": dict(sorted(per_q.items()))}


# ----------------------------------------------------------------------------- job inventory / overlap
def scan_qe_jobs(roots: list, now: float,
                 exclude_parts=("/tmp/", "/save/", "/out/", "/_ph0/", "/src/")) -> list:
    """Every QE output under the roots (reference outputs under src/ are not runs and are skipped)."""
    jobs = []
    for root in roots:
        if not root.exists():
            continue
        for dp, dn, fn in os.walk(root, followlinks=False):
            dn[:] = [d for d in dn if d not in ("tmp", "save", "out", "_ph0", "elph_dir", "structures")]
            for f in fn:
                if not re.search(r"\.out(\.part\d+)?$", f):
                    continue
                p = Path(dp) / f
                sp = str(p)
                if any(x in sp + "/" for x in exclude_parts):
                    continue
                txt = read_text(p, limit=4000)
                if not txt or "starts on" not in txt:
                    continue
                q = parse_qe_out(p, now)
                if not q or q.get("start_ts") is None:
                    continue
                end = q["end_ts"] or (now if q["status"] == "running" else mtime(p))
                jobs.append({"id": tpath(p), "program": q.get("program"), "start_ts": q["start_ts"],
                             "end_ts": end, "status": q["status"], "ranks": q.get("mpi_procs") or 1,
                             "threads": q.get("threads_per_mpi") or 1,
                             "program_wall_s": q.get("program_wall_s"), "program_cpu_s": q.get("program_cpu_s")})
    return jobs


def concurrency(target_start, target_end, jobs: list, exclude_ids=(), self_ranks=0) -> dict | None:
    """Time-weighted number of OTHER busy ranks (ranks*threads) while [start,end] ran, plus overlaps."""
    if target_start is None or target_end is None or target_end <= target_start:
        return None
    T = target_end - target_start
    ev = []
    overl = []
    for j in jobs:
        if j["id"] in exclude_ids or j["start_ts"] is None or j["end_ts"] is None:
            continue
        a, b = max(target_start, j["start_ts"]), min(target_end, j["end_ts"])
        if b <= a:
            continue
        w = (j.get("ranks") or 1) * (j.get("threads") or 1)
        ev.append((a, w))
        ev.append((b, -w))
        overl.append({"job": j["id"], "program": j.get("program"), "ranks_x_threads": w,
                      "overlap_s": round(b - a, 1), "overlap_frac": round((b - a) / T, 3)})
    ev.sort()
    cur, last, area, peak = 0, target_start, 0.0, 0
    for t, d in ev:
        area += cur * (t - last)
        cur += d
        last = t
        peak = max(peak, cur)
    area += cur * (target_end - last)
    overl.sort(key=lambda o: -o["overlap_s"] * o["ranks_x_threads"])
    other = area / T
    return {"interval_s": round(T, 1), "other_busy_ranks_time_avg": round(other, 2),
            "other_busy_ranks_peak": peak,
            "total_ranks_time_avg_incl_self": round(other + self_ranks, 2),
            "oversubscription_vs_8_cores": round((other + self_ranks) / 8.0, 2),
            "oversubscription_vs_16_threads": round((other + self_ranks) / 16.0, 2),
            "overlapping_jobs": overl[:25], "n_overlapping_jobs": len(overl)}


def proc_snapshot() -> dict:
    snap = {}
    try:
        snap["loadavg"] = read_text(Path("/proc/loadavg")).split()[:3]
    except Exception:
        pass
    snap["nproc"] = os.cpu_count()
    try:
        mi = read_text(Path("/proc/meminfo")) or ""
        m = re.search(r"MemTotal:\s+(\d+)", mi)
        a = re.search(r"MemAvailable:\s+(\d+)", mi)
        snap["mem_total_GB"] = round(int(m.group(1)) / 1e6, 2) if m else None
        snap["mem_available_GB"] = round(int(a.group(1)) / 1e6, 2) if a else None
    except Exception:
        pass
    procs = {}
    states = {}
    for pid in os.listdir("/proc"):
        if not pid.isdigit():
            continue
        st = read_text(Path(f"/proc/{pid}/status"))
        if not st:
            continue
        name = re.search(r"^Name:\s+(\S+)", st, re.M)
        state = re.search(r"^State:\s+(\S)", st, re.M)
        thr = re.search(r"^Threads:\s+(\d+)", st, re.M)
        nm = name.group(1) if name else "?"
        if nm not in ("pw.x", "ph.x", "epw.x", "q2r.x", "matdyn.x", "alpha2f.x", "wannier90.x", "python",
                      "python3", "mpirun", "prterun", "orted"):
            continue
        d = procs.setdefault(nm, {"count": 0, "threads": 0, "states": {}})
        d["count"] += 1
        d["threads"] += int(thr.group(1)) if thr else 0
        s = state.group(1) if state else "?"
        d["states"][s] = d["states"].get(s, 0) + 1
        states[s] = states.get(s, 0) + 1
    snap["relevant_processes"] = procs
    snap["note"] = ("state R=running, S=sleeping, D=uninterruptible (usually disk I/O wait; D tasks add to "
                    "the load average without using CPU)")
    return snap


# ----------------------------------------------------------------------------- per-material bench
def bench_material(M: str, work: Path, queue: list, now: float, manifests: list) -> dict:
    bdir = work / "bench" / M
    res: dict = {"material": M, "bench_dir": tpath(bdir)}
    if not bdir.is_dir():
        res["status"] = "missing"
        return res
    cfg = load_toml(bdir / "config.toml")
    res["config"] = config_summary(cfg)
    # run logs
    logs = {}
    for p in sorted(bdir.glob("run*.log")):
        txt = read_text(p) or ""
        logs[p.name] = {"bytes": len(txt), "mtime": iso(mtime(p)), **parse_time_block(txt)}
    res["run_logs"] = logs
    # queue
    qev = [e for e in queue if e.get("material") == M]
    res["queue_events"] = [{k: e[k] for k in ("time", "event", "text")} for e in qev]
    runs = sorted((bdir / "run").glob("*/dfpt")) if (bdir / "run").is_dir() else []
    if not runs:
        res["status"] = "not_started"
        return res
    dfpt = runs[0]
    res["run_id"] = dfpt.parent.name
    verdict = load_json(dfpt.parent / "verdict.json")
    if verdict:
        res["verdict"] = {"status": verdict.get("status"),
                          "reasons_head": [str(r)[:240] for r in (verdict.get("reasons") or [])[:2]],
                          "rule_ids": verdict.get("rule_ids")}
    state = load_json(dfpt / "dfpt_state.json") or {}
    steps_state = state.get("steps", {}) if isinstance(state, dict) else {}
    # QE outputs
    qe = {}
    for s in STEPS:
        p = dfpt / f"{s}.out"
        if p.exists():
            qe[s] = parse_qe_out(p, now)
    steps = {}
    for s in STEPS:
        st = steps_state.get(s, {}) or {}
        rec = st.get("record", {}) or {}
        q = qe.get(s) or {}
        d = {"state_status": st.get("status"), "recover": st.get("recover"),
             "elapsed_s (rtsc wall)": rnd(rec.get("elapsed_s"), 2),
             "started": iso(rec.get("started_at")),
             "returncode": rec.get("returncode"), "timed_out": rec.get("timed_out"),
             "qe_status": q.get("status"), "qe_program_wall_s": rnd(q.get("program_wall_s"), 2),
             "qe_program_cpu_s_rank0": rnd(q.get("program_cpu_s"), 2),
             "cpu_wall_ratio_rank0": q.get("cpu_wall_ratio_rank0"),
             "qe_time_resolution_s": q.get("program_time_resolution_s"),
             "mpi_procs": q.get("mpi_procs"), "threads_per_mpi": q.get("threads_per_mpi"), "npool": q.get("npool"),
             "n_k_points": q.get("n_k_points"), "n_ks_states": q.get("n_ks_states"),
             "fft_dims": q.get("fft_dims"), "dense_G_vectors": q.get("dense_G_vectors"),
             "scf_iterations_total": q.get("scf_iterations_total"), "bfgs": q.get("bfgs"),
             "est_ram_per_process": q.get("est_ram_per_process"),
             "qe_start": iso(q.get("start_ts")), "qe_end_or_last_write": iso(q.get("end_or_last_write_ts"))}
        if rec.get("started_at") and rec.get("elapsed_s") is not None:
            d["_interval"] = (rec["started_at"], rec["started_at"] + rec["elapsed_s"])
        elif q.get("start_ts"):
            d["_interval"] = (q["start_ts"], q.get("end_or_last_write_ts") or now)
        if q.get("routines"):
            d["top_routines_by_wall"] = top_routines(q["routines"], n=12)
        steps[s] = {k: v for k, v in d.items() if v is not None}
    res["steps"] = steps
    # ph.x details
    restart_ts = [e["ts"] for e in qev if e["event"] == "RESUME" and e.get("ts")]
    marks = completion_marks(dfpt / "elph_dir", ["elph.inp_lambda.{i}"]) or completion_marks(dfpt, ["dyn{i}"])
    if "ph" in qe:
        res["ph"] = {"breakdown": ph_breakdown(qe["ph"]),
                     "top_routines_by_wall": top_routines(qe["ph"]["routines"], n=25),
                     "q": ph_qblocks(qe["ph"], marks)}
        res["ph"]["per_q_summary"] = per_q_summary(res["ph"]["q"]["per_q"])
    res["q_completion_timeline (all runs incl. pre-restart)"] = marks_timeline(marks, restart_ts)
    dyns = dyn_timeline(dfpt)
    res["dyn_files"] = {"n": len(dyns), "n_nonempty": sum(1 for d in dyns if d["bytes"] > 0)}
    # restarts / crashes
    rs = []
    for e in qev:
        if e["event"] == "RESUME":
            rs.append({"time": e["time"], "source": "queue.log", "text": e["text"]})
    if (steps_state.get("ph") or {}).get("recover"):
        rs.append({"source": "dfpt_state.json", "text": "ph step marked recover=true (ph.x restarted with recover)"})
    if "ph" in res and res["ph"]["q"].get("restarted_output (recover)"):
        rs.append({"source": "ph.out", "text": "ph.out is a recover run: 'Restart in Phonon calculation'; "
                                               "it covers only the resumed part"})
    for name, lg in logs.items():
        if lg.get("bytes", 0) == 0:
            rs.append({"source": name, "text": "empty log (process killed or still running; `time` output absent)"})
    res["restart_events"] = rs
    # errors
    errs = {}
    for s in STEPS:
        t = read_text(dfpt / f"{s}.err")
        if t:
            lines = [l for l in t.splitlines() if l.strip()]
            errs[s] = {"lines": len(lines), "distinct": sorted(set(lines))[:5]}
    res["stderr"] = errs
    # end-to-end
    e2e = {}
    for name in ("run.log",):
        lg = logs.get(name, {})
        if "real_s" in lg:
            e2e = {"source": f"{name} `time`", "real_s": lg["real_s"], "user_s": lg.get("user_s"),
                   "sys_s": lg.get("sys_s"),
                   "avg_busy_cores (user+sys)/real": ratio((lg.get("user_s") or 0) + (lg.get("sys_s") or 0),
                                                         lg["real_s"], 2)}
    starts = [v["_interval"][0] for v in steps.values() if "_interval" in v]
    ends = [v["_interval"][1] for v in steps.values() if "_interval" in v]
    if starts:
        e2e["first_step_start"] = iso(min(starts))
        e2e["last_step_end"] = iso(max(ends))
        e2e["steps_span_s"] = rnd(max(ends) - min(starts), 1)
        done = [v.get("elapsed_s (rtsc wall)") for v in steps.values() if v.get("elapsed_s (rtsc wall)")]
        e2e["sum_step_elapsed_s"] = rnd(sum(done), 1)
        if "real_s" in e2e:
            e2e["rtsc_overhead_s (real - sum steps)"] = rnd(e2e["real_s"] - sum(done), 1)
    qs = [e for e in qev if e["event"] == "START"]
    qe_end = [e for e in qev if e["event"] == "END"]
    if qs:
        e2e["queue_start"] = qs[0]["time"]
        if qe_end:
            e2e["queue_end"] = qe_end[-1]["time"]
            e2e["queue_span_s"] = rnd(qe_end[-1]["ts"] - qs[0]["ts"], 1)
        else:
            e2e["elapsed_so_far_s (now - queue start; includes crash downtime)"] = rnd(now - qs[0]["ts"], 1)
    res["end_to_end"] = e2e
    # status
    phs = (steps_state.get("ph") or {}).get("status")
    if verdict:
        res["status"] = "finished"
    elif phs == "running" or any((qe.get(s) or {}).get("status") == "running" for s in STEPS):
        res["status"] = "running"
    elif steps_state:
        res["status"] = "partial"
    else:
        res["status"] = "not_started"
    # disk
    res["disk"] = {"current_dfpt_dir": walk_sizes(dfpt),
                   "deleted_by_pruning (from manifests)": deleted_bytes_under(manifests, tpath(dfpt) + "/")}
    cur = res["disk"]["current_dfpt_dir"].get("bytes", 0) or 0
    dele = res["disk"]["deleted_by_pruning (from manifests)"]["bytes"]
    res["disk"]["total_written_lower_bound_GB"] = round((cur + dele) / 1e9, 3)
    res["disk"]["note"] = ("lower bound = files present now + files deleted by pruning (manifest bytes). "
                           "Files overwritten in place or deleted by QE itself are not counted.")
    pq_bytes = {}
    for src in (res["disk"]["current_dfpt_dir"].get("per_q_dir_bytes", {}),
                res["disk"]["deleted_by_pruning (from manifests)"]["per_q_dir_bytes"]):
        for k, v in src.items():
            pq_bytes[k] = pq_bytes.get(k, 0) + v
    if pq_bytes:
        vals = sorted(pq_bytes.values())
        res["disk"]["per_q_scratch"] = {"n_q_dirs_seen": len(vals), "GB_min": round(vals[0] / 1e9, 3),
                                        "GB_max": round(vals[-1] / 1e9, 3),
                                        "GB_mean": round(sum(vals) / len(vals) / 1e9, 3)}
    return res


# ----------------------------------------------------------------------------- validation
def validation_runs(work: Path, now: float) -> dict:
    vroot = work / "validate"
    out = {"chain_log": parse_chain_log(vroot / "chain.log"), "runs": []}
    for ev in out["chain_log"]:
        ev.pop("ts", None)
    if not vroot.is_dir():
        return out
    for mdir in sorted(p for p in vroot.iterdir() if p.is_dir()):
        for tdir in sorted(p for p in mdir.iterdir() if p.is_dir()):
            r = {"material": mdir.name, "tag": tdir.name, "dir": tpath(tdir)}
            txt = read_text(tdir / "run.log") or ""
            r.update(parse_time_block(txt))
            r["config"] = (config_summary(load_toml(tdir / "config.toml")) or {}).get("certify")
            cstart = mtime(tdir / "config.toml")
            rl = mtime(tdir / "run.log")
            r["status"] = "done" if "wall_s" in r else ("running" if rl and now - rl < 6 * 3600 else "incomplete")
            r["_interval"] = (cstart, (rl if "wall_s" in r else now))
            recs = list((tdir / "run").glob("*/record.json")) if (tdir / "run").is_dir() else []
            if recs:
                rec = load_json(recs[0]) or {}
                c = ((rec.get("stages") or {}).get("certify") or {})
                cert = ((rec.get("verdict") or {}).get("certificate") or {})
                v = rec.get("verdict") or {}
                r["verdict"] = {"status": v.get("status"), "T_low_K": v.get("T_low_K"), "T_up_K": v.get("T_up_K"),
                                "rule_ids": v.get("rule_ids")}
                r["certify_started"] = c.get("started_utc")
                r["certify_finished"] = c.get("finished_utc")
                r["certificate_seconds"] = rnd(cert.get("seconds"), 2)
                r["certify_params"] = cert.get("params")
                alts = []
                for a in cert.get("alternatives") or []:
                    su, sl = a.get("search_upper") or {}, a.get("search_lower") or {}
                    alts.append({"label": a.get("label"), "seconds": rnd(a.get("seconds"), 2),
                                 "search_upper_s": rnd(su.get("seconds"), 2), "search_upper_attempts": su.get("attempts"),
                                 "search_lower_s": rnd(sl.get("seconds"), 2), "search_lower_attempts": sl.get("attempts"),
                                 "n_cells": len(((a.get("cells") or {}).get("cells")) or []),
                                 "cell_check_s_sum": rnd(sum(((x.get("check") or {}).get("seconds") or 0)
                                                             for x in ((a.get("cells") or {}).get("cells") or [])), 3),
                                 "classify": a.get("classifyAlt")})
                r["alternatives"] = alts
                if alts:
                    r["alternatives_seconds_sum"] = rnd(sum(a["seconds"] or 0 for a in alts), 2)
                    r["search_seconds_sum"] = rnd(sum((a["search_upper_s"] or 0) + (a["search_lower_s"] or 0)
                                                      for a in alts), 2)
                st = load_json(recs[0].parent / "dfpt" / "dfpt_state.json") or {}
                r["dfpt_steps_adopted"] = {k: ("adopted" if (s or {}).get("adopted") else (s or {}).get("status"))
                                           for k, s in (st.get("steps") or {}).items()}
            out["runs"].append(r)
    return out


# ----------------------------------------------------------------------------- tetra and EPW
def tetra_runs(work: Path, now: float) -> dict:
    root = work / "tetra_Al"
    res = {"dir": tpath(root)}
    if not root.is_dir():
        res["status"] = "missing"
        return res
    runs = {}
    for rel in ["q6s_opt/scf.out", "q6s_opt/ph.out.part1", "q6s_opt/ph.out", "q6s_opt/elph_k16/elph.out",
                "q6s_opt/elph_k16/a2F.out", "q6s_opt/elph_k24_crashed_20261027/elph.out",
                "q6s_opt/elph_k24/elph.out", "q6s_lin/scf.out", "q6s_lin/ph.out"]:
        p = root / rel
        if not p.exists():
            continue
        q = parse_qe_out(p, now)
        if not q:
            continue
        d = {"program": q.get("program"), "status": q["status"], "start": iso(q.get("start_ts")),
             "end_or_last_write": iso(q.get("end_or_last_write_ts")),
             "mpi_procs": q.get("mpi_procs"), "npool": q.get("npool"), "n_k_points": q.get("n_k_points"),
             "program_wall_s": rnd(q.get("program_wall_s"), 2), "program_cpu_s_rank0": rnd(q.get("program_cpu_s"), 2),
             "cpu_wall_ratio_rank0": q.get("cpu_wall_ratio_rank0"),
             "time_resolution_s": q.get("program_time_resolution_s")}
        if q.get("program") == "PHONON":
            if "elph" in p.name:
                mdir = p.parent if p.parent.name != "elph_k24" else p.parent.parent   # k24 runs in place
                marks = completion_marks(mdir, ["al.dyn{i}.elph.{i}"])
            else:
                marks = completion_marks(p.parent, ["al.dyn{i}"])
            qb = ph_qblocks(q, marks)
            d["n_irreducible_q"] = qb["n_irreducible_q"]
            d["q_blocks"] = qb["q_blocks_in_this_output"]
            d["q_complete"] = sum(1 for b in qb["per_q"] if b["complete"])
            d["per_q"] = [{k: b.get(k) for k in ("q_index", "wall_s", "nscf_wall_s", "dfpt_elph_wall_s",
                                                  "window_cpu_wall_ratio_rank0", "nscf_k_points",
                                                  "dfpt_scf_iterations", "n_irreps", "complete")}
                          for b in qb["per_q"]]
            d["per_q_summary"] = per_q_summary(qb["per_q"])
            if q.get("routines"):
                d["breakdown"] = ph_breakdown(q)
                d["top_routines_by_wall"] = top_routines(q["routines"], n=15)
        runs[rel] = d
    res["runs"] = runs
    tim = {}
    for rel in ["q6s_opt/timing.txt", "q6s_opt/elph_k16/timing.txt"]:
        t = read_text(root / rel)
        if t:
            tim[rel] = [l.strip() for l in t.splitlines() if l.strip()]
    res["timing_txt"] = tim
    res["disk"] = {"q6s_opt_current": walk_sizes(root / "q6s_opt")}
    return res


def epw_runs(work: Path, now: float) -> dict:
    root = work / "epw_Al"
    res = {"dir": tpath(root)}
    if not root.is_dir():
        res["status"] = "missing"
        return res
    jobs = []
    for p in sorted(root.rglob("*.out")):
        sp = str(p)
        if any(x in sp for x in ("/tmp/", "/save/", "/out/")):
            continue
        head = read_text(p, limit=3000) or ""
        if "starts on" not in head:
            continue
        q = parse_qe_out(p, now)
        if not q:
            continue
        jobs.append({"path": tpath(p.relative_to(root)), "program": q.get("program"), "status": q["status"],
                     "start": iso(q.get("start_ts")), "end_or_last_write": iso(q.get("end_or_last_write_ts")),
                     "mpi_procs": q.get("mpi_procs"), "npool": q.get("npool"),
                     "program_wall_s": rnd(q.get("program_wall_s"), 2),
                     "program_cpu_s_rank0": rnd(q.get("program_cpu_s"), 2),
                     "cpu_wall_ratio_rank0": q.get("cpu_wall_ratio_rank0"),
                     "top_routines_by_wall": top_routines(q.get("routines") or {}, n=6)})
    res["qe_jobs"] = jobs
    logs = {}
    for p in sorted(list(root.glob("*.log")) + list(root.glob("*.out"))):
        t = read_text(p) or ""
        walls = [float(x) for x in re.findall(r"^wall\s+([\d.]+)\s*s", t, re.M)]
        rcs = re.findall(r"rc=(\d+), (\d+)s\)", t)
        if walls or rcs or p.suffix == ".log":
            logs[p.name] = {"bytes": len(t), "mtime": iso(mtime(p)), "wall_s": walls,
                            "failed_commands (rc, s)": [[int(a), int(b)] for a, b in rcs],
                            "has_traceback": "Traceback" in t}
    res["driver_logs"] = logs
    res["crash_files"] = [tpath(p.relative_to(root)) for p in root.rglob("CRASH")
                          if "/tmp/" not in str(p)]
    res["disk"] = {"epw_Al_current": walk_sizes(root)}
    return res


# ----------------------------------------------------------------------------- main
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--work", default=os.path.join(HOME, "rtsc_work"))
    ap.add_argument("--repo", default="/mnt/g/ClaudeWorkspace/RTSC_Project")
    ap.add_argument("--out", default=None)
    ap.add_argument("--no-proc", action="store_true", help="skip the /proc snapshot")
    a = ap.parse_args(argv)
    work, repo = Path(os.path.expanduser(a.work)), Path(a.repo)
    out = Path(a.out) if a.out else repo / "bench_tools" / "perf" / "perf_metrics.json"
    now = time.time()
    t0 = time.time()

    queue = parse_queue_log(work / "bench" / "queue.log")
    manifests = parse_manifests(repo / "reports" / "data" / "pruned_manifests")

    result = {"schema": SCHEMA, "generated": iso(now), "work_root": tpath(work),
              "conventions": {
                  "qe_cpu": "QE CPU = rank-0 CPU seconds; cpu_wall_ratio_rank0 = CPU/WALL of rank 0",
                  "time_user": "`time` user/sys in run.log = summed over all processes",
                  "per_q_wall": "from cumulative PHONON lines (resolution 0.01 s < 1 h, 60 s above) and dyn<i> mtimes",
                  "disk": "bytes = present now + deleted by pruning (manifests); a lower bound of bytes written",
                  "times": "ISO 8601 local time of the run machine"},
              "machine": {"logical_cpus": os.cpu_count(), "notes": "8 physical cores / 16 logical (SMT); WSL2"}}
    if not a.no_proc:
        try:
            result["machine"]["snapshot_now"] = proc_snapshot()
        except Exception as e:     # never fail on the snapshot
            result["machine"]["snapshot_now"] = {"error": repr(e)}

    result["queue_log"] = [{k: e[k] for k in ("time", "event", "text", "material")} for e in queue]
    bench = {M: bench_material(M, work, queue, now, manifests) for M in MATERIALS}
    validation = validation_runs(work, now)
    tetra = tetra_runs(work, now)
    epw = epw_runs(work, now)

    # job inventory and concurrency
    jobs = scan_qe_jobs([work / "bench", work / "tetra_Al", work / "epw_Al"], now)
    for r in validation["runs"]:
        iv = r.get("_interval")
        if iv and iv[0]:
            jobs.append({"id": f"validate/{r['material']}/{r['tag']} (rtsc certify, serial python)",
                         "program": "rtsc-certify", "start_ts": iv[0], "end_ts": iv[1],
                         "status": r["status"], "ranks": 1, "threads": 1})
    for M, b in bench.items():
        for s, d in (b.get("steps") or {}).items():
            iv = d.pop("_interval", None)
            if iv:
                self_id = tpath(work / "bench" / M / "run" / b.get("run_id", "?") / "dfpt" / f"{s}.out")
                d["concurrency"] = concurrency(iv[0], iv[1], jobs, exclude_ids={self_id},
                                               self_ranks=(d.get("mpi_procs") or 1) * (d.get("threads_per_mpi") or 1))
    # contention per q (time-averaged number of other busy ranks while each q ran)
    for M, b in bench.items():
        own = tpath(work / "bench" / M / "run" / b.get("run_id", "?") / "dfpt" / "ph.out")
        tl = b.get("q_completion_timeline (all runs incl. pre-restart)") or []
        prev = None
        for ent in tl:
            t = dt.datetime.fromisoformat(ent["finished"]).timestamp()
            if prev is not None and not ent.get("spans_restart"):
                c = concurrency(prev, t, jobs, exclude_ids={own}, self_ranks=8)
                if c:
                    ent["other_busy_ranks_time_avg"] = c["other_busy_ranks_time_avg"]
                    ent["overlapping"] = sorted({o["job"].split("rtsc_work/")[-1].split("/")[0] + ":" +
                                                 (o["program"] or "") for o in c["overlapping_jobs"]})
            prev = t
    for r in validation["runs"]:
        iv = r.pop("_interval", None)
        if iv and iv[0]:
            r["concurrency"] = concurrency(iv[0], iv[1], jobs,
                                           exclude_ids={f"validate/{r['material']}/{r['tag']} (rtsc certify, serial python)"},
                                           self_ranks=1)
            r["interval"] = [iso(iv[0]), iso(iv[1])]
    result["bench"] = bench
    result["validation"] = validation
    result["tetra_Al"] = tetra
    result["epw_Al"] = epw
    result["job_inventory"] = sorted(
        [{"id": j["id"], "program": j["program"], "start": iso(j["start_ts"]), "end_or_now": iso(j["end_ts"]),
          "wall_s": rnd((j["end_ts"] or now) - j["start_ts"], 1), "status": j["status"],
          "ranks": j["ranks"], "threads": j["threads"]} for j in jobs], key=lambda d: d["start"] or "")
    result["pruned_manifests"] = [{k: v for k, v in m.items() if k != "_rows"} for m in manifests]

    # compact cross-material summary
    summ = {}
    for M, b in bench.items():
        s = {"status": b.get("status")}
        e = b.get("end_to_end") or {}
        s["end_to_end_s"] = e.get("real_s") or e.get("queue_span_s") or e.get(
            "elapsed_so_far_s (now - queue start; includes crash downtime)")
        s["steps_wall_s"] = {k: v.get("elapsed_s (rtsc wall)") for k, v in (b.get("steps") or {}).items()
                             if v.get("elapsed_s (rtsc wall)") is not None}
        run_steps = {k: {"qe_status": v.get("qe_status"), "wall_s_of_current_qe_output": v.get("qe_program_wall_s")}
                     for k, v in (b.get("steps") or {}).items()
                     if v.get("elapsed_s (rtsc wall)") is None and v.get("qe_status")}
        if run_steps:
            s["steps_not_finished"] = run_steps
        tl = b.get("q_completion_timeline (all runs incl. pre-restart)") or []
        if tl:
            s["q_finished"] = len(tl)
        s["steps_cpu_wall_ratio_rank0"] = {k: v.get("cpu_wall_ratio_rank0") for k, v in (b.get("steps") or {}).items()}
        ph = b.get("ph") or {}
        if ph.get("top_routines_by_wall"):
            s["ph_top5_routines"] = [(t["routine"], t["wall_s"]) for t in ph["top_routines_by_wall"][:5]]
        if ph.get("q"):
            s["n_irreducible_q"] = ph["q"].get("n_irreducible_q")
        cfg = b.get("config") or {}
        s["layout"] = cfg.get("launcher")
        s["grids"] = {k: (cfg.get("settings") or {}).get(k) for k in ("k_scf", "k_fit", "q", "k_matdyn")}
        s["disk_written_lower_bound_GB"] = (b.get("disk") or {}).get("total_written_lower_bound_GB")
        s["restart_events"] = len(b.get("restart_events") or [])
        summ[M] = s
    result["summary"] = summ
    result["extractor_runtime_s"] = round(time.time() - t0, 2)

    result = sanitize(result)
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(".json.tmp")
    with open(tmp, "w") as f:
        json.dump(result, f, indent=1, default=str)
    os.replace(tmp, out)
    print(f"wrote {tpath(out)}  ({os.path.getsize(out)} bytes, {result['extractor_runtime_s']} s)")
    for M, s in summ.items():
        print(f"  {M:5s} {s['status']:12s} e2e={s['end_to_end_s']}  steps={s['steps_wall_s']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
