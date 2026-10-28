"""Executors for the DFT stage of a discovery round (DISCOVERY.md section 4).

``run_pipeline(cfg, candidates, on_submit=None, existing_job=None)`` runs the EXISTING pipeline
(``rtsc.pipeline.run`` / ``collect``) on one goal's pipeline config and returns
{material_id: record dict or None}.  None means "no finished record" (job failed, killed or still
running); the loop records that as a retryable ERROR and retries it on ``continue``.

* ``local``: ``pipeline.run`` in this process (process pool of ``[discovery] workers``).
* ``slurm``: ``pipeline.run(plan_only=True)``, then ONE sbatch job array whose task K runs
  ``python -m rtsc.discovery.executors shard <job.json> K N`` (= ``pipeline.run(shard="K/N")``),
  submitted with ``--hold`` (its job id is recorded in the ledger) and then ``scontrol release``d,
  polling ``squeue`` until the array has left the queue and reading final states from ``sacct``;
  then ``pipeline.collect``.  On ``continue`` a recorded job id that is still queued is waited for
  instead of being resubmitted.  Keys of ``[discovery.slurm]``: partition, account, time, nodes,
  ntasks, cpus_per_task, gpus, extra (list of raw ``#SBATCH`` options), setup (list of shell lines,
  e.g. ``module load qe``), max_array (64), poll_seconds (30), timeout_hours (none), python,
  sbatch / squeue / sacct / scontrol (command names), gpu_partition (used for generation by adapters that
  submit their own jobs; recorded only).
"""
from __future__ import annotations

import json
import shlex
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Callable

from .. import pipeline
from ..sources.base import Candidate


class LocalExecutor:
    name = "local"

    def run_pipeline(self, cfg: dict[str, Any], candidates: list[Candidate], *,
                     on_submit: Callable[[str], None] | None = None,
                     existing_job: str | None = None) -> dict[str, dict | None]:
        pipeline.run(cfg, candidates)
        wd = cfg["run"]["workdir"]
        return {c.material_id: pipeline.read_record(wd, c.material_id) for c in candidates}


class SlurmError(RuntimeError):
    pass


class SlurmExecutor:
    name = "slurm"

    def __init__(self, scfg: dict[str, Any] | None = None):
        s = dict(scfg or {})
        self.cfg = s
        self.max_array = int(s.get("max_array", 64))
        self.poll = float(s.get("poll_seconds", 30))
        self.timeout_h = s.get("timeout_hours")
        self.python = s.get("python", sys.executable)
        self.sbatch = s.get("sbatch", "sbatch")
        self.squeue = s.get("squeue", "squeue")
        self.sacct = s.get("sacct", "sacct")
        self.scontrol = s.get("scontrol", "scontrol")
        self.last_states: dict[str, str] = {}

    # -------------------------------------------------------------- slurm plumbing
    def _run(self, cmd: list[str]) -> str:
        p = subprocess.run(cmd, capture_output=True, text=True)
        if p.returncode != 0:
            raise SlurmError(f"{' '.join(cmd)} failed ({p.returncode}): {p.stderr.strip()}")
        return p.stdout

    def script(self, job_file: Path, n: int, wd: Path) -> str:
        s = self.cfg
        lines = ["#!/bin/bash", f"#SBATCH --job-name=rtsc-{wd.parent.name}-{wd.name}",
                 f"#SBATCH --output={wd}/slurm-%A_%a.out"]
        for key, opt in (("partition", "partition"), ("account", "account"), ("time", "time"),
                         ("nodes", "nodes"), ("ntasks", "ntasks"), ("cpus_per_task", "cpus-per-task"),
                         ("gpus", "gpus")):
            if s.get(key) is not None:
                lines.append(f"#SBATCH --{opt}={s[key]}")
        lines += [f"#SBATCH {x}" for x in s.get("extra", [])]
        lines += list(s.get("setup", []))
        lines.append(f"exec {shlex.quote(self.python)} -m rtsc.discovery.executors shard "
                     f"{shlex.quote(str(job_file))} \"$SLURM_ARRAY_TASK_ID\" {n}")
        return "\n".join(lines) + "\n"

    def submit(self, script_path: Path, n: int) -> str:
        out = self._run([self.sbatch, "--parsable", "--hold", f"--array=0-{n - 1}",
                         str(script_path)]).strip()
        tok = out.split()[-1] if out else ""
        jid = tok.split(";")[0]
        if not jid.isdigit():
            raise SlurmError(f"cannot parse sbatch output {out!r}")
        return jid

    def alive(self, jid: str) -> bool:
        try:
            out = self._run([self.squeue, "-h", "-j", str(jid), "-o", "%i %T"])
        except SlurmError:
            return False          # squeue errors for unknown (finished) job ids
        return bool(out.strip())

    def wait(self, jid: str) -> dict[str, str]:
        t0 = time.monotonic()
        while self.alive(jid):
            if self.timeout_h is not None and time.monotonic() - t0 > float(self.timeout_h) * 3600:
                raise SlurmError(f"job {jid} still queued after {self.timeout_h} h")
            time.sleep(self.poll)
        states: dict[str, str] = {}
        try:
            for line in self._run([self.sacct, "-j", str(jid), "-n", "-P", "-o", "JobID,State"]).splitlines():
                parts = line.strip().split("|")
                if len(parts) >= 2 and parts[0]:
                    states[parts[0]] = parts[1]
        except SlurmError:
            pass
        self.last_states = states
        return states

    def failure_note(self) -> str:
        bad = {k: v for k, v in self.last_states.items() if not str(v).startswith("COMPLETED")}
        return f"sacct states of the last array: {bad}" if bad else ""

    # -------------------------------------------------------------- executor API
    def run_pipeline(self, cfg: dict[str, Any], candidates: list[Candidate], *,
                     on_submit: Callable[[str], None] | None = None,
                     existing_job: str | None = None) -> dict[str, dict | None]:
        wd = Path(cfg["run"]["workdir"])
        wd.mkdir(parents=True, exist_ok=True)
        if existing_job and self.alive(existing_job):
            self.wait(existing_job)
        pipeline.run(cfg, candidates, plan_only=True, write_reports=False)
        todo = [c.material_id for c in candidates
                if (pipeline.read_record(wd, c.material_id) or {}).get("verdict") is None]
        if todo:
            job_file = wd / "discovery_job.json"
            job_file.write_text(json.dumps({"cfg": cfg, "candidates": [c.to_dict() for c in candidates]},
                                           default=str))
            n = max(1, min(self.max_array, len(candidates)))
            sp = wd / "array.sbatch"
            sp.write_text(self.script(job_file, n, wd))
            jid = self.submit(sp, n)          # submitted HELD: the id is in the ledger before it runs
            if on_submit:
                on_submit(jid)
            self._run([self.scontrol, "release", jid])
            self.wait(jid)
        pipeline.collect(cfg)
        return {c.material_id: _finished(pipeline.read_record(wd, c.material_id)) for c in candidates}


def _finished(rec: dict | None) -> dict | None:
    return rec if rec and rec.get("verdict") else None


def make_executor(dcfg: dict[str, Any]):
    if dcfg.get("executor", "local") == "slurm":
        return SlurmExecutor(dcfg.get("slurm"))
    return LocalExecutor()


def _shard_main(job_file: str, k: str, n: str) -> int:
    job = json.loads(Path(job_file).read_text())
    cands = [Candidate.from_dict(d) for d in job["candidates"]]
    pipeline.run(job["cfg"], cands, shard=f"{int(k)}/{int(n)}", write_reports=False)
    return 0


if __name__ == "__main__":   # pragma: no cover - exercised through a fake sbatch in tests
    if len(sys.argv) == 5 and sys.argv[1] == "shard":
        sys.exit(_shard_main(*sys.argv[2:]))
    print("usage: python -m rtsc.discovery.executors shard JOB.json K N", file=sys.stderr)
    sys.exit(2)
