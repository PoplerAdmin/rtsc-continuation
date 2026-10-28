"""Discovery loop: fingerprints, policy determinism, config, slurm executor (fake sbatch), and
the end-to-end ``rtsc discover`` CLI on the SYNTHETIC demo config (DISCOVERY.md)."""
from __future__ import annotations

import csv
import json
import os
import stat
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from discovery_mocks import EXAMPLE_DIR, ROOT, make_dc  # noqa: E402

from rtsc.cli import main  # noqa: E402
from rtsc.discovery.fingerprint import fingerprint  # noqa: E402
from rtsc.discovery.interfaces import (GeneratedStructure, GenerationRequest, Generator,  # noqa: E402
                                       Prescreener, PrescreenResult)
from rtsc.discovery.generators import (DirectoryGenerator, NullPrescreener,  # noqa: E402
                                       PrototypeSubstitutionGenerator)
from rtsc.discovery.loop import DiscoveryLoop  # noqa: E402
from rtsc.discovery.policy import ArmStats, BanditPolicy, build_arms  # noqa: E402

NaCl = {"lattice": [[4.0, 0, 0], [0, 4.0, 0], [0, 0, 4.0]], "species": ["La", "H"],
        "frac_coords": [[0, 0, 0], [0.5, 0.5, 0.5]]}


# ------------------------------------------------------------------------------- fingerprints
def test_fingerprint_invariances():
    fp = fingerprint(NaCl, 0.0)
    shifted = {**NaCl, "frac_coords": [[0.5, 0.5, 0.5], [0.0, 0.0, 0.0]], "species": ["H", "La"]}
    assert fingerprint(shifted, 0.0) == fp                     # permutation + origin shift
    sc = {"lattice": [[8.0, 0, 0], [0, 4.0, 0], [0, 0, 4.0]], "species": ["La", "H", "La", "H"],
          "frac_coords": [[0, 0, 0], [0.25, 0.5, 0.5], [0.5, 0, 0], [0.75, 0.5, 0.5]]}
    assert fingerprint(sc, 0.0) == fp                          # 2x1x1 supercell
    scaled = {**NaCl, "lattice": [[4.4, 0, 0], [0, 4.4, 0], [0, 0, 4.4]]}
    assert fingerprint(scaled, 0.0) == fp                      # isotropic scale
    assert fingerprint(NaCl, 100.0) != fp                      # other pressure
    assert fingerprint({**NaCl, "species": ["Y", "H"]}, 0.0) != fp
    distorted = {**NaCl, "frac_coords": [[0, 0, 0], [0.5, 0.5, 0.3]]}
    assert fingerprint(distorted, 0.0) != fp


def test_interfaces_and_builtin_generators(tmp_path):
    g = PrototypeSubstitutionGenerator(str(EXAMPLE_DIR / "prototypes"))
    assert isinstance(g, Generator) and isinstance(NullPrescreener(), Prescreener)
    req = GenerationRequest(("H", "La"), None, 100.0, "G-RT", 50, 1, {})
    outs = g.generate(req)
    assert outs and all(set(o.structure["species"]) == {"H", "La"} for o in outs)
    assert all(o.request_id == req.request_id for o in outs)
    # S-H prototype keeps H on H sites (keep_common)
    h3la = [o for o in outs if o.provenance["prototype"] == "H3S_Im-3m"]
    assert h3la and h3la[0].provenance["mapping"] == {"H": "H", "S": "La"}
    res = NullPrescreener().prescreen(outs[:2], 100.0)
    assert len(res) == 2 and all(isinstance(r, PrescreenResult) and r.ok for r in res)
    d = DirectoryGenerator(str(EXAMPLE_DIR / "prototypes"))
    got = d.generate(GenerationRequest(("H", "S"), None, 0.0, "G-RT", 5, 0, {}))
    assert len(got) == 1 and got[0].provenance["source_id"] == "H3S_Im-3m"
    from rtsc.discovery.generators import make_generator, make_prescreener
    rg = make_generator("rg", {"type": "python", "class": "discovery_mocks:RandomGenerator", "seed": 1})
    assert rg.name == "randgen" and rg.generate(req)
    assert make_prescreener({"type": "null"}).name == "null"
    try:
        import rtsc.generate.mattergen  # noqa: F401
        assert make_generator("mg", {"type": "mattergen"}).name
    except ImportError:
        pass


# ------------------------------------------------------------------------------- policy
def test_policy_deterministic_with_exploration_floor(tmp_path):
    dc = make_dc(tmp_path / "w")
    arms = build_arms(dc.enabled_goals, ["prototypes"])
    assert len(arms) == 3 * 2 + 3 * 1        # G-RT: 3 chemsys x 2 pressures, G-AMB: 3 x 1
    pol = BanditPolicy(arms, seed=5, exploration_floor=0.5)
    stats = {a.key: ArmStats(pulls=3, S=2.0, F=1.0) for a in arms[:-2]}
    a1 = [a.key for a in pol.propose(4, 6, stats)]
    a2 = [a.key for a in BanditPolicy(arms, seed=5, exploration_floor=0.5).propose(4, 6, stats)]
    assert a1 == a2
    assert set(a1[:2]) <= {a.key for a in arms[-2:]}        # floor: unexplored arms first
    ucb = BanditPolicy(arms, seed=5, kind="ucb", exploration_floor=0.0)
    full = {a.key: ArmStats(pulls=10, S=1.0, F=9.0) for a in arms}
    full[arms[3].key] = ArmStats(pulls=10, S=9.0, F=1.0)
    assert ucb.propose(1, 1, full)[0] == arms[3]
    th = BanditPolicy(arms, seed=5, exploration_floor=0.0)
    picks = [th.propose(r, 1, full)[0] for r in range(40)]
    assert picks.count(arms[3]) > 20                         # Thompson favours the good arm


def test_same_seed_same_ledger(tmp_path):
    snaps = []
    for name in ("a", "b"):
        loop = DiscoveryLoop(make_dc(tmp_path / name, stop={"max_rounds": 2, "certified_per_goal": 100}))
        loop.run()
        snaps.append([(m["material_id"], m["fingerprint"], m["arm"]) for m in loop.ledger.materials()]
                     + [json.dumps(r["requests"]) for r in loop.ledger.rounds()])
    assert snaps[0] == snaps[1]


def test_stop_file_and_dft_budget(tmp_path):
    dc = make_dc(tmp_path / "w", stop={"max_rounds": 5, "certified_per_goal": 100, "max_dft_jobs": 4})
    loop = DiscoveryLoop(dc)
    res = loop.run()
    assert "max_dft_jobs" in res["stop_reason"] and len(loop.ledger.jobs()) == 4
    (dc.workdir / "STOP").write_text("")
    res = DiscoveryLoop(dc).run("continue")
    assert "STOP" in res["stop_reason"] and res["rounds_this_session"] == 0


# ------------------------------------------------------------------------------- slurm executor
FAKE_SBATCH = """#!/usr/bin/env python3
import os, re, subprocess, sys
args = sys.argv[1:]
arr = [a for a in args if a.startswith("--array=")][0]
lo, hi = map(int, re.match(r"--array=(\\d+)-(\\d+)", arr).groups())
script = args[-1]
log = os.environ["FAKE_SLURM_LOG"]
with open(log, "a") as fh:
    fh.write("sbatch " + " ".join(args) + "\\n")
for k in range(lo, hi + 1):
    env = dict(os.environ, SLURM_ARRAY_TASK_ID=str(k))
    subprocess.run(["bash", script], env=env, check=True)
print("4242")
"""
FAKE_SQUEUE = """#!/usr/bin/env python3
import os
p = os.environ["FAKE_SLURM_LOG"] + ".squeue"
n = int(open(p).read()) if os.path.exists(p) else 0
open(p, "w").write(str(n + 1))
if n % 3 < 2:
    print("4242_0 RUNNING")
"""
FAKE_SCONTROL = """#!/usr/bin/env python3
import os, sys
with open(os.environ["FAKE_SLURM_LOG"], "a") as fh:
    fh.write("scontrol " + " ".join(sys.argv[1:]) + "\\n")
"""
FAKE_SACCT = """#!/usr/bin/env python3
print("4242_0|COMPLETED")
print("4242_1|COMPLETED")
"""


def test_slurm_executor_with_fake_sbatch(tmp_path, monkeypatch):
    bin_ = tmp_path / "bin"
    bin_.mkdir()
    for name, body in (("sbatch", FAKE_SBATCH), ("squeue", FAKE_SQUEUE), ("sacct", FAKE_SACCT),
                       ("scontrol", FAKE_SCONTROL)):
        p = bin_ / name
        p.write_text(body)
        p.chmod(p.stat().st_mode | stat.S_IEXEC)
    monkeypatch.setenv("PATH", f"{bin_}{os.pathsep}{os.environ['PATH']}")
    monkeypatch.setenv("PYTHONPATH", str(ROOT))
    monkeypatch.setenv("FAKE_SLURM_LOG", str(tmp_path / "slurm.log"))
    dc = make_dc(tmp_path / "w", executor="slurm",
                 slurm={"poll_seconds": 0.01, "max_array": 2, "partition": "cpu", "setup": ["true"]},
                 stop={"max_rounds": 1, "certified_per_goal": 100})
    loop = DiscoveryLoop(dc)
    loop.run()
    log = (tmp_path / "slurm.log").read_text()
    assert "--array=0-1" in log and "--hold" in log and "scontrol release 4242" in log
    jobs = loop.ledger.jobs()
    assert jobs and all(j["external_id"] == "4242" and j["state"] == "DONE" for j in jobs)
    rows = [r for r in loop.ledger.current("G-RT") if r["stage"] == "pipeline"]
    assert len(rows) == len(jobs) and all(r["status"] != "ERROR" for r in rows)
    script = next((dc.workdir / "rounds").rglob("array.sbatch")).read_text()
    assert "#SBATCH --partition=cpu" in script and "rtsc.discovery.executors shard" in script
    # results equal a local run of the same config
    ref = DiscoveryLoop(make_dc(tmp_path / "ref", stop={"max_rounds": 1, "certified_per_goal": 100}))
    ref.run()
    key = lambda lp: sorted((r["material_id"], r["goal_id"], r["status"], r["T_low_K"], r["T_up_K"])
                            for g in ("G-RT", "G-AMB") for r in lp.ledger.current(g))
    assert key(loop) == key(ref)


# ------------------------------------------------------------------------------- (5) end to end
def test_cli_discover_run_continue_status_list_export(tmp_path, capsys):
    wd = tmp_path / "demo"
    cfg = str(EXAMPLE_DIR / "config.toml")
    assert main(["discover", "run", "-c", cfg, "-w", str(wd)]) == 0
    out = capsys.readouterr().out
    assert "certified_per_goal=1 reached" in out
    rep = (wd / "discovery_report.md").read_text()
    assert "SYNTHETIC DATA" in rep and "Model assumptions" in rep
    for g in ("G-RT", "G-AMB"):
        rows = list(csv.DictReader(open(wd / f"candidates_{g}.csv")))
        assert any(r["status"] == "CERTIFIED_CANDIDATE" for r in rows)
        assert all(r["material_id"].startswith("SYNTH-") for r in rows)
        cert = [r for r in rows if r["status"] == "CERTIFIED_CANDIDATE"]
        assert all("SYNTHETIC_DATA" in r["flags"] and "HYPOTHETICAL_STRUCTURE" in r["flags"] for r in cert)
        assert all(Path(r["certificate_path"]).is_file() for r in cert)
    assert main(["discover", "run", "-c", cfg, "-w", str(wd)]) == 2      # refuses to overwrite
    capsys.readouterr()
    from rtsc.discovery.ledger import Ledger
    led = Ledger(wd / "ledger.sqlite")
    before = {g: {r["material_id"] for r in led.current(g) if r["status"] == "CERTIFIED_CANDIDATE"}
              for g in ("G-RT", "G-AMB")}
    led.close()
    assert main(["discover", "continue", "-w", str(wd)]) == 0
    out = capsys.readouterr().out
    assert "certified_per_goal=1 reached" in out
    led = Ledger(wd / "ledger.sqlite")
    for g in ("G-RT", "G-AMB"):
        after = {r["material_id"] for r in led.current(g) if r["status"] == "CERTIFIED_CANDIDATE"}
        assert len(after - before[g]) >= 1, g                 # continue found more
        assert before[g] <= after                             # nothing lost
    led.close()
    assert main(["discover", "status", "-w", str(wd)]) == 0
    assert "G-AMB" in capsys.readouterr().out
    assert main(["discover", "list", "-w", str(wd), "--goal", "G-RT",
                 "--status", "CERTIFIED_CANDIDATE"]) == 0
    assert "CERTIFIED_CANDIDATE" in capsys.readouterr().out
    assert main(["discover", "export", "-w", str(wd)]) == 0
    data = json.loads((wd / "ledger_export.json").read_text())
    assert data["accounting"]["materials"] == len(data["materials"])
    assert main(["discover", "export", "-w", str(wd), "--finetune"]) == 0
    assert " 0 labelled" in capsys.readouterr().out                       # synthetic excluded
    assert main(["discover", "export", "-w", str(wd), "--finetune", "--include-synthetic",
                 "--out", str(tmp_path / "ft")]) == 0
    ft = tmp_path / "ft"
    tr = list(ft.rglob("train.csv"))
    assert tr
    rows = list(csv.DictReader(open(tr[0]))) + list(csv.DictReader(open(tr[0].parent / "val.csv")))
    assert rows and all(r["cif"].startswith("data_") for r in rows)
    # package fallback exporter
    from rtsc.discovery.report import export_finetune
    led = Ledger(wd / "ledger.sqlite")
    from rtsc.discovery.goals import load_discovery_config
    res = export_finetune(led, load_discovery_config(cfg).enabled_goals, tmp_path / "ft2",
                          include_synthetic=True)
    rows = list(csv.DictReader(open(tmp_path / "ft2" / "train.csv")))
    assert res["n"] > 0 and all(float(r["A_K_max"]) > 0 for r in rows)
