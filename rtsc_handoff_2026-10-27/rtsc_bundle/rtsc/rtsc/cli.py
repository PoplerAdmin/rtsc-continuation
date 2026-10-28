"""Command line interface.

    rtsc run -c config.toml [--shard K/N] [--plan-only] [--workers N] [--no-resume]
    rtsc collect -c config.toml                  # merge records of a sharded run, write reports
    rtsc certify a2F.dat --units meV --mu-hi 0.15 --target 300 [--envelope 0.1] [--out v.json]
    rtsc demo [--workdir DIR]                    # SYNTHETIC spectra, no DFT
    rtsc recheck verdict.json                    # re-verify a certificate (rtsc.certify)
    rtsc discover run|continue|status|list|export -c discovery.toml   # discovery loop

Exit codes: 0 success, 1 failure (exception), 2 usage error, 3 recheck mismatch.
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Any

from . import pipeline


def _print_summary(summ) -> None:
    c = summ.counts()
    width = max(len(k) for k in c)
    print(f"inputs: {len(summ.inputs)}  (config {summ.config_hash})")
    for k, v in c.items():
        print(f"  {k:<{width}}  {v}")
    print(f"  {'total':<{width}}  {sum(c.values())}")


def cmd_run(a: argparse.Namespace) -> int:
    cfg = pipeline.load_config(a.config)
    if a.workers is not None:
        cfg["run"]["workers"] = a.workers
    if a.no_resume:
        cfg["run"]["resume"] = False
    if a.workdir:
        cfg["run"]["workdir"] = str(Path(a.workdir).resolve())
    summ = pipeline.run(cfg, shard=a.shard, plan_only=a.plan_only)
    _print_summary(summ)
    if not a.shard and not a.plan_only:
        out = cfg["report"].get("outdir") or cfg["run"]["workdir"]
        print(f"reports: {out}/candidates.csv, summary.md, certificates/")
    return 0


def cmd_collect(a: argparse.Namespace) -> int:
    cfg = pipeline.load_config(a.config)
    summ = pipeline.collect(cfg)
    _print_summary(summ)
    return 0


def _certify_cfg(a: argparse.Namespace) -> dict[str, Any]:
    c: dict[str, Any] = {"T_target_K": a.target, "mu_hi": a.mu_hi, "cutoff_rule": "block",
                         "mode": a.mode, "rejection_model": a.rejection_model,
                         "allow_unverified_rejection": a.allow_unverified_rejection}
    if a.omega_c is not None:
        c["omega_c_K"] = a.omega_c
    else:
        c["omega_c_factor"] = a.omega_c_factor
    return c


def cmd_certify(a: argparse.Namespace) -> int:
    spec = {"path": a.file, "format": a.format, "units": a.units, "label": Path(a.file).name,
            "columns": [int(x) for x in a.columns.split(",")], "rule": a.rule,
            "negative_policy": a.negative_policy, "expected_lambda": a.expected_lambda,
            "expected_A_K": a.expected_A}
    if a.format == "two_column" and not a.units:
        print("error: --units is required for a two-column file", file=sys.stderr)
        return 2
    # same validation as a pipeline run (audit F4): bad T*, mu*, omega_c -> usage error
    cfg = pipeline.normalize_config({"certify": _certify_cfg(a),
                                     "spectra": {"rel_envelope": a.envelope,
                                                 "two_column_rule": a.rule,
                                                 "negative_policy": a.negative_policy}})
    spectra = pipeline._read_spectrum_file(spec, cfg["spectra"])
    spectra = pipeline._with_envelope(spectra, float(cfg["spectra"]["rel_envelope"]))
    mid = a.material_id or Path(a.file).stem
    bcs = []
    if a.band_coupling:
        from .dfpt import BandCoupling
        bcs = [BandCoupling.from_file(x) for x in a.band_coupling]
    v = pipeline.certify_spectra(mid, spectra, cfg["certify"], band_couplings=bcs)
    for s in spectra:
        for w in s.provenance.get("warnings", []):
            print(f"input warning: {w}", file=sys.stderr)
    d = v.to_dict()
    text = json.dumps(d, indent=1, default=str)
    if a.out:
        Path(a.out).write_text(text)
        print(f"verdict written to {a.out}")
    print(f"{mid}: {v.status.value}  T_low={v.T_low_K}  T_up={v.T_up_K}  rules={v.rule_ids}")
    for r in v.reasons:
        print(f"  - {r}")
    print("  assumptions: " + "; ".join(v.assumptions))
    if not a.out:
        print(text)
    return 0


def cmd_demo(a: argparse.Namespace) -> int:
    print("rtsc demo: SYNTHETIC spectra (invented; not data for any material), no DFT.\n")
    raw = {"run": {"name": "rtsc-demo (SYNTHETIC)", "workdir": str(Path(a.workdir).resolve()),
                   "workers": a.workers, "resume": not a.no_resume},
           "sources": [{"type": "synthetic_demo"}],
           "budget": {"max_materials": 4},
           "spectra": {"rel_envelope": 0.05},
           "certify": {"T_target_K": 300.0, "mu_hi": 0.15, "omega_c_factor": 10.0,
                       "cutoff_rule": "block", "mode": "decide"}}
    cfg = pipeline.normalize_config(raw)
    summ = pipeline.run(cfg)
    _print_summary(summ)
    print()
    rows = json.loads((Path(cfg["run"]["workdir"]) / "candidates.json").read_text())
    for r in rows:
        lo = "-" if r["T_low_K"] is None else repr(r["T_low_K"])
        up = "-" if r["T_up_K"] is None else repr(r["T_up_K"])
        print(f"  {r['material_id']:<38} {r['status']:<20} T_low={lo:<20} T_up={up:<20} "
              f"{','.join(r['rule_ids'])}")
    print(f"\nreports in {cfg['run']['workdir']} (summary.md, candidates.csv, certificates/)")
    return 0


def _load_recheck_inputs(path: Path, spectra_path: str | None):
    """Return (verdict dict, spectra dicts, band-coupling file entries) from a verdict.json, a
    certificates/<id>.json or a record.json.  For a bare verdict.json the spectra (and band
    coupling files) come from --spectra or from the record.json next to it."""
    d = json.loads(path.read_text())
    spectra = None
    bcf: list[dict] = []
    if "verdict" in d and "status" not in d:          # certificates/<id>.json or record.json
        spectra = d.get("spectra")
        bcf = list((d.get("candidate") or {}).get("band_coupling_files") or [])
        d = d["verdict"]
    if spectra_path:
        s = json.loads(Path(spectra_path).read_text())
        spectra = s.get("spectra", s) if isinstance(s, dict) else s
    if spectra is None and (path.parent / "record.json").is_file():
        rec = json.loads((path.parent / "record.json").read_text())
        if rec.get("material_id") == d.get("material_id"):
            spectra = rec.get("spectra")
            bcf = list((rec.get("candidate") or {}).get("band_coupling_files") or [])
    if not spectra:
        raise FileNotFoundError("no spectra found: pass --spectra FILE (list of Spectrum dicts, "
                                "a record.json or a certificates/<id>.json)")
    return d, spectra, bcf


def cmd_recheck(a: argparse.Namespace) -> int:
    from .certify import recheck
    from .spectrum import Spectrum
    from .status import Verdict
    vd, sd, bcf = _load_recheck_inputs(Path(a.verdict), a.spectra)
    v = Verdict.from_dict(vd)
    paths = [(x, None) for x in (a.band_coupling or [])] or \
        [(x["path"], x.get("table")) for x in bcf]
    kw = {}
    if paths:
        from .dfpt import BandCoupling
        kw["band_couplings"] = [BandCoupling.from_file(p, t) for p, t in paths]
    if getattr(a, "config", None):
        from .pipeline import _certify_params, load_config
        kw["expected_params"] = _certify_params(load_config(a.config).get("certify", {}))
    else:
        print("warning: no --config given; run parameters (T*, mu_hi, omega_c, row_check, "
              "rejection_model) are taken from the certificate itself", file=sys.stderr)
    res = recheck(v, [Spectrum.from_dict(x) for x in sd], **kw)
    ok = bool(res.get("ok"))
    if a.verbose_json:
        print(json.dumps(res, indent=1, default=str))
    print(f"{v.material_id}: stored {v.status.value} [{v.T_low_K!r}, {v.T_up_K!r}] -> rechecked "
          f"{res.get('status')} [{res.get('T_low_K')!r}, {res.get('T_up_K')!r}]")
    for e in res.get("errors", []):
        print(f"  error: {e}")
    print("RECHECK OK" if ok else "RECHECK FAILED")
    return 0 if ok else 3


def _discovery_config(a: argparse.Namespace):
    from .discovery.goals import load_discovery_config, normalize_discovery
    from .discovery.ledger import Ledger
    if getattr(a, "config", None):
        dc = load_discovery_config(a.config)
        if getattr(a, "workdir", None):
            dc.d["workdir"] = str(Path(a.workdir).resolve())
        return dc
    if not getattr(a, "workdir", None):
        raise pipeline.ConfigError("give -c CONFIG or -w WORKDIR")
    lp = Path(a.workdir) / "ledger.sqlite"
    if not lp.is_file():
        raise FileNotFoundError(f"no ledger at {lp}")
    led = Ledger(lp)
    stored = led.get_meta("config")
    led.close()
    dc = normalize_discovery(stored["raw"], stored.get("base_dir"))
    dc.d["workdir"] = str(Path(a.workdir).resolve())
    return dc


def _limits(a: argparse.Namespace) -> dict[str, Any]:
    return {"certified_per_goal": a.certified_per_goal, "max_rounds": a.max_rounds,
            "max_wall_hours": a.max_wall_hours, "max_dft_jobs": a.max_dft_jobs}


def _print_discovery(res: dict[str, Any]) -> None:
    print(f"stop: {res['stop_reason']}")
    print(f"rounds this session: {res['rounds_this_session']}")
    for g, ids in res["new_certified"].items():
        print(f"  {g}: {len(ids)} newly certified this session, {res['certified_total'][g]} in total")
    acc = res["accounting"]
    print(f"ledger: {acc['generated']} generated = {acc['materials']} materials + "
          f"{acc['duplicates']} duplicates (S4 checked)")


def cmd_discover_run(a: argparse.Namespace) -> int:
    from .discovery.loop import DiscoveryLoop
    dc = _discovery_config(a)
    lp = dc.workdir / "ledger.sqlite"
    if lp.is_file():
        from .discovery.ledger import Ledger
        led = Ledger(lp)
        n = len(led.rounds())
        led.close()
        if n:
            print(f"error: {lp} already has {n} rounds; use `rtsc discover continue`", file=sys.stderr)
            return 2
    res = DiscoveryLoop(dc).run("run", **_limits(a))
    _print_discovery(res)
    print(f"report: {dc.workdir / 'discovery_report.md'}")
    return 0


def cmd_discover_continue(a: argparse.Namespace) -> int:
    from .discovery.loop import DiscoveryLoop
    dc = _discovery_config(a)
    if not (dc.workdir / "ledger.sqlite").is_file():
        raise FileNotFoundError(f"no ledger in {dc.workdir}: use `rtsc discover run`")
    res = DiscoveryLoop(dc).run("continue", **_limits(a))
    _print_discovery(res)
    print(f"report: {dc.workdir / 'discovery_report.md'}")
    return 0


def _open_ledger(a: argparse.Namespace):
    from .discovery.ledger import Ledger
    dc = _discovery_config(a)
    lp = dc.workdir / "ledger.sqlite"
    if not lp.is_file():
        raise FileNotFoundError(f"no ledger at {lp}")
    return dc, Ledger(lp)


def cmd_discover_status(a: argparse.Namespace) -> int:
    from .status import Status
    dc, led = _open_ledger(a)
    goals = [g for g in dc.enabled_goals if g.id in (led.get_meta("goal_ids") or [])]
    acc = led.accounting([g.id for g in goals])
    rounds = led.rounds()
    print(f"workdir: {dc.workdir}")
    print(f"rounds: {len(rounds)} (last phase: {rounds[-1]['phase'] if rounds else '-'})")
    print(f"generated: {acc['generated']} = {acc['materials']} materials + {acc['duplicates']} duplicates")
    from .discovery.report import NOT_COUNTING, goal_rows
    for g in goals:
        rows = [r for r in goal_rows(led, g.id, dc) if r["in_scope"]]
        c = {s.value: sum(1 for r in rows if r["status_shown"] == s.value) for s in Status}
        c["CERTIFIED (does not count)"] = sum(1 for r in rows if r["status_shown"] == NOT_COUNTING)
        print(f"  {g.id} (T*={g.T_star_K:g} K, P<={g.P_max_GPa:g} GPa): "
              + "  ".join(f"{k}={v}" for k, v in c.items()))
    s = led.q("SELECT * FROM sessions ORDER BY id DESC LIMIT 1")
    if s:
        print(f"last session: {s[0]['command']} started {s[0]['started_utc']}, "
              f"stop: {s[0]['stop_reason'] or '(unfinished: crashed or running)'}")
    return 0


def cmd_discover_list(a: argparse.Namespace) -> int:
    from .discovery.report import goal_rows
    dc, led = _open_ledger(a)
    for g in dc.enabled_goals:
        if a.goal and g.id != a.goal:
            continue
        rows = [r for r in goal_rows(led, g.id, dc) if not a.status or r["status_shown"] == a.status]
        print(f"== {g.id} (T*={g.T_star_K:g} K)")
        for r in rows[:a.limit]:
            lo = "-" if r["T_low_K"] is None else f"{r['T_low_K']:.4g}"
            up = "-" if r["T_up_K"] is None else f"{r['T_up_K']:.4g}"
            print(f"  {r['material_id']:<28} {r['pressure_GPa']:>7g} GPa  {r['status_shown']:<20} "
                  f"T_low={lo:<10} T_up={up:<10} {r['stage']}")
    return 0


def cmd_discover_export(a: argparse.Namespace) -> int:
    from .discovery.report import export_finetune, export_json
    dc, led = _open_ledger(a)
    goals = dc.enabled_goals
    if a.finetune:
        out = Path(a.out or dc.workdir / "finetune")
        try:   # MatterGen-format exporter of rtsc.generate (verified against the MatterGen repo)
            from .generate.finetune_export import export_finetune_dataset, records_from_ledger
        except ImportError:
            export_finetune_dataset = None
        if export_finetune_dataset is not None:
            synth = {r["material_id"] for g in goals for r in led.current(g.id)
                     if "SYNTHETIC_DATA" in r["flags"]}
            recs = [r for r in records_from_ledger(led.path)
                    if a.include_synthetic or r["material_id"] not in synth]
            export_finetune_dataset(recs, out)
            print(f"finetune export: {len(recs)} labelled structures in {out} "
                  f"(MatterGen CSV layout; {'including' if a.include_synthetic else 'excluding'} "
                  "SYNTHETIC)")
        else:
            res = export_finetune(led, goals, out, include_synthetic=a.include_synthetic)
            print(f"finetune export: {res['n']} labelled structures ({res['train']} train, "
                  f"{res['val']} val) in {out}")
    else:
        out = Path(a.out or dc.workdir / "ledger_export.json")
        export_json(led, goals, out)
        print(f"ledger exported to {out}")
    return 0


def _add_discover(sub) -> None:
    dp = sub.add_parser("discover", help="autonomous discovery loop (DISCOVERY.md)")
    ds = dp.add_subparsers(dest="discover_cmd", required=True)

    def common(p, need_config: bool):
        p.add_argument("-c", "--config", required=need_config, help="discovery TOML config")
        p.add_argument("-w", "--workdir", help="discovery workdir (overrides [discovery] workdir)")

    def limits(p):
        p.add_argument("--certified-per-goal", type=int)
        p.add_argument("--max-rounds", type=int)
        p.add_argument("--max-wall-hours", type=float)
        p.add_argument("--max-dft-jobs", type=int)

    r = ds.add_parser("run", help="start a new discovery run")
    common(r, True)
    limits(r)
    r.set_defaults(fn=cmd_discover_run)
    c = ds.add_parser("continue", help="resume from the ledger with fresh stop limits")
    common(c, False)
    limits(c)
    c.set_defaults(fn=cmd_discover_continue)
    s = ds.add_parser("status", help="progress summary from the ledger")
    common(s, False)
    s.set_defaults(fn=cmd_discover_status)
    li = ds.add_parser("list", help="materials and current statuses per goal")
    common(li, False)
    li.add_argument("--goal")
    li.add_argument("--status")
    li.add_argument("--limit", type=int, default=50)
    li.set_defaults(fn=cmd_discover_list)
    e = ds.add_parser("export", help="export the ledger (JSON) or fine-tuning labels (--finetune)")
    common(e, False)
    e.add_argument("--finetune", action="store_true")
    e.add_argument("--include-synthetic", action="store_true")
    e.add_argument("--out")
    e.set_defaults(fn=cmd_discover_export)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="rtsc", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-v", "--verbose", action="store_true")
    sub = p.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("run", help="run the pipeline from a TOML config")
    r.add_argument("-c", "--config", required=True)
    r.add_argument("--shard", help="K/N: compute only every N-th material (job arrays)")
    r.add_argument("--plan-only", action="store_true", help="stages 1-5 only (sources..budget)")
    r.add_argument("--workers", type=int)
    r.add_argument("--workdir")
    r.add_argument("--no-resume", action="store_true")
    r.set_defaults(fn=cmd_run)

    c = sub.add_parser("collect", help="merge records of a sharded run and write reports")
    c.add_argument("-c", "--config", required=True)
    c.set_defaults(fn=cmd_collect)

    k = sub.add_parser("certify", help="classify one alpha^2F file")
    k.add_argument("file")
    k.add_argument("--units", help="unit of the omega column: K meV eV THz cm-1 Ry mRy")
    k.add_argument("--format", default="two_column", choices=["two_column", "qe_a2f_dos", "epw"])
    k.add_argument("--columns", default="0,1")
    k.add_argument("--rule", default="bins", choices=["bins", "trapezoid", "midpoint"])
    k.add_argument("--negative-policy", default="error", choices=["error", "clip"])
    k.add_argument("--mu-hi", type=float, default=0.15)
    k.add_argument("--target", type=float, default=300.0, help="T* in kelvin")
    k.add_argument("--envelope", type=float, default=0.10, help="relative weight envelope e")
    k.add_argument("--omega-c", type=float, help="fixed Coulomb cutoff omega_c in kelvin")
    k.add_argument("--omega-c-factor", type=float, default=10.0,
                   help="omega_c = factor * max omega (if --omega-c not given)")
    k.add_argument("--mode", default="bracket", choices=["bracket", "decide"],
                   help="bracket: tightest certified bracket; decide: stop once decided")
    k.add_argument("--expected-lambda", type=float, help="known lambda of the full file (completeness check)")
    k.add_argument("--expected-A", type=float, help="known A in kelvin of the full file (completeness check)")
    k.add_argument("--allow-unverified-rejection", action="store_true",
                   help="keep REJECTED although the file's completeness is unverified")
    k.add_argument("--rejection-model", default="isotropic", choices=["isotropic", "multiband"])
    k.add_argument("--band-coupling", action="append",
                   help="rtsc.dfpt.BandCoupling .json/.toml (repeatable; for multiband)")
    k.add_argument("--material-id")
    k.add_argument("--out")
    k.set_defaults(fn=cmd_certify)

    d = sub.add_parser("demo", help="end-to-end run on SYNTHETIC spectra (no DFT)")
    d.add_argument("--workdir", default="rtsc_demo")
    d.add_argument("--workers", type=int, default=1)
    d.add_argument("--no-resume", action="store_true")
    d.set_defaults(fn=cmd_demo)

    v = sub.add_parser("recheck", help="re-verify a verdict/certificate JSON")
    v.add_argument("verdict", help="verdict.json, certificates/<id>.json or record.json")
    v.add_argument("--spectra", help="JSON list of Spectrum dicts (if not next to the verdict)")
    v.add_argument("--band-coupling", action="append",
                   help="BandCoupling file(s) if not recorded with the verdict")
    v.add_argument("--json", dest="verbose_json", action="store_true", help="print full result")
    v.add_argument("-c", "--config", help="run config TOML: bind the expected [certify] parameters")
    v.set_defaults(fn=cmd_recheck)

    _add_discover(sub)
    return p


def main(argv: list[str] | None = None) -> int:
    a = build_parser().parse_args(argv)
    logging.basicConfig(level=logging.INFO if a.verbose else logging.WARNING,
                        format="%(levelname)s %(name)s: %(message)s")
    try:
        return a.fn(a)
    except (pipeline.ConfigError, FileNotFoundError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    except Exception as exc:  # noqa: BLE001
        if a.verbose:
            raise
        print(f"error: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
