"""ORCH: end-to-end orchestration with tiny synthetic inputs.

S4 accounting, S1 origin guard, ARCH 4b F1 (zero DFPT spectra) and F3 (cutoff recorded),
resume, sharding + collect, process pool, reports.  Uses the TEST-ONLY certify stub unless the
real rtsc.certify.classify exists (tests that need specific outcomes force the stub).
"""
from __future__ import annotations

import csv
import json
import tomllib
from pathlib import Path

import pytest

import rtsc.pipeline as P
from rtsc.sources import Candidate
from rtsc.sources.synthetic import demo_candidates, gaussian_peak_spectrum
from rtsc.status import Status, Verdict

import test_pipeline_stubs as stubs


def _cfg(tmp_path, **over):
    raw = {"run": {"workdir": str(tmp_path / "work"), "resume": True},
           "budget": {"max_materials": 4},
           "spectra": {"rel_envelope": 0.05}}
    for k, v in over.items():
        raw[k] = v
    return P.normalize_config(raw)


def _mixed_inputs():
    good = demo_candidates()
    dup_a = Candidate("dupA", "NbH", spacegroup="Fm-3m", pressure_GPa=0, enthalpy_eV_per_atom=-1.0,
                      properties={"is_metal": True})
    dup_b = Candidate("dupB", "NbH", spacegroup="Fm-3m", pressure_GPa=0, enthalpy_eV_per_atom=-2.0,
                      properties={"is_metal": True})
    broken = Candidate("broken", {}, properties={"load_error": "ValueError: garbage"})
    return good + [dup_a, dup_b, broken]


def test_every_input_has_exactly_one_status(tmp_path, monkeypatch):
    stubs.install(monkeypatch, force=True)
    cfg = _cfg(tmp_path, budget={"max_materials": 5})
    cands = _mixed_inputs()
    ids = [c.material_id for c in cands]
    summ = P.run(cfg, cands)
    assert sorted(summ.inputs) == sorted(ids)
    assert sum(summ.counts().values()) == len(ids)
    by = {m: v.status for m, v in summ.verdicts.items()}
    assert by["broken"] is Status.ERROR
    assert by["dupA"] is Status.DEFERRED                      # duplicate of the lower-enthalpy dupB
    assert summ.duplicate_of == {"dupA": "dupB"}
    assert by["SYNTHETIC-above-hull"] is Status.EXCLUDED_ASSUMPTION
    assert by["SYNTHETIC-insulator"] is Status.EXCLUDED_ASSUMPTION
    assert by["SYNTHETIC-weak-einstein"] is Status.REJECTED_PROVEN
    assert by["SYNTHETIC-strong-coupling"] is Status.CERTIFIED_CANDIDATE
    # budget 5 among SYNTHETIC weak/high-freq/strong/borderline/low-priority + dupB
    deferred = [m for m, s in by.items() if s is Status.DEFERRED]
    assert "SYNTHETIC-low-priority" in deferred
    # dupB has no spectrum and DFPT is off -> DEFERRED (not computed), never a proven status
    assert by["dupB"] is Status.DEFERRED
    # every synthetic verdict carries the SYNTHETIC_DATA flag
    for m, v in summ.verdicts.items():
        if m.startswith("SYNTHETIC-"):
            assert "SYNTHETIC_DATA" in v.flags, m
    # reports: one row per input, sorted by status then T_low
    work = Path(cfg["run"]["workdir"])
    rows = list(csv.DictReader(open(work / "candidates.csv")))
    assert len(rows) == len(ids)
    order = [r["status"] for r in rows]
    rank = {s: i for i, s in enumerate(["CERTIFIED_CANDIDATE", "UNDECIDED", "REJECTED_PROVEN",
                                        "EXCLUDED_ASSUMPTION", "DEFERRED", "ERROR"])}
    assert order == sorted(order, key=rank.get)
    md = (work / "summary.md").read_text()
    assert "SYNTHETIC DATA" in md and "isotropic Migdal-Eliashberg" in md
    assert f"| **total** | **{len(ids)}** |" in md
    for m in ids:
        assert (work / "certificates" / f"{P.safe_id(m)}.json").is_file()
    v = json.loads((work / P.safe_id("SYNTHETIC-strong-coupling") / "verdict.json").read_text())
    assert Verdict.from_dict(v).status is Status.CERTIFIED_CANDIDATE


def test_accounting_detects_violations():
    v = Verdict("a", Status.DEFERRED)
    with pytest.raises(P.AccountingError):
        P.account(["a", "b"], {"a": v})                         # missing
    with pytest.raises(P.AccountingError):
        P.account(["a"], {"a": v, "z": Verdict("z", Status.ERROR)})   # unknown id
    with pytest.raises(P.AccountingError):
        P.account(["a", "a"], {"a": v})                         # duplicate input
    with pytest.raises(P.AccountingError):                      # S1: proven status, wrong origin
        P.account(["a"], {"a": Verdict("a", Status.REJECTED_PROVEN)}, {"a": "pipeline.budget"})
    P.account(["a"], {"a": v}, {"a": "pipeline.budget"})


def test_run_raises_if_a_stage_drops_a_material(tmp_path, monkeypatch):
    stubs.install(monkeypatch, force=True)
    real = P.compute_material

    def lossy(cand, rec, cfg, workdir):
        out = real(cand, rec, cfg, workdir)
        if cand["material_id"] == "SYNTHETIC-strong-coupling":
            out["verdict"] = None                               # simulate a lost result
        return out
    monkeypatch.setattr(P, "compute_material", lossy)
    with pytest.raises(P.AccountingError):
        P.run(_cfg(tmp_path), demo_candidates(), write_reports=False)


def test_F1_zero_dfpt_spectra_is_error_never_classified(tmp_path, monkeypatch):
    """ARCH 4b F1: DFPT returning zero alternatives -> ERROR; classify is never called."""
    called = []

    def spy(spectra, params, **kw):
        called.append(len(spectra))
        return stubs.stub_classify(spectra, params)
    stubs.install(monkeypatch, force=True)
    import rtsc.certify as C
    monkeypatch.setattr(C, "classify", spy)
    monkeypatch.setattr(P, "run_dfpt", lambda c, d, cfg: ([], {"note": "all smearings failed"}))
    c = Candidate("x", "NbH", lattice=[[3, 0, 0], [0, 3, 0], [0, 0, 3]], species=["Nb", "H"],
                  frac_coords=[[0, 0, 0], [.5, .5, .5]], properties={"is_metal": True})
    cfg = _cfg(tmp_path, dfpt={"enabled": True})
    summ = P.run(cfg, [c], write_reports=False)
    assert summ.verdicts["x"].status is Status.ERROR
    assert "F1" in summ.verdicts["x"].reasons[0]
    assert called == []
    with pytest.raises(ValueError):
        P.certify_spectra("x", [], cfg["certify"])          # hard guard before classify


def test_F3_cutoff_in_assumptions_and_config(tmp_path, monkeypatch):
    stubs.install(monkeypatch, force=True)
    cfg = _cfg(tmp_path, certify={"T_target_K": 300.0, "mu_hi": 0.15, "omega_c_K": 50000.0,
                                  "cutoff_rule": "block"})
    summ = P.run(cfg, [demo_candidates()[2]])
    v = summ.verdicts["SYNTHETIC-strong-coupling"]
    cut = [a for a in v.assumptions if "omega_c" in a]
    assert cut and "50000 K" in cut[0] and "block" in cut[0]
    rows = json.loads((Path(cfg["run"]["workdir"]) / "candidates.json").read_text())
    assert "omega_c" in rows[0]["cutoff"] and any("omega_c" in a for a in rows[0]["assumptions"])
    assert "omega_c" in (Path(cfg["run"]["workdir"]) / "summary.md").read_text()
    with pytest.raises(P.ConfigError):
        _cfg(tmp_path, certify={"cutoff_rule": "smooth"})
    d = P.cutoff_description({"omega_c_factor": 10.0}, [P.Spectrum((100.0,), (1.0,))])
    assert "10*max(omega)" in d and "1000" in d


def test_resume_skips_finished_materials(tmp_path, monkeypatch):
    stubs.install(monkeypatch, force=True)
    cfg = _cfg(tmp_path)
    s1 = P.run(cfg, demo_candidates())
    calls = []
    import rtsc.certify as C
    orig = C.classify
    monkeypatch.setattr(C, "classify", lambda s, p, **kw: calls.append(1) or orig(s, p))
    s2 = P.run(cfg, demo_candidates())
    assert calls == []
    assert {m: v.status for m, v in s1.verdicts.items()} == {m: v.status for m, v in s2.verdicts.items()}
    cfg2 = _cfg(tmp_path, certify={"T_target_K": 250.0, "mu_hi": 0.15})   # new hash -> recompute
    P.run(cfg2, demo_candidates())
    assert calls


def test_shards_and_collect(tmp_path, monkeypatch):
    stubs.install(monkeypatch, force=True)
    cfg = _cfg(tmp_path)
    ids = [c.material_id for c in demo_candidates()]
    P.run(cfg, demo_candidates(), shard="0/2")
    with_missing = P.collect(cfg, write=False)                 # shard 1 not run: explicit ERRORs
    assert sum(with_missing.counts().values()) == len(ids)
    assert with_missing.counts()["ERROR"] >= 1
    P.run(cfg, demo_candidates(), shard="1/2")
    full = P.collect(cfg)
    assert full.counts()["ERROR"] == 0 and sorted(full.inputs) == sorted(ids)


def test_process_pool(tmp_path, monkeypatch):
    import multiprocessing as mp
    if mp.get_start_method() != "fork":
        pytest.skip("the test stub reaches workers only with the fork start method")
    stubs.install(monkeypatch, force=True)
    cfg = _cfg(tmp_path, run={"workdir": str(tmp_path / "pool"), "workers": 2, "resume": False})
    summ = P.run(cfg, demo_candidates())
    assert sum(summ.counts().values()) == len(demo_candidates())


def test_worker_exception_becomes_error(tmp_path, monkeypatch):
    stubs.install(monkeypatch, force=True)
    import rtsc.certify as C

    def boom(spectra, params, **kw):
        raise RuntimeError("engine failure")
    monkeypatch.setattr(C, "classify", boom)
    summ = P.run(_cfg(tmp_path), demo_candidates()[:1])
    v = summ.verdicts["SYNTHETIC-weak-einstein"]
    assert v.status is Status.ERROR and "engine failure" in v.reasons[0]


def test_attached_two_column_file(tmp_path, monkeypatch):
    stubs.install(monkeypatch, force=True)
    f = tmp_path / "a2f.dat"
    # SYNTHETIC two-column (omega meV, a2F): tiny weight -> stub ceiling rejects
    f.write_text("\n".join(f"{w} {0.01}" for w in range(1, 40)) + "\n")
    c = Candidate("SYNTHETIC-file", {}, formula="(synthetic)", synthetic=True,
                  spectra_files=[{"path": str(f), "units": "meV"}], properties={"is_metal": True})
    summ = P.run(_cfg(tmp_path), [c])
    v = summ.verdicts["SYNTHETIC-file"]
    # no completeness check for a plain two-column file: the rejection is blocked (policy)
    assert v.status is Status.UNDECIDED and "blocked" in v.reasons[-1]
    assert any(a.startswith("A-FILE-COMPLETE(unverified)") for a in v.assumptions)
    summ = P.run(_cfg(tmp_path, certify={"T_target_K": 300.0, "mu_hi": 0.15,
                                         "allow_unverified_rejection": True}), [c])
    v = summ.verdicts["SYNTHETIC-file"]
    assert v.status is Status.REJECTED_PROVEN
    assert any(a.startswith("TRUST: python-flint") for a in v.assumptions)
    rec = P.read_record(_cfg(tmp_path)["run"]["workdir"], "SYNTHETIC-file")
    assert rec["spectra"][0]["rel_envelope"] == 0.05      # config envelope applied (widened)
    assert rec["spectra"][0]["provenance"]["synthetic"] is True


def test_config_validation(tmp_path):
    with pytest.raises(P.ConfigError):
        P.normalize_config({"bogus": {}})
    with pytest.raises(P.ConfigError):
        P.normalize_config({"run": {"wrkdir": "x"}})
    with pytest.raises(ValueError):
        P.normalize_config({"validity": {"nonadiabatic": 1}})


def test_example_configs_parse():
    root = Path(__file__).resolve().parents[1] / "examples"
    files = sorted(root.rglob("config.toml"))
    assert len(files) >= 3
    for f in files:
        if "discovery" in tomllib.loads(f.read_text()):
            # discovery-loop configs: [discovery] + pipeline sections, no [[sources]]
            from rtsc.discovery.goals import load_discovery_config
            dc = load_discovery_config(f)
            assert dc.enabled_goals and dc.d["generators"]
            continue
        cfg = P.load_config(f)
        assert cfg["sources"]


def test_F1_dfpt_zero_even_with_attached_spectra(tmp_path, monkeypatch):
    stubs.install(monkeypatch, force=True)
    monkeypatch.setattr(P, "run_dfpt", lambda c, d, cfg: ([], {}))
    c = demo_candidates()[2]                       # has attached SYNTHETIC spectra
    cfg = _cfg(tmp_path, dfpt={"enabled": True, "also_when_attached": True})
    summ = P.run(cfg, [c], write_reports=False)
    assert summ.verdicts[c.material_id].status is Status.ERROR


class _FakeResult:
    def __init__(self, spectra, flags):
        self.spectra, self.flags = spectra, flags

    def to_dict(self):
        return {"spectra": [s.to_dict() for s in self.spectra], "flags": self.flags,
                "warnings": [], "provenance": {"synthetic": True}}


@pytest.mark.parametrize("min_cm1,check_flag,expect", [
    (-50.0, "DFPT_CHECK_FAIL:lambda_vs_lambda_file:sigma1", Status.EXCLUDED_ASSUMPTION),
    (-1.0, "DFPT_CHECK_FAIL:lambda_vs_lambda_file:sigma1", Status.CERTIFIED_CANDIDATE),
    (-1.0, "DFPT_CHECK_FAIL:lambda_vs_a2F_footer:sigma1", Status.ERROR)])
def test_dfpt_adapter_with_fake_workflow(tmp_path, monkeypatch, min_cm1, check_flag, expect):
    """run_dfpt drives rtsc.dfpt.EPhWorkflow (faked here, SYNTHETIC spectra); the matdyn-grid
    minimum frequency feeds the dynamic-stability filter.  A convergence-type check failure
    (lambda_vs_lambda_file) is only a flag; a parse/units failure (lambda_vs_a2F_footer) is
    ERROR before classify (audit F8).  use_smearings gives the SMEARINGS_SUBSET flag."""
    dfpt = pytest.importorskip("rtsc.dfpt")
    stubs.install(monkeypatch, force=True)
    seen = {}

    class FakeWF:
        def __init__(self, material_id, structure, workdir, settings, pseudos, launcher=None, masses=None):
            seen.update(mid=material_id, structure=structure, workdir=str(workdir),
                        pressure=settings.pressure_GPa, nprocs=launcher.nprocs)

        def run(self):
            from rtsc.spectrum import Spectrum
            sp = [Spectrum.from_dict(gaussian_peak_spectrum(2000.0, 100.0, 3000.0 + 50 * k,
                                                            label=f"SYNTHETIC sigma{k}"))
                  for k in range(3)]
            for s in sp:
                s.provenance["min_freq_matdyn_grid_cm1"] = min_cm1
            return _FakeResult(sp, [check_flag])
    monkeypatch.setattr(dfpt, "EPhWorkflow", FakeWF)
    c = Candidate("SYNTHETIC-dfpt", "NbH", lattice=[[3, 0, 0], [0, 3, 0], [0, 0, 3]],
                  species=["Nb", "H"], frac_coords=[[0, 0, 0], [.5, .5, .5]], pressure_GPa=100.0,
                  properties={"is_metal": True}, synthetic=True)
    cfg = _cfg(tmp_path, dfpt={"enabled": True, "launcher": {"kind": "mpirun", "nprocs": 4},
                               "pseudos": {"dir": str(tmp_path)},
                               "settings": {"do_relax": True, "use_smearings": [1, 2, 3]}},
               filters={"dynamic": {"imaginary_tol": 10.0, "imaginary_tol_unit": "cm-1"}})
    summ = P.run(cfg, [c], write_reports=False)
    v = summ.verdicts["SYNTHETIC-dfpt"]
    assert seen["pressure"] == 100.0 and seen["nprocs"] == 4 and seen["mid"] == "SYNTHETIC-dfpt"
    st = seen["structure"]
    assert hasattr(st, "lattice") or st["masses"]["Nb"] > 92   # pymatgen Structure or dict
    assert v.status is expect
    if expect is Status.CERTIFIED_CANDIDATE:
        assert check_flag in v.flags and "SMEARINGS_SUBSET" in v.flags
        rec = P.read_record(cfg["run"]["workdir"], "SYNTHETIC-dfpt")
        assert len(rec["spectra"]) == 3                  # all alternatives went to classify
    elif expect is Status.ERROR:
        assert "lambda_vs_a2F_footer" in v.reasons[0]
    else:
        assert v.rule_ids == ["X-DYNAMIC"]


def test_safe_id_is_injective_and_read_record_checks_id(tmp_path):
    ids = ["X/1", "X_1", "X 1", "a.b", "a_b", "..", "", "H3S#0", "H3S_0"]
    paths = {P.record_path(tmp_path, m) for m in ids}
    assert len(paths) == len(ids)
    assert P.safe_id("LaH10") == "LaH10"                 # plain ids unchanged
    p = P.record_path(tmp_path, "A")
    p.parent.mkdir(parents=True)
    p.write_text(json.dumps({"material_id": "B", "verdict": {"status": "REJECTED_PROVEN"}}))
    assert P.read_record(tmp_path, "A") is None           # belongs to another material
    with pytest.raises(P.ConfigError):
        P.check_path_collisions(["abc", "ABC"])          # case-insensitive file systems


def test_input_digest_tracks_file_contents(tmp_path):
    f = tmp_path / "s.vasp"
    f.write_text("x")
    c = Candidate("m", "NbH", provenance={"file": str(f)})
    d1 = P.input_digest(c)
    f.write_text("y")
    assert P.input_digest(c) != d1


def test_annotations_cutoff_discretisation_singleband_coulomb():
    sp = P.Spectrum((100.0, 200.0), (1.0, 1.0), 0.1, "t",
                    {"rule": "bins: atom at w_i", "source_file": "/x/a2F.dat"})
    v = Verdict("m", Status.CERTIFIED_CANDIDATE, T_low_K=400.0, T_up_K=500.0,
                certificate={"params": {"omega_c": 1000.0}})
    out = P.annotate_model(v, {"omega_c_factor": 5.0}, [sp])
    assert "DISCRETISATION(bins, /x/a2F.dat)" in out.assumptions
    assert "COULOMB_INACTIVE" in out.flags                  # pi*400 > 1000
    r = P.annotate_model(Verdict("m", Status.REJECTED_PROVEN, T_up_K=10.0,
                                 rule_ids=["R-UNIV"]), {"omega_c_factor": 10.0}, [sp])
    assert any(a.startswith("M-ISO-SINGLEBAND") for a in r.assumptions)
    assert r.status is Status.REJECTED_PROVEN and r.T_up_K == 10.0


def test_band_coupling_files_reach_classify(tmp_path, monkeypatch):
    pytest.importorskip("rtsc.dfpt")
    stubs.install(monkeypatch, force=True)
    bc = tmp_path / "bc.json"
    bc.write_text(json.dumps({"labels": ["a", "b"], "unit": "K", "kind": "computed",
                              "A": [[100.0, 20.0], [20.0, 50.0]], "N": [1.0, 1.0],
                              "provenance": {"note": "SYNTHETIC test data"}}))
    got = {}
    import rtsc.certify as C

    def spy(spectra, params, material_id="", band_couplings=None, **kw):
        got["bcs"] = band_couplings
        return stubs.stub_classify(spectra, params, material_id)
    monkeypatch.setattr(C, "classify", spy)
    d = tmp_path / "structs"
    d.mkdir()
    (d / "m.vasp").write_text("t\n1.0\n3 0 0\n0 3 0\n0 0 3\nNb H\n1 1\nDirect\n0 0 0\n.5 .5 .5\n")
    (d / "a2f.dat").write_text("\n".join(f"{w} 0.01" for w in range(1, 30)) + "\n")
    (d / "manifest.toml").write_text('[m]\nband_coupling = "../bc.json"\nsynthetic = true\n'
                                     'properties = { is_metal = true }\n'
                                     'spectra_files = [ { path = "a2f.dat", units = "meV" } ]\n')
    cfg = _cfg(tmp_path, sources=[{"type": "local", "path": str(d), "backend": "builtin"}])
    summ = P.run(cfg)
    assert summ.verdicts["m"].status is not Status.ERROR, summ.verdicts["m"].reasons
    assert got["bcs"] and got["bcs"][0].labels == ("a", "b")


def test_collect_refuses_other_config_and_report_uses_certificate_params(tmp_path, monkeypatch):
    stubs.install(monkeypatch, force=True)
    cfg = _cfg(tmp_path)
    P.run(cfg, demo_candidates()[:1])
    with pytest.raises(P.ConfigError):
        P.collect(_cfg(tmp_path, certify={"T_target_K": 30.0, "mu_hi": 0.15}))
    from rtsc.report import row_for
    rec = P.read_record(cfg["run"]["workdir"], "SYNTHETIC-weak-einstein")
    rec["verdict"]["certificate"] = {"params": {"T_star": 300.0, "mu_hi": 0.15}}
    row = row_for(rec, _cfg(tmp_path, certify={"T_target_K": 30.0, "mu_hi": 0.15}), tmp_path)
    assert row["T_target_K"] == 300.0 and "PARAMS_MISMATCH" in row["flags"]


@pytest.mark.parametrize("tag,status,allow,expect", [
    ("A-FILE-COMPLETE(unverified): no end marker", Status.REJECTED_PROVEN, False, Status.UNDECIDED),
    ("A-DROP-NONPOS(12.5K): weight at omega<=0", Status.REJECTED_PROVEN, False, Status.UNDECIDED),
    ("A-DROP-NONPOS(0K): none", Status.REJECTED_PROVEN, False, Status.REJECTED_PROVEN),
    ("A-FILE-INCOMPLETE: no footer", Status.CERTIFIED_CANDIDATE, False, Status.UNDECIDED),
    ("A-NEG-CLIP(3K): clipped", Status.CERTIFIED_CANDIDATE, False, Status.UNDECIDED),
    ("A-NEG-CLIP(3K): clipped", Status.REJECTED_PROVEN, False, Status.REJECTED_PROVEN),
    ("A-FILE-COMPLETE(unverified): no end marker", Status.REJECTED_PROVEN, True, Status.REJECTED_PROVEN),
    ("A-NEG-CLIP(3K): clipped", Status.CERTIFIED_CANDIDATE, True, Status.CERTIFIED_CANDIDATE)])
def test_input_quality_policy_only_downgrades(tag, status, allow, expect):
    sp = [P.Spectrum((100.0,), (1.0,), 0.1, "alt0"),
          P.Spectrum((100.0,), (1.0,), 0.1, "alt1", {"assumptions": [tag]})]
    v = Verdict("m", status, T_low_K=400.0 if status is Status.CERTIFIED_CANDIDATE else None,
                T_up_K=50.0, certificate={"c": 1})
    key = ("allow_unverified_rejection" if status is Status.REJECTED_PROVEN
           else "allow_unverified_certification")
    out = P.apply_input_policy(P.annotate_model(v, {"omega_c_factor": 10.0}, sp), sp, {key: allow})
    assert out.status is expect and out.certificate == {"c": 1}
    assert tag in out.assumptions                         # always propagated
    if expect is Status.UNDECIDED:
        assert "blocked" in out.reasons[-1] and tag in out.reasons[-1]
    if allow:
        assert any(a.startswith("UNVERIFIED-INPUT-ACCEPTED") for a in out.assumptions)
