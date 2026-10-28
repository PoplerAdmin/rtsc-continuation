"""Adversarial audit of the orchestration layer (pipeline / report / cli / parsers).

Every test here asserts the behaviour the soundness contract (ARCH section 3) requires.  At the
time of the audit each test FAILS, demonstrating a concrete way the contract can be violated.
They use the real rtsc.certify.classify only on spectra that are decided instantly (R-UNIV), so
they are fast; the orchestration bugs do not depend on the engine.
"""
from __future__ import annotations

import csv
import json
import math
import os
import sys
import time
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # tests/ (for the stub helpers)

import rtsc.certify as C
import rtsc.pipeline as P
from rtsc.certify import CertifyParams, classify, spectrum_digest
from rtsc.sources import Candidate
from rtsc.sources.synthetic import gaussian_peak_spectrum
from rtsc.spectrum import Spectrum
from rtsc.status import Status, Verdict

A0 = 4.33148145
WEAK = dict(center_K=500.0, width_K=30.0, A_K=200.0)      # A_hi ~ 210 K << a0*300 K: R-UNIV
STRONG = dict(center_K=2000.0, width_K=120.0, A_K=3000.0)  # A_hi > a0*300: never R-UNIV;
#                                                           certified (T_low ~ 413 K) by the engine


def _cfg(tmp_path, **over):
    raw = {"run": {"workdir": str(tmp_path / "work"), "resume": True},
           "dedupe": {"method": "none"},
           "spectra": {"rel_envelope": 0.05},
           "certify": {"mode": "decide"}}
    raw.update(over)
    return P.normalize_config(raw)


def _cand(mid, peak, **kw):
    return Candidate(mid, {}, formula="(synthetic)", synthetic=True, source="audit",
                     properties={"is_metal": True},
                     spectra=[gaussian_peak_spectrum(peak["center_K"], peak["width_K"], peak["A_K"],
                                                     label=f"SYNTHETIC {mid}")], **kw)


def _write_two_column(path: Path, peak) -> None:
    sp = gaussian_peak_spectrum(peak["center_K"], peak["width_K"], peak["A_K"], n=60)
    xs = np.asarray(sp["omega_K"])
    h = xs[1] - xs[0]
    a2f = np.asarray(sp["weight_K"]) / h          # density per kelvin
    path.write_text("# omega[K] a2F\n" + "\n".join(f"{x:.10g} {y:.10g}" for x, y in zip(xs, a2f)))


# ------------------------------------------------------------------------------------------
# FINDING 1 (critical): resume reuses a finished verdict although the INPUT changed.
# config_hash() covers only the config; the record's candidate / spectra are never compared
# with the current input.  A REJECTED_PROVEN about the old alpha^2F is reported for the new one.
# ------------------------------------------------------------------------------------------
def test_resume_reuses_stale_verdict_after_attached_file_changed(tmp_path):
    f = tmp_path / "a2F.dat"
    _write_two_column(f, WEAK)
    cand = lambda: Candidate("M1", {}, formula="(synthetic)", synthetic=True,
                             properties={"is_metal": True},
                             spectra_files=[{"path": str(f), "units": "K"}])
    # (coordinator) two-column files are tagged A-FILE-COMPLETE(unverified) and rejections on
    # them are downgraded by policy; opt in so this test still exercises stale-resume.
    cfg = _cfg(tmp_path, certify={"mode": "decide", "allow_unverified_rejection": True})
    s1 = P.run(cfg, [cand()])
    assert s1.verdicts["M1"].status is Status.REJECTED_PROVEN
    _write_two_column(f, STRONG)                     # new DFPT result for the same material
    s2 = P.run(cfg, [cand()])
    current = P._with_envelope(P.attached_spectra(cand(), cfg["spectra"]), 0.05)
    used = [a["spectrum_sha256"] for a in s2.verdicts["M1"].certificate["alternatives"]]
    # contract: the verdict must be about the CURRENT input spectrum
    assert used == [spectrum_digest(s) for s in current], (
        f"stale verdict {s2.verdicts['M1'].status.value} reused for a changed alpha^2F file")


def test_resume_reuses_stale_verdict_for_a_different_candidate_with_same_id(tmp_path):
    """Same id, different material (e.g. 'H3S' / 'H3S~2' swap when a file is added to a local
    directory, or a seed list is edited): the old REJECTED_PROVEN is silently reused."""
    cfg = _cfg(tmp_path)
    assert P.run(cfg, [_cand("X", WEAK)]).verdicts["X"].status is Status.REJECTED_PROVEN
    s2 = P.run(cfg, [_cand("X", STRONG)])
    # STRONG has A_hi > a0*T*, and the engine certifies it: REJECTED_PROVEN is false for it
    assert s2.verdicts["X"].status is not Status.REJECTED_PROVEN


# ------------------------------------------------------------------------------------------
# FINDING 2 (critical): `collect` merges records regardless of their config hash and labels
# them with the CURRENT config's T*, mu_hi in candidates.csv / summary.md.
# ------------------------------------------------------------------------------------------
def test_collect_relabels_verdicts_made_under_another_config(tmp_path):
    cfg300 = _cfg(tmp_path)
    s = P.run(cfg300, [_cand("W", WEAK)])
    v = s.verdicts["W"]
    assert v.status is Status.REJECTED_PROVEN and v.T_up_K > 40.0     # T_up ~ 48 K
    cfg30 = _cfg(tmp_path, certify={"mode": "decide", "T_target_K": 30.0})
    assert P.config_hash(cfg30) != P.config_hash(cfg300)
    try:
        P.collect(cfg30)                               # e.g. after editing config.toml
    except (P.ConfigError, P.AccountingError):
        return                                         # acceptable fix: refuse to merge
    rows = list(csv.DictReader(open(tmp_path / "work" / "candidates.csv")))
    row = rows[0]
    # The report must not claim "Tc_model < 30 K proved" when only Tc_model < 300 K was proved
    # (and the certified upper bound is ~48 K > 30 K).
    assert not (row["status"] == "REJECTED_PROVEN" and float(row["T_target_K"]) == 30.0), row


# ------------------------------------------------------------------------------------------
# FINDING 3 (high): two different material ids map to the same <workdir>/<safe_id>/ directory
# (record.json, verdict.json, certificates/<id>.json, and the DFPT run directory).
# After an interrupted run, one material is classified with the OTHER one's spectra.
# ------------------------------------------------------------------------------------------
def test_safe_id_collision_distinct_record_paths(tmp_path):
    for a, b in [("X/1", "X_1"), ("H3S#0", "H3S_0"), ("H3S~2", "H3S_2"), ("LaH10 Fm-3m", "LaH10_Fm-3m")]:
        assert P.record_path(tmp_path, a) != P.record_path(tmp_path, b), (a, b)


def test_safe_id_collision_gives_wrong_rejection_after_interrupt(tmp_path, monkeypatch):
    cfg = _cfg(tmp_path)
    strong, weak = "LaH10 Fm-3m", "LaH10_Fm-3m"         # different materials, same safe_id
    real = C.classify

    def interrupted(spectra, params, material_id="", **kw):
        if material_id == weak:
            raise KeyboardInterrupt("node preempted while certifying")   # job killed
        return Verdict(material_id, Status.UNDECIDED, reasons=["(fast placeholder)"])
    monkeypatch.setattr(C, "classify", interrupted)
    with pytest.raises(KeyboardInterrupt):
        P.run(cfg, [_cand(strong, STRONG), _cand(weak, WEAK)])
    monkeypatch.setattr(C, "classify", real)
    s = P.run(cfg, [_cand(strong, STRONG), _cand(weak, WEAK)])      # resume after preemption
    # STRONG (A_hi = 3150 K > a0*300 K, certified by the engine) must never be REJECTED_PROVEN
    assert s.verdicts[strong].status is not Status.REJECTED_PROVEN, (
        "classified with the other material's spectra: " + "; ".join(s.verdicts[strong].reasons))


# ------------------------------------------------------------------------------------------
# FINDING 4 (high): [certify] T_target_K and mu_hi are not validated.
# ------------------------------------------------------------------------------------------
@pytest.mark.parametrize("bad", [{"T_target_K": -5.0}, {"T_target_K": 0.0},
                                 {"T_target_K": float("nan")}, {"T_target_K": float("inf")},
                                 {"mu_hi": -0.1}, {"mu_hi": float("nan")}, {"mu_hi": True},
                                 {"T_target_K": "300"}])
def test_invalid_certify_parameters_rejected(bad):
    with pytest.raises(P.ConfigError):
        P.normalize_config({"certify": bad})


def test_negative_target_gives_vacuous_certified():
    sp = [Spectrum.from_dict(gaussian_peak_spectrum(**{"center_K": 500.0, "width_K": 30.0,
                                                        "A_K": 200.0})).with_envelope(0.05)]
    # (coordinator, after fix): a negative target must be refused both by normalize_config and
    # by certify_spectra when normalize_config is bypassed; never a vacuous CERTIFIED verdict.
    with pytest.raises(P.ConfigError):
        P.normalize_config({"certify": {"T_target_K": -5.0, "mode": "decide"}})
    with pytest.raises(P.ConfigError):
        P.certify_spectra("W", sp, {"T_target_K": -5.0, "mode": "decide"})


# ------------------------------------------------------------------------------------------
# FINDING 5 (high): a truncated QE a2F.dos (no footer; e.g. copied while matdyn was writing, or
# disk full) is accepted silently.  The missing high-frequency bins make A too SMALL, and the
# material is REJECTED_PROVEN with no flag.  Attached QE files also get no footer-lambda check.
# ------------------------------------------------------------------------------------------
def _qe_a2f_dos(path: Path, truncate_at_Ry: float | None = None) -> None:
    d = 1.0e-4                                         # Delta in Ry
    E = d * (np.arange(200) + 0.5)
    a2f = 7.6 * np.exp(-0.5 * ((E - 0.012) / 0.001) ** 2)
    lam = float(np.sum(2.0 * a2f / E * d))
    lines = ["", " # Eliashberg function a2F (per both spin)", " #  frequencies in Rydberg  ",
             " # DOS normalized to E in Rydberg: a2F_total, a2F(mode) ", ""]
    for e, a in zip(E, a2f):
        if truncate_at_Ry is not None and e > truncate_at_Ry:
            break
        lines.append(f"   {e:16.6E}{a:16.6E}{a / 2:16.6E}{a / 2:16.6E}")
    else:
        lines.append(f"  lambda =   {lam:.8f}   Delta =   {d:.8E}")
    path.write_text("\n".join(lines) + "\n")


def test_truncated_qe_a2f_dos_is_rejected_not_classified(tmp_path):
    full, trunc = tmp_path / "a2F.dos1", tmp_path / "trunc" / "a2F.dos1"
    trunc.parent.mkdir()
    _qe_a2f_dos(full)
    _qe_a2f_dos(trunc, truncate_at_Ry=0.0075)
    scfg = P.DEFAULTS["spectra"]
    A_full = P._read_spectrum_file({"path": str(full), "format": "qe_a2f_dos"}, scfg)[0].A()
    assert A_full * 1.05 > A0 * 300.0                 # the complete spectrum: no R-UNIV
    cand = Candidate("Q", {}, properties={"is_metal": True},
                     spectra_files=[{"path": str(trunc), "format": "qe_a2f_dos"}])
    s = P.run(_cfg(tmp_path), [cand])
    v = s.verdicts["Q"]
    # contract: an incomplete input must not yield a proof (ERROR, or at least never REJECTED)
    assert v.status is not Status.REJECTED_PROVEN, (v.status, v.flags)


# ------------------------------------------------------------------------------------------
# FINDING 6 (medium): one malformed property aborts the whole run (no reports at all) instead
# of an ERROR for that input (S4: "unreadable inputs end with an explicit status").
# ------------------------------------------------------------------------------------------
def test_malformed_property_becomes_error_not_abort(tmp_path):
    bad = _cand("bad", WEAK)
    bad.properties["band_gap"] = "n/a"
    s = P.run(_cfg(tmp_path), [bad, _cand("ok", WEAK)])
    assert s.verdicts["bad"].status is Status.ERROR
    assert s.verdicts["ok"].status is Status.REJECTED_PROVEN


# ------------------------------------------------------------------------------------------
# FINDING 7 (medium): a single crashed worker (OOM-kill, segfault in QE wrapper) breaks the
# process pool; EVERY pending material is persisted as ERROR "worker failed" and, being
# 'finished', is never retried on resume.
# ------------------------------------------------------------------------------------------
_ORIG_COMPUTE = P.compute_material


def _crashy_compute(cand_dict, rec, cfg, workdir):
    if cand_dict["material_id"] == "boom":
        os._exit(1)                                   # e.g. the OOM killer
    time.sleep(1.5)
    return _ORIG_COMPUTE(cand_dict, rec, cfg, workdir)


def test_broken_pool_errors_are_retried_on_resume(tmp_path, monkeypatch):
    import multiprocessing as mp
    if mp.get_start_method() != "fork":
        pytest.skip("needs fork")
    cfg = _cfg(tmp_path, run={"workdir": str(tmp_path / "work"), "workers": 2, "resume": True})
    cands = lambda: [_cand("boom", STRONG), _cand("a", WEAK), _cand("b", WEAK)]
    monkeypatch.setattr(P, "compute_material", _crashy_compute)
    s1 = P.run(cfg, cands())
    assert s1.verdicts["boom"].status is Status.ERROR
    monkeypatch.setattr(P, "compute_material", _ORIG_COMPUTE)
    s2 = P.run(cfg, cands())                           # resume: a and b were never computed
    assert s2.verdicts["a"].status is not Status.ERROR
    assert s2.verdicts["b"].status is not Status.ERROR


# ------------------------------------------------------------------------------------------
# FINDING 8 (medium): DFPT consistency-check failures that indicate a PARSE/UNITS error
# (lambda_vs_a2F_footer) only add a flag; the material can still be REJECTED_PROVEN.
# ------------------------------------------------------------------------------------------
class _FakeResult:
    def __init__(self, spectra, flags, warnings=()):
        self.spectra, self.flags, self.warnings = spectra, flags, list(warnings)

    def to_dict(self):
        return {"spectra": [s.to_dict() for s in self.spectra], "flags": self.flags,
                "warnings": self.warnings}


def test_dfpt_units_check_failure_blocks_proof(tmp_path, monkeypatch):
    dfpt = pytest.importorskip("rtsc.dfpt")

    class FakeWF:
        def __init__(self, *a, **k):
            pass

        def run(self):
            sp = Spectrum.from_dict(gaussian_peak_spectrum(500.0, 30.0, 200.0, label="sigma1"))
            return _FakeResult([sp], ["DFPT_CHECK_FAIL:lambda_vs_a2F_footer:sigma1"],
                               ["only smearing indices [1] of 1..10 kept (explicit config)"])
    monkeypatch.setattr(dfpt, "EPhWorkflow", FakeWF)
    c = Candidate("D", "NbH", lattice=[[3, 0, 0], [0, 3, 0], [0, 0, 3]], species=["Nb", "H"],
                  frac_coords=[[0, 0, 0], [.5, .5, .5]], properties={"is_metal": True})
    cfg = _cfg(tmp_path, dfpt={"enabled": True, "pseudos": {"dir": str(tmp_path)}})
    v = P.run(cfg, [c], write_reports=False).verdicts["D"]
    assert v.status is not Status.REJECTED_PROVEN, v.flags


# ------------------------------------------------------------------------------------------
# FINDING 9 (medium): a small omega_c_factor (any value > 0 is accepted) switches the Coulomb
# term off (N_c(T) = 0 for T >= omega_c/pi), so CERTIFIED no longer depends on mu*, while the
# verdict and summary.md still advertise "U-MU[0,0.15]" / "mu* in [0, 0.15]".
# ------------------------------------------------------------------------------------------
def test_tiny_omega_c_factor_rejected_or_flagged():
    from rtsc.sources.synthetic import demo_candidates
    try:
        cfg = P.normalize_config({"certify": {"omega_c_factor": 0.01, "mode": "decide"}})
    except P.ConfigError:
        return                                          # acceptable fix
    c = {x.material_id: x for x in demo_candidates()}["SYNTHETIC-borderline"]
    sp = [Spectrum.from_dict(d).with_envelope(0.05) for d in c.spectra]
    v = P.certify_spectra("b", sp, cfg["certify"])   # UNDECIDED with the default factor 10
    assert not (v.status is Status.CERTIFIED_CANDIDATE
                and any(a.startswith("U-MU[0,0.15]") for a in v.assumptions)
                and not any("COULOMB" in f for f in v.flags)), (v.status, v.assumptions)
