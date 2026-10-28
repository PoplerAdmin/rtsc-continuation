"""ORCH x CERT integration with the REAL rtsc.certify.classify (skipped if it is absent).

Only cheap SYNTHETIC cases (decided by R-UNIV / a few R-CELLS cells in a few seconds).
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

import rtsc.pipeline as P
from rtsc.sources import Candidate
from rtsc.sources.synthetic import gaussian_peak_spectrum
from rtsc.status import Status

import test_pipeline_stubs as stubs

pytestmark = pytest.mark.skipif(not stubs.real_classify_available(),
                                reason="rtsc.certify.classify not available")


def _cand(name, center, width, A):
    return Candidate(f"SYNTHETIC-{name}", {}, formula="(synthetic)", source="test", synthetic=True,
                     spectra=[gaussian_peak_spectrum(center, width, A, n=8, label=f"SYNTHETIC {name}")],
                     properties={"is_metal": True})


def test_real_classify_through_pipeline_and_recheck(tmp_path, capsys):
    cfg = P.normalize_config({"run": {"workdir": str(tmp_path / "w"), "resume": False},
                              "spectra": {"rel_envelope": 0.05},
                              "certify": {"T_target_K": 300.0, "mu_hi": 0.15, "mode": "decide",
                                          "omega_c_factor": 10.0}})
    cands = [_cand("weak", 500.0, 30.0, 200.0), _cand("mid", 6000.0, 300.0, 1400.0)]
    summ = P.run(cfg, cands)
    weak, mid = summ.verdicts["SYNTHETIC-weak"], summ.verdicts["SYNTHETIC-mid"]
    assert weak.status is Status.REJECTED_PROVEN and "R-UNIV" in weak.rule_ids
    assert mid.status is Status.REJECTED_PROVEN and "R-CELLS" in mid.rule_ids
    for v in (weak, mid):
        assert any("omega_c" in a for a in v.assumptions)          # F3
        assert "SYNTHETIC_DATA" in v.flags
        assert summ.origins[v.material_id] == P.CLASSIFY_ORIGIN    # S1
    from rtsc.cli import main
    work = Path(cfg["run"]["workdir"])
    assert main(["recheck", str(work / "SYNTHETIC-mid" / "verdict.json")]) == 0
    assert "RECHECK OK" in capsys.readouterr().out
    cert = work / "certificates" / "SYNTHETIC-weak.json"
    assert main(["recheck", str(cert)]) == 0
    # tampering is detected
    d = json.loads(cert.read_text())
    d["verdict"]["status"] = "CERTIFIED_CANDIDATE"
    bad = tmp_path / "bad.json"
    bad.write_text(json.dumps(d))
    assert main(["recheck", str(bad)]) == 3


def test_cli_certify_two_column(tmp_path, capsys):
    f = tmp_path / "synthetic_a2f.dat"
    f.write_text("# SYNTHETIC test data (omega meV, a2F)\n" +
                 "\n".join(f"{w} 0.01" for w in range(1, 40)) + "\n")
    from rtsc.cli import main
    out = tmp_path / "v.json"
    rc = main(["certify", str(f), "--units", "meV", "--mu-hi", "0.15", "--target", "300",
               "--mode", "decide", "--expected-A", "4.525762067404532", "--out", str(out)])
    assert rc == 0
    v = json.loads(out.read_text())
    assert v["status"] == "REJECTED_PROVEN"
    assert any("omega_c" in a for a in v["assumptions"])
