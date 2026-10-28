"""Adversarial audit of the discovery loop (DISCOVERY.md promises a-e).

Every test here encodes a promise of DISCOVERY.md / README "Discovery loop" / ARCH S1-S6 and is
expected to FAIL against the audited code: each failure is a demonstrated violation.  No package
code is modified; DFT is faked by monkeypatching rtsc.dfpt.EPhWorkflow (test-only).
"""
from __future__ import annotations

import copy
import csv
import json
import math
import sqlite3
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
sys.path.insert(0, str(HERE.parent))

from discovery_mocks import FixedGenerator, base_raw, make_dc  # noqa: E402

import rtsc.discovery.loop as L  # noqa: E402
import rtsc.dfpt as D  # noqa: E402
import rtsc.pipeline as P  # noqa: E402
from rtsc.discovery.fingerprint import fingerprint  # noqa: E402
from rtsc.discovery.goals import normalize_discovery  # noqa: E402
from rtsc.discovery.loop import DiscoveryLoop  # noqa: E402
from rtsc.sources.synthetic import gaussian_peak_spectrum  # noqa: E402
from rtsc.spectrum import Spectrum  # noqa: E402

ROOT = HERE.parents[1]
CUBE = {"lattice": [[3.0, 0, 0], [0, 3.0, 0], [0, 0, 3.0]], "species": ["La", "H"],
        "frac_coords": [[0, 0, 0], [0.5, 0.5, 0.5]]}
ONE = {"chemical_systems": [["H", "La"]], "requests_per_round": 1, "structures_per_request": 1}


def _spec(center, A):
    def f(fp, chemsys, p, scfg):
        return gaussian_peak_spectrum(center, center * 0.06, A, label="SYNTHETIC audit")
    return f


def _hex_layers(stacking: str, a: float = 2.5, h: float = 7.0) -> dict:
    """Triangular B layers, spacing h; 'AA' (eclipsed) or 'AB' (staggered) stacking."""
    lat = [[a, 0, 0], [-a / 2, a * math.sqrt(3) / 2, 0], [0, 0, 2 * h]]
    second = [0, 0, 0.5] if stacking == "AA" else [1 / 3, 2 / 3, 0.5]
    return {"lattice": lat, "species": ["B", "B"], "frac_coords": [[0, 0, 0], second]}


def _interlayer_nn(s: dict) -> float:
    """Shortest distance between atom 0 and atom 1 (incl. periodic images): a NON scale-free
    geometric invariant that differs between the two stackings."""
    import itertools
    import numpy as np
    Lm = np.asarray(s["lattice"])
    f0, f1 = np.asarray(s["frac_coords"][0]), np.asarray(s["frac_coords"][1])
    return min(float(np.linalg.norm((f1 - f0 + np.array(t)) @ Lm))
               for t in itertools.product((-1, 0, 1), repeat=3))


# --------------------------------------------------------------------------------------------
# F1  fingerprint collision merges two DIFFERENT materials (promise c: never loses a material;
#     the duplicate is "linked to its existing entry", i.e. reported as being that material)
# --------------------------------------------------------------------------------------------
def test_F1_fingerprint_merges_distinct_stacking_polytypes(tmp_path, monkeypatch):
    AA, AB = _hex_layers("AA"), _hex_layers("AB")
    assert abs(_interlayer_nn(AA) - _interlayer_nn(AB)) > 0.1      # genuinely different crystals
    monkeypatch.setattr(L, "synthetic_spectrum", _spec(800.0, 600.0))
    dc = make_dc(tmp_path / "w", goals={"G-RT": {"pressures_GPa": [0.0]}},
                 search={**ONE, "chemical_systems": [["B"]]}, stop={"max_rounds": 1})
    loop = DiscoveryLoop(dc, generators={"fixed": FixedGenerator([AA, AB])})
    loop.run()
    mats = loop.ledger.materials()
    dups = loop.ledger.q("SELECT * FROM duplicates")
    assert len(mats) == 2, (
        f"two distinct polytypes collapsed into {len(mats)} material(s); {len(dups)} linked as "
        f"duplicate via {[d['how'] for d in dups]} -> the AB structure is never computed and is "
        f"reported as {[d['material_id'] for d in dups]} (fingerprints equal: "
        f"{fingerprint(AA, 0.0) == fingerprint(AB, 0.0)})")


# --------------------------------------------------------------------------------------------
# F2  pressure_tol_GPa (1 GPa) bins P=0 with P=0.4/0.5: a structure first generated for G-RT at
#     0.4 GPa swallows the same structure requested for G-AMB at 0 GPa; G-AMB never sees it.
# --------------------------------------------------------------------------------------------
def test_F2_pressure_bin_merge_loses_ambient_material(tmp_path, monkeypatch):
    monkeypatch.setattr(L, "synthetic_spectrum", _spec(800.0, 600.0))
    dc = make_dc(tmp_path / "w", goals={"G-RT": {"pressures_GPa": [0.4]}},
                 search=ONE, stop={"max_rounds": 1})
    loop = DiscoveryLoop(dc, generators={"fixed": FixedGenerator([CUBE])})
    sid = loop.ledger.start_session("run", {})
    base = {"generator": "fixed", "chemical_system": ["H", "La"], "space_group": None, "n": 1,
            "seed": 1, "extra": {}}
    reqs = [{**base, "arm": "G-RT|H-La|any|P0.4|fixed", "pressure_GPa": 0.4, "goal_id": "G-RT"},
            {**base, "arm": "G-AMB|H-La|any|P0|fixed", "pressure_GPa": 0.0, "goal_id": "G-AMB"}]
    loop.ledger.start_round(1, sid, reqs)
    loop._run_round(1, None)
    ps = sorted(m["pressure_GPa"] for m in loop.ledger.materials())
    amb = dc.goal("G-AMB")
    assert any(amb.in_scope(p) for p in ps), (
        f"materials only at {ps} GPa; the P=0 request for G-AMB was linked as a duplicate of the "
        f"0.4 GPa material (out of G-AMB scope) -> never evaluated for the ambient goal. "
        f"fingerprint(P=0)==fingerprint(P=0.5): {fingerprint(CUBE, 0.0) == fingerprint(CUBE, 0.5)}")


# --------------------------------------------------------------------------------------------
# F3  [dfpt.settings] pressure_GPa overrides every material's pressure: spectra computed at
#     150 GPa are certified and counted for the AMBIENT goal (promise b).
# --------------------------------------------------------------------------------------------
class _FakeRes:
    def __init__(self, spectra):
        self.spectra = spectra

    def to_dict(self):
        return {}


def test_F3_dft_pressure_override_counts_for_ambient_goal(tmp_path, monkeypatch):
    seen: list[float] = []

    class FakeWF:
        def __init__(self, mid, structure, wd, settings, pseudos, launcher=None, masses=None):
            seen.append(float(settings.pressure_GPa))

        def run(self, *a, **k):
            return _FakeRes([Spectrum.from_dict(gaussian_peak_spectrum(800.0, 48.0, 700.0, label="fake DFPT"))])

    monkeypatch.setattr(D, "EPhWorkflow", FakeWF)
    raw = base_raw(tmp_path / "w", goals={"G-RT": {"pressures_GPa": [0.0]}}, search=ONE,
                   stop={"max_rounds": 1}, synthetic_spectra={"enabled": False})
    raw["dfpt"] = {"enabled": True, "pseudos": {"dir": str(tmp_path)},
                   "settings": {"pressure_GPa": 150.0}}
    try:
        dc = normalize_discovery(raw, ROOT)
        loop = DiscoveryLoop(dc, generators={"fixed": FixedGenerator([CUBE])})
        loop.run()
    except P.ConfigError:
        return                                     # acceptable fix: config rejected
    (m,) = loop.ledger.materials()
    cert = loop.certified_ids(dc.goal("G-AMB"))
    assert not (m["material_id"] in cert and any(p != m["pressure_GPa"] for p in seen)), (
        f"{m['material_id']} counted CERTIFIED for G-AMB (P<=1e-4 GPa, ledger pressure "
        f"{m['pressure_GPa']}) but its spectra were computed with DFT pressure(s) {seen} GPa")


# --------------------------------------------------------------------------------------------
# F4  config drift on `continue`: a certificate made under the OLD [certify] config keeps
#     counting under the NEW one although classify on the same spectra no longer certifies
#     (promise c: no reuse of a verdict whose config differs).
# --------------------------------------------------------------------------------------------
def test_F4_stale_certificate_counted_after_certify_config_change(tmp_path, monkeypatch):
    monkeypatch.setattr(L, "synthetic_spectrum", _spec(800.0, 600.0))
    wd = tmp_path / "w"
    dc1 = make_dc(wd, goals={"G-RT": {"pressures_GPa": [0.0]}}, search=ONE, stop={"max_rounds": 1})
    l1 = DiscoveryLoop(dc1, generators={"fixed": FixedGenerator([CUBE])})
    l1.run()
    (mid,) = l1.certified_ids(dc1.goal("G-AMB"))
    l1.ledger.close()
    raw2 = base_raw(wd, goals={"G-RT": {"pressures_GPa": [0.0]}}, search=ONE, stop={"max_rounds": 1})
    raw2["certify"]["mu_hi"] = 0.30                       # stricter model range
    dc2 = normalize_discovery(raw2, ROOT)
    l2 = DiscoveryLoop(dc2, generators={"fixed": FixedGenerator([CUBE])})
    res = l2.run("continue", max_rounds=0)
    g = dc2.goal("G-AMB")
    cur = l2.ledger.current_one(mid, g.id)
    rec = P.read_record(Path(cur["record_path"]).parent.parent, mid)
    cfg2 = dc2.goal_pipeline_config(g, tmp_path / "scratch")
    v_new = P.certify_spectra(mid, [Spectrum.from_dict(d) for d in rec["spectra"]], cfg2["certify"])
    assert not (mid in l2.certified_ids(g) and v_new.status.value != "CERTIFIED_CANDIDATE"), (
        f"under the continued config (mu_hi=0.30) classify gives {v_new.status.value} for {mid}, "
        f"yet it is counted certified for G-AMB (certified_total={res['certified_total']}); the "
        "stored config was never compared and the stale verdict is reused")


# --------------------------------------------------------------------------------------------
# F5  report overstates after a T* change: discovery_report.md / candidates_<goal>.csv list a
#     77 K certificate under "G-AMB: certified (Tc_model >= 150 K proved)"; certified_ids says 0.
# --------------------------------------------------------------------------------------------
def test_F5_report_lists_certificate_at_old_T_star(tmp_path, monkeypatch):
    monkeypatch.setattr(L, "synthetic_spectrum", _spec(800.0, 600.0))
    wd = tmp_path / "w"
    dc1 = make_dc(wd, goals={"G-RT": {"pressures_GPa": [0.0]}}, search=ONE, stop={"max_rounds": 1})
    l1 = DiscoveryLoop(dc1, generators={"fixed": FixedGenerator([CUBE])})
    l1.run()
    (mid,) = l1.certified_ids(dc1.goal("G-AMB"))
    l1.ledger.close()
    dc2 = make_dc(wd, goals={"G-RT": {"pressures_GPa": [0.0]}, "G-AMB": {"T_star_K": 150.0}},
                  search=ONE, stop={"max_rounds": 1})
    l2 = DiscoveryLoop(dc2, generators={"fixed": FixedGenerator([CUBE])})
    l2.run("continue", max_rounds=0)
    assert l2.certified_ids(dc2.goal("G-AMB")) == set()          # the counter is right ...
    rep = (wd / "discovery_report.md").read_text()
    sec = rep.split("## G-AMB: certified candidates")[1].split("## ")[0]
    rows = list(csv.DictReader(open(wd / "candidates_G-AMB.csv")))
    cert_csv = [r["material_id"] for r in rows if r["status"] == "CERTIFIED_CANDIDATE"]
    assert mid not in sec and not cert_csv, (                    # ... the report is not
        f"report section header:{sec.splitlines()[0]!r} lists {mid} (certified at T*=77 K only); "
        f"CSV CERTIFIED rows: {cert_csv}")


# --------------------------------------------------------------------------------------------
# F6  a goal enabled on `continue` is never evaluated for materials already computed: status
#     stays "DEFERRED ... awaiting selection" forever and the loop declares itself exhausted.
# --------------------------------------------------------------------------------------------
def test_F6_goal_added_on_continue_never_evaluated(tmp_path, monkeypatch):
    monkeypatch.setattr(L, "synthetic_spectrum", _spec(800.0, 600.0))
    wd = tmp_path / "w"
    dc1 = make_dc(wd, goals={"G-RT": {"pressures_GPa": [0.0]}, "G-AMB": {"enabled": False}},
                  search=ONE, stop={"max_rounds": 1})
    l1 = DiscoveryLoop(dc1, generators={"fixed": FixedGenerator([CUBE])})
    l1.run()
    (mid,) = [m["material_id"] for m in l1.ledger.materials()]
    l1.ledger.close()
    dc2 = make_dc(wd, goals={"G-RT": {"pressures_GPa": [0.0]}}, search=ONE, stop={"max_rounds": 3})
    l2 = DiscoveryLoop(dc2, generators={"fixed": FixedGenerator([CUBE])})
    res = l2.run("continue")
    cur = l2.ledger.current_one(mid, "G-AMB")
    assert cur["stage"] == "pipeline", (
        f"{mid} (P=0, in scope) is {cur['status']}/{cur['stage']} for G-AMB: {cur['reasons']}; "
        f"loop stopped with {res['stop_reason']!r}")


# --------------------------------------------------------------------------------------------
# F7  SQL-level bypasses of "append-only" and of the S1/S2/S3 guard (promise a, e)
# --------------------------------------------------------------------------------------------
def test_F7a_insert_or_replace_rewrites_verdict_history(tmp_path, monkeypatch):
    monkeypatch.setattr(L, "synthetic_spectrum", _spec(800.0, 300.0))    # REJECTED at 77 K
    dc = make_dc(tmp_path / "w", goals={"G-RT": {"pressures_GPa": [0.0]}}, search=ONE,
                 stop={"max_rounds": 1})
    loop = DiscoveryLoop(dc, generators={"fixed": FixedGenerator([CUBE])})
    loop.run()
    led = loop.ledger
    v = [h for h in led.history() if h["status"] == "REJECTED_PROVEN"][0]
    row = dict(led.q("SELECT * FROM verdicts WHERE id=?", (v["id"],))[0])
    row.update(status="CERTIFIED_CANDIDATE", T_low_K=999.0)
    cols = ",".join(row)
    with pytest.raises(sqlite3.DatabaseError):
        led.con.execute(f"INSERT OR REPLACE INTO verdicts({cols}) VALUES ({','.join('?' * len(row))})",
                        tuple(row.values()))


def test_F7b_raw_insert_bypasses_stage_guard(tmp_path):
    from rtsc.discovery.ledger import Ledger
    led = Ledger(tmp_path / "l.sqlite")
    led.add_material(material_id="X1", seq=1, fingerprint="fp", coarse_key="k", formula="H",
                     chemsys=("H",), pressure_GPa=0.0, structure={}, generator="g", request_id="r",
                     request={}, arm="G-AMB|x", round_=1, provenance={})
    with pytest.raises(sqlite3.DatabaseError):
        led.con.execute(
            "INSERT INTO verdicts(material_id, goal_id, round, stage, origin, status, T_star_K, "
            "pressure_GPa, rule_ids, assumptions, flags, reasons, created_utc) "
            "VALUES ('X1','G-AMB',1,'prescreen','discovery.prescreen:evil','CERTIFIED_CANDIDATE',77,0,"
            "'[]','[]','[]','[]','now')")
