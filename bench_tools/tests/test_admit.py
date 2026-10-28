"""P3 admission / G4 anchor / P5 envelope (ephproto.admit) on SYNTHETIC spectra, incl. §V sabotage."""
import math
import re

import numpy as np
import pytest

from ephproto import admit as AD
from ephproto.admit import (AdmitError, admit_k, admit_q, admit_sequence, decide, envelope,
                            envelope_assumption, envelope_provenance, error_record, observables,
                            richardson_item, tetra_anchor_gate)
from ephproto.tags import EPH_ERROR_TAG_RE, EPH_PIN_VIOLATION, unconverged_eph_tag

OMAX = 40.0                                    # meV, coarse-DFPT max omega (windows of 2 meV)
BASE_OM = np.array([3.3, 7.1, 11.6, 15.2, 19.4, 23.1, 27.1, 31.5, 35.3])
BASE_W = np.array([0.30, 0.80, 1.20, 1.50, 1.10, 0.90, 1.40, 0.70, 0.40])
PATTERN = np.array([1.0, -0.5, 0.8, -1.0, 0.6, 0.3, -0.7, 1.0, -0.4])


def spec(n, *, c=0.2, p=2.0, cw=0.01, extra=None):
    """Grid n of a synthetic sequence converging with order p (weights and positions)."""
    r = 2.0 ** (-p * n)
    om = BASE_OM * (1.0 + cw * r)
    w = BASE_W * (1.0 + c * r * PATTERN)
    if extra is not None:
        w = w * (1.0 + extra)
    return om, w


def seq(tag="k", sigma=10.0, *, efs=(11.10, 11.1001, 11.10012), grids=None, **kw):
    out = []
    for i, n in enumerate((1, 2, 3)):
        om, w = grids[i] if grids is not None else spec(n, **kw)
        out.append({"omega": om, "weight": w, "hash": f"{tag}-{sigma:g}-G{i + 1}",
                    "ef_eV": None if efs is None else efs[i]})
    return out


def run_seq(direction="k", sigma=10.0, **kw):
    return admit_sequence(seq(direction, sigma, **kw), direction=direction, sigma_meV=sigma,
                          omega_max=OMAX)


def k_scan(bad=(), make_bad=None):
    make_bad = make_bad or (lambda s: run_seq("k", s, p=-1.0))
    return {s: (make_bad(s) if s in bad else run_seq("k", s)) for s in AD.SIGMA_MEV}


ANCHOR_SEQ = [{"lambda": 1.0 + a, "A": 10.0 + 10 * a, "omega_log": 20.0 * (1 + a)}
              for a in (0.004, 0.001, 0.00025)]


def protocol_values(kres, smin):
    ob = kres[smin]["observables"][-1]
    return {x: ob[x] for x in AD.SCALARS}


def full(kres, qres=None, *, anchor_shift=0.0, anchor_k=ANCHOR_SEQ, anchor_q=ANCHOR_SEQ, g6=0.03,
         extended=False):
    kg = admit_k(kres, extended=extended)
    if kg["status"] != "PASS":
        return decide(kg, None, None, None), kg
    adm = kg["admitted_sigma_meV"]
    qres = qres if qres is not None else {s: run_seq("q", s, efs=None) for s in adm}
    qg = admit_q(qres, adm)
    pv = protocol_values(kres, kg["sigma_min_adm_meV"])
    anchor = {x: v * (1 + anchor_shift) for x, v in pv.items()}
    tg = tetra_anchor_gate(anchor, pv, anchor_k_values=anchor_k, anchor_q_values=anchor_q)
    env = envelope(kres, qres, adm, g6_rel_dgamma=g6, g4_gaps=tg["gaps"]) \
        if qg["status"] == "PASS" else None
    return decide(kg, qg, tg, env), kg


# ----------------------------------------------------------------------------- observables
def test_observables_definitions():
    ob = observables([2.0, 4.0], [1.0, 2.0], 8.0, n_windows=4)
    assert ob["A"] == pytest.approx(3.0)
    assert ob["lambda"] == pytest.approx(2 * (0.5 + 0.5))
    assert ob["omega_log"] == pytest.approx(math.exp((0.5 * math.log(2) + 0.5 * math.log(4)) / 1.0))
    assert ob["W"] == pytest.approx([1.0, 2.0, 0.0, 0.0])   # right-closed windows (0,2], (2,4]
    with pytest.raises(AdmitError):
        observables([1.0], [1.0, 2.0], 8.0)


def test_observables_flags_atoms_outside_windows():
    ob = observables([2.0, 9.0, -1.0], [1.0, 1.0, 1.0], 8.0)
    assert ob["outside"]["n_above_omega_max"] == 1
    assert ob["outside"]["n_nonpositive_or_negative_weight"] == 1


# ----------------------------------------------------------------------------- one item
def test_item_exact_richardson_p2():
    a = 0.64
    it = richardson_item([1 + a / 4, 1 + a / 16, 1 + a / 64], 0.02)
    assert it["branch"] == "rate" and it["p"] == pytest.approx(2.0)
    assert it["R"] == pytest.approx((a * 3 / 64) / 3)
    assert it["pass"]


def test_item_p_half_factor_and_cap():
    d1, d2 = 1e-2, 1e-2 / 2 ** 0.50001
    it = richardson_item([1.0, 1.0 + d1, 1.0 + d1 + d2], 0.02)
    assert it["p"] == pytest.approx(0.5, abs=1e-4) and it["pass_rate"]
    assert it["R"] == pytest.approx(d2 / (math.sqrt(2) - 1), rel=1e-4)       # factor 2.41
    it = richardson_item([1.0, 1.1, 1.1 + 0.1 / 64], 0.02)                    # p = 6 -> capped at 4
    assert it["p"] == pytest.approx(6.0) and it["R"] == pytest.approx(0.1 / 64 / 15)


def test_item_converged_branch_R_is_max_step_never_zero():
    it = richardson_item([1.0, 1.001, 1.0005], 0.02)         # tiny, opposite signs
    assert it["branch"] == "converged" and it["pass"]
    assert it["R"] == pytest.approx(0.001)
    it = richardson_item([2.0, 2.0, 2.0], 0.02)
    assert it["pass"] and it["R"] == pytest.approx(5e-5 * 2.0) and it["R"] > 0


@pytest.mark.parametrize("x,why", [
    ([1.0, 1.1, 1.05], "oscillating"),
    ([1.0, 1.01, 1.05], "diverging"),
    ([1.0, 1.1, 1.18], "too slow"),                   # p = 0.32
    ([1.0, 1.1, 1.1], "p undefined"),                 # d_2 = 0 outside converged branch
    ([1.0, 0.0, 1.0], "X_ref = 0"),
    ([1.0, float("nan"), 1.0], "non-finite"),
])
def test_item_fails(x, why):
    it = richardson_item(x, 0.02)
    assert not it["pass"] and why in it["why"], it


def test_item_remainder_limit():
    it = richardson_item([1.0, 1.08, 1.12], 0.02)            # p = 1, R = 0.04
    assert it["pass_rate"] and not it["pass_remainder"] and not it["pass"]


# ----------------------------------------------------------------------------- one sequence
def test_converging_sequence_admitted():
    r = run_seq("k", 10.0)
    assert r["admitted"], r["failures"]
    assert set(r["tracked_windows"]) == {1, 3, 5, 7, 9, 11, 13, 15, 17}
    assert all(r["items"][f"W{b}"]["pass"] for b in r["tracked_windows"])
    assert 0 < r["max_R_rel"] < 0.02


def test_window_tracked_by_lambda_share_only():
    om = np.array([0.5, 20.0])                 # low window: < 1 % of A but > 1 % of lambda
    w = np.array([0.004, 1.0])
    grids = [(om, w * (1 + 0.1 * 4.0 ** -n * np.array([1, 0]))) for n in (1, 2, 3)]
    r = admit_sequence(seq(grids=grids), direction="k", sigma_meV=10.0, omega_max=OMAX)
    assert 0 in r["tracked_windows"]


def test_oscillating_single_window_fails():
    extra = [np.zeros(9), np.zeros(9), np.zeros(9)]
    extra[1][4] = 0.05                          # window of 19.4 meV goes up then back down
    grids = [spec(n, extra=extra[n - 1]) for n in (1, 2, 3)]
    r = admit_sequence(seq(grids=grids), direction="k", sigma_meV=10.0, omega_max=OMAX)
    assert not r["admitted"] and any(f.startswith("W9:") and "oscillating" in f for f in r["failures"])


def test_atom_above_omega_max_fails():
    grids = [spec(n) for n in (1, 2, 3)]
    grids[2] = (np.append(grids[2][0], 41.0), np.append(grids[2][1], 1e-3))
    r = admit_sequence(seq(grids=grids), direction="k", sigma_meV=10.0, omega_max=OMAX)
    assert not r["admitted"] and any("above omega_max" in f for f in r["failures"])


def test_ef_shift_gate_k_only():
    r = run_seq("k", 10.0, efs=(11.1, 11.1, 11.1026))        # 2.6 meV > 10/4
    assert not r["admitted"] and any("E_F shift" in f for f in r["failures"])
    r = run_seq("k", 10.0, efs=None)
    assert not r["admitted"] and any("E_F of G2/G3" in f for f in r["failures"])
    assert run_seq("q", 10.0, efs=None)["admitted"]


def test_spectrum_like_input():
    class Sp:
        def __init__(self, n):
            self.omega_K, self.weight_K = spec(n)
            self.label = f"g{n}"
            self.provenance = {"sources": {"epw_in": {"sha256": f"h{n}"}},
                               "stdout_checks": {"ef_eV": 11.1}}
    r = admit_sequence([Sp(1), Sp(2), Sp(3)], direction="k", sigma_meV=10.0, omega_max=OMAX)
    assert r["admitted"] and r["labels"] == ["g1", "g2", "g3"]


# ----------------------------------------------------------------------------- whole protocol
def test_full_protocol_admitted_and_envelope():
    kres = k_scan()
    out, kg = full(kres)
    assert out["status"] == "ADMITTED", out
    assert out["admitted_sigma_meV"] == list(AD.SIGMA_MEV) and out["sigma_min_adm_meV"] == 5.0
    assert out["e"] == pytest.approx(0.10)                   # floor dominates
    assert out["assumption"] == "E-ENV(0.1, estimated by eph_protocol v3.1)"


def test_envelope_terms():
    kres = k_scan()
    adm = list(AD.SIGMA_MEV)
    qres = {s: run_seq("q", s, efs=None) for s in adm}
    env = envelope(kres, qres, adm, g6_rel_dgamma=0.0, g4_gaps={"lambda": 0.9})
    assert env["status"] == "PASS" and env["estimated"]
    assert env["e"] == pytest.approx(max(0.10, 3 * env["max_R_rel"]))
    assert env["g4_gaps_reported"] == {"lambda": 0.9}          # recorded, not used (A6)
    qres[20.0] = run_seq("q", 20.0, efs=None, c=2.0)            # large window remainder
    assert qres[20.0]["admitted"], qres[20.0]["failures"]
    env = envelope(kres, qres, adm, g6_rel_dgamma=0.0)
    assert env["max_R_rel_term"].startswith("q:20:") and env["e"] == pytest.approx(3 * env["max_R_rel"])
    assert env["e"] > 0.10
    env = envelope(kres, qres, adm, g6_rel_dgamma=0.42)
    assert env["e"] == pytest.approx(0.42)                     # G6 enters without factor 3
    with pytest.raises(AdmitError):
        envelope(kres, qres, adm, g6_rel_dgamma=float("nan"))


def test_envelope_provenance():
    p = envelope_provenance(0.12, "K8xQ8", config_floor=0.15)
    assert p["rel_envelope"] == 0.15 and p["envelope_status"] == "estimated"
    assert p["discretisation"] == "eph_protocol_v3.1:K8xQ8"
    assert p["assumptions"] == [envelope_assumption(0.12)] == ["E-ENV(0.12, estimated by eph_protocol v3.1)"]


def test_k4_extension_and_sequence():
    assert AD.k_multiples() == (2, 4, 8) and AD.k_multiples(True) == (4, 8, 16)
    kres = k_scan(bad=AD.SIGMA_MEV)
    assert admit_k(kres)["k4_extension_allowed"] is True
    assert admit_k(kres, extended=True)["k4_extension_allowed"] is False


# ----------------------------------------------------------------------------- §V sabotage: ERROR, never admission
SABOTAGE = {
    # name: (kwargs for full(), expected primary reason)
    "k_nonconvergent_all_sigma": (dict(kres=lambda: k_scan(bad=AD.SIGMA_MEV)), "UNCONVERGED_EPH(k)"),
    "k_oscillating_all_sigma": (dict(kres=lambda: k_scan(bad=AD.SIGMA_MEV, make_bad=lambda s: run_seq(
        "k", s, grids=[spec(1), spec(3), spec(2)]))), "UNCONVERGED_EPH(k)"),
    "k_wrong_order_all_sigma": (dict(kres=lambda: k_scan(bad=AD.SIGMA_MEV, make_bad=lambda s: run_seq(
        "k", s, p=0.3, c=0.02))), "UNCONVERGED_EPH(k)"),
    "k_sigma_gap": (dict(kres=lambda: k_scan(bad=(20.0,))), "UNCONVERGED_EPH(k)"),
    "sigma_min_adm_too_large": (dict(kres=lambda: k_scan(bad=(5.0, 10.0, 20.0))),
                                "UNCONVERGED_EPH(sigma)"),
    "E_F_shift": (dict(kres=lambda: k_scan(bad=AD.SIGMA_MEV, make_bad=lambda s: run_seq(
        "k", s, efs=(11.1, 11.1, 11.4))), ), "UNCONVERGED_EPH(k)"),
    "K1_eq_K2_eq_K3": (dict(kres=lambda: k_scan(bad=(10.0,), make_bad=lambda s: admit_sequence(
        [{"omega": BASE_OM, "weight": BASE_W, "hash": "same", "ef_eV": 11.1}] * 3,
        direction="k", sigma_meV=s, omega_max=OMAX))), EPH_PIN_VIOLATION),
    "identical_data_distinct_hashes": (dict(kres=lambda: k_scan(bad=(10.0,), make_bad=lambda s: run_seq(
        "k", s, grids=[spec(2), spec(2), spec(3)]))), EPH_PIN_VIOLATION),
    "q_unconverged": (dict(qres=lambda: {s: run_seq("q", s, efs=None, p=-0.5) for s in AD.SIGMA_MEV}),
                      "UNCONVERGED_EPH(q)"),
    "q_oscillating_one_sigma": (dict(qres=lambda: {s: (run_seq("q", s, efs=None,
                                                               grids=[spec(1), spec(3), spec(2)])
                                                       if s == 40.0 else run_seq("q", s, efs=None))
                                                   for s in AD.SIGMA_MEV}), "UNCONVERGED_EPH(q)"),
    "q_missing_sigma": (dict(qres=lambda: {s: run_seq("q", s, efs=None) for s in AD.SIGMA_MEV[:-1]}),
                        "UNCONVERGED_EPH(q)"),
    "anchor_disagreement": (dict(anchor_shift=0.08), "UNCONVERGED_EPH(tetra)"),
    "anchor_disagreement_low": (dict(anchor_shift=-0.06), "UNCONVERGED_EPH(tetra)"),
    "anchor_unconverged": (dict(anchor_k=[{"lambda": v, "A": 10.0, "omega_log": 20.0}
                                          for v in (1.0, 1.1, 1.05)]), "UNCONVERGED_EPH(tetra)"),
    "anchor_convergence_missing": (dict(anchor_q=None), "UNCONVERGED_EPH(tetra)"),
    "envelope_too_large": (dict(g6=0.6), "UNCONVERGED_EPH(envelope)"),
}


@pytest.mark.parametrize("name", sorted(SABOTAGE))
def test_sabotage_is_error_never_admission(name):
    kw, expect = SABOTAGE[name]
    kw = dict(kw)
    kres = kw.pop("kres", k_scan)()
    if "qres" in kw:
        kw["qres"] = kw["qres"]()
    out, _ = full(kres, **kw)
    assert out["status"] == "ERROR" and out["reason"] == expect, (name, out)
    rec = out["error"]
    assert rec["stage"] == "pipeline" and rec["origin"] == "eph_protocol" and rec["retryable"] == 0
    assert all(EPH_ERROR_TAG_RE.match(t) for t in out["all_reasons"])


def test_anchor_within_tolerance_passes_and_gaps():
    kres = k_scan()
    pv = protocol_values(kres, 5.0)
    g = tetra_anchor_gate({x: v * 1.04 for x, v in pv.items()}, pv,
                          anchor_k_values=ANCHOR_SEQ, anchor_q_values=ANCHOR_SEQ)
    assert g["status"] == "PASS" and g["max_gap"] == pytest.approx(0.04) and not g["used_in_envelope"]
    g = tetra_anchor_gate({x: v / 1.04 for x, v in pv.items()}, pv,   # gap vs min(|a|,|p|)
                          anchor_k_values=ANCHOR_SEQ, anchor_q_values=ANCHOR_SEQ)
    assert g["max_gap"] == pytest.approx(0.04)
    g = tetra_anchor_gate({"lambda": pv["lambda"]}, pv, anchor_k_values=ANCHOR_SEQ,
                          anchor_q_values=ANCHOR_SEQ)
    assert g["status"] == "ERROR" and g["reason"] == "UNCONVERGED_EPH(tetra)"


def test_missing_or_extra_sigma_is_pin_violation():
    kres = k_scan()
    del kres[160.0]
    assert admit_k(kres)["reason"] == EPH_PIN_VIOLATION
    kres = k_scan()
    kres[7.0] = run_seq("k", 7.0)
    assert admit_k(kres)["reason"] == EPH_PIN_VIOLATION


def test_missing_hash_is_pin_violation():
    s = seq()
    s[1]["hash"] = None
    r = admit_sequence(s, direction="k", sigma_meV=10.0, omega_max=OMAX)
    assert not r["admitted"] and r["pin_violation"]


def test_decide_requires_every_gate():
    kres = k_scan()
    kg = admit_k(kres)
    out = decide(kg, None, None, None)
    assert out["status"] == "ERROR" and out["reason"] == "UNCONVERGED_EPH(q)"


def test_multiple_failures_ordered():
    kres = k_scan()
    out, _ = full(kres, anchor_shift=0.2, qres={s: run_seq("q", s, efs=None, p=-1.0)
                                                for s in AD.SIGMA_MEV})
    assert out["all_reasons"] == ["UNCONVERGED_EPH(q)", "UNCONVERGED_EPH(tetra)"]
    assert out["reason"] == "UNCONVERGED_EPH(q)"


# ----------------------------------------------------------------------------- tags
def test_error_tag_spelling():
    for g in ("k", "q", "sigma", "envelope", "tetra", "imaginary", "wannier"):
        t = unconverged_eph_tag(g)
        assert t == f"UNCONVERGED_EPH({g})" and not re.search(r"\s", t) and EPH_ERROR_TAG_RE.match(t)
    assert EPH_ERROR_TAG_RE.match(EPH_PIN_VIOLATION) and " " not in EPH_PIN_VIOLATION
    for bad in ("K", "k ", "grid", "envelope)"):
        with pytest.raises(ValueError):
            unconverged_eph_tag(bad)
    assert not EPH_ERROR_TAG_RE.match("UNCONVERGED_EPH( k)")
    with pytest.raises(AdmitError):
        error_record("UNCONVERGED_EPH(k) ", [])


def test_error_tags_are_inert_in_rtsc_blocking():
    P = pytest.importorskip("rtsc.pipeline")
    for t in [unconverged_eph_tag(g) for g in ("k", "q", "envelope")] + [EPH_PIN_VIOLATION]:
        sp = P.Spectrum((100.0,), (1.0,), 0.1, "alt", {"assumptions": [t]})
        r, c = P.blocking_tags([sp])
        assert not r and not c
