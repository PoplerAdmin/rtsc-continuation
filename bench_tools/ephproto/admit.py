"""P3 admission, G4 tetrahedron anchor and P5 envelope of the alpha^2F protocol (v3 + v3.1 + v3.2).

Source: reports/EPH_CONVERGENCE_PROTOCOL.md (P2, P3, P4 G4, P5, P7, §A A1/A2/A3/A4/A6, §D2, §V).
Everything here is pure arithmetic on already-parsed spectra (``phself.build_phself_spectrum``
output, or synthetic ``{'omega', 'weight', 'hash', 'ef_eV'}`` samples).  No DFT is run.

Rules implemented (per sigma, per direction; three nested ratio-2 grids G_1, G_2, G_3)::

    observables  X in {lambda, omega_log, A, W_b}; W_b = alpha^2F weight in window b of 20 equal
                 windows over (0, omega_max], omega_max = coarse-DFPT max omega (A2, fixed);
                 window b tracked if >= 1 % of A or >= 1 % of lambda
    steps        d_1 = X(G_2) - X(G_1),  d_2 = X(G_3) - X(G_2)
    converged    |d_1| <= 0.2 % X and |d_2| <= 0.2 % X  -> pass, R = max(|d_1|, |d_2|)      (A1)
    otherwise    same sign, p = log2(|d_1|/|d_2|) >= 0.5 -> pass,
                 R = |d_2| / (2^min(p, 4) - 1)
    admitted     every item passes; R/X <= 2 % (lambda, omega_log, A), <= 10 % (tracked W_b);
                 k direction: |E_F(G_3) - E_F(G_2)| <= sigma/4
    k            admitted sigma = passing sigma; contiguous block incl. the smallest passing one;
                 sigma_min,adm <= 20 meV                        else UNCONVERGED_EPH(k) / (sigma)
    q            every admitted sigma passes at the densest admitted k  else UNCONVERGED_EPH(q)
    G4           tetra anchor (itself P3-converged on two k and two q refinements) vs the protocol
                 at sigma_min,adm: lambda, A, omega_log within 5 % (A6)  else UNCONVERGED_EPH(tetra)
    P5           e = max(0.10, 3 max R/X over admitted sigma, k and q, G6 max rel. dgamma)  (A2)
                 G4 not in e (A6); e > 0.5 -> UNCONVERGED_EPH(envelope); recorded as
                 ``E-ENV(e, estimated by eph_protocol v3.1)``; DISCRETISATION stays an assumption
    pins         identical input hashes or bit-identical spectra in a sequence, missing or extra
                 sigma -> EPH_PIN_VIOLATION (A1, A3)

ERROR records follow P7: stage ``pipeline``, origin ``eph_protocol``, ``retryable = 0``.

Spec interpretations
--------------------
Where the text admits more than one reading, the more conservative one (more ERROR, larger
remainders) is implemented:

I1  "% of X" uses X_ref = min |X(G_i)| over the three grids, not X(G_3).
I2  R is floored at the print quantum: R >= 5e-5 * X_ref (E15.5 lambda rounding, phself), in both
    branches, so R is never 0 (A1 "never 0"; review 3 M1).
I3  Outside the converged branch a zero step (d_1 = 0 or d_2 = 0) leaves p undefined -> fail
    (review 3 M1: d_2 = 0 would otherwise give p = inf, R = 0).
I4  The converged branch waives p and the sign test only (spec: "regardless of p"; tiny
    oscillation is within R = max(|d_1|, |d_2|)).
I5  A window is tracked if it meets the 1 % share on ANY of the three grids; a tracked window
    with zero weight on some grid fails (X_ref = 0).
I6  An atom above omega_max, or omega <= 0 / negative weight, fails the sequence (no atom may sit
    outside the fixed windows).
I7  W_b is the alpha^2F weight (A share) of window b; lambda share is used only for tracking.
I8  Admitted sigma with a gap (e.g. 5, 10 pass, 20 fails, 40 passes) -> ERROR(k), rather than
    silently taking the lower block.
I9  Any pair of bit-identical spectra in one sequence is treated like identical input hashes:
    EPH_PIN_VIOLATION.  A missing or extra sigma in the k scan is also EPH_PIN_VIOLATION (it
    departs from the pre-registered set).  Grid hashes are mandatory inputs.
I10 Missing E_F on G_2 or G_3 in the k direction fails the E_F test.  If E_F is given in the q
    direction it is checked the same way.
I11 G4 remains a hard gate (ERROR on failure) in v3.1 although A6 says it is "reported"; the
    tolerance is A6's 5 %, the gap is relative to min(|anchor|, |protocol|), and the anchor's
    own P3 convergence (k and q, lambda/A/omega_log, R/X <= 2 %) is a mandatory input.
I12 G4 gaps are not in e (A6, explicit); they are recorded alongside it.  The G6 term is a
    required input: e cannot be computed without it.
I13 The E-ENV assumption is the exact A4 string (an informational assumption, not a status
    reason); only the ERROR reasons are guaranteed whitespace-free.
I14 When several gates fail, every failing reason is listed; the primary one is the first in the
    order pin, k, sigma, q, tetra, envelope.
"""
from __future__ import annotations

import math
from typing import Any, Mapping, Sequence

import numpy as np

from .tags import EPH_ERROR_TAG_RE, EPH_PIN_VIOLATION, unconverged_eph_tag


class AdmitError(ValueError):
    """Malformed input to the admission arithmetic (a programming error, not a gate result)."""


PROTOCOL_VERSION = "v3.2-draft"
ENVELOPE_LABEL = "eph_protocol v3.1"     # A4 assumption string
SIGMA_MEV = (5.0, 10.0, 20.0, 40.0, 80.0, 160.0)
K_MULTIPLES = (2, 4, 8)                  # K_1..K_3 = 2, 4, 8 K_c
K_MULTIPLES_EXTENDED = (4, 8, 16)        # after the K_4 extension: (K_2, K_3, K_4)  (A2)
Q_MULTIPLES = (2, 4, 8)

THRESHOLDS = {
    "p_min": 0.5,
    "p_cap": 4.0,
    "converged_rel": 0.002,
    "r_rel_max_scalar": 0.02,
    "r_rel_max_window": 0.10,
    "n_windows": 20,
    "window_track_share": 0.01,
    "ef_shift_sigma_frac": 0.25,
    "rel_resolution": 5e-5,
    "sigma_min_adm_max_meV": 20.0,
    "g4_tol": 0.05,
    "envelope_floor": 0.10,
    "envelope_factor": 3.0,
    "envelope_max": 0.5,
}
SCALARS = ("lambda", "omega_log", "A")


# ----------------------------------------------------------------------------- tags / records
def error_record(tag: str, reasons: Sequence[str]) -> dict[str, Any]:
    """P7 ERROR record for a gate failure or pin violation."""
    if not EPH_ERROR_TAG_RE.match(tag):
        raise AdmitError(f"not a protocol ERROR reason: {tag!r}")
    return {"status": "ERROR", "reason": tag, "stage": "pipeline", "origin": "eph_protocol",
            "retryable": 0, "details": list(reasons)}


def k_multiples(extended: bool = False) -> tuple[int, int, int]:
    """Multiples of K_c forming the k sequence (A2: after K_4 the sequence is K_2, K_3, K_4)."""
    return K_MULTIPLES_EXTENDED if extended else K_MULTIPLES


# ----------------------------------------------------------------------------- observables
def observables(omega, weight, omega_max: float, *, n_windows: int = 20) -> dict[str, Any]:
    """A, lambda, omega_log and the window weights of an atomic spectrum.

    ``omega`` and ``omega_max`` in one unit (any); ``weight`` = alpha^2F atom weights in that unit.
    Windows b = 0..n-1 are (b*dw, (b+1)*dw], dw = omega_max/n.  Atoms outside (0, omega_max] are
    counted in ``outside`` (they fail the sequence, I6)."""
    om = np.asarray(omega, dtype=float).ravel()
    w = np.asarray(weight, dtype=float).ravel()
    if om.shape != w.shape or om.size == 0:
        raise AdmitError(f"omega/weight shapes {om.shape} / {w.shape} (non-empty, equal)")
    if not (math.isfinite(omega_max) and omega_max > 0):
        raise AdmitError(f"omega_max must be finite > 0, got {omega_max!r}")
    if not (np.all(np.isfinite(om)) and np.all(np.isfinite(w))):
        raise AdmitError("non-finite omega or weight")
    bad = (om <= 0) | (w < 0)
    over = om > omega_max
    ok = ~(bad | over)
    omk, wk = om[ok], w[ok]
    A = math.fsum(wk.tolist())
    lam = 2.0 * math.fsum((wk / omk).tolist())
    wlog = (math.exp(2.0 * math.fsum((wk / omk * np.log(omk)).tolist()) / lam)
            if lam > 0 else float("nan"))
    dw = omega_max / n_windows
    idx = np.clip(np.ceil(omk / dw).astype(int) - 1, 0, n_windows - 1)
    W = np.zeros(n_windows)
    WL = np.zeros(n_windows)
    np.add.at(W, idx, wk)
    np.add.at(WL, idx, 2.0 * wk / omk)
    return {"A": A, "lambda": lam, "omega_log": wlog, "W": W.tolist(), "W_lambda": WL.tolist(),
            "outside": {"n_nonpositive_or_negative_weight": int(np.count_nonzero(bad)),
                        "n_above_omega_max": int(np.count_nonzero(over & ~bad)),
                        "weight_above_omega_max": math.fsum(w[over & ~bad].tolist())}}


def _sample(s: Any, i: int) -> dict[str, Any]:
    """Normalise one grid sample: a mapping {'omega','weight','hash'[, 'ef_eV']} or an rtsc
    Spectrum (omega_K / weight_K, hash = provenance sources epw_in sha256, ef from stdout_checks)."""
    if isinstance(s, Mapping):
        out = {"omega": s["omega"], "weight": s["weight"], "hash": s.get("hash"),
               "ef_eV": s.get("ef_eV"), "label": s.get("label", f"G{i + 1}")}
    elif hasattr(s, "omega_K") and hasattr(s, "weight_K"):
        prov = getattr(s, "provenance", {}) or {}
        h = ((prov.get("sources") or {}).get("epw_in") or {}).get("sha256")
        ef = (prov.get("stdout_checks") or {}).get("ef_eV")
        out = {"omega": s.omega_K, "weight": s.weight_K, "hash": h, "ef_eV": ef,
               "label": getattr(s, "label", f"G{i + 1}")}
    else:
        raise AdmitError(f"grid sample {i + 1}: expected a mapping or a Spectrum, got {type(s)}")
    return out


# ----------------------------------------------------------------------------- one item
def richardson_item(x: Sequence[float], rel_limit: float, *,
                    th: Mapping[str, float] = THRESHOLDS) -> dict[str, Any]:
    """P3 rule for one observable on three grids.  Returns the full record; ``pass`` is the
    conjunction of the rate test (``pass_rate``) and the remainder limit (``pass_remainder``)."""
    if len(x) != 3:
        raise AdmitError(f"need X on exactly 3 grids, got {len(x)}")
    x1, x2, x3 = (float(v) for v in x)
    rec: dict[str, Any] = {"x": [x1, x2, x3], "d1": None, "d2": None, "p": None, "branch": None,
                           "X_ref": None, "R": None, "R_rel": None, "rel_limit": rel_limit,
                           "pass_rate": False, "pass_remainder": False, "pass": False, "why": None}
    if not all(math.isfinite(v) for v in (x1, x2, x3)):
        rec["why"] = "non-finite value"
        return rec
    X = min(abs(x1), abs(x2), abs(x3))
    d1, d2 = x2 - x1, x3 - x2
    rec.update(d1=d1, d2=d2, X_ref=X)
    if X <= 0:
        rec["why"] = "X_ref = 0 (observable vanishes on some grid)"
        return rec
    floor = th["rel_resolution"] * X
    if d1 != 0 and d2 != 0:
        rec["p"] = math.log2(abs(d1) / abs(d2))
    c = th["converged_rel"] * X
    if abs(d1) <= c and abs(d2) <= c:
        rec.update(branch="converged", pass_rate=True, R=max(abs(d1), abs(d2), floor))
    else:
        rec["branch"] = "rate"
        if d1 == 0 or d2 == 0:
            rec["why"] = "zero step outside the converged branch: p undefined"
            return rec
        if (d1 > 0) != (d2 > 0):
            rec["why"] = "oscillating (d_1, d_2 of opposite sign)"
            return rec
        p = rec["p"]
        if p < th["p_min"]:
            rec["why"] = (f"diverging (p = {p:.4g} < 0)" if p < 0 else
                          f"rate too slow (p = {p:.4g} < {th['p_min']})")
            return rec
        pc = min(p, th["p_cap"])
        rec.update(pass_rate=True, R=max(abs(d2) / (2.0 ** pc - 1.0), floor))
    rec["R_rel"] = rec["R"] / X
    rec["pass_remainder"] = rec["R_rel"] <= rel_limit
    rec["pass"] = rec["pass_rate"] and rec["pass_remainder"]
    if not rec["pass_remainder"]:
        rec["why"] = f"R/X = {rec['R_rel']:.4g} > {rel_limit}"
    return rec


# ----------------------------------------------------------------------------- one sequence
def admit_sequence(samples: Sequence[Any], *, direction: str, sigma_meV: float, omega_max: float,
                   th: Mapping[str, float] = THRESHOLDS) -> dict[str, Any]:
    """P3 on one nested ratio-2 sequence (k at Q_2, or q at the densest admitted k).

    ``samples``: three grid samples, coarse to dense.  Returns {'admitted', 'pin_violation'
    (list), 'failures' (list), 'items', 'tracked_windows', 'R_rel', 'max_R_rel', 'ef_shift_meV'}."""
    if direction not in ("k", "q"):
        raise AdmitError(f"direction must be 'k' or 'q', got {direction!r}")
    if len(samples) != 3:
        raise AdmitError(f"a P3 sequence has exactly 3 grids, got {len(samples)}")
    S = [_sample(s, i) for i, s in enumerate(samples)]
    pin: list[str] = []
    fail: list[str] = []
    hashes = [s["hash"] for s in S]
    if not all(isinstance(h, str) and h for h in hashes):
        pin.append(f"input hash missing for grid(s) {[i + 1 for i, h in enumerate(hashes) if not h]}")
    elif len(set(hashes)) != 3:
        pin.append(f"grids share an input hash: {hashes}")
    arrs = [(np.asarray(s["omega"], dtype=float).ravel(), np.asarray(s["weight"], dtype=float).ravel())
            for s in S]
    for i in range(3):
        for j in range(i + 1, 3):
            if arrs[i][0].shape == arrs[j][0].shape and np.array_equal(arrs[i][0], arrs[j][0]) \
                    and np.array_equal(arrs[i][1], arrs[j][1]):
                pin.append(f"grids G{i + 1} and G{j + 1} hold bit-identical spectra")

    n = int(th["n_windows"])
    obs = [observables(o, w, omega_max, n_windows=n) for o, w in arrs]
    for i, ob in enumerate(obs):
        out = ob["outside"]
        if out["n_nonpositive_or_negative_weight"] or out["n_above_omega_max"]:
            fail.append(f"G{i + 1}: {out['n_nonpositive_or_negative_weight']} atom(s) with omega <= 0 "
                        f"or weight < 0, {out['n_above_omega_max']} above omega_max = {omega_max}")

    items: dict[str, dict[str, Any]] = {}
    for name in SCALARS:
        items[name] = richardson_item([ob[name] for ob in obs], th["r_rel_max_scalar"], th=th)
    share = th["window_track_share"]
    tracked = [b for b in range(n)
               if any(ob["A"] > 0 and ob["W"][b] >= share * ob["A"] for ob in obs)
               or any(ob["lambda"] > 0 and ob["W_lambda"][b] >= share * ob["lambda"] for ob in obs)]
    for b in tracked:
        items[f"W{b}"] = richardson_item([ob["W"][b] for ob in obs], th["r_rel_max_window"], th=th)
    for name, it in items.items():
        if not it["pass"]:
            fail.append(f"{name}: {it['why']}")

    ef_shift = None
    efs = [s["ef_eV"] for s in S]
    if direction == "k" or efs[1] is not None or efs[2] is not None:
        if efs[1] is None or efs[2] is None:
            fail.append("E_F of G2/G3 not supplied")
        else:
            ef_shift = abs(float(efs[2]) - float(efs[1])) * 1e3
            lim = th["ef_shift_sigma_frac"] * sigma_meV
            if not (ef_shift <= lim):
                fail.append(f"E_F shift G2->G3 = {ef_shift:.6g} meV > sigma/4 = {lim:.6g} meV")

    R_rel = {k: v["R_rel"] for k, v in items.items()}
    admitted = not pin and not fail
    return {"direction": direction, "sigma_meV": float(sigma_meV), "labels": [s["label"] for s in S],
            "hashes": hashes, "admitted": admitted, "pin_violation": pin, "failures": fail,
            "items": items, "tracked_windows": tracked, "observables": obs, "R_rel": R_rel,
            "max_R_rel": max(R_rel.values()) if admitted else None, "ef_shift_meV": ef_shift}


# ----------------------------------------------------------------------------- k, sigma, q
def _result(gate_tag: str | None, reasons: list[str], **kw) -> dict[str, Any]:
    out = {"status": "PASS" if gate_tag is None else "ERROR", "reason": gate_tag,
           "error": error_record(gate_tag, reasons) if gate_tag else None, "details": reasons}
    out.update(kw)
    return out


def admit_k(k_results: Mapping[float, Mapping[str, Any]], *, extended: bool = False,
            sigmas: Sequence[float] = SIGMA_MEV,
            th: Mapping[str, float] = THRESHOLDS) -> dict[str, Any]:
    """k admission over the pre-registered sigma set (results of ``admit_sequence(direction='k')``
    keyed by sigma in meV).  ``k4_extension_allowed`` is True only when no sigma passed and the
    K_4 extension has not been used yet."""
    keys = sorted(float(s) for s in k_results)
    want = sorted(float(s) for s in sigmas)
    if keys != want:
        return _result(EPH_PIN_VIOLATION, [f"k scan sigma set {keys} != pre-registered {want}"],
                       admitted_sigma_meV=[], sigma_min_adm_meV=None, k4_extension_allowed=False)
    res = {float(s): r for s, r in k_results.items()}
    pin = [f"sigma={s:g}: {m}" for s in want for m in res[s]["pin_violation"]]
    if pin:
        return _result(EPH_PIN_VIOLATION, pin, admitted_sigma_meV=[], sigma_min_adm_meV=None,
                       k4_extension_allowed=False)
    for s in want:
        if res[s]["direction"] != "k" or float(res[s]["sigma_meV"]) != s:
            raise AdmitError(f"k result keyed {s} is direction={res[s]['direction']!r}, "
                             f"sigma={res[s]['sigma_meV']!r}")
    passing = [s for s in want if res[s]["admitted"]]
    detail = {s: res[s]["failures"][:5] for s in want if not res[s]["admitted"]}
    if not passing:
        return _result(unconverged_eph_tag("k"), [f"no sigma admitted in k: {detail}"],
                       admitted_sigma_meV=[], sigma_min_adm_meV=None,
                       k4_extension_allowed=not extended)
    ix = [want.index(s) for s in passing]
    if ix != list(range(ix[0], ix[0] + len(ix))):
        return _result(unconverged_eph_tag("k"),
                       [f"admitted sigma {passing} are not a contiguous block of {want}: {detail}"],
                       admitted_sigma_meV=[], sigma_min_adm_meV=None, k4_extension_allowed=False)
    smin = passing[0]
    if smin > th["sigma_min_adm_max_meV"]:
        return _result(unconverged_eph_tag("sigma"),
                       [f"sigma_min,adm = {smin:g} meV > {th['sigma_min_adm_max_meV']:g} meV"],
                       admitted_sigma_meV=passing, sigma_min_adm_meV=smin, k4_extension_allowed=False)
    return _result(None, [], admitted_sigma_meV=passing, sigma_min_adm_meV=smin,
                   k4_extension_allowed=False)


def admit_q(q_results: Mapping[float, Mapping[str, Any]],
            admitted_sigma_meV: Sequence[float]) -> dict[str, Any]:
    """q admission: every admitted sigma must pass at the densest admitted k."""
    if not admitted_sigma_meV:
        raise AdmitError("admit_q needs a non-empty admitted sigma set (run admit_k first)")
    res = {float(s): r for s, r in q_results.items()}
    pin = [f"sigma={s:g}: {m}" for s in admitted_sigma_meV for m in
           (res.get(float(s)) or {}).get("pin_violation", [])]
    if pin:
        return _result(EPH_PIN_VIOLATION, pin)
    bad: list[str] = []
    for s in admitted_sigma_meV:
        r = res.get(float(s))
        if r is None:
            bad.append(f"sigma={s:g}: no q sequence")
            continue
        if r["direction"] != "q" or float(r["sigma_meV"]) != float(s):
            raise AdmitError(f"q result keyed {s} is direction={r['direction']!r}, "
                             f"sigma={r['sigma_meV']!r}")
        if not r["admitted"]:
            bad.append(f"sigma={s:g}: {r['failures'][:5]}")
    if bad:
        return _result(unconverged_eph_tag("q"), bad)
    return _result(None, [])


# ----------------------------------------------------------------------------- G4 tetra anchor
def tetra_anchor_gate(anchor: Mapping[str, float], protocol: Mapping[str, float], *,
                      anchor_k_values: Sequence[Mapping[str, float]] | None,
                      anchor_q_values: Sequence[Mapping[str, float]] | None,
                      th: Mapping[str, float] = THRESHOLDS) -> dict[str, Any]:
    """Gate G4 (mandatory for NC inputs).

    ``anchor``: tetrahedron lambda, A, omega_log at its densest grid; ``protocol``: the same at
    sigma_min,adm and (K_densest, Q_3).  ``anchor_k_values`` / ``anchor_q_values``: the anchor's
    lambda/A/omega_log on its three k and three q grids (two refinements each), which must pass
    P3 with R/X <= 2 %.  None -> fail (convergence of the anchor not shown).
    gap_X = |X_anchor - X_protocol| / min(|X_anchor|, |X_protocol|) <= 5 % (A6).  Not in e (A6)."""
    fails: list[str] = []
    conv: dict[str, Any] = {}
    for dname, seq in (("k", anchor_k_values), ("q", anchor_q_values)):
        if seq is None or len(seq) != 3:
            fails.append(f"anchor {dname} convergence not demonstrated (need 3 grids)")
            continue
        for x in SCALARS:
            try:
                vals = [float(g[x]) for g in seq]
            except (KeyError, TypeError, ValueError):
                fails.append(f"anchor {dname}: {x} missing on some grid")
                continue
            it = richardson_item(vals, th["r_rel_max_scalar"], th=th)
            conv[f"{dname}:{x}"] = it
            if not it["pass"]:
                fails.append(f"anchor {dname}:{x}: {it['why']}")
    gaps: dict[str, float | None] = {}
    for x in SCALARS:
        a, p = anchor.get(x), protocol.get(x)
        if a is None or p is None or not (math.isfinite(a) and math.isfinite(p)) \
                or min(abs(a), abs(p)) <= 0:
            gaps[x] = None
            fails.append(f"G4 {x}: anchor {a!r} / protocol {p!r} not finite and nonzero")
            continue
        gaps[x] = abs(a - p) / min(abs(a), abs(p))
        if not gaps[x] <= th["g4_tol"]:
            fails.append(f"G4 {x}: gap {gaps[x]:.4g} > {th['g4_tol']}")
    finite = [g for g in gaps.values() if g is not None]
    return _result(unconverged_eph_tag("tetra") if fails else None, fails, gaps=gaps,
                   max_gap=max(finite) if finite else None, anchor_convergence=conv,
                   tol=th["g4_tol"], used_in_envelope=False)


# ----------------------------------------------------------------------------- P5 envelope
def envelope(k_results: Mapping[float, Mapping[str, Any]],
             q_results: Mapping[float, Mapping[str, Any]],
             admitted_sigma_meV: Sequence[float], *, g6_rel_dgamma: float,
             g4_gaps: Mapping[str, float | None] | None = None,
             th: Mapping[str, float] = THRESHOLDS) -> dict[str, Any]:
    """P5 / A2: e = max(0.10, 3 * max R/X, G6 max relative dgamma).

    The max runs over the admitted sigma, both directions, and every item of the sequence
    (lambda, omega_log, A, every tracked W_b).  G4 gaps are recorded, not used (A6).  The result is
    an ESTIMATE; discretisation stays an assumption (section 0)."""
    if not admitted_sigma_meV:
        raise AdmitError("envelope needs a non-empty admitted sigma set")
    g6 = float(g6_rel_dgamma)
    if not (math.isfinite(g6) and g6 >= 0):
        raise AdmitError(f"G6 relative dgamma must be finite >= 0, got {g6_rel_dgamma!r}")
    terms: dict[str, float] = {}
    for dname, res in (("k", k_results), ("q", q_results)):
        res = {float(s): r for s, r in res.items()}
        for s in admitted_sigma_meV:
            r = res.get(float(s))
            if r is None or not r["admitted"]:
                raise AdmitError(f"envelope: sigma={s:g} has no admitted {dname} sequence")
            for name, v in r["R_rel"].items():
                terms[f"{dname}:{float(s):g}:{name}"] = v
    rmax_key = max(terms, key=terms.get)
    rmax = terms[rmax_key]
    e = max(th["envelope_floor"], th["envelope_factor"] * rmax, g6)
    tag = unconverged_eph_tag("envelope") if e > th["envelope_max"] else None
    reasons = [f"e = {e:.6g} > {th['envelope_max']}"] if tag else []
    return _result(tag, reasons, e=e, estimated=True, max_R_rel=rmax, max_R_rel_term=rmax_key,
                   g6_rel_dgamma=g6, g4_gaps_reported=dict(g4_gaps or {}), terms=terms,
                   rule="e = max(0.10, 3*max R/X (admitted sigma, k and q), G6 rel dgamma); "
                        "G4 not used (A6)",
                   assumption=envelope_assumption(e))


def envelope_assumption(e: float) -> str:
    """The A4 assumption string ``E-ENV(e, estimated by eph_protocol v3.1)``."""
    return f"E-ENV({e:.6g}, estimated by {ENVELOPE_LABEL})"


def envelope_provenance(e: float, grid_label: str, *, config_floor: float = 0.0) -> dict[str, Any]:
    """Provenance additions for one forwarded protocol Spectrum (A4): ``rel_envelope`` (config
    floor still applies), the ``discretisation`` key that makes rtsc emit DISCRETISATION, and the
    E-ENV assumption.  The caller merges ``assumptions`` into the spectrum's list."""
    if not (math.isfinite(e) and e > 0):
        raise AdmitError(f"envelope must be finite > 0, got {e!r}")
    return {"rel_envelope": max(float(e), float(config_floor)),
            "envelope_status": "estimated",
            "discretisation": f"eph_protocol_{ENVELOPE_LABEL.split()[-1]}:{grid_label}",
            "assumptions": [envelope_assumption(e)]}


# ----------------------------------------------------------------------------- whole decision
_ORDER = (EPH_PIN_VIOLATION, "UNCONVERGED_EPH(k)", "UNCONVERGED_EPH(sigma)", "UNCONVERGED_EPH(q)",
          "UNCONVERGED_EPH(tetra)", "UNCONVERGED_EPH(envelope)")


def decide(k_gate: Mapping[str, Any], q_gate: Mapping[str, Any] | None,
           tetra_gate: Mapping[str, Any] | None, env: Mapping[str, Any] | None) -> dict[str, Any]:
    """Combine the gate results into the protocol outcome.  Any missing gate (None) or failed gate
    is ERROR; admission requires every gate to have run and passed."""
    tags: list[str] = []
    details: list[str] = []
    for name, g in (("k", k_gate), ("q", q_gate), ("tetra", tetra_gate), ("envelope", env)):
        if g is None:
            details.append(f"gate {name} not evaluated")
            continue
        if g["status"] != "PASS":
            tags.append(g["reason"])
            details += [f"{g['reason']}: {d}" for d in g["details"]]
    missing = [d for d in details if d.startswith("gate ")]
    if missing and not tags:
        # a gate that never ran cannot be admitted; attribute it to the first missing gate
        first = missing[0].split()[1]
        tags.append(unconverged_eph_tag({"k": "k", "q": "q", "tetra": "tetra",
                                         "envelope": "envelope"}[first]))
    tags = sorted(set(tags), key=_ORDER.index)
    if tags:
        rec = error_record(tags[0], details)
        rec["all_reasons"] = tags
        return {"status": "ERROR", "reason": tags[0], "all_reasons": tags, "error": rec,
                "protocol_version": PROTOCOL_VERSION}
    return {"status": "ADMITTED", "reason": None, "all_reasons": [], "error": None,
            "draft": True, "protocol_version": PROTOCOL_VERSION,
            "admitted_sigma_meV": list(k_gate["admitted_sigma_meV"]),
            "sigma_min_adm_meV": k_gate["sigma_min_adm_meV"], "e": env["e"],
            "assumption": env["assumption"]}
