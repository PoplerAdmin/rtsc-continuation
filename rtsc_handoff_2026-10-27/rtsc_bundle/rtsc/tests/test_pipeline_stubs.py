"""TEST-ONLY STUB of rtsc.certify.classify (team CERT owns the real one).

!!! This stub is NOT rigorous and is used only by the ORCH tests when the real
!!! rtsc.certify.classify is not importable yet.  It never ships in the package.

Its rules (float arithmetic, labelled "STUB-*"):
  * zero spectra -> UNDECIDED (mirrors ARCH 4b F1);
  * max (1+e)A < a0 T*            -> REJECTED_PROVEN  (rule "STUB-R-UNIV");
  * min lambda >= 2.5 over spectra -> CERTIFIED_CANDIDATE ("STUB-CERT");
  * else UNDECIDED.
"""
from __future__ import annotations

import importlib

from rtsc.status import Status, Verdict

A0 = 4.33148145


def stub_classify(spectra, params, material_id=""):
    p = params if isinstance(params, dict) else vars(params)
    T = float(p.get("T_target_K", 300.0))
    base = dict(material_id="", T_low_K=None, T_up_K=None, flags=[],
                assumptions=["STUB-CERTIFY (not rigorous; tests only)"],
                certificate={"stub": True})
    if not spectra:
        return Verdict(status=Status.UNDECIDED, reasons=["no spectra"], rule_ids=[], **base)
    A_hi = max((1 + s.rel_envelope) * s.A() for s in spectra)
    if A_hi < A0 * T:
        b = dict(base, T_up_K=A_hi / A0)
        return Verdict(status=Status.REJECTED_PROVEN, reasons=["stub universal ceiling"],
                       rule_ids=["STUB-R-UNIV"], **b)
    if min(s.lambda_() for s in spectra) >= 2.5:
        b = dict(base, T_low_K=T + 1.0)
        return Verdict(status=Status.CERTIFIED_CANDIDATE, reasons=["stub"], rule_ids=["STUB-CERT"], **b)
    return Verdict(status=Status.UNDECIDED, reasons=["stub"], rule_ids=[], **base)


def real_classify_available() -> bool:
    try:
        mod = importlib.import_module("rtsc.certify")
    except ImportError:
        return False
    return hasattr(mod, "classify")


def install(monkeypatch, force: bool = False):
    """Use the real classify when present (unless ``force``), else patch in the stub."""
    mod = importlib.import_module("rtsc.certify")
    if force or not hasattr(mod, "classify"):
        import rtsc.pipeline as P
        monkeypatch.setattr(mod, "classify", stub_classify, raising=False)
        monkeypatch.setattr(P, "_certify_params", lambda ccfg: dict(ccfg))
        return "stub"
    return "real"
