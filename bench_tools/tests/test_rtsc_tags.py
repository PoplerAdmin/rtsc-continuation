"""How rtsc AS-IS treats provenance tags on protocol spectra (no rtsc code is modified).

Source of truth: rtsc/pipeline.py ``blocking_tags`` / ``apply_input_policy`` (called at the end of
``certify_spectra``) and ``annotate_model`` (copies every spectrum assumption into the verdict).
"""
import re

import pytest

import synthetic_epw as syn
from ephproto import tags as T
from ephproto.phself import build_phself_spectrum

P = pytest.importorskip("rtsc.pipeline")
from rtsc.status import Status  # noqa: E402

LF, WF, QF = syn.LAMBDA_FILE, syn.LINEWIDTH_FILE, syn.QFILE


def test_regex_copy_matches_rtsc():
    assert P._TAG_MASS.pattern == T.RTSC_TAG_MASS_RE.pattern


def sp_with(tag):
    return P.Spectrum((100.0,), (1.0,), 0.1, "alt", {"assumptions": [tag]})


# (tag, blocks REJECTED, blocks CERTIFIED) -- rtsc as-is
MATRIX = [
    ("A-DROP-NONPOS(12.5K): excluded weight", True, False),
    ("A-DROP-NONPOS(1e-30K): tiny but positive", True, False),
    ("A-DROP-NONPOS(1e+400K): overflows to inf", True, False),
    ("A-DROP-NONPOS(0K): none", False, False),
    ("A-DROP-NONPOS(-3K): negative mass", False, False),
    ("A-DROP-NONPOS(5 K): space before K", False, False),          # protocol text spelling!
    ("A-DROP-NONPOS(m K)", False, False),
    ("A-DROP-NONPOS(5K)", True, False),                            # no description needed
    ("EPH: A-DROP-NONPOS(5K): not at start", False, False),
    ("a-drop-nonpos(5K): lower case", False, False),
    ("A-DROP-EXCLUDED(3 modes)", False, False),                    # new names are inert
    ("A-NEG-CLIP(3K): clipped", False, True),
    ("A-NEG-CLIP(0K): none", False, False),
    ("A-FILE-COMPLETE(unverified): unknown dropped weight", True, False),
    ("A-FILE-COMPLETE(unverified)", True, False),
    ("A-FILE-COMPLETE: verified", False, False),
    ("A-FILE-INCOMPLETE: truncated", True, True),
    ("A-FILE-INCOMPLETENESS-anything", True, True),                # prefix match
    (T.EPH_PROTOCOL_TAG, False, False),
    ("UNCONVERGED_EPH(envelope)", False, False),
]


@pytest.mark.parametrize("tag,rej,cert", MATRIX)
def test_blocking_matrix_as_is(tag, rej, cert):
    r, c = P.blocking_tags([sp_with(tag)])
    assert (bool(r), bool(c)) == (rej, cert), (tag, r, c)


def test_mass_without_digits_crashes_blocking_tags():
    # '[-+0-9.eE]+' admits 'e', '+', '.', so float() raises inside the certify path
    with pytest.raises(ValueError):
        P.blocking_tags([sp_with("A-DROP-NONPOS(eK): bogus")])


def test_tag_builders_produce_blocking_strings():
    assert P.blocking_tags([sp_with(T.drop_nonpos_tag(2.5, "x"))])[0]
    assert P.blocking_tags([sp_with(T.drop_nonpos_tag(1e-300, "x"))])[0]
    assert P.blocking_tags([sp_with(T.completeness_unverified_tag("x"))]) [0]
    r, c = P.blocking_tags([sp_with(T.incomplete_tag("x"))])
    assert r and c
    for bad in (0.0, -1.0, float("nan"), float("inf")):
        with pytest.raises(ValueError):
            T.drop_nonpos_tag(bad, "x")


# ----------------------------------------------------------------------------- end to end
@pytest.fixture(scope="module")
def protocol_spectrum():
    d = syn.fixture_dir("normal")
    sp = build_phself_spectrum(d / LF, d / WF, d / QF, d / "epw.in", d / "epw.out")
    return sp.with_envelope(0.10)


def ccfg(**kw):
    c = {"T_target_K": 300.0, "mu_hi": 0.15, "mode": "decide"}
    c.update(kw)
    return P.normalize_config({"certify": c})["certify"]


def with_tags(sp, *extra):
    prov = dict(sp.provenance)
    prov["assumptions"] = list(prov.get("assumptions", [])) + list(extra)
    return P.Spectrum(sp.omega_K, sp.weight_K, sp.rel_envelope, sp.label, prov)


def test_protocol_spectrum_e2e_rejected_and_tag_propagated(protocol_spectrum):
    v = P.certify_spectra("SYNTHETIC-phself", [protocol_spectrum], ccfg())
    assert v.status is Status.REJECTED_PROVEN, (v.status, v.reasons)
    assert T.EPH_PROTOCOL_TAG in v.assumptions           # informational tag reaches the verdict


@pytest.mark.parametrize("tag,expect", [
    ("A-DROP-NONPOS(5K): excluded modes", Status.UNDECIDED),
    ("A-DROP-NONPOS(5 K): excluded modes", Status.REJECTED_PROVEN),     # spelling gotcha
    ("A-DROP-NONPOS(0K): excluded modes", Status.REJECTED_PROVEN),
    ("A-DROP-EXCLUDED(3 modes)", Status.REJECTED_PROVEN),
    ("A-FILE-COMPLETE(unverified): excluded modes, mass unknown", Status.UNDECIDED),
])
def test_drop_tag_blocks_rejection_e2e(protocol_spectrum, tag, expect):
    alts = [protocol_spectrum, with_tags(protocol_spectrum, tag)]
    v = P.certify_spectra("SYNTHETIC-phself", alts, ccfg())
    assert v.status is expect, (tag, v.status, v.reasons)
    assert tag in v.assumptions
    if expect is Status.UNDECIDED:
        assert v.T_up_K is None and re.search("rejection blocked", v.reasons[-1])


def test_opt_in_unblocks_with_recorded_assumption(protocol_spectrum):
    tag = "A-DROP-NONPOS(5K): excluded modes"
    v = P.certify_spectra("SYNTHETIC-phself", [with_tags(protocol_spectrum, tag)],
                          ccfg(allow_unverified_rejection=True))
    assert v.status is Status.REJECTED_PROVEN
    assert any(a.startswith("UNVERIFIED-INPUT-ACCEPTED(allow_unverified_rejection)") for a in v.assumptions)
