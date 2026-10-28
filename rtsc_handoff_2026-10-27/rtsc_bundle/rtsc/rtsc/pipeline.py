"""Stage orchestration (ARCH section 6).

Stages and the only statuses each may produce:

    1 source            ERROR (unreadable input file)
    2 dedupe            DEFERRED  (duplicate of another input; not computed separately)
    3 existence filters EXCLUDED_ASSUMPTION  (rtsc.stability; criterion/value/threshold recorded)
    4 priority+budget   DEFERRED  (rtsc.descriptors; heuristic, S2)
    5 ML pre-check      DEFERRED  (optional; heuristic, S2)
    6 spectra / DFPT    DEFERRED (no spectrum and DFPT disabled), EXCLUDED_ASSUMPTION (imaginary
                        first-principles modes), ERROR
    7 certify.classify  REJECTED_PROVEN | CERTIFIED_CANDIDATE | UNDECIDED  (only here, S1)
    8 validity flags    no status change (S6)
    9 report

Per-material records: ``<workdir>/<safe id>/record.json`` (atomic writes; ``safe_id`` is
injective and plan() refuses ids that would still share a directory).  With ``resume = true`` a
record is reused only when its material_id, configuration hash AND input digest (candidate +
contents of its structure, a2F and band-coupling files) match; stored spectra are reused when
only the certify stage is missing; ERRORs from crashed worker processes are retried.

S4: ``account()`` checks that every input id has exactly one verdict and that no verdict has an id
that is not an input; ``run()`` raises ``AccountingError`` otherwise.  S1: every
REJECTED_PROVEN / CERTIFIED_CANDIDATE verdict must carry origin ``rtsc.certify.classify``.
"""
from __future__ import annotations

import concurrent.futures as cf
import copy
import dataclasses
import datetime as _dt
import hashlib
import json
import logging
import math
import os
import re
import tomllib
import traceback
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Sequence

from .descriptors import PriorityWeights, plan_budget, priority
from .descriptors.ml import get_ml_hook
from .sources import Candidate, dedupe, load_sources, safe_id
from .spectrum import Spectrum
from .stability import FilterSet
from .status import PROVEN_STATUSES, Status, Verdict
from .validity import ValidityConfig, apply_flags, compute_flags

log = logging.getLogger("rtsc.pipeline")

RECORD_SCHEMA = "rtsc.record/1"
CLASSIFY_ORIGIN = "rtsc.certify.classify"

# ----------------------------------------------------------------------------- configuration
DEFAULTS: dict[str, Any] = {
    "run": {"name": "rtsc-run", "workdir": "rtsc_run", "workers": 1, "resume": True},
    "sources": [],
    "dedupe": {"method": "auto", "ltol": 0.2, "stol": 0.3, "angle_tol": 5.0, "symprec": 0.1,
               "pressure_tol_GPa": 1.0},
    "filters": {"hull": {}, "metallic": {}, "magnetic": {}, "dynamic": {}},
    "priority": {},
    "budget": {},                       # max_materials, min_score
    "ml": {"enabled": False},
    "spectra": {"use_attached": True, "rel_envelope": 0.10, "two_column_rule": "bins",
                "negative_policy": "error"},
    "dfpt": {"enabled": False},
    # omega_c_K (fixed, kelvin) overrides omega_c_factor * max(omega); cutoff_rule: "block"
    "certify": {"T_target_K": 300.0, "mu_hi": 0.15, "omega_c_factor": 10.0, "cutoff_rule": "block"},
    "validity": {},
    "report": {"top_n": 25},
}
_SECTIONS = set(DEFAULTS)
_RUN_KEYS = {"name", "workdir", "workers", "resume"}
_SPECTRA_KEYS = {"use_attached", "rel_envelope", "two_column_rule", "negative_policy"}
_BUDGET_KEYS = {"max_materials", "min_score"}
_DEDUPE_KEYS = set(DEFAULTS["dedupe"])
_REPORT_KEYS = {"top_n", "outdir"}


class ConfigError(ValueError):
    pass


class AccountingError(AssertionError):
    """S4 violated: some input has no status, several, or an unknown id appeared."""


def _merge(base: dict, over: dict) -> dict:
    out = copy.deepcopy(base)
    for k, v in over.items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _merge(out[k], v)
        else:
            out[k] = copy.deepcopy(v)
    return out


def _real(cc: dict[str, Any], key: str) -> float:
    v = cc.get(key)
    if isinstance(v, bool) or not isinstance(v, (int, float)):
        raise ConfigError(f"certify.{key} must be a real number, got {v!r}")
    v = float(v)
    if not math.isfinite(v):
        raise ConfigError(f"certify.{key} must be finite, got {v!r}")
    return v


def validate_certify_config(cc: dict[str, Any]) -> None:
    """[certify] values that change what a verdict claims are checked strictly (audit F4/F9)."""
    if cc.get("cutoff_rule", "block") != "block":
        raise ConfigError("certify.cutoff_rule: only 'block' (ARCH section 4 N_c rule) is defined")
    T = _real(cc, "T_target_K")
    if not T > 0:
        raise ConfigError(f"certify.T_target_K must be > 0, got {T}")
    mu = _real(cc, "mu_hi")
    if not mu >= 0:
        raise ConfigError(f"certify.mu_hi must be >= 0, got {mu}")
    if cc.get("omega_c_K") is not None:
        oc = _real(cc, "omega_c_K")
        if not oc >= math.pi * T:
            raise ConfigError(f"certify.omega_c_K = {oc} K < pi*T* = {math.pi * T:.6g} K: the "
                              "Coulomb block would be empty at T* (N_c(T*) = 0)")
    for k in ("allow_unverified_rejection", "allow_unverified_certification"):
        if k in cc and not isinstance(cc[k], bool):
            raise ConfigError(f"certify.{k} must be true or false")
    if "omega_c_factor" in cc:
        f = _real(cc, "omega_c_factor")
        if not f >= 1.0:
            raise ConfigError(f"certify.omega_c_factor must be >= 1 (omega_c >= max phonon "
                              f"frequency), got {f}")


def normalize_config(raw: dict[str, Any], base_dir: str | os.PathLike | None = None) -> dict[str, Any]:
    unknown = set(raw) - _SECTIONS
    if unknown:
        raise ConfigError(f"unknown config sections: {sorted(unknown)}")
    cfg = _merge(DEFAULTS, raw)
    for sec, keys in (("run", _RUN_KEYS), ("spectra", _SPECTRA_KEYS), ("budget", _BUDGET_KEYS),
                      ("dedupe", _DEDUPE_KEYS), ("report", _REPORT_KEYS)):
        bad = set(cfg[sec]) - keys
        if bad:
            raise ConfigError(f"unknown keys in [{sec}]: {sorted(bad)}")
    if not isinstance(cfg["sources"], list):
        raise ConfigError("[[sources]] must be an array of tables")
    # validate the sub-configs we own by constructing them
    FilterSet.from_config(cfg["filters"])
    PriorityWeights.from_config(cfg["priority"])
    ValidityConfig.from_config(cfg["validity"])
    cc = cfg["certify"]
    validate_certify_config(cc)
    try:   # validate [certify] keys against team CERT's CertifyParams when it is importable
        _certify_params(cc)
    except ImportError:
        pass
    env = float(cfg["spectra"]["rel_envelope"])
    if not (0.0 <= env < 1.0):
        raise ConfigError("spectra.rel_envelope must be in [0, 1)")
    if base_dir is not None:   # resolve relative paths against the config file's directory
        b = Path(base_dir)
        for s in cfg["sources"]:
            for key in ("path", "manifest"):
                if key in s and not Path(s[key]).is_absolute():
                    s[key] = str((b / s[key]).resolve())
        h = cfg["filters"].get("hull", {})
        if h.get("reference_file") and not Path(h["reference_file"]).is_absolute():
            h["reference_file"] = str((b / h["reference_file"]).resolve())
        for sec, key in (("run", "workdir"), ("report", "outdir")):
            if cfg[sec].get(key) and not Path(cfg[sec][key]).is_absolute():
                cfg[sec][key] = str((b / cfg[sec][key]).resolve())
    return cfg


def load_config(path: str | os.PathLike) -> dict[str, Any]:
    p = Path(path)
    with open(p, "rb") as fh:
        raw = tomllib.load(fh)
    return normalize_config(raw, p.parent)


def config_hash(cfg: dict[str, Any]) -> str:
    """Hash of everything that can change a verdict (run.workers/resume/report excluded)."""
    c = copy.deepcopy(cfg)
    c.pop("report", None)
    for k in ("workers", "resume", "name"):
        c.get("run", {}).pop(k, None)
    return hashlib.sha256(json.dumps(c, sort_keys=True, default=str).encode()).hexdigest()[:16]


# ----------------------------------------------------------------------------- records
def _now() -> str:
    return _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds")


def _atomic_write_json(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + f".tmp{os.getpid()}")
    tmp.write_text(json.dumps(data, indent=1, default=str, allow_nan=True))
    os.replace(tmp, path)


def record_path(workdir: str | os.PathLike, material_id: str) -> Path:
    return Path(workdir) / safe_id(material_id) / "record.json"


def read_record(workdir, material_id) -> dict[str, Any] | None:
    """The record of ``material_id``, or None.  A record whose stored material_id differs (a
    path collision, or a hand-copied directory) is never returned (audit F3)."""
    p = record_path(workdir, material_id)
    if not p.is_file():
        return None
    try:
        rec = json.loads(p.read_text())
    except (OSError, json.JSONDecodeError):
        return None
    if not isinstance(rec, dict) or rec.get("material_id") != material_id:
        log.warning("ignoring %s: it belongs to %r, not %r", p,
                    rec.get("material_id") if isinstance(rec, dict) else None, material_id)
        return None
    return rec


def _file_digest(h, path: str | None, tag: str) -> None:
    h.update(f"\0{tag}\0{path}\0".encode())
    if not path:
        return
    try:
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
    except OSError as exc:
        h.update(f"<unreadable: {type(exc).__name__}>".encode())


def input_digest(c: Candidate) -> str:
    """sha256 of everything a verdict depends on besides the config: the candidate record and
    the CONTENTS of its attached a2F files, band-coupling files and structure file (audit F1)."""
    h = hashlib.sha256()
    h.update(json.dumps(c.to_dict(), sort_keys=True, default=str).encode())
    for sf in c.spectra_files:
        _file_digest(h, sf.get("path"), "spectrum")
        if sf.get("lambda_file"):
            _file_digest(h, sf["lambda_file"], "lambda")
    for bc in c.band_coupling_files:
        _file_digest(h, bc.get("path"), "band_coupling")
    _file_digest(h, c.provenance.get("file"), "structure")
    return h.hexdigest()


def new_record(c: Candidate, chash: str, digest: str | None = None) -> dict[str, Any]:
    return {"schema": RECORD_SCHEMA, "material_id": c.material_id, "config_hash": chash,
            "input_digest": digest or input_digest(c),
            "candidate": c.to_dict(), "stages": {}, "spectra": [], "verdict": None,
            "verdict_origin": None, "created_utc": _now()}


def record_reusable(old: dict[str, Any] | None, c: Candidate, chash: str, digest: str) -> bool:
    """Resume only a record made for this very input under this very configuration."""
    return bool(old and old.get("schema") == RECORD_SCHEMA and old.get("material_id") == c.material_id
                and old.get("config_hash") == chash and old.get("input_digest") == digest
                and not old.get("retryable"))


def _pipeline_verdict(stage: str, material_id: str, status: Status, reasons: list[str], *,
                      rule_ids: list[str] | None = None, certificate: dict | None = None) -> Verdict:
    """Verdicts made outside rtsc.certify: only DEFERRED / EXCLUDED / ERROR (S1, S2, S3)."""
    allowed = {"dedupe": {Status.DEFERRED}, "budget": {Status.DEFERRED},
               "ml": {Status.DEFERRED}, "spectra": {Status.DEFERRED},
               "error": {Status.ERROR}}[stage]
    if status not in allowed:
        raise AssertionError(f"stage {stage!r} may not produce {status}")
    return Verdict(material_id=material_id, status=status, T_low_K=None, T_up_K=None,
                   reasons=list(reasons), rule_ids=list(rule_ids or []), assumptions=[], flags=[],
                   certificate=dict(certificate or {"kind": stage, "proved": False}))


def _finish(rec: dict[str, Any], verdict: Verdict, origin: str) -> dict[str, Any]:
    if verdict.status in PROVEN_STATUSES and origin != CLASSIFY_ORIGIN:
        raise AssertionError(f"S1 violated: {verdict.status} from {origin}")
    rec["verdict"] = verdict.to_dict()
    rec["verdict_origin"] = origin
    rec["finished_utc"] = _now()
    return rec


# ----------------------------------------------------------------------------- spectra
#: DFT consistency checks whose failure means a PARSE / UNITS problem: fatal (ERROR)
FATAL_DFPT_CHECKS = ("lambda_vs_a2F_footer", "lambda_vs_epw")


def _fatal_checks(checks) -> list[str]:
    return [f"{c.name}: ours {c.ours:.8g} vs reference {c.reference:.8g} (rel {c.rel_diff:.3g})"
            for c in checks if not c.ok]


def _read_spectrum_file(spec: dict[str, Any], scfg: dict[str, Any]) -> list[Spectrum]:
    """Read an attached alpha^2F file through team DFT's ``rtsc.dfpt.load_attached_spectra``,
    which parses AND verifies it (QE footer / lambda file, EPW lambda, optional expected lambda
    or A for two-column files) and raises on incomplete or inconsistent files (audit F5).
    No fallback reader: without rtsc.dfpt the file is not read.

    ``spectra_files`` keys: path, format (two_column | qe_a2f_dos | epw), units (two_column
    only), columns, rule, label, negative_policy, lambda_file, expected_lambda, expected_A_K.
    """
    try:
        from .dfpt import load_attached_spectra
    except ImportError as exc:
        raise ImportError(f"reading {spec.get('path')} needs rtsc.dfpt: {exc}") from exc
    fmt = spec.get("format", "two_column")
    unit = spec.get("units") or spec.get("unit")
    if fmt == "two_column" and not unit:
        raise ValueError(f"spectra_files entry {spec['path']}: 'units' is required")
    kw: dict[str, Any] = {
        "label": spec.get("label", Path(spec["path"]).name),
        "rule": spec.get("rule", scfg.get("two_column_rule", "bins")),
        "negative_policy": spec.get("negative_policy", scfg.get("negative_policy", "error")),
        "columns": tuple(spec.get("columns", (0, 1))),
    }
    for k in ("lambda_file", "expected_lambda", "expected_A_K"):
        if spec.get(k) is not None:
            kw[k] = spec[k]
    return list(load_attached_spectra(spec["path"], fmt, unit if fmt == "two_column" else None,
                                      **kw))


def load_band_couplings(c: Candidate) -> list[Any]:
    """rtsc.dfpt.BandCoupling objects attached via the manifest key ``band_coupling``."""
    if not c.band_coupling_files:
        return []
    from .dfpt import BandCoupling
    return [BandCoupling.from_file(bc["path"], bc.get("table")) for bc in c.band_coupling_files]


def attached_spectra(c: Candidate, scfg: dict[str, Any]) -> list[Spectrum]:
    out = [Spectrum.from_dict(d) for d in c.spectra]
    for sf in c.spectra_files:
        out.extend(_read_spectrum_file(sf, scfg))
    return out


def _mark_synthetic(spectra: list[Spectrum]) -> list[Spectrum]:
    out = []
    for s in spectra:
        prov = dict(s.provenance)
        prov["synthetic"] = True
        lab = s.label if "SYNTHETIC" in s.label.upper() else f"SYNTHETIC {s.label}".strip()
        out.append(Spectrum(s.omega_K, s.weight_K, s.rel_envelope, lab, prov))
    return out


def _with_envelope(spectra: list[Spectrum], e: float) -> list[Spectrum]:
    """Conservative: the envelope used is max(spectrum's own, config) (never narrowed)."""
    return [s if s.rel_envelope >= e else s.with_envelope(e) for s in spectra]


def dfpt_structure(c: Candidate) -> Any:
    """Structure handed to team DFT: pymatgen Structure if available, else a plain dict with
    explicit masses (IUPAC standard atomic weights; recorded)."""
    structure = c.pymatgen_structure()
    if structure is not None:
        return structure
    if c.lattice is None or not c.species or c.frac_coords is None:
        raise ValueError("no ordered crystal structure for DFPT")
    from .sources.composition import ATOMIC_MASS
    return {"lattice": c.lattice, "species": c.species, "frac_coords": c.frac_coords,
            "masses": {e: ATOMIC_MASS[e] for e in set(c.species)}}


def run_dfpt(c: Candidate, mat_dir: Path, dcfg: dict[str, Any]) -> tuple[list[Spectrum], dict[str, Any]]:
    """Delegate one material to team DFT (``rtsc.dfpt.EPhWorkflow``; its own MPI launch via the
    configured ``Launcher``).  Returns (alternative spectra, info).  ``[dfpt]`` keys:
    ``launcher`` (rtsc.dfpt.Launcher), ``pseudos`` (rtsc.dfpt.PseudoTable), ``settings``
    (rtsc.dfpt.EPhSettings; ``pressure_GPa`` defaults to the candidate's target pressure)."""
    from .dfpt import EPhSettings, EPhWorkflow, Launcher, PseudoTable
    pressure = float(c.properties.get("target_pressure_GPa", c.pressure_GPa or 0.0))
    scfg = dict(dcfg.get("settings", {}))
    scfg.setdefault("pressure_GPa", pressure)
    settings = EPhSettings.from_config(scfg)
    if "pseudos" not in dcfg:
        raise ConfigError("[dfpt.pseudos] is required when DFPT is enabled")
    pseudos = PseudoTable.from_config(dcfg["pseudos"])
    launcher = Launcher.from_config(dcfg.get("launcher"))
    structure = dfpt_structure(c)
    wf = EPhWorkflow(c.material_id, structure, mat_dir / "dfpt", settings, pseudos, launcher)
    res = wf.run()
    if res is None:
        return [], {"note": "workflow returned no result"}
    info = _jsonable_info(res.to_dict() if hasattr(res, "to_dict") else {})
    info["structure_source"] = "pymatgen" if not isinstance(structure, dict) else \
        "candidate arrays + IUPAC standard atomic weights"
    return list(res.spectra), info


def _jsonable_info(info: Any) -> dict[str, Any]:
    from .sources.base import _jsonable
    d = _jsonable(info)
    if isinstance(d, dict):
        d.pop("spectra", None)
        return d
    return {"info": d}


def _phonon_from_dfpt(info: dict[str, Any], spectra: list[Spectrum]) -> dict[str, Any] | None:
    """First-principles phonon data for the dynamic filter, if the DFT stage provided it."""
    for key in ("phonon", "phonons"):
        if isinstance(info.get(key), dict):
            return info[key]
    for key, unit in (("min_frequency_THz", "THz"), ("min_frequency_cm1", "cm-1")):
        if info.get(key) is not None:
            return {"min_frequency": float(info[key]), "unit": unit, "source": "dfpt"}
    vals = [s.provenance.get("min_freq_matdyn_grid_cm1") for s in spectra]
    vals = [float(v) for v in vals if v is not None]
    if vals:   # rtsc.dfpt.collect_qe_eph: minimum over the matdyn freq.dos grid (cm^-1)
        return {"min_frequency": min(vals), "unit": "cm-1",
                "source": "dfpt matdyn grid (min_freq_matdyn_grid_cm1)"}
    return None


# ----------------------------------------------------------------------------- certify
#: [certify] config key -> rtsc.certify.CertifyParams field (other keys pass through unchanged)
_CERTIFY_KEYMAP = {"T_target_K": "T_star", "omega_c_K": "omega_c"}
_CERTIFY_ORCH_ONLY = {"cutoff_rule",           # validated here ("block" is the only rule)
                      "allow_unverified_rejection", "allow_unverified_certification"}


def _certify_params(ccfg: dict[str, Any]):
    """Build ``rtsc.certify.CertifyParams`` from the ``[certify]`` table."""
    from .certify import CertifyParams
    names = {f.name for f in dataclasses.fields(CertifyParams)}
    kw = {}
    for k, v in ccfg.items():
        if k in _CERTIFY_ORCH_ONLY:
            continue
        name = _CERTIFY_KEYMAP.get(k, k)
        if name not in names:
            raise ConfigError(f"unknown [certify] key {k!r} (CertifyParams fields: {sorted(names)})")
        kw[name] = v
    return CertifyParams(**kw)


def certify_spectra(material_id: str, spectra: Sequence[Spectrum], ccfg: dict[str, Any],
                    band_couplings: Sequence[Any] | None = None) -> Verdict:
    """The single call into rtsc.certify.classify (S1).  All alternatives are passed together;
    the S5 conjunction is certify's job.  Afterwards only assumptions and flags are added
    (cutoff rule, discretisation, single-band model, Coulomb-inactive); status and bounds are
    never touched."""
    validate_certify_config(ccfg)
    spectra = list(spectra)
    if not spectra:   # ARCH 4b F1: the S5 conjunction over zero alternatives is vacuous
        raise ValueError("refusing to classify zero alternative spectra (ARCH 4b F1)")
    from .certify import classify
    kw: dict[str, Any] = {"material_id": material_id}
    if band_couplings:
        kw["band_couplings"] = list(band_couplings)
    v = classify(spectra, _certify_params(ccfg), **kw)
    if not isinstance(v, Verdict):
        raise TypeError(f"certify.classify returned {type(v).__name__}, expected Verdict")
    if v.material_id != material_id:
        v = dataclasses.replace(v, material_id=material_id)
    v = ensure_cutoff_assumption(v, ccfg, spectra)
    v = annotate_model(v, ccfg, spectra)
    v = dataclasses.replace(v, assumptions=list(v.assumptions) + [trust_base()])
    return apply_input_policy(v, spectra, ccfg)


#: python-flint version the trusted base was verified with (correspondence audit G5)
VERIFIED_PYTHON_FLINT = "0.9.0"


def trust_base() -> str:
    """'TRUST: python-flint X / FLINT Y' for Verdict.assumptions; warns on an unverified version."""
    import warnings
    try:
        import flint
        pv, fv = getattr(flint, "__version__", "?"), getattr(flint, "__FLINT_VERSION__", "?")
    except ImportError:
        pv, fv = "missing", "?"
    if pv != VERIFIED_PYTHON_FLINT:
        warnings.warn(f"python-flint {pv} installed; the trusted base was verified with "
                      f"{VERIFIED_PYTHON_FLINT}", RuntimeWarning, stacklevel=2)
    return f"TRUST: python-flint {pv} / FLINT {fv}"


_TAG_MASS = re.compile(r"^A-(DROP-NONPOS|NEG-CLIP)\(([-+0-9.eE]+)K\)")


_EXCL_ITEM = re.compile(r"sigma(\d+):lam=([-+0-9.eE]+):A=([-+0-9.eE]+)K")


def _sigma_exclusion_directions(tag: str, spectra: Sequence[Spectrum]) -> tuple[bool, bool]:
    """(blocks REJECTED, blocks CERTIFIED) for an A-SIGMA-EXCLUDED tag.

    The excluded smearings violate alpha^2F_nu >= 0, so they are not admissible alternatives; the
    question is only whether they point outside what the admitted alternatives already cover. Each
    admitted alternative s is treated by the engine as the box (1 +- e_s) w_j, so its lambda and A
    range over [(1-e_s) X_s, (1+e_s) X_s]. An excluded smearing whose lambda or A lies ABOVE the
    union of these ranges blocks REJECTED (it might be harder to reject); BELOW blocks CERTIFIED.
    This is a moment-level comparison (lambda, A), not a pointwise one: inside the range the
    exclusion is recorded but not blocking. A malformed tag blocks both."""
    items = _EXCL_ITEM.findall(tag)
    if not items or not spectra:
        return True, True
    rng = {}
    for name, fn in (("lam", lambda s: s.lambda_()), ("A", lambda s: s.A())):
        rng[name] = (min((1 - s.rel_envelope) * fn(s) for s in spectra),
                     max((1 + s.rel_envelope) * fn(s) for s in spectra))
    above = below = False
    for _, lam, A in items:
        for name, x in (("lam", float(lam)), ("A", float(A))):
            above |= x > rng[name][1]
            below |= x < rng[name][0]
    return above, below


def blocking_tags(spectra: Sequence[Spectrum]) -> tuple[list[str], list[str]]:
    """Input-quality tags (team DFT, Spectrum.provenance['assumptions']) that make one
    direction of proof unsafe.  Returns (tags blocking REJECTED, tags blocking CERTIFIED):
      REJECTED : A-FILE-COMPLETE(unverified), A-FILE-INCOMPLETE..., A-DROP-NONPOS(m K) with m > 0
                 (missing weight lowers A, which is anti-conservative for an upper bound);
      CERTIFIED: A-NEG-CLIP(m K) with m > 0 (clipping raises A), A-FILE-INCOMPLETE...
      either   : A-SIGMA-EXCLUDED(...) (smearings excluded by smearing_admission='mode_positive'):
                 see _sigma_exclusion_directions
    """
    rej, cert = [], []
    seen_excl: set[str] = set()
    for s in spectra:
        for a in s.provenance.get("assumptions", []) or []:
            a = str(a)
            m = _TAG_MASS.match(a)
            mass = float(m.group(2)) if m else 0.0
            tag = f"{a} [{s.label}]"
            if a.startswith("A-FILE-INCOMPLETE"):
                rej.append(tag)
                cert.append(tag)
            elif a.startswith("A-FILE-COMPLETE(unverified)"):
                rej.append(tag)
            elif m and m.group(1) == "DROP-NONPOS" and mass > 0:
                rej.append(tag)
            elif m and m.group(1) == "NEG-CLIP" and mass > 0:
                cert.append(tag)
            elif a.startswith("A-SIGMA-EXCLUDED(") and a not in seen_excl:
                seen_excl.add(a)
                r, c = _sigma_exclusion_directions(a, spectra)
                if r:
                    rej.append(tag)
                if c:
                    cert.append(tag)
    return rej, cert


def apply_input_policy(v: Verdict, spectra: Sequence[Spectrum], ccfg: dict[str, Any]) -> Verdict:
    """Coordinator policy (S1-compatible: it can only DOWNGRADE).  A proven status whose
    direction is blocked by an input-quality tag on any alternative becomes UNDECIDED, unless
    [certify] allow_unverified_rejection / allow_unverified_certification is true (then the
    tag is kept as an assumption).  The certificate is kept unchanged; the untrusted bound of the
    blocked direction is dropped from the bracket."""
    rej, cert = blocking_tags(spectra)
    reasons, assumptions = list(v.reasons), list(v.assumptions)
    status, T_low, T_up = v.status, v.T_low_K, v.T_up_K
    for kind, tags, key in (("rejection", rej, "allow_unverified_rejection"),
                            ("certification", cert, "allow_unverified_certification")):
        hit = (status is Status.REJECTED_PROVEN and kind == "rejection") or \
              (status is Status.CERTIFIED_CANDIDATE and kind == "certification")
        if not (hit and tags):
            continue
        if ccfg.get(key, False):
            assumptions.extend(f"UNVERIFIED-INPUT-ACCEPTED({key}): {t}" for t in tags)
            continue
        reasons.append(f"certify.classify returned {status.value}; {kind} blocked: "
                       + "; ".join(tags))
        if kind == "rejection":
            T_up = None
        else:
            T_low = None
        status = Status.UNDECIDED
    if status is v.status and reasons == v.reasons and assumptions == v.assumptions:
        return v
    out = dataclasses.replace(v, status=status, T_low_K=T_low, T_up_K=T_up, reasons=reasons,
                              assumptions=assumptions)
    assert out.status is v.status or out.status is Status.UNDECIDED   # downgrade only
    assert out.certificate == v.certificate
    return out


def annotate_model(v: Verdict, ccfg: dict[str, Any], spectra: Sequence[Spectrum]) -> Verdict:
    """Add the model assumptions the report must state (audit low items, F9).  Only
    ``assumptions`` and ``flags`` change."""
    assumptions = list(v.assumptions)
    flags = list(v.flags)
    seen = set()
    for s in spectra:
        rule = s.provenance.get("rule") or s.provenance.get("discretisation")
        src = s.provenance.get("source_file") or s.provenance.get("attached_file") or \
            ("inline/synthetic" if s.provenance.get("synthetic") else "inline")
        if rule:
            a = f"DISCRETISATION({str(rule).split(':')[0]}, {src})"
            if a not in seen:
                seen.add(a)
                assumptions.append(a)
    for s in spectra:                       # e.g. team DFT's A-FILE-COMPLETE(unverified)
        for a in s.provenance.get("assumptions", []) or []:
            if a not in assumptions:
                assumptions.append(a)
    has_sb = any(a.startswith("M-ISO-SINGLEBAND") for a in assumptions)
    if v.status is Status.REJECTED_PROVEN and "R-MULTI" not in v.rule_ids and not has_sb:
        assumptions.append("M-ISO-SINGLEBAND (rejection uses the DOS-averaged single-band A)")
    elif v.status in (Status.UNDECIDED, Status.CERTIFIED_CANDIDATE) and not has_sb and \
            (v.T_low_K is not None or v.T_up_K is not None):
        # Lean note MB-N1: the bracket [T_low, T_up] is the single-band (DOS-averaged) one,
        # also when rejection_model = "multiband"
        assumptions.append("M-ISO-SINGLEBAND (the bracket [T_low, T_up] is the single-band "
                           "DOS-averaged one)")
    # F9: is the Coulomb block empty at the certified lower bound?  N_c(T) = 0 iff pi*T > omega_c
    oc = (v.certificate or {}).get("params", {}).get("omega_c")
    if oc is None:
        oc = (float(ccfg["omega_c_K"]) if ccfg.get("omega_c_K") is not None else
              float(ccfg.get("omega_c_factor", 10.0)) * max(s.omega_max() for s in spectra))
    if v.T_low_K is not None and math.pi * float(v.T_low_K) > float(oc) and \
            "COULOMB_INACTIVE" not in flags:
        flags.append("COULOMB_INACTIVE")
        assumptions.append(f"COULOMB_INACTIVE: N_c(T_low={v.T_low_K} K) = 0 for omega_c = {oc} K; "
                           "the lower bound does not depend on mu*")
    out = dataclasses.replace(v, assumptions=assumptions, flags=flags)
    assert (out.status, out.T_low_K, out.T_up_K, out.certificate) == \
        (v.status, v.T_low_K, v.T_up_K, v.certificate)
    return out


def cutoff_description(ccfg: dict[str, Any], spectra: Sequence[Spectrum] = ()) -> str:
    """ARCH 4b F3: the Coulomb block-cutoff rule and omega_c used (config [certify])."""
    rule = ccfg.get("cutoff_rule", "block")
    if ccfg.get("omega_c_K") is not None:
        oc = f"omega_c={float(ccfg['omega_c_K']):g} K (fixed)"
    else:
        f = float(ccfg.get("omega_c_factor", 10.0))
        vals = sorted({round(f * s.omega_max(), 6) for s in spectra})
        oc = f"omega_c={f:g}*max(omega) per spectrum" + (f" = {vals} K" if vals else "")
    return f"CUTOFF({rule}: N_c(T)=#{{n>=0 : (2n+1)*pi*T <= omega_c}}; {oc})"


def ensure_cutoff_assumption(v: Verdict, ccfg: dict[str, Any], spectra: Sequence[Spectrum]) -> Verdict:
    """Make sure Verdict.assumptions names the cutoff rule and omega_c (ARCH 4b F3).  Only the
    assumptions list is extended; status, bounds and certificate are untouched."""
    if any("omega_c" in a for a in v.assumptions):
        return v
    out = dataclasses.replace(v, assumptions=list(v.assumptions) + [cutoff_description(ccfg, spectra)])
    assert (out.status, out.T_low_K, out.T_up_K) == (v.status, v.T_low_K, v.T_up_K)
    return out


# ----------------------------------------------------------------------------- per-material work
def compute_material(cand_dict: dict[str, Any], rec: dict[str, Any], cfg: dict[str, Any],
                     workdir: str) -> dict[str, Any]:
    """Stages 6-8 for one material.  Runs in a worker process; returns the finished record."""
    c = Candidate.from_dict(cand_dict)
    mat_dir = Path(workdir) / safe_id(c.material_id)
    try:
        return _compute_material(c, rec, cfg, mat_dir)
    except Exception as exc:  # noqa: BLE001 - every failure becomes an explicit ERROR (S4)
        tb = traceback.format_exc()
        rec.setdefault("stages", {})["error"] = {"exception": f"{type(exc).__name__}: {exc}",
                                                 "traceback": tb}
        v = _pipeline_verdict("error", c.material_id, Status.ERROR,
                              [f"{type(exc).__name__}: {exc}"],
                              certificate={"kind": "error", "proved": False})
        _finish(rec, v, "pipeline.error")
        _atomic_write_json(mat_dir / "record.json", rec)
        return rec


def _compute_material(c: Candidate, rec: dict[str, Any], cfg: dict[str, Any], mat_dir: Path):
    scfg = cfg["spectra"]
    stages = rec.setdefault("stages", {})
    spectra: list[Spectrum] = []
    if rec.get("spectra"):
        spectra = [Spectrum.from_dict(d) for d in rec["spectra"]]
        stages.setdefault("spectra", {})["resumed"] = True
    else:
        info: dict[str, Any] = {"sources": []}
        if scfg["use_attached"]:
            att = attached_spectra(c, scfg)
            if att:
                spectra.extend(att)
                info["sources"].append(f"attached ({len(att)})")
        dcfg = dict(cfg["dfpt"])
        if dcfg.get("enabled") and (not spectra or dcfg.get("also_when_attached", False)):
            dsp, dinfo = run_dfpt(c, mat_dir, {k: v for k, v in dcfg.items()
                                              if k not in ("enabled", "also_when_attached")})
            stages["dfpt"] = dinfo
            info["sources"].append(f"dfpt ({len(dsp)})")
            if not dsp:
                # ARCH 4b F1: every DFPT alternative failed -> ERROR; never reaches classify
                stages["spectra"] = info
                v = _pipeline_verdict("error", c.material_id, Status.ERROR,
                                      ["DFPT produced zero alternative spectra (all failed); "
                                       "no verdict about the model is possible (ARCH 4b F1)"],
                                      certificate={"kind": "error", "proved": False,
                                                   "errata": "F1"})
                _finish(rec, v, "pipeline.dfpt")
                _atomic_write_json(mat_dir / "record.json", rec)
                return rec
            fatal = [f for f in dinfo.get("flags", []) if str(f).startswith("DFPT_CHECK_FAIL:")
                     and str(f).split(":")[1] in FATAL_DFPT_CHECKS]
            if fatal:
                # audit F8: a parse/units consistency failure must never reach a proof
                stages["spectra"] = info
                v = _pipeline_verdict("error", c.material_id, Status.ERROR,
                                      ["DFT consistency check indicates a parse/units error: "
                                       + ", ".join(map(str, fatal))],
                                      certificate={"kind": "error", "proved": False,
                                                   "dfpt_flags": list(map(str, fatal))})
                _finish(rec, v, "pipeline.dfpt")
                _atomic_write_json(mat_dir / "record.json", rec)
                return rec
            subset = (dcfg.get("settings", {}) or {}).get("use_smearings") is not None or any(
                "only smearing indices" in str(w) for w in dinfo.get("warnings", []))
            if subset:
                dinfo["smearings_subset"] = True
            spectra.extend(dsp)
            phonon = _phonon_from_dfpt(dinfo, dsp)
            if phonon is not None:
                fs = FilterSet.from_config(cfg["filters"])
                v_ex, fr = fs.evaluate_dynamic(c, phonon)
                if fr is not None:
                    stages.setdefault("filters", []).append(fr.to_dict())
                if v_ex is not None:
                    rec["spectra"] = [s.to_dict() for s in spectra]
                    _finish(rec, v_ex, "rtsc.stability.dynamic")
                    _atomic_write_json(mat_dir / "record.json", rec)
                    return rec
        if not spectra:      # DFPT disabled and nothing attached: simply not computed
            stages["spectra"] = info
            v = _pipeline_verdict("spectra", c.material_id, Status.DEFERRED,
                                  ["not computed: no alpha^2F attached and DFPT disabled"])
            _finish(rec, v, "pipeline.spectra")
            _atomic_write_json(mat_dir / "record.json", rec)
            return rec
        if c.synthetic:
            spectra = _mark_synthetic(spectra)
        spectra = _with_envelope(spectra, float(scfg["rel_envelope"]))
        info["n_spectra"] = len(spectra)
        info["rel_envelope_applied"] = float(scfg["rel_envelope"])
        info["summary"] = [s.summary() for s in spectra]
        stages["spectra"] = info
        rec["spectra"] = [s.to_dict() for s in spectra]
        _atomic_write_json(mat_dir / "record.json", rec)     # checkpoint before certify

    # stage 7: the only place REJECTED_PROVEN / CERTIFIED_CANDIDATE can arise
    ccfg = dict(cfg["certify"])
    stages["certify"] = {"params": ccfg, "started_utc": _now()}
    bcs = load_band_couplings(c)
    if bcs:
        stages["certify"]["band_coupling_files"] = list(c.band_coupling_files)
    verdict = certify_spectra(c.material_id, spectra, ccfg, band_couplings=bcs)
    stages["certify"]["finished_utc"] = _now()
    origin = CLASSIFY_ORIGIN

    # stage 8: validity flags (never change status, S6)
    vcfg = ValidityConfig.from_config(cfg["validity"])
    pmin = None
    dyn = [f for f in stages.get("filters", []) if f.get("criterion") == "X-DYNAMIC"
           and f.get("details", {}).get("lowest_frequency") is not None]
    if dyn:
        d = dyn[-1]["details"]
        from .spectrum import unit_factor
        pmin = d["lowest_frequency"] * unit_factor(d.get("unit", "THz")) / unit_factor("THz")
    flags, unchecked = compute_flags(c, spectra, vcfg, mu_hi=ccfg.get("mu_hi"), phonon_min_THz=pmin)
    from .validity import Flag
    for f in (stages.get("dfpt") or {}).get("flags", []):   # team DFT's diagnostics
        flags.append(Flag(str(f), "reported by rtsc.dfpt"))
    if (stages.get("dfpt") or {}).get("smearings_subset"):
        flags.append(Flag("SMEARINGS_SUBSET", "only a subset of the DFPT smearings was used as "
                          "alternative spectra (explicit use_smearings or smearing_admission)"))
    warn = [w for s in spectra for w in s.provenance.get("warnings", [])]
    if warn:
        flags.append(Flag("SPECTRUM_WARNINGS", f"{len(warn)} input warnings (see record)",
                          {"warnings": warn[:20]}))
    verdict = apply_flags(verdict, flags)
    stages["validity"] = {"flags": [f.to_dict() for f in flags], "unchecked": unchecked}
    _finish(rec, verdict, origin)
    _atomic_write_json(mat_dir / "record.json", rec)
    return rec


# ----------------------------------------------------------------------------- run summary (S4)
@dataclass
class RunSummary:
    inputs: list[str]
    verdicts: dict[str, Verdict]
    origins: dict[str, str] = field(default_factory=dict)
    duplicate_of: dict[str, str] = field(default_factory=dict)
    config_hash: str = ""
    started_utc: str = ""
    finished_utc: str = ""
    shard: str | None = None

    def by_status(self) -> dict[str, list[str]]:
        out: dict[str, list[str]] = {s.value: [] for s in Status}
        for mid in self.inputs:
            out[self.verdicts[mid].status.value].append(mid)
        return out

    def counts(self) -> dict[str, int]:
        return {k: len(v) for k, v in self.by_status().items()}

    def to_dict(self) -> dict[str, Any]:
        return {"n_inputs": len(self.inputs), "counts": self.counts(),
                "by_status": self.by_status(), "duplicate_of": dict(self.duplicate_of),
                "config_hash": self.config_hash, "started_utc": self.started_utc,
                "finished_utc": self.finished_utc, "shard": self.shard,
                "accounting": "every input has exactly one status (S4 checked)"}


def account(inputs: Sequence[str], verdicts: dict[str, Verdict],
            origins: dict[str, str] | None = None) -> None:
    """Raise AccountingError unless every input has exactly one verdict (S4) and proven
    statuses come only from certify.classify (S1)."""
    seen: set[str] = set()
    dup: list[str] = []
    for m in inputs:
        if m in seen:
            dup.append(m)
        seen.add(m)
    if dup:
        raise AccountingError(f"duplicate input ids: {dup[:10]}")
    missing = [m for m in inputs if m not in verdicts]
    extra = [m for m in verdicts if m not in seen]
    if missing or extra:
        raise AccountingError(f"S4 violated: {len(missing)} inputs without status {missing[:10]}, "
                              f"{len(extra)} unknown ids {extra[:10]}")
    for m, v in verdicts.items():
        if not isinstance(v.status, Status):
            raise AccountingError(f"{m}: status {v.status!r} is not a Status")
        if v.material_id != m:
            raise AccountingError(f"{m}: verdict carries material_id {v.material_id!r}")
        if origins is not None and v.status in PROVEN_STATUSES and origins.get(m) != CLASSIFY_ORIGIN:
            raise AccountingError(f"S1 violated: {m} is {v.status.value} from {origins.get(m)}")
    total = sum(1 for _ in inputs)
    if total != len(verdicts):
        raise AccountingError("S4 violated: counts do not add up")


# ----------------------------------------------------------------------------- the run
@dataclass
class Plan:
    """Output of stages 1-5: records for every input, and the ordered ids to compute."""
    candidates: dict[str, Candidate]
    inputs: list[str]
    records: dict[str, dict[str, Any]]
    to_compute: list[str]
    duplicate_of: dict[str, str]
    config_hash: str


def plan(cfg: dict[str, Any], candidates: list[Candidate] | None = None, *,
         ml_hook: Callable[[Candidate], dict[str, Any]] | None = None,
         rester_factory=None) -> Plan:
    chash = config_hash(cfg)
    workdir = cfg["run"]["workdir"]
    if candidates is None:
        candidates = load_sources(cfg["sources"], rester_factory=rester_factory)
    else:
        from .sources import ensure_unique_ids
        ensure_unique_ids(candidates)
    inputs = [c.material_id for c in candidates]
    cmap = {c.material_id: c for c in candidates}
    check_path_collisions(inputs)
    records: dict[str, dict[str, Any]] = {}
    for c in candidates:
        digest = input_digest(c)
        old = read_record(workdir, c.material_id) if cfg["run"]["resume"] else None
        if record_reusable(old, c, chash, digest):
            records[c.material_id] = old
        else:
            records[c.material_id] = new_record(c, chash, digest)

    def done(mid: str) -> bool:
        return records[mid].get("verdict") is not None

    def input_error(mid: str, stage: str, exc: BaseException) -> None:
        """Audit F6: a malformed property ends THIS input as ERROR, never the whole run."""
        records[mid]["stages"].setdefault("errors", []).append(
            {"stage": stage, "exception": f"{type(exc).__name__}: {exc}",
             "traceback": traceback.format_exc()})
        if not done(mid):
            v = _pipeline_verdict("error", mid, Status.ERROR,
                                  [f"malformed input at stage {stage}: {type(exc).__name__}: {exc}"])
            _finish(records[mid], v, f"pipeline.{stage}")

    # stage 1: unreadable inputs -> ERROR
    for c in candidates:
        if done(c.material_id):
            continue
        if c.properties.get("load_error"):
            v = _pipeline_verdict("error", c.material_id, Status.ERROR,
                                  [f"input could not be read: {c.properties['load_error']}"])
            _finish(records[c.material_id], v, "pipeline.source")

    # stage 2: dedupe -> DEFERRED for duplicates
    live = [c for c in candidates if not c.properties.get("load_error")]
    dcfg = dict(cfg["dedupe"])
    dres = dedupe(live, **dcfg)
    for dup, rep in dres.duplicate_of.items():
        records[dup]["stages"]["dedupe"] = {"duplicate_of": rep, "settings": dcfg}
        if not done(dup):
            v = _pipeline_verdict("dedupe", dup, Status.DEFERRED,
                                  [f"duplicate of {rep} (dedupe method "
                                   f"{dres.method_used.get(rep, '?')}); not computed separately"],
                                  certificate={"kind": "dedupe", "proved": False, "duplicate_of": rep,
                                               "settings": dcfg})
            _finish(records[dup], v, "pipeline.dedupe")
    # Stages 3-5 are evaluated for every surviving input, including resumed ones, so that the
    # plan (order, budget, shard assignment) does not depend on what has already finished.
    survivors = [c for c in dres.kept if not c.properties.get("load_error")]

    # stage 3: existence filters -> EXCLUDED_ASSUMPTION
    fs = FilterSet.from_config(cfg["filters"])
    fs.prepare(candidates)
    passed: list[Candidate] = []
    for c in survivors:
        try:
            v_ex, results = fs.evaluate_static(c)
        except Exception as exc:  # noqa: BLE001
            input_error(c.material_id, "filters", exc)
            continue
        if not done(c.material_id):
            records[c.material_id]["stages"]["filters"] = [r.to_dict() for r in results]
            if v_ex is not None:
                _finish(records[c.material_id], v_ex, "rtsc.stability")
        if v_ex is None:
            passed.append(c)

    # stage 5 (computed before ordering so it can inform it): optional ML pre-check
    mlcfg = dict(cfg["ml"])
    hook = ml_hook or get_ml_hook(mlcfg)
    ml_results: dict[str, dict[str, Any]] = {}
    if hook is not None:
        for c in passed:
            cached = records[c.material_id]["stages"].get("ml")
            try:
                ml_results[c.material_id] = cached if cached else hook(c)
            except Exception as exc:  # noqa: BLE001 - heuristic: losing it only loses ordering
                ml_results[c.material_id] = {"ok": False, "error": f"{type(exc).__name__}: {exc}"}
            records[c.material_id]["stages"]["ml"] = ml_results[c.material_id]

    # stage 4: priority + budget -> DEFERRED
    weights = PriorityWeights.from_config({**{"T_target_K": float(cfg["certify"].get("T_target_K", 300.0))},
                                           **cfg["priority"]})
    scores = []
    for c in passed:
        try:
            scores.append(priority(c, weights, ml_results.get(c.material_id)))
        except Exception as exc:  # noqa: BLE001
            input_error(c.material_id, "priority", exc)
    bcfg = cfg["budget"]
    bplan = plan_budget(scores, max_materials=bcfg.get("max_materials"),
                        min_score=bcfg.get("min_score"), ml_results=ml_results,
                        defer_if_ml_unstable=bool(mlcfg.get("defer_if_ml_unstable", False)))
    for s in scores:
        if not done(s.material_id):
            records[s.material_id]["stages"]["priority"] = s.to_dict()
    for mid, why in bplan.deferred.items():
        if done(mid):
            continue
        stage = "ml" if "ML pre-check" in why else "budget"
        v = _pipeline_verdict(stage, mid, Status.DEFERRED, [why],
                              certificate={"kind": stage, "proved": False,
                                           "budget": dict(bcfg),
                                           "priority": bplan.scores[mid].to_dict()})
        _finish(records[mid], v, f"pipeline.{stage}")
    for rank, mid in enumerate(bplan.order):
        if not done(mid):
            records[mid]["stages"]["priority"]["rank"] = rank
    for mid in inputs:   # synthetic inputs carry the SYNTHETIC_DATA flag whatever their status
        r = records[mid]
        if cmap[mid].synthetic and r.get("verdict") and "SYNTHETIC_DATA" not in r["verdict"]["flags"]:
            r["verdict"]["flags"].append("SYNTHETIC_DATA")
    return Plan(cmap, inputs, records, list(bplan.order), dres.duplicate_of, chash)


def check_path_collisions(ids: Sequence[str]) -> None:
    """Two inputs must never share <workdir>/<safe_id>/ (audit F3).  Compared case-insensitively
    because macOS/Windows file systems are."""
    seen: dict[str, str] = {}
    clashes = []
    for m in ids:
        k = safe_id(m).lower()
        if k in seen and seen[k] != m:
            clashes.append((seen[k], m))
        seen.setdefault(k, m)
    if clashes:
        raise ConfigError(f"material ids map to the same record directory: {clashes[:10]}; "
                          "rename them")


def _inputs_file(workdir: Path, sh: tuple[int, int] | None) -> Path:
    return workdir / ("inputs.json" if sh is None else f"inputs_shard_{sh[0]}_of_{sh[1]}.json")


def _parse_shard(shard: str | None) -> tuple[int, int] | None:
    if not shard:
        return None
    k, n = (int(x) for x in shard.split("/"))
    if not (n >= 1 and 0 <= k < n):
        raise ConfigError(f"bad shard {shard!r}; use K/N with 0 <= K < N")
    return k, n


def run(cfg: dict[str, Any], candidates: list[Candidate] | None = None, *,
        ml_hook: Callable[[Candidate], dict[str, Any]] | None = None, rester_factory=None,
        shard: str | None = None, plan_only: bool = False, write_reports: bool = True
        ) -> RunSummary:
    """Run all stages.  ``shard='K/N'`` computes only every N-th material of the global order
    (cluster job arrays); ``rtsc collect`` / ``collect()`` then merges all records."""
    started = _now()
    workdir = Path(cfg["run"]["workdir"])
    workdir.mkdir(parents=True, exist_ok=True)
    sh = _parse_shard(shard)
    p = plan(cfg, candidates, ml_hook=ml_hook, rester_factory=rester_factory)
    for f in workdir.glob("inputs*.json"):   # metadata of runs under another configuration
        try:
            if json.loads(f.read_text()).get("config_hash") != p.config_hash:
                f.unlink()
        except (OSError, json.JSONDecodeError):
            pass
    _atomic_write_json(_inputs_file(workdir, sh),
                       {"config_hash": p.config_hash, "inputs": p.inputs, "to_compute": p.to_compute,
                        "duplicate_of": p.duplicate_of, "shard": shard})
    order_all = p.to_compute + [m for m in p.inputs if m not in set(p.to_compute)]
    mine = set(order_all if sh is None else order_all[sh[0]::sh[1]])
    for mid in p.inputs:
        rec = p.records[mid]
        if mid in mine and rec.get("verdict") is not None:
            _atomic_write_json(record_path(workdir, mid), rec)
    if plan_only:
        for mid in p.to_compute:
            if p.records[mid].get("verdict") is None:
                _atomic_write_json(record_path(workdir, mid), p.records[mid])
        return _summary(p, {m: p.records[m] for m in p.inputs if p.records[m].get("verdict")},
                        started, shard, check=False)
    todo = [m for m in p.to_compute if m in mine and p.records[m].get("verdict") is None]
    workers = max(1, int(cfg["run"]["workers"]))
    results: dict[str, dict[str, Any]] = {}
    args = [(p.candidates[m].to_dict(), p.records[m], cfg, str(workdir)) for m in todo]
    if workers == 1 or len(todo) <= 1:
        for a in args:
            mid = a[0]["material_id"]
            r = compute_material(*a)
            if r.get("material_id") != mid:
                raise AccountingError(f"record for {mid!r} came back as {r.get('material_id')!r}")
            results[mid] = r
    else:
        with cf.ProcessPoolExecutor(max_workers=workers) as ex:
            futs = {ex.submit(compute_material, *a): a[0]["material_id"] for a in args}
            for fut in cf.as_completed(futs):
                mid = futs[fut]
                try:
                    r = fut.result()
                    if r.get("material_id") != mid:
                        raise AccountingError(f"record for {mid!r} came back as "
                                              f"{r.get('material_id')!r}")
                except AccountingError:
                    raise
                except Exception as exc:  # noqa: BLE001 - worker crash -> ERROR (S4), retryable
                    # audit F7: a dead worker (OOM kill, segfault) breaks the whole pool; the
                    # pending materials did not fail themselves.  Their ERROR is reported for
                    # this run but marked retryable, so a resumed run computes them again.
                    rec = p.records[mid]
                    v = _pipeline_verdict("error", mid, Status.ERROR,
                                          [f"worker process failed ({type(exc).__name__}: {exc}); "
                                           "not a property of the material; retried on resume"],
                                          certificate={"kind": "error", "proved": False,
                                                       "retryable": True})
                    r = _finish(rec, v, "pipeline.worker")
                    r["retryable"] = True
                    _atomic_write_json(record_path(workdir, mid), r)
                results[mid] = r
    for mid, r in results.items():
        p.records[mid] = r
    scope = [m for m in p.inputs if m in mine]
    summ = _summary(p, {m: p.records[m] for m in scope}, started, shard, check=True, scope=scope)
    if write_reports and sh is None:
        from .report import write_reports as _wr
        _wr(cfg, summ, [p.records[m] for m in p.inputs])
    elif sh is not None:
        _atomic_write_json(workdir / f"shard_{sh[0]}_of_{sh[1]}.summary.json", summ.to_dict())
    return summ


def _summary(p: Plan, recs: dict[str, dict[str, Any]], started: str, shard: str | None, *,
             check: bool, scope: list[str] | None = None) -> RunSummary:
    verdicts = {m: Verdict.from_dict(r["verdict"]) for m, r in recs.items() if r.get("verdict")}
    origins = {m: r.get("verdict_origin") for m, r in recs.items()}
    scope = scope if scope is not None else [m for m in p.inputs if m in verdicts]
    s = RunSummary(scope, verdicts, origins, {k: v for k, v in p.duplicate_of.items() if k in scope},
                   p.config_hash, started, _now(), shard)
    if check:
        account(scope, verdicts, origins)
    return s


def collect(cfg: dict[str, Any], *, write: bool = True) -> RunSummary:
    """Merge all per-material records of a (sharded) run.

    Refuses (ConfigError) when the run's input lists were written under another configuration
    (audit F2).  The inputs are the union of ``inputs.json`` and every
    ``inputs_shard_K_of_N.json``.  A missing or unfinished record, or a record made under another
    configuration, becomes an explicit ERROR, so S4 still holds."""
    workdir = Path(cfg["run"]["workdir"])
    chash = config_hash(cfg)
    metas = []
    for f in sorted(workdir.glob("inputs*.json")):
        m = json.loads(f.read_text())
        if m.get("config_hash") != chash:
            raise ConfigError(f"{f} was written under configuration {m.get('config_hash')}, the "
                              f"current one is {chash}: refusing to merge (re-run instead)")
        metas.append(m)
    if not metas:
        raise FileNotFoundError(f"no inputs*.json in {workdir}: nothing to collect")
    inputs: list[str] = []
    dup_of: dict[str, str] = {}
    for m in metas:
        for mid in m["inputs"]:
            if mid not in inputs:
                inputs.append(mid)
        dup_of.update(m.get("duplicate_of", {}))
    records, verdicts, origins = [], {}, {}
    for mid in inputs:
        rec = read_record(workdir, mid)
        why = None
        if rec is None or not rec.get("verdict"):
            why = "no finished record: shard not run, still running, or crashed"
        elif rec.get("config_hash") != chash:
            why = (f"record from another configuration (config hash {rec.get('config_hash')} "
                   f"!= {chash}); not merged")
        if why:
            rec = {"schema": RECORD_SCHEMA, "material_id": mid, "stages": {}, "config_hash": chash,
                   "candidate": (rec or {}).get("candidate", {"material_id": mid, "composition": {}})}
            _finish(rec, _pipeline_verdict("error", mid, Status.ERROR, [why]), "pipeline.collect")
        records.append(rec)
        verdicts[mid] = Verdict.from_dict(rec["verdict"])
        origins[mid] = rec.get("verdict_origin")
    s = RunSummary(list(inputs), verdicts, origins, dup_of, chash, "", _now(), "collected")
    account(inputs, verdicts, origins)
    if write:
        from .report import write_reports as _wr
        _wr(cfg, s, records)
    return s
