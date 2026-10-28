"""MatterGen generator adapter (``rtsc.discovery.interfaces.Generator``), CLI subprocess only.

MatterGen runs in its OWN environment (Linux, CUDA GPU, Python 3.10 per the official README; its
PyPI metadata pins ``torch==2.2.1+cu118`` and ``numpy<2``), which cannot share an interpreter with
MatterSim (Python >= 3.12, numpy >= 2) or with rtsc's main py3.12 environment.  So this adapter
never imports ``mattergen``: it calls the ``mattergen-generate`` executable through a configurable
command (absolute path, ``conda run -n <env>``, ``srun -p gpu ...`` prefix, ``ssh host`` + remote
shell, ...) and parses the files it writes.

Verified against the official sources (github.com/microsoft/mattergen @ 92423660, 2026-08-27;
README, ``mattergen/scripts/generate.py``, ``checkpoints/*/config.yaml``,
``mattergen/common/utils/globals.py``; PyPI ``mattergen`` 1.0.3, MIT licence):

* CLI ``mattergen-generate OUTPUT_PATH`` (python-fire over ``generate.main``) with
  ``--pretrained_name=<name>`` XOR ``--model_path=<dir>``, ``--batch_size`` (default 64),
  ``--num_batches`` (default 1), ``--properties_to_condition_on={...}`` (python dict literal),
  ``--diffusion_guidance_factor`` (classifier-free guidance gamma; 0 = unconditional; README uses
  2.0), ``--seed``, ``--record_trajectories`` (default True), ``--checkpoint_epoch`` (last|best|int),
  ``--sampling_config_overrides``, ``--target_compositions`` (CSP-mode models only).
* Pretrained checkpoints (Hugging Face ``microsoft/mattergen``, files
  ``checkpoints/<name>/checkpoints/last.ckpt`` + ``config.yaml``; also Git LFS in the repo) and the
  properties each one is conditioned on (``property_embeddings_adapt`` in its config.yaml):
  see :data:`CHECKPOINT_PROPERTIES`.
* Outputs in OUTPUT_PATH: ``generated_crystals.extxyz`` (one frame per structure),
  ``generated_crystals_cif.zip`` (``gen_<i>.cif``), and ``generated_trajectories.zip`` when
  trajectories are recorded.
* Elements MatterGen can generate: ``SELECTED_ATOMIC_NUMBERS`` (76 elements; no noble gases,
  Tc, Pm or Z >= 84) -- see :data:`MATTERGEN_ELEMENTS`.
* Chemical systems are passed as ``"El1-El2-..."`` strings; space groups as int 1..230;
  ``energy_above_hull`` (eV/atom), ``dft_band_gap`` (eV), ``ml_bulk_modulus`` (GPa),
  ``dft_mag_density`` (mu_B/A^3), ``hhi_score`` are floats.

PRESSURE: MatterGen has no pressure conditioning.  Its training data (MP-20 / Alex-MP-20) are
ambient-pressure DFT relaxations, so generated structures are ambient-like proposals; the
prescreener (``rtsc.mlff.mattersim``) relaxes them at the request's target pressure.  Recorded
in every structure's provenance.

Every generated structure is HYPOTHETICAL_STRUCTURE (provenance flag).
"""
from __future__ import annotations

import hashlib
import io
import json
import math
import os
import re
import shlex
import subprocess
import time
import zipfile
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Sequence

from ..discovery.interfaces import GeneratedStructure, GenerationRequest

GENERATOR_NAME = "mattergen"
EXTXYZ_NAME = "generated_crystals.extxyz"
CIF_ZIP_NAME = "generated_crystals_cif.zip"
PRESSURE_NOTE = ("MatterGen has no pressure conditioning: structure sampled from an ambient-pressure "
                 "training distribution; relaxed at the target pressure by the ML prescreener")

# property_embeddings_adapt of each official checkpoint (checkpoints/<name>/config.yaml @ 92423660)
CHECKPOINT_PROPERTIES: dict[str, frozenset[str]] = {
    "mattergen_base": frozenset(),
    "mp_20_base": frozenset(),
    "chemical_system": frozenset({"chemical_system"}),
    "space_group": frozenset({"space_group"}),
    "dft_mag_density": frozenset({"dft_mag_density"}),
    "dft_band_gap": frozenset({"dft_band_gap"}),
    "ml_bulk_modulus": frozenset({"ml_bulk_modulus"}),
    "dft_mag_density_hhi_score": frozenset({"dft_mag_density", "hhi_score"}),
    "chemical_system_energy_above_hull": frozenset({"chemical_system", "energy_above_hull"}),
}

# PROPERTY_SOURCE_IDS in mattergen/common/utils/globals.py (names usable in datasets / fine-tuning)
PROPERTY_SOURCE_IDS = ("dft_mag_density", "dft_bulk_modulus", "dft_shear_modulus", "energy_above_hull",
                       "formation_energy_per_atom", "space_group", "hhi_score", "ml_bulk_modulus",
                       "chemical_system", "dft_band_gap")

_SYMBOLS = ("X H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge "
            "As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm "
            "Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi").split()
_SELECTED_Z = [1, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 19, 20, 21, 22, 23, 24, 25, 26, 27,
               28, 29, 30, 31, 32, 33, 34, 35, 37, 38, 39, 40, 41, 42, 44, 45, 46, 47, 48, 49, 50, 51,
               52, 53, 55, 56, 57, 58, 59, 60, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75,
               76, 77, 78, 79, 80, 81, 82, 83]
MATTERGEN_ELEMENTS = frozenset(_SYMBOLS[z] for z in _SELECTED_Z)


class UnsupportedConditioningError(ValueError):
    """The request asks for conditioning the configured checkpoint was not trained for."""


class MatterGenRunError(RuntimeError):
    pass


@dataclass
class MatterGenConfig:
    pretrained_name: str | None = "chemical_system"     # one of CHECKPOINT_PROPERTIES, or None
    model_path: str | None = None                        # fine-tuned / own model directory instead
    # properties supported by ``model_path``; None => read property_embeddings_adapt from
    # <model_path>/config.yaml (needs PyYAML) and refuse if that is impossible
    model_properties: Sequence[str] | None = None
    checkpoint_epoch: str = "last"
    executable: Sequence[str] = ("mattergen-generate",)  # e.g. ["conda","run","-n","mattergen","mattergen-generate"]
    command_prefix: Sequence[str] = ()                    # e.g. ["srun","-p","gpu","--gres=gpu:1","-t","60"]
    shell: bool = False             # True: join argv with shlex.quote into ONE string (ssh-style prefixes)
    env: dict[str, str] = field(default_factory=dict)
    batch_size: int = 16
    max_batches: int = 64
    guidance_factor: float = 2.0                          # used only when conditioning is non-empty
    record_trajectories: bool = False
    sampling_config_overrides: Sequence[str] = ()
    timeout_s: float = 6 * 3600
    workdir: str = "mattergen_runs"
    reuse_outputs: bool = True                            # same invocation + finished outputs => no rerun
    # request fields the checkpoint cannot condition on that the user ALLOWS to be dropped
    # (e.g. ["space_group"]).  Dropped fields are recorded; chemical_system is then enforced by
    # post-filtering (elements subset).  Anything else unsupported raises.
    allow_unconditioned: Sequence[str] = ()
    postfilter_chemical_system: bool = True

    @classmethod
    def from_dict(cls, d: dict[str, Any] | None) -> "MatterGenConfig":
        return cls(**dict(d or {}))


def _model_properties(cfg: MatterGenConfig) -> frozenset[str]:
    if (cfg.pretrained_name is None) == (cfg.model_path is None):
        raise ValueError("MatterGen: set exactly one of pretrained_name or model_path")
    if cfg.pretrained_name is not None:
        if cfg.pretrained_name not in CHECKPOINT_PROPERTIES:
            raise ValueError(f"unknown MatterGen pretrained_name {cfg.pretrained_name!r}; "
                             f"official names: {sorted(CHECKPOINT_PROPERTIES)}")
        return CHECKPOINT_PROPERTIES[cfg.pretrained_name]
    if cfg.model_properties is not None:
        return frozenset(cfg.model_properties)
    cpath = Path(cfg.model_path) / "config.yaml"
    try:
        import yaml
        d = yaml.safe_load(cpath.read_text())
        adapt = d["lightning_module"]["diffusion_module"]["model"].get("property_embeddings_adapt") or {}
        return frozenset(adapt)
    except Exception as exc:  # noqa: BLE001
        raise UnsupportedConditioningError(
            f"cannot determine the conditioning properties of model_path={cfg.model_path} "
            f"({type(exc).__name__}: {exc}); set model_properties explicitly") from exc


def conditioning_for(request: GenerationRequest, cfg: MatterGenConfig) -> tuple[dict[str, Any], list[str]]:
    """Map a request onto ``properties_to_condition_on``.  Returns (conditioning, dropped fields).

    Raises UnsupportedConditioningError for anything the checkpoint does not support unless the
    field is listed in ``allow_unconditioned``; never silently ignores a request field."""
    supported = _model_properties(cfg)
    bad_el = sorted(set(request.chemical_system) - MATTERGEN_ELEMENTS)
    if bad_el:
        raise UnsupportedConditioningError(f"MatterGen cannot generate element(s) {bad_el} "
                                           "(not in SELECTED_ATOMIC_NUMBERS)")
    wanted: dict[str, Any] = {}
    if request.chemical_system:
        wanted["chemical_system"] = "-".join(sorted(request.chemical_system))
    if request.space_group is not None:
        sg = int(request.space_group)
        if not 1 <= sg <= 230:
            raise UnsupportedConditioningError(f"space group {sg} not in 1..230")
        wanted["space_group"] = sg
    extra = dict(request.extra or {}).get("mattergen_conditioning", {})
    for k, v in extra.items():
        if k not in PROPERTY_SOURCE_IDS and k not in supported:
            raise UnsupportedConditioningError(f"{k!r} is not a MatterGen property id")
        wanted[k] = v
    cond, dropped = {}, []
    for k, v in wanted.items():
        if k in supported:
            cond[k] = v
        elif k in cfg.allow_unconditioned:
            dropped.append(k)
        else:
            model = cfg.pretrained_name or cfg.model_path
            raise UnsupportedConditioningError(
                f"MatterGen checkpoint {model!r} cannot condition on {k!r} (supports "
                f"{sorted(supported) or 'nothing: unconditional'}); choose a checkpoint that does "
                f"or list {k!r} in allow_unconditioned")
    return cond, dropped


def _literal(d: dict[str, Any]) -> str:
    """python-fire dict literal without whitespace (README: no whitespace between key and value)."""
    parts = []
    for k, v in d.items():
        parts.append(f"'{k}':" + (f"'{v}'" if isinstance(v, str) else repr(v)))
    return "{" + ",".join(parts) + "}"


def build_command(request: GenerationRequest, cfg: MatterGenConfig, outdir: str | Path,
                  cond: dict[str, Any], batch_size: int, num_batches: int) -> list[str]:
    argv = list(cfg.command_prefix) + list(cfg.executable) + [str(outdir)]
    if cfg.pretrained_name is not None:
        argv.append(f"--pretrained_name={cfg.pretrained_name}")
    else:
        argv += [f"--model_path={cfg.model_path}", f"--checkpoint_epoch={cfg.checkpoint_epoch}"]
    argv += [f"--batch_size={batch_size}", f"--num_batches={num_batches}", f"--seed={int(request.seed)}",
             f"--record_trajectories={bool(cfg.record_trajectories)}"]
    if cond:
        argv += [f"--properties_to_condition_on={_literal(cond)}",
                 f"--diffusion_guidance_factor={float(cfg.guidance_factor)}"]
    if cfg.sampling_config_overrides:
        argv.append("--sampling_config_overrides=[" +
                    ",".join(repr(s) for s in cfg.sampling_config_overrides) + "]")
    return argv


def _safe(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.=-]+", "_", s)[:120]


def parse_outputs(outdir: str | Path) -> tuple[list[dict[str, Any]], str]:
    """Structures written by mattergen-generate: extxyz preferred, CIF zip as fallback."""
    from ..sources.fileformats import read_cif, read_extxyz
    outdir = Path(outdir)
    xyz = outdir / EXTXYZ_NAME
    raws = []
    if xyz.exists():
        raws, src = read_extxyz(xyz.read_text()), EXTXYZ_NAME
    elif (outdir / CIF_ZIP_NAME).exists():
        src = CIF_ZIP_NAME
        with zipfile.ZipFile(outdir / CIF_ZIP_NAME) as z:
            names = sorted(z.namelist(), key=lambda n: int(re.sub(r"\D", "", n) or 0))
            for n in names:
                if n.endswith(".cif"):
                    raws.extend(read_cif(io.TextIOWrapper(z.open(n), encoding="utf-8").read()))
    else:
        raise MatterGenRunError(f"no {EXTXYZ_NAME} or {CIF_ZIP_NAME} in {outdir}")
    return [{"lattice": r.lattice, "species": r.species,
             "frac_coords": [[x % 1.0 for x in fc] for fc in r.frac_coords]} for r in raws], src


class MatterGenGenerator:
    """Generator protocol implementation; see module docstring."""

    name = GENERATOR_NAME

    def __init__(self, config: MatterGenConfig | dict | None = None):
        self.cfg = config if isinstance(config, MatterGenConfig) else MatterGenConfig.from_dict(config)

    def plan(self, request: GenerationRequest) -> dict[str, Any]:
        """Conditioning, batching and argv for a request (raises on unsupported conditioning)."""
        cond, dropped = conditioning_for(request, self.cfg)
        n = max(1, int(request.n))
        bs = max(1, min(self.cfg.batch_size, n))
        nb = min(self.cfg.max_batches, int(math.ceil(n / bs)))
        key = hashlib.sha256(json.dumps({"r": request.request_id, "c": cond,
                                         "cfg": asdict(self.cfg)}, sort_keys=True, default=str)
                             .encode()).hexdigest()[:12]
        outdir = Path(self.cfg.workdir) / f"{_safe(request.request_id)}_{key}"
        argv = build_command(request, self.cfg, outdir.resolve(), cond, bs, nb)
        return {"conditioning": cond, "dropped": dropped, "batch_size": bs, "num_batches": nb,
                "outdir": outdir, "argv": argv}

    def _run(self, argv: list[str], outdir: Path) -> dict[str, Any]:
        env = dict(os.environ)
        env.update(self.cfg.env)
        cmd: Any = " ".join(shlex.quote(a) for a in argv) if self.cfg.shell else argv
        t0 = time.time()
        try:
            cp = subprocess.run(cmd, shell=self.cfg.shell, env=env, capture_output=True, text=True,
                                timeout=self.cfg.timeout_s)
        except FileNotFoundError as exc:
            raise MatterGenRunError(f"MatterGen executable not found: {argv[0]!r} ({exc}). MatterGen "
                                    "must be installed in its own environment; set executable / "
                                    "command_prefix") from exc
        except subprocess.TimeoutExpired as exc:
            raise MatterGenRunError(f"mattergen-generate timed out after {self.cfg.timeout_s} s") from exc
        (outdir / "rtsc_mattergen.log").write_text(
            f"$ {cmd if isinstance(cmd, str) else shlex.join(cmd)}\n--- stdout\n{cp.stdout}\n--- stderr\n{cp.stderr}")
        if cp.returncode != 0:
            tail = (cp.stderr or cp.stdout or "")[-2000:]
            raise MatterGenRunError(f"mattergen-generate exited with {cp.returncode}: {tail}")
        return {"wall_s": round(time.time() - t0, 2), "returncode": cp.returncode}

    def generate(self, request: GenerationRequest) -> list[GeneratedStructure]:
        p = self.plan(request)
        outdir: Path = p["outdir"]
        outdir.mkdir(parents=True, exist_ok=True)
        stamp = outdir / "rtsc_invocation.json"
        argv_s = shlex.join(p["argv"])
        reused = False
        run_info: dict[str, Any] = {}
        if (self.cfg.reuse_outputs and stamp.exists()
                and json.loads(stamp.read_text()).get("argv") == argv_s
                and json.loads(stamp.read_text()).get("finished")):
            reused = True
        else:
            for f in (EXTXYZ_NAME, CIF_ZIP_NAME):
                (outdir / f).unlink(missing_ok=True)
            stamp.write_text(json.dumps({"argv": argv_s, "finished": False}))
            run_info = self._run(p["argv"], outdir)
            stamp.write_text(json.dumps({"argv": argv_s, "finished": True, **run_info}))
        structs, src = parse_outputs(outdir)
        allowed = set(request.chemical_system)
        out: list[GeneratedStructure] = []
        n_off = 0
        for i, s in enumerate(structs):
            if (self.cfg.postfilter_chemical_system and allowed and not set(s["species"]) <= allowed):
                n_off += 1
                continue
            prov = {
                "flags": ["HYPOTHETICAL_STRUCTURE"],
                "model": "MatterGen",
                "checkpoint": self.cfg.pretrained_name or str(self.cfg.model_path),
                "checkpoint_source": ("huggingface:microsoft/mattergen" if self.cfg.pretrained_name
                                      else f"model_path (epoch={self.cfg.checkpoint_epoch})"),
                "conditioning": p["conditioning"],
                "not_conditioned": p["dropped"],
                "guidance_factor": self.cfg.guidance_factor if p["conditioning"] else 0.0,
                "seed": int(request.seed),
                "batch_size": p["batch_size"], "num_batches": p["num_batches"],
                "cmdline": argv_s,
                "output_file": str(outdir / src), "index_in_output": i,
                "reused_previous_output": reused,
                "pressure_note": PRESSURE_NOTE,
                "target_pressure_GPa": request.pressure_GPa,
                **run_info,
            }
            out.append(GeneratedStructure(structure=s, generator=self.name,
                                          request_id=request.request_id, provenance=prov))
            if len(out) >= request.n:
                break
        if n_off:
            for g in out:
                g.provenance["postfiltered_off_system"] = n_off
        return out


__all__ = ["MatterGenGenerator", "MatterGenConfig", "UnsupportedConditioningError", "MatterGenRunError",
           "CHECKPOINT_PROPERTIES", "MATTERGEN_ELEMENTS", "PROPERTY_SOURCE_IDS", "conditioning_for",
           "build_command", "parse_outputs"]
