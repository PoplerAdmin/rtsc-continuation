"""Rebuild the Spectrum list exactly as `rtsc certify` does (cli.cmd_certify) and dump it for `rtsc recheck --spectra`."""
import sys, json
from pathlib import Path
from rtsc import pipeline
path, unit, rule, pol, env, out = sys.argv[1:7]
expl = float(sys.argv[7]) if len(sys.argv) > 7 and sys.argv[7] != "none" else None
spec = {"path": path, "format": "two_column", "units": unit, "label": Path(path).name,
        "columns": [0, 1], "rule": rule, "negative_policy": pol, "expected_lambda": expl, "expected_A_K": None}
cfg = pipeline.normalize_config({"certify": {"T_target_K": 300.0, "mu_hi": 0.15, "cutoff_rule": "block", "mode": "bracket",
                                             "omega_c_factor": 10.0},
                                 "spectra": {"rel_envelope": float(env), "two_column_rule": rule, "negative_policy": pol}})
sp = pipeline._read_spectrum_file(spec, cfg["spectra"])
sp = pipeline._with_envelope(sp, float(cfg["spectra"]["rel_envelope"]))
Path(out).write_text(json.dumps([s.to_dict() for s in sp], default=str))
print("wrote", out, len(sp))
