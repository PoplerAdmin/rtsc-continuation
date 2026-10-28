"""Export DFT-level ledger labels as a MatterGen fine-tuning dataset (documented; not run here).

Format verified against github.com/microsoft/mattergen @ 92423660 (README "Fine-tune on your own
property data", ``mattergen/scripts/csv_to_dataset.py``, ``CrystalDatasetBuilder.from_csv`` in
``mattergen/common/data/dataset.py``, ``data-release/mp-20/mp_20.zip``):

* a folder of CSV files (``train.csv``, ``val.csv``; every ``*.csv`` becomes a split);
* required columns ``material_id`` and ``cif`` (CIF text, parsed with pymatgen ``CifParser``,
  primitive cell); every other column whose name is in ``PROPERTY_SOURCE_IDS`` becomes a property;
  others are ignored (so the ``rtsc_*`` columns below must be added to ``PROPERTY_SOURCE_IDS``);
* missing values are allowed (``filter_sparse_properties`` drops rows lacking a selected property);
* then ``csv-to-dataset --csv-folder <dir> --dataset-name <name> --cache-folder datasets/cache``,
  a ``property_embeddings/<prop>.yaml`` per new float property (copy of ``dft_mag_density.yaml``),
  a ``data_module/<name>.yaml`` pointing at the cache, and ``mattergen-finetune ...``.

Columns written: ``material_id, cif, chemical_system`` (an official MatterGen property id), and the
float labels ``rtsc_A_K`` (A = integral alpha^2F d omega, kelvin), ``rtsc_lambda``,
``rtsc_omega_log_K``, ``rtsc_Tc_low_K``, ``rtsc_Tc_up_K`` (certified bracket on Tc_model),
``rtsc_pressure_GPa``; metadata ``rtsc_status``, ``rtsc_goals``, ``rtsc_formula`` (ignored by
MatterGen).  Only DFT-stage verdicts (REJECTED_PROVEN / CERTIFIED_CANDIDATE / UNDECIDED) are
exported by default.  Labels are statements about the isotropic Migdal-Eliashberg MODEL (ARCH S6).

NOTE: these structures may be high-pressure phases while MatterGen's training data are ambient;
a pressure-aware generator needs ``rtsc_pressure_GPa`` as a conditioning property too.
"""
from __future__ import annotations

import csv
import json
import math
import random
import sqlite3
from pathlib import Path
from typing import Any, Iterable

LABEL_STATUSES = ("REJECTED_PROVEN", "CERTIFIED_CANDIDATE", "UNDECIDED")
FLOAT_PROPERTIES = ("rtsc_A_K", "rtsc_lambda", "rtsc_omega_log_K", "rtsc_Tc_low_K", "rtsc_Tc_up_K",
                    "rtsc_pressure_GPa")
COLUMNS = ("material_id", "cif", "chemical_system") + FLOAT_PROPERTIES + (
    "rtsc_status", "rtsc_goals", "rtsc_formula")

# verbatim structure of mattergen/conf/lightning_module/diffusion_module/model/property_embeddings/dft_mag_density.yaml
EMBEDDING_YAML = """_target_: mattergen.property_embeddings.PropertyEmbedding
name: {name}
unconditional_embedding_module:
  _target_: mattergen.property_embeddings.EmbeddingVector
  hidden_dim: ${{lightning_module.diffusion_module.model.hidden_dim}}
conditional_embedding_module:
  _target_: mattergen.diffusion.model_utils.NoiseLevelEncoding
  d_model: ${{lightning_module.diffusion_module.model.hidden_dim}}
scaler:
  _target_: mattergen.common.utils.data_utils.StandardScalerTorch
"""

# mattergen/conf/data_module/mp_20.yaml with root_dir / average_density replaced
DATA_MODULE_YAML = """_target_: mattergen.common.data.datamodule.CrystDataModule
_recursive_: true
properties: []

transforms:
- _target_: mattergen.common.data.transform.symmetrize_lattice
  _partial_: true
- _target_: mattergen.common.data.transform.set_chemical_system_string
  _partial_: true

dataset_transforms:
  - _target_: mattergen.common.data.dataset_transform.filter_sparse_properties
    _partial_: true

average_density: {density} # atoms/Angstrom**3 (mean over the exported structures)
root_dir: ${{oc.env:PROJECT_ROOT}}/../datasets/cache/{name}

train_dataset:
  _target_: mattergen.common.data.dataset.CrystalDataset.from_cache_path
  cache_path: ${{data_module.root_dir}}/train
  properties: ${{data_module.properties}}
  transforms: ${{data_module.transforms}}
  dataset_transforms: ${{data_module.dataset_transforms}}

val_dataset:
  _target_: mattergen.common.data.dataset.CrystalDataset.from_cache_path
  cache_path: ${{data_module.root_dir}}/val
  properties: ${{data_module.properties}}
  transforms: ${{data_module.transforms}}
  dataset_transforms: ${{data_module.dataset_transforms}}

num_workers:
  train: 0
  val: 0

batch_size:
  train: ${{eval:'(512 // ${{trainer.accumulate_grad_batches}}) // (${{trainer.devices}} * ${{trainer.num_nodes}})'}}
  val: ${{eval:'(64 // (${{trainer.devices}} * ${{trainer.num_nodes}}))'}}
  test: ${{eval:'(64 // (${{trainer.devices}} * ${{trainer.num_nodes}}))'}}

max_epochs: 900
"""


def structure_to_cif(s: dict[str, Any], name: str = "rtsc") -> str:
    """Minimal P1 CIF (cell parameters + fractional coordinates), readable by pymatgen CifParser."""
    import numpy as np
    lat = np.asarray(s["lattice"], dtype=float)
    a, b, c = (float(np.linalg.norm(v)) for v in lat)

    def ang(u, v):
        return math.degrees(math.acos(float(np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v)))))
    al, be, ga = ang(lat[1], lat[2]), ang(lat[0], lat[2]), ang(lat[0], lat[1])
    lines = [f"data_{name}", "_symmetry_space_group_name_H-M   'P 1'", "_symmetry_Int_Tables_number   1",
             f"_cell_length_a   {a:.8f}", f"_cell_length_b   {b:.8f}", f"_cell_length_c   {c:.8f}",
             f"_cell_angle_alpha   {al:.8f}", f"_cell_angle_beta   {be:.8f}", f"_cell_angle_gamma   {ga:.8f}",
             "loop_", " _symmetry_equiv_pos_site_id", " _symmetry_equiv_pos_as_xyz", "  1  'x, y, z'",
             "loop_", " _atom_site_type_symbol", " _atom_site_label", " _atom_site_fract_x",
             " _atom_site_fract_y", " _atom_site_fract_z", " _atom_site_occupancy"]
    for k, (el, f) in enumerate(zip(s["species"], s["frac_coords"])):
        lines.append(f"  {el}  {el}{k}  {f[0] % 1.0:.8f}  {f[1] % 1.0:.8f}  {f[2] % 1.0:.8f}  1")
    return "\n".join(lines) + "\n"


def _find_number(obj: Any, keys: tuple[str, ...]) -> float | None:
    """First finite number stored under any of ``keys`` (recursive; best effort over records)."""
    if isinstance(obj, dict):
        for k in keys:
            v = obj.get(k)
            if isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v):
                return float(v)
        for v in obj.values():
            r = _find_number(v, keys)
            if r is not None:
                return r
    elif isinstance(obj, list):
        for v in obj:
            r = _find_number(v, keys)
            if r is not None:
                return r
    return None


def _spectral_labels(record_path: str | None) -> dict[str, float | None]:
    out: dict[str, float | None] = {"rtsc_A_K": None, "rtsc_lambda": None, "rtsc_omega_log_K": None}
    if not record_path or not Path(record_path).exists():
        return out
    try:
        rec = json.loads(Path(record_path).read_text())
    except Exception:  # noqa: BLE001
        return out
    out["rtsc_A_K"] = _find_number(rec, ("A_K", "A"))
    out["rtsc_lambda"] = _find_number(rec, ("lambda", "lambda_", "lam"))
    out["rtsc_omega_log_K"] = _find_number(rec, ("omega_log_K", "omega_log"))
    return out


def records_from_ledger(ledger_path: str | Path, statuses: Iterable[str] = LABEL_STATUSES) -> list[dict]:
    """One label record per material from the ledger's current verdicts (read-only connection).

    Brackets from several goals at the material's own pressure are intersected
    (T_low = max, T_up = min): each is a valid bound on the same Tc_model."""
    uri = f"file:{Path(ledger_path).resolve()}?mode=ro"
    con = sqlite3.connect(uri, uri=True)
    con.row_factory = sqlite3.Row
    rows = con.execute(
        "SELECT v.material_id, v.goal_id, v.status, v.T_low_K, v.T_up_K, v.pressure_GPa AS vp, "
        "v.record_path, m.structure, m.formula, m.chemsys, m.pressure_GPa AS mp, m.dup_of "
        "FROM current_status c JOIN verdicts v ON v.id=c.verdict_id "
        "JOIN materials m ON m.material_id=c.material_id WHERE v.stage='pipeline' ORDER BY m.seq").fetchall()
    con.close()
    want = set(statuses)
    by: dict[str, dict] = {}
    for r in rows:
        if r["status"] not in want or r["dup_of"]:
            continue
        if r["vp"] is not None and r["mp"] is not None and abs(float(r["vp"]) - float(r["mp"])) > 1e-9:
            continue
        d = by.setdefault(r["material_id"], {
            "material_id": r["material_id"], "structure": json.loads(r["structure"]),
            "formula": r["formula"], "chemsys": r["chemsys"], "pressure_GPa": r["mp"],
            "T_low_K": None, "T_up_K": None, "statuses": {}, "record_path": r["record_path"]})
        d["statuses"][r["goal_id"]] = r["status"]
        if r["T_low_K"] is not None:
            d["T_low_K"] = r["T_low_K"] if d["T_low_K"] is None else max(d["T_low_K"], r["T_low_K"])
        if r["T_up_K"] is not None:
            d["T_up_K"] = r["T_up_K"] if d["T_up_K"] is None else min(d["T_up_K"], r["T_up_K"])
        d["record_path"] = d["record_path"] or r["record_path"]
    for d in by.values():                 # label the geometry the DFT stage used, not the generated one
        s = _dft_structure(d.get("record_path"))
        if s is not None:
            d["structure"] = s
    return list(by.values())


def _dft_structure(record_path: str | None) -> dict[str, Any] | None:
    """Structure from the pipeline record: the DFPT final cell if it carries lattice/species/
    frac_coords, else the candidate structure handed to the DFT stage."""
    if not record_path or not Path(record_path).exists():
        return None
    try:
        rec = json.loads(Path(record_path).read_text())
    except Exception:  # noqa: BLE001
        return None
    cell = (((rec.get("stages") or {}).get("dfpt") or {}).get("provenance") or {}).get("cell") or {}
    for src in (cell, rec.get("candidate") or {}):
        if all(src.get(k) for k in ("lattice", "species", "frac_coords")):
            return {"lattice": src["lattice"], "species": list(src["species"]),
                    "frac_coords": src["frac_coords"]}
    return None


def to_row(rec: dict[str, Any]) -> dict[str, Any]:
    s = rec["structure"]
    els = sorted(set(s["species"]))
    row = {"material_id": rec["material_id"], "cif": structure_to_cif(s, rec["material_id"]),
           "chemical_system": "-".join(els),
           "rtsc_Tc_low_K": rec.get("T_low_K"), "rtsc_Tc_up_K": rec.get("T_up_K"),
           "rtsc_pressure_GPa": rec.get("pressure_GPa"),
           "rtsc_status": json.dumps(rec.get("statuses", {}), sort_keys=True),
           "rtsc_goals": ",".join(sorted(rec.get("statuses", {}))),
           "rtsc_formula": rec.get("formula", "")}
    lab = _spectral_labels(rec.get("record_path"))
    for k in ("rtsc_A_K", "rtsc_lambda", "rtsc_omega_log_K"):
        row[k] = rec.get(k, lab[k])
    return row


def export_finetune_dataset(records: list[dict[str, Any]], outdir: str | Path, *,
                            dataset_name: str = "rtsc_ledger", val_fraction: float = 0.1,
                            seed: int = 0) -> dict[str, Any]:
    """Write ``<outdir>/<dataset_name>/{train,val}.csv`` + config templates + instructions."""
    import numpy as np
    out = Path(outdir)
    ds = out / dataset_name
    ds.mkdir(parents=True, exist_ok=True)
    rows = [to_row(r) for r in records]
    idx = list(range(len(rows)))
    random.Random(seed).shuffle(idx)
    n_val = int(round(val_fraction * len(rows))) if len(rows) > 1 else 0
    if len(rows) > 1:
        n_val = min(max(n_val, 1), len(rows) - 1)
    val = set(idx[:n_val])
    for split, sel in (("train", [i for i in range(len(rows)) if i not in val]),
                       ("val", [i for i in range(len(rows)) if i in val])):
        with open(ds / f"{split}.csv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(COLUMNS))
            w.writeheader()
            for i in sel:
                w.writerow({k: ("" if rows[i][k] is None else rows[i][k]) for k in COLUMNS})
    dens = [len(r["structure"]["species"]) / abs(float(np.linalg.det(np.asarray(r["structure"]["lattice"]))))
            for r in records] or [0.05771451654022283]
    conf = out / "mattergen_conf"
    (conf / "property_embeddings").mkdir(parents=True, exist_ok=True)
    (conf / "data_module").mkdir(parents=True, exist_ok=True)
    for p in FLOAT_PROPERTIES:
        (conf / "property_embeddings" / f"{p}.yaml").write_text(EMBEDDING_YAML.format(name=p))
    (conf / "data_module" / f"{dataset_name}.yaml").write_text(
        DATA_MODULE_YAML.format(name=dataset_name, density=float(np.mean(dens))))
    prop = "rtsc_A_K"
    instr = f"""MatterGen fine-tuning on rtsc ledger labels (generated by rtsc.generate.finetune_export)
Verified against github.com/microsoft/mattergen README "Fine-tune on your own property data".
Run inside the MatterGen environment (Python 3.10, CUDA GPU), from the mattergen repo root:

1. Add these names to PROPERTY_SOURCE_IDS in mattergen/common/utils/globals.py:
   {", ".join(FLOAT_PROPERTIES)}
2. Copy mattergen_conf/property_embeddings/*.yaml to
   mattergen/conf/lightning_module/diffusion_module/model/property_embeddings/
   and mattergen_conf/data_module/{dataset_name}.yaml to mattergen/conf/data_module/
3. csv-to-dataset --csv-folder {ds.resolve()} --dataset-name {dataset_name} --cache-folder datasets/cache
   (the data_module yaml expects the cache at datasets/cache/{dataset_name})
4. mattergen-finetune adapter.pretrained_name=mattergen_base data_module={dataset_name} \\
     +lightning_module/diffusion_module/model/property_embeddings@adapter.adapter.property_embeddings_adapt.{prop}={prop} \\
     ~trainer.logger data_module.properties=["{prop}"]
   (add chemical_system / rtsc_pressure_GPa the same way for joint conditioning)
5. Generate with --model_path=<output dir> and set MatterGenConfig.model_path; the adapter reads
   the supported properties from <model_path>/config.yaml.

Small ledgers overfit: consider appending these rows to alex_mp_20 train/val CSVs (columns
absent there stay empty; filter_sparse_properties keeps only labelled rows for the selected
properties).  Labels refer to the isotropic Migdal-Eliashberg MODEL (ARCH S6); structures may be
high-pressure phases (rtsc_pressure_GPa), unlike MatterGen's ambient training data.
rows: {len(rows)} (train {len(rows) - n_val}, val {n_val})
"""
    (out / "FINETUNE_INSTRUCTIONS.txt").write_text(instr)
    return {"dataset_dir": str(ds), "n_rows": len(rows), "n_train": len(rows) - n_val, "n_val": n_val,
            "properties": list(FLOAT_PROPERTIES)}


def export_from_ledger(ledger_path: str | Path, outdir: str | Path, **kw: Any) -> dict[str, Any]:
    return export_finetune_dataset(records_from_ledger(ledger_path), outdir, **kw)


__all__ = ["export_finetune_dataset", "export_from_ledger", "records_from_ledger", "structure_to_cif",
           "COLUMNS", "FLOAT_PROPERTIES"]
