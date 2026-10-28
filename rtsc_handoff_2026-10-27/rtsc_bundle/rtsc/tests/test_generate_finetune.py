"""Fine-tuning export: CSV layout verified against MatterGen's csv-to-dataset (material_id, cif,
property columns).  Labels below are test fixtures, not results for real compounds."""
from __future__ import annotations

import csv
import json

import pytest

from rtsc.generate.finetune_export import (COLUMNS, FLOAT_PROPERTIES, export_finetune_dataset,
                                           records_from_ledger, structure_to_cif)
from rtsc.sources.fileformats import read_cif

S = {"lattice": [[3.07, 0, 0], [-1.535, 2.6587, 0], [0, 0, 3.52]], "species": ["Mg", "B", "B"],
     "frac_coords": [[0, 0, 0], [1 / 3, 2 / 3, 0.5], [2 / 3, 1 / 3, 0.5]]}


def test_cif_roundtrip():
    r = read_cif(structure_to_cif(S, "x"))[0]
    assert r.species == ["Mg", "B", "B"]
    for a, b in zip(r.lattice, S["lattice"]):
        assert all(abs(x - y) < 1e-6 for x, y in zip(a, b))


def test_export_layout(tmp_path):
    csv.field_size_limit(10**8)
    recs = [{"material_id": f"t{i}", "structure": S, "formula": "MgB2", "pressure_GPa": 0.0,
             "T_low_K": 10.0 + i, "T_up_K": 50.0 + i, "statuses": {"G-AMB": "UNDECIDED"},
             "rtsc_A_K": 100.0 + i} for i in range(5)]
    info = export_finetune_dataset(recs, tmp_path, dataset_name="ds", val_fraction=0.2)
    assert info["n_train"] == 4 and info["n_val"] == 1
    rows = list(csv.DictReader(open(tmp_path / "ds" / "train.csv")))
    assert tuple(rows[0].keys()) == COLUMNS and rows[0]["chemical_system"] == "B-Mg"
    assert rows[0]["cif"].startswith("data_t") and rows[0]["rtsc_lambda"] == ""
    for p in FLOAT_PROPERTIES:
        assert (tmp_path / "mattergen_conf" / "property_embeddings" / f"{p}.yaml").exists()
    dm = (tmp_path / "mattergen_conf" / "data_module" / "ds.yaml").read_text()
    assert "datasets/cache/ds" in dm and "average_density" in dm
    assert "PROPERTY_SOURCE_IDS" in (tmp_path / "FINETUNE_INSTRUCTIONS.txt").read_text()
    yaml = pytest.importorskip("yaml")
    assert yaml.safe_load((tmp_path / "mattergen_conf/property_embeddings/rtsc_A_K.yaml").read_text())[
        "name"] == "rtsc_A_K"


def test_records_from_ledger(tmp_path):
    from rtsc.discovery.ledger import Ledger
    from rtsc.pipeline import CLASSIFY_ORIGIN
    from rtsc.status import Status
    rec = tmp_path / "record.json"
    rec.write_text(json.dumps({"spectra": [{"A_K": 321.0, "lambda": 1.2, "omega_log_K": 900.0}]}))
    led = Ledger(tmp_path / "ledger.sqlite")
    led.add_material(material_id="m1", seq=1, fingerprint="f1", coarse_key="c", formula="MgB2",
                     chemsys=("B", "Mg"), pressure_GPa=0.0, structure=S, generator="test",
                     request_id="r", request={}, arm="a", round_=1, provenance={})
    for goal, lo, up in (("G-AMB", 20.0, 60.0), ("G-RT", 25.0, 70.0)):
        led.append_verdict(material_id="m1", goal_id=goal, round_=1, stage="pipeline",
                           origin=CLASSIFY_ORIGIN, status=Status.UNDECIDED, T_star_K=77.0,
                           pressure_GPa=0.0, T_low_K=lo, T_up_K=up, record_path=str(rec))
    led.close()
    recs = records_from_ledger(tmp_path / "ledger.sqlite")
    assert len(recs) == 1 and recs[0]["T_low_K"] == 25.0 and recs[0]["T_up_K"] == 60.0
    info = export_finetune_dataset(recs, tmp_path / "out", val_fraction=0.0)
    row = next(csv.DictReader(open(tmp_path / "out" / "rtsc_ledger" / "train.csv")))
    assert float(row["rtsc_A_K"]) == 321.0 and float(row["rtsc_lambda"]) == 1.2
    assert info["n_rows"] == 1
