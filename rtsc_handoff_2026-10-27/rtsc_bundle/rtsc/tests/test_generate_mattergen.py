"""MatterGen adapter tests.  MatterGen itself is NOT run here: a FAKE ``mattergen-generate`` script
(clearly labelled, written into tmp_path) emits the documented output files so the adapter's
command construction, parsing, provenance and reuse logic can be tested.  Its structures are
placeholders, not MatterGen samples."""
from __future__ import annotations

import json
import os
import stat
import sys
import textwrap
import zipfile
from pathlib import Path

import pytest

from rtsc.discovery.interfaces import GenerationRequest
from rtsc.generate.mattergen import (CHECKPOINT_PROPERTIES, MATTERGEN_ELEMENTS, MatterGenConfig,
                                     MatterGenGenerator, MatterGenRunError, UnsupportedConditioningError,
                                     build_command, conditioning_for, parse_outputs)

FAKE = textwrap.dedent('''\
    #!{py}
    # FAKE mattergen-generate for rtsc tests: writes the documented output files with placeholder
    # structures (NOT MatterGen output).  Records its argv in fake_calls.jsonl next to itself.
    import json, os, sys, zipfile, io, ast
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "fake_calls.jsonl"), "a") as fh:
        fh.write(json.dumps(sys.argv[1:]) + "\\n")
    out = sys.argv[1]
    opts = dict(a[2:].split("=", 1) for a in sys.argv[2:] if a.startswith("--"))
    n = int(opts["batch_size"]) * int(opts["num_batches"])
    cond = ast.literal_eval(opts.get("properties_to_condition_on", "{{}}"))
    els = cond.get("chemical_system", "Mg-B").split("-")
    os.makedirs(out, exist_ok=True)
    frames, cifs = [], []
    for i in range(n):
        a = 3.0 + 0.01 * i
        sp = [els[0]] + [els[-1]] * 2
        pos = [(0, 0, 0), (a / 2, a / 2, a / 2), (0, a / 2, a / 2)]
        frames.append("%d\\nLattice=\\"%g 0 0 0 %g 0 0 0 %g\\" Properties=species:S:1:pos:R:3 pbc=\\"T T T\\"\\n" % (3, a, a, a)
                      + "".join("%s %g %g %g\\n" % (s, *p) for s, p in zip(sp, pos)))
        cifs.append("data_gen%d\\n_cell_length_a %g\\n_cell_length_b %g\\n_cell_length_c %g\\n"
                    "_cell_angle_alpha 90\\n_cell_angle_beta 90\\n_cell_angle_gamma 90\\n"
                    "_symmetry_space_group_name_H-M 'P 1'\\nloop_\\n_atom_site_type_symbol\\n"
                    "_atom_site_fract_x\\n_atom_site_fract_y\\n_atom_site_fract_z\\n" % (i, a, a, a)
                    + "".join("%s %g %g %g\\n" % (s, p[0] / a, p[1] / a, p[2] / a) for s, p in zip(sp, pos)))
    if os.environ.get("FAKE_MG_MODE") != "zip_only":
        open(os.path.join(out, "generated_crystals.extxyz"), "w").write("".join(frames))
    with zipfile.ZipFile(os.path.join(out, "generated_crystals_cif.zip"), "w") as z:
        for i, c in enumerate(cifs):
            z.writestr("gen_%d.cif" % i, c)
    sys.exit(int(os.environ.get("FAKE_MG_EXIT", "0")))
''')


@pytest.fixture
def fake_cli(tmp_path):
    b = tmp_path / "bin"
    b.mkdir()
    exe = b / "mattergen-generate"
    exe.write_text(FAKE.format(py=sys.executable))
    exe.chmod(exe.stat().st_mode | stat.S_IEXEC)
    return exe


def req(**kw):
    d = dict(chemical_system=("B", "Mg"), space_group=None, pressure_GPa=0.0, goal_id="G-AMB", n=3, seed=7)
    d.update(kw)
    return GenerationRequest(**d)


def test_checkpoint_table_matches_official_names():
    assert set(CHECKPOINT_PROPERTIES) == {
        "mattergen_base", "mp_20_base", "chemical_system", "space_group", "dft_mag_density",
        "dft_band_gap", "ml_bulk_modulus", "dft_mag_density_hhi_score", "chemical_system_energy_above_hull"}
    assert len(MATTERGEN_ELEMENTS) == 76 and "H" in MATTERGEN_ELEMENTS and "He" not in MATTERGEN_ELEMENTS


def test_conditioning_refuses_unsupported():
    cfg = MatterGenConfig(pretrained_name="chemical_system")
    assert conditioning_for(req(), cfg) == ({"chemical_system": "B-Mg"}, [])
    with pytest.raises(UnsupportedConditioningError, match="space_group"):
        conditioning_for(req(space_group=191), cfg)
    cfg2 = MatterGenConfig(pretrained_name="chemical_system", allow_unconditioned=("space_group",))
    assert conditioning_for(req(space_group=191), cfg2) == ({"chemical_system": "B-Mg"}, ["space_group"])
    with pytest.raises(UnsupportedConditioningError, match="unconditional"):
        conditioning_for(req(), MatterGenConfig(pretrained_name="mattergen_base"))
    with pytest.raises(UnsupportedConditioningError, match="He"):
        conditioning_for(req(chemical_system=("He", "Mg")), cfg)
    with pytest.raises(UnsupportedConditioningError, match="not a MatterGen property"):
        conditioning_for(req(extra={"mattergen_conditioning": {"pressure": 100}}), cfg)
    j = MatterGenConfig(pretrained_name="chemical_system_energy_above_hull")
    c, _ = conditioning_for(req(extra={"mattergen_conditioning": {"energy_above_hull": 0.05}}), j)
    assert c == {"chemical_system": "B-Mg", "energy_above_hull": 0.05}
    with pytest.raises(ValueError):
        conditioning_for(req(), MatterGenConfig(pretrained_name="nope"))


def test_model_path_reads_config(tmp_path):
    pytest.importorskip("yaml")
    (tmp_path / "config.yaml").write_text(
        "lightning_module:\n  diffusion_module:\n    model:\n      property_embeddings_adapt:\n"
        "        space_group: {}\n")
    cfg = MatterGenConfig(pretrained_name=None, model_path=str(tmp_path))
    assert conditioning_for(req(chemical_system=(), space_group=225), cfg)[0] == {"space_group": 225}
    with pytest.raises(UnsupportedConditioningError):
        conditioning_for(req(), cfg)


def test_build_command():
    cfg = MatterGenConfig(command_prefix=("srun", "-p", "gpu"), executable=("conda", "run", "-n", "mg",
                                                                          "mattergen-generate"))
    argv = build_command(req(), cfg, "/o", {"chemical_system": "B-Mg"}, 3, 1)
    assert argv[:7] == ["srun", "-p", "gpu", "conda", "run", "-n", "mg"]
    assert "--pretrained_name=chemical_system" in argv and "--seed=7" in argv
    assert "--properties_to_condition_on={'chemical_system':'B-Mg'}" in argv
    assert "--diffusion_guidance_factor=2.0" in argv and "--record_trajectories=False" in argv
    assert not any(" " in a for a in argv)


def test_generate_with_fake_cli(fake_cli, tmp_path):
    cfg = MatterGenConfig(executable=(str(fake_cli),), workdir=str(tmp_path / "runs"), batch_size=2)
    g = MatterGenGenerator(cfg)
    out = g.generate(req(n=3))
    assert len(out) == 3
    s = out[0]
    assert s.generator == "mattergen" and s.request_id == req(n=3).request_id
    assert set(s.structure["species"]) <= {"B", "Mg"} and len(s.structure["frac_coords"]) == 3
    p = s.provenance
    assert p["flags"] == ["HYPOTHETICAL_STRUCTURE"] and p["conditioning"] == {"chemical_system": "B-Mg"}
    assert p["seed"] == 7 and "no pressure conditioning" in p["pressure_note"]
    assert p["num_batches"] == 2 and p["output_file"].endswith("generated_crystals.extxyz")
    calls = (fake_cli.parent / "fake_calls.jsonl").read_text().splitlines()
    assert len(calls) == 1
    out2 = g.generate(req(n=3))                     # reuse finished outputs, no second run
    assert len((fake_cli.parent / "fake_calls.jsonl").read_text().splitlines()) == 1
    assert out2[0].provenance["reused_previous_output"] is True
    assert out2[0].structure == out[0].structure


def test_cif_zip_fallback_and_errors(fake_cli, tmp_path, monkeypatch):
    monkeypatch.setenv("FAKE_MG_MODE", "zip_only")
    g = MatterGenGenerator(MatterGenConfig(executable=(str(fake_cli),), workdir=str(tmp_path / "r")))
    out = g.generate(req(n=2))
    assert len(out) == 2 and out[0].provenance["output_file"].endswith("generated_crystals_cif.zip")
    assert abs(out[1].structure["lattice"][0][0] - 3.01) < 1e-6
    monkeypatch.setenv("FAKE_MG_EXIT", "3")
    with pytest.raises(MatterGenRunError, match="exited with 3"):
        MatterGenGenerator(MatterGenConfig(executable=(str(fake_cli),), workdir=str(tmp_path / "e"))
                           ).generate(req(n=1))
    with pytest.raises(MatterGenRunError, match="not found"):
        MatterGenGenerator(MatterGenConfig(executable=(str(tmp_path / "missing"),),
                                           workdir=str(tmp_path / "m"))).generate(req(n=1))
    with pytest.raises(MatterGenRunError):
        parse_outputs(tmp_path)


def test_no_mattergen_import():
    src = Path(__file__).resolve().parents[1] / "rtsc" / "generate" / "mattergen.py"
    text = src.read_text()
    assert "import mattergen" not in text and "from mattergen" not in text


def test_parse_ase_written_outputs(tmp_path):
    """MatterGen writes its outputs with ase.io.write (eval_utils.save_structures); parse those."""
    ase_io = pytest.importorskip("ase.io")
    import io as _io
    from ase.build import bulk
    frames = [bulk("Mg", "hcp", a=3.2, c=5.2), bulk("Al", "fcc", a=4.05, cubic=True)]
    ase_io.write(tmp_path / "generated_crystals.extxyz", frames)
    s, src = parse_outputs(tmp_path)
    assert src == "generated_crystals.extxyz" and [len(x["species"]) for x in s] == [2, 4]
    (tmp_path / "generated_crystals.extxyz").unlink()
    with zipfile.ZipFile(tmp_path / "generated_crystals_cif.zip", "w") as z:
        for i, a in enumerate(frames):
            b = _io.BytesIO()
            ase_io.write(b, a, format="cif")
            z.writestr(f"gen_{i}.cif", b.getvalue())
    s2, src2 = parse_outputs(tmp_path)
    assert src2 == "generated_crystals_cif.zip" and [len(x["species"]) for x in s2] == [2, 4]
