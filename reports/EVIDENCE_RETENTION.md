# Evidence retention and outside verification

Adopted 2026-10-27, after the WSL disk image filled the host drive. Tool: `bench_tools/jobs/prune.py`.

**SHA-256 checksums are the method for outside verification of results in this project.** Every file a run produced is either kept, or deleted only after its SHA-256 is recorded in a manifest. Kept files are listed with their checksums too. A verifier checks the evidence against these manifests and does not need to trust our copies of it.

## What is kept and what is deleted

| Category | Examples | Policy |
|---|---|---|
| Inputs | structures, QE/EPW inputs, pseudopotentials, `config.toml` | Kept permanently |
| Physics outputs | `dynN`, `*.dvscf*`, `elph_dir/`, `a2F.dos*`, `al.dyn*.elph.*`, EPW `epmatwp`, logs | Kept permanently |
| Results | `record.json`, `verdict.json`, certificates, `RUN_MANIFEST.tsv` | Kept permanently |
| scf restart data | everything inside `<prefix>.save/` | Kept |
| **Per-q wavefunctions** | `_ph0/<prefix>.q_N/<prefix>.wfcK` | **Hashed, recorded, then deleted** |
| Copies and dead scratch | duplicated `tmp/` trees, scratch of runs killed by the crash | Hashed, recorded, then deleted |

Only wavefunction files are deleted from live runs. They take 90–95 % of a ph.x scratch folder, nothing downstream reads them once their q point is finished, and the same inputs regenerate them.

### When a file may be deleted

- **Live ph.x run (`qdone`).** q point N is finished only when `dynN` is complete **and** `dyn(N+1)` exists. ph.x creates every q folder up front, so folder presence is not used.
- **`lambda_tetra` run (`qmark`).** q point N is finished when the per-q outputs for N and N+1 both exist.
- **Never:** the last q point of a live run, or a run on hold. Pb is on hold until its `ERROR` verdict is resolved.
- **Finished run (`wfc`):** all remaining wavefunctions outside `.save/`.

`prune_daemon.sh` applies these rules every 120 s.

## Manifests

| File | Contents |
|---|---|
| `<root>/PRUNED_MANIFEST.tsv` | One row per deleted file: UTC time, machine label, path, bytes, SHA-256, file mtime, reason |
| `<root>/RUN_MANIFEST.tsv` | One row per kept file of a finished run (from `prune.py inventory`) |
| `reports/data/pruned_manifests/` | Copies of both, outside the WSL disk image and committed to this repository |

Order of operations for each file:
1. Hash it.
2. Append the row to both manifests and `fsync` them.
3. Only then delete it.

A crash at any point leaves the file either still present or recorded.

Git commits timestamp the central copies, so a manifest cannot be changed later without that showing in history.

Paths are recorded relative to the run account's home directory (`~/rtsc_work/...`), and the machine is recorded by a neutral label (`run-machine-1`). Hostnames and usernames are never recorded.

## How an outside party verifies a result

1. **Certificate.** Run `rtsc recheck <verdict.json>`. It re-derives the verdict from the α²F spectrum and the certificate parameters, without the DFT scratch.
2. **Kept evidence.** Run `python prune.py verify <RUN_MANIFEST.tsv> <base>`, where `<base>` is the directory that stands for `~` on the verifier's copy. Every kept file is re-hashed, and any file that is `MISMATCH` or `MISSING` makes the command exit 1. The spectra, dyn files, dvscf files and logs are therefore exactly the ones the result was computed from.
3. **Deleted intermediates.** Regenerate them by rerunning the recorded inputs, then run `prune.py verify <PRUNED_MANIFEST.tsv> <base>` on the regenerated files.

### Limits of step 3

A bitwise SHA-256 match on regenerated wavefunctions needs the same QE build, libraries, MPI ranks and pools, and in practice the same hardware.

Wavefunctions carry an arbitrary phase per state, and parallel floating-point reductions are not order-stable. A rerun elsewhere can therefore produce different bytes and still be equally correct.

In that case:
- the manifest still proves what existed and when;
- the comparison moves to the physics outputs (dyn frequencies, dvscf, α²F, λ, ω_log) within stated tolerances.

Whether QE reproduces wavefunction bytes on this machine has not yet been tested; see *Open items*.

## Candidates

For any material that rtsc reports as a candidate (`CERTIFIED_CANDIDATE`, or `UNDECIDED` with an upper bound above the target):
- its full scratch is archived to a second drive before any pruning;
- it is rerun independently with a different code or settings.

An independent rerun, not the archive, is the standard of evidence for a claimed superconductor.

## Open items

- Test whether a rerun with identical layout reproduces wavefunction SHA-256 values on this machine. Record the result here.
- Pb: the verdict is `ERROR` because rtsc refuses to clip small negative α²F mass at σ = 1 and 2. Its scratch is on hold until that is decided.
