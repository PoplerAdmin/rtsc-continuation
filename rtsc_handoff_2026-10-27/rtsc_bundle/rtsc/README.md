# rtsc: rigorous screening of phonon-mediated superconductors against a target Tc

`rtsc` takes candidate compounds and gives each one a verdict about whether its phonon-mediated
critical temperature **in a stated model** can reach a target T\* (default 300 K). Verdicts that
claim something (REJECTED_PROVEN, CERTIFIED_CANDIDATE) rely on rules proved in Lean 4 and on
rigorous ball-arithmetic checks (python-flint / Arb). The other verdicts make no such claim and
are marked as such.

Design and soundness contract: [`ARCH.md`](ARCH.md). Lean proofs: `/home/claude/tc/lean`
(library `TcR1` for the universal bound, library **`PipelineLogic`** for the classifier rules).
Python-to-Lean correspondence: [`docs/CORRESPONDENCE_LEAN.md`](docs/CORRESPONDENCE_LEAN.md) for
the `PipelineLogic` library and `docs/CORRESPONDENCE.md` (to be written) for the full mapping of
code to theorems.

---

## 1. What the verdicts mean, and what they do not mean

**Every verdict is a statement about a MODEL, not a prediction of a measured Tc.** The model:

* isotropic Migdal–Eliashberg theory with the linearised gap equation. Tc_model is the supremum of
  the temperatures where the linearised problem is unstable. It is at least the usual
  highest-solution Tc;
* the given Eliashberg function α²F, taken from DFPT or from an attached file. The true atom
  weights are assumed to lie in a relative envelope `[(1-e)w, (1+e)w]`. Uncertainty in the
  *frequencies* enters only through alternative spectra (for example each DFPT smearing), and
  all alternatives must agree (S5);
* a Coulomb pseudopotential μ\* in `[0, mu_hi]` with the block-form cutoff
  `N_c(T) = #{n >= 0 : (2n+1)·π·T <= omega_c}`.

| status | meaning |
|---|---|
| `CERTIFIED_CANDIDATE` | Proved: Tc_model ≥ T\* for every spectrum in the envelope, every alternative spectrum and every μ\* ≤ mu_hi, **under the configured cutoff rule and omega_c** (ARCH §4b F3). |
| `REJECTED_PROVEN` | Proved: Tc_model < T\* for every spectrum in the envelope, every alternative, every μ\* ≥ 0 and every cutoff. |
| `UNDECIDED` | The certified bracket `[T_low, T_up]` contains T\*, or the alternatives disagree. |
| `EXCLUDED_ASSUMPTION` | Removed by a physical criterion that is **not proved**: convex hull, metallicity, magnetism, or imaginary first-principles phonon modes. The criterion, value, threshold and data source are recorded. |
| `DEFERRED` | Not computed. Reasons: budget, heuristic ordering, duplicate of another input, or no spectrum while DFPT is disabled. |
| `ERROR` | A stage failed. Causes: an unreadable or malformed input; a failed DFT run; zero DFPT spectra (F1); a DFT parse/units check failure (`lambda_vs_a2F_footer`, `lambda_vs_epw`); an incomplete or inconsistent attached a²F file; an engine exception; or a crashed worker process (retried on resume). See the record. |

**What the verdicts do not cover.** The model leaves out anharmonicity and quantum nuclear motion
(both large for hydrogen), non-adiabatic effects (ω/E_F not small), anisotropic or multiband gaps,
spin fluctuations, vertex corrections, the accuracy of DFT and DFPT, and whether the material can
be synthesised. The pipeline marks these risks with **flags** (`ANHARMONIC_H`, `NONADIABATIC`,
`MAGNETIC_PROXIMITY`, `MULTIBAND_ANISOTROPY`, `MU_STAR_UNCERTAIN`, `SPECTRUM_SPREAD`, `SOFT_MODES`,
`STRONG_COUPLING`, `SPECTRUM_WARNINGS`, `SYNTHETIC_DATA`, `MULTIBAND_RISK`, `COULOMB_INACTIVE`,
`SMEARINGS_SUBSET`, `PARAMS_MISMATCH`, and the DFT flags such as `DFPT_CHECK_FAIL:...`).
**A flag never changes a status
or a bound** (S6). For example, "REJECTED_PROVEN + ANHARMONIC_H" means that the material is
rejected within the harmonic model and that the harmonic model may describe it poorly.

**Single-band versus multiband rejections.** With the default `rejection_model = "isotropic"`,
every rejection uses the DOS-averaged single-band α²F. Such verdicts carry the assumption
`M-ISO-SINGLEBAND` and the flag `MULTIBAND_RISK`, because the averaged A can underestimate Tc in a
multiband material (ARCH §4c). `rejection_model = "multiband"` requires band-resolved coupling
data (rule R-MULTI). You attach it per material with the manifest key `band_coupling`, which
points to a `rtsc.dfpt.BandCoupling` `.json`/`.toml` file, or with `rtsc certify --band-coupling`.
Without band data, multiband mode cannot reject: the result is UNDECIDED.

**`COULOMB_INACTIVE`.** This flag means π·T_low > omega_c. The Coulomb block is then empty at the
certified lower bound, so that bound does not depend on μ\*.

Soundness contract, in brief (see ARCH §3):

* S1: only `rtsc.certify.classify` can produce REJECTED_PROVEN or CERTIFIED_CANDIDATE. The
  pipeline records each verdict's origin and refuses (raises) any proven status with another
  origin.
* S2: heuristics (priority scores, budgets, ML potentials) only reorder candidates or DEFER them.
  `tests/test_pipeline_s2.py` enforces this.
* S3: each existence filter yields EXCLUDED_ASSUMPTION with its exact criterion, and each one can
  be disabled.
* S4: every input ends with exactly one status. The run checks this and raises `AccountingError`
  otherwise. Unreadable files, duplicates and crashed workers all end with an explicit status.
* S5/F1: all alternative spectra go to `classify` together. If DFPT yields zero alternatives,
  the material is ERROR and is never classified.

## 2. Installation

```sh
pip install -e .                 # core: numpy, scipy, python-flint==0.9.0 (certify + demo)
pip install -e '.[mp]'           # Materials Project source (mp-api, pymatgen)
pip install -e '.[dft]'          # pymatgen, ASE, phonopy, spglib, PyYAML for structures/phonons
pip install -e '.[ml]'           # MACE pre-check (or '.[ml-chgnet]'); heuristic ordering only
pip install -e '.[test]' && pytest
```

`python-flint` is pinned to 0.9.0, the version the trusted base was verified with. If another
version is installed, the pipeline warns, and every verdict records
`TRUST: python-flint X / FLINT Y` in its assumptions.

Quantum ESPRESSO (`pw.x`, `ph.x`, `q2r.x`, `matdyn.x`) is installed separately. Its path is set by
`[dfpt.launcher] bin_dir` or by `PATH`.

## 3. Commands

```sh
rtsc run -c config.toml [--workers N] [--workdir DIR] [--no-resume]
rtsc run -c config.toml --plan-only          # stages 1-5 only (sources ... budget, ML)
rtsc run -c config.toml --shard K/N          # compute every N-th material (job arrays)
rtsc collect -c config.toml                  # merge records of a sharded run, write reports
rtsc certify a2F.dat --units meV --mu-hi 0.15 --target 300 [--envelope 0.10]
             [--rejection-model isotropic|multiband] [--band-coupling bc.toml ...]
             [--format two_column|qe_a2f_dos|epw] [--columns 0,1] [--rule bins]
             [--omega-c K | --omega-c-factor 10] [--mode bracket|decide] [--out verdict.json]
rtsc demo [--workdir rtsc_demo]              # SYNTHETIC spectra end to end, no DFT
rtsc recheck runs/x/<id>/verdict.json        # re-verify a certificate (rtsc.certify); spectra and
                                             # band couplings come from the record next to it
```

`rtsc demo` uses **invented** spectra labelled `SYNTHETIC-*`. Their verdicts test the software
and say nothing about any real material. The reports print a SYNTHETIC banner.

## 4. Outputs

In `[report] outdir`, which defaults to `[run] workdir`:

* `candidates.csv` / `candidates.json`: one row per input. Rows are sorted by status
  (CERTIFIED, UNDECIDED, REJECTED, EXCLUDED, DEFERRED, ERROR), then by certified T_low
  (descending). Columns: `material_id formula source synthetic status T_low_K T_up_K T_target_K
  bracket_model mu_hi rel_envelope cutoff n_spectra A_max_K lambda_min lambda_max dropped_nonpos_mass_K
  pressure_GPa priority
  proof_basis rule_ids assumptions flags reasons verdict_origin record`.
  - `T_target_K` and `mu_hi` are the values certify actually used, read from the verdict's
    certificate. If they differ from the config, the row gets the flag `PARAMS_MISMATCH`.
  - `bracket_model` says which model the bracket refers to. It is single-band
    (DOS-averaged) except for R-MULTI rejections, also in multiband mode (Lean note MB-N1).
    Such UNDECIDED and CERTIFIED rows carry `M-ISO-SINGLEBAND`.
  - `proof_basis` is "Arb + Lean", or "Arb + hand proof (NOT all in Lean)" when the verdict's
    assumptions say so (e.g. `row_check = "cform"`). `summary.md` lists such verdicts
    separately.
  - Bounds are printed at full precision (299.9999999997 is never shown as 300.0).
* `certificates/<id>.json`: the verdict and its origin, the spectra, and the pipeline context
  (filter results with thresholds, priority, validity details, config hash).
* `summary.md`: the counts (they add up to the number of inputs), the model assumptions in plain
  words, the top candidates, exclusions by criterion, deferrals, errors and flags.
* `run_summary.json`; `inputs.json` (or `inputs_shard_K_of_N.json` for a shard). These are
  the lists of inputs that `collect` merges.
* `<workdir>/<id>/record.json`: the resumable per-material record. `<workdir>/<id>/verdict.json`
  holds the bare Verdict that `rtsc recheck` takes as input. `<workdir>/<id>/dfpt/` holds the QE
  run directory.

## 5. Configuration reference (TOML)

Relative paths are resolved against the directory of the config file. Unknown keys are errors in
every section the orchestrator owns. The `[dfpt.*]` and `[certify]` keys are validated by team
DFT's and team CERT's code. Examples: `examples/local_structures/config.toml`,
`examples/mp_hydrides/config.toml`, `examples/demo/config.toml`.

### `[run]`
| key | default | meaning |
|---|---|---|
| `name` | `"rtsc-run"` | label used in reports |
| `workdir` | `"rtsc_run"` | records, DFT directories, reports |
| `workers` | `1` | materials processed at the same time (process pool). Each DFPT job starts its own MPI launch, so the cores in use are `workers × [dfpt.launcher] nprocs` |
| `resume` | `true` | reuse a finished record only if three things all match: the material id, the configuration hash, and the **input digest**. The digest is a sha256 of the candidate and of the contents of its structure file, attached a²F files and band-coupling files, so a changed input is always recomputed. Spectra are reused if only certify is missing. ERROR records from a crashed worker are retried |

### `[[sources]]` (array; every entry needs `type`)
| type | keys |
|---|---|
| `local` | `path`; `recursive` (false); `pressure_GPa`; `backend` (`auto` = pymatgen→ASE→built-in, or `pymatgen`/`ase`/`builtin`); `manifest` (default `manifest.toml`/`.json` in the directory: per-stem `material_id`, `pressure_GPa`, `enthalpy_eV_per_atom`, `spacegroup`, `properties`, `spectra_files`, `band_coupling`, `synthetic`); `id_prefix`. Files: `*.cif`, `*.vasp`, `*.poscar`, `POSCAR*`, `CONTCAR*`, `*.xyz`/`*.extxyz`, `*.res` |
| `materials_project` | `elements`, `chemsys`, `exclude_elements`, `energy_above_hull_max` (0.1), `is_metal`, `magnetic_ordering` (`NM`/`FM`/`AFM`/`FiM`), `num_elements` ([min,max]), `material_ids`, `fields`, `extra_fields`, `limit`, `api_key_env` (`"MP_API_KEY"`; the key itself is never read from the config), `pressure_GPa` (target pressure for DFT). The query, the fields, the database version and the client version are stored in provenance |
| `airss` | `path` (`*.res`; TITL enthalpy is eV per cell, divided by n_ions), `recursive`, `pressure_GPa`, `id_prefix` |
| `uspex` | `path`, `poscars` (`gatheredPOSCARS`), `individuals` (`Individuals`), `enthalpy_unit` (`eV/cell`, the documented assumption; check it for your version), `pressure_GPa` (else `ExternalPressure` from `INPUT.txt`), `id_prefix` |
| `calypso` | `path`, `table` (`Analysis_Output.dat`), `enthalpy_unit` (`eV/atom`; best effort, check it), `pressure_GPa`, `id_prefix` |
| `seed_list` | `path` (CSV with `material_id, formula, pressure_GPa, ...`; extra columns become properties. JSON or TOML lists are also accepted), `id_prefix` |
| `prototype_substitution` | `path` (directory of prototypes), `substitutions` (`{La = ["Y","Sc"]}`), `max_per_prototype`, `include_original`, `pressure_GPa`. The output structures are hypothetical and unrelaxed |
| `synthetic_demo` | none. SYNTHETIC demo candidates |

Recognised candidate `properties`: `energy_above_hull`, `energy_above_hull_pressure_GPa`,
`is_metal`, `band_gap`, `ordering`, `total_magnetization`, `magnetization_per_atom`, `nsites`,
`dos_at_fermi` (states/eV/atom), `fermi_energy_eV`, `electronic_bandwidth_eV`, `n_fermi_sheets`,
`phonon` (`{min_frequency, unit}` or `{frequencies, qpoints, unit}` or `{band_yaml}`),
`target_pressure_GPa`. Attached spectra: `spectra_files = [{path, units, format = "two_column" |
"qe_a2f_dos" | "epw", columns, rule, label, negative_policy, lambda_file, expected_lambda,
expected_A_K}]`. Attached files are read and verified by team DFT's
`rtsc.dfpt.load_attached_spectra`; there is no fallback reader. It checks:

- QE: matdyn's footer must be present and consistent, plus the optional `lambda` file.
- EPW: the reported λ.
- Two-column files: the optional `expected_lambda` / `expected_A_K`. Without them, the verdict
  carries the assumption `A-FILE-COMPLETE(unverified)`.

An incomplete or inconsistent file makes that material ERROR.

**Input-quality tags.** Team DFT marks problems in a spectrum's `provenance.assumptions`, and all
of these tags are copied into `Verdict.assumptions`. Some tags make one direction of proof
unsafe:

| tag | blocks |
|---|---|
| `A-FILE-COMPLETE(unverified)` | REJECTED |
| `A-FILE-INCOMPLETE...` | REJECTED and CERTIFIED |
| `A-DROP-NONPOS(m K)` with m > 0 | REJECTED |
| `A-NEG-CLIP(m K)` with m > 0 | CERTIFIED |

If any alternative spectrum carries a tag that blocks the status certify returned, the pipeline
downgrades the verdict to UNDECIDED. It keeps the certificate and adds the reason
"... blocked: <tag>". To keep the proven status anyway, set `[certify]
allow_unverified_rejection = true` or `allow_unverified_certification = true`; the tag is then
recorded as an assumption. The dropped mass appears in the `dropped_nonpos_mass_K` column of
`candidates.csv`.
- The discretisation rule and the file are recorded as the assumption
  `DISCRETISATION(rule, file)`.

Material ids are mapped to directory names injectively: a hash is appended whenever an id
contains unsafe characters. Ids that would still collide (for example on a case-insensitive
file system) are refused.

### `[dedupe]`
`method` (`auto` | `structure_matcher` | `composition_spacegroup` | `none`), `ltol` (0.2),
`stol` (0.3), `angle_tol` (5.0), `symprec` (0.1, spglib), `pressure_tol_GPa` (1.0).
Duplicates become DEFERRED ("duplicate of X"). The representative is the member with the lowest
enthalpy.

### `[filters.hull]`, `[filters.metallic]`, `[filters.magnetic]`, `[filters.dynamic]` (S3)
Every filter has `enabled` (true). When a filter has no applicable data it reports "not
evaluated" and the material proceeds.

* hull (`X-HULL`): `max_e_above_hull_eV_per_atom` (0.1), `pressure_tol_GPa` (1.0),
  `include_run_candidates` (true: inputs with enthalpies at the same pressure form the hull),
  `reference_file` (CSV/JSON with `formula, enthalpy_eV_per_atom, pressure_GPa`), `backend`
  (`auto` = pymatgen PhaseDiagram, else a scipy LP), `use_stored_value` (true: use a stored
  `energy_above_hull` only when it refers to the candidate's target pressure; MP values refer to
  0 GPa). Elemental references are required. Do not mix energies from different functionals or
  codes.
* metallic (`X-METAL`): `max_band_gap_eV` (0.0). Uses `band_gap`, else `is_metal`.
* magnetic (`X-MAGNETIC`): `max_abs_magnetization_muB_per_atom` (0.05), `exclude_orderings`
  (`["FM","AFM","FiM"]`).
* dynamic (`X-DYNAMIC`): `imaginary_tol` (0.25), `imaginary_tol_unit` (`THz`; also `cm-1`,
  `meV`, `K`), `ignore_acoustic_gamma` (true), `gamma_q_tol` (0.05). Applied to supplied
  first-principles phonons and, after DFPT, to the minimum over the matdyn grid. ML phonons are
  never used here.

### `[priority]` (ordering only, S2)
Weights `hydrogen` (1.0), `light` (0.5), `phonon_scale` (1.0), `dos` (0.5), `spectrum` (1.0),
`ml` (0.5); `light_elements` (`["H","Li","Be","B","C","N"]`); `dos_scale` (0.5); `a0`.
The score is the weighted mean of the components that are available.

### `[budget]` (DEFERRED only, S2)
`max_materials` (unlimited), `min_score` (none).

### `[ml]` (heuristic, S2)
`enabled` (false), `model` (`mace` | `chgnet`), `defer_if_ml_unstable` (false), `relax` (true),
`fmax` (0.05), `steps` (300), `supercell` (2), `mesh` (4), `imaginary_tol_THz` (0.3),
`model_path`, `device`.

### `[spectra]`
`use_attached` (true), `rel_envelope` (0.10; the envelope used is max(the spectrum's own, this)),
`two_column_rule` (`bins` | `trapezoid` | `midpoint`), `negative_policy` (`error` | `clip`).

### `[dfpt]` (team DFT, `rtsc.dfpt.EPhWorkflow`)
`enabled` (false), `also_when_attached` (false).
`[dfpt.launcher]`: `kind` (`plain` | `mpirun` | `srun`), `nprocs`, `npool`, `mpirun`, `srun`,
`mpi_args`, `bin_dir`, `nice`, `env`, `omp_threads`, `serial_codes`.
`[dfpt.pseudos]`: `dir`, `sssp_json`, `files`, `cutoffs_Ry`.
`[dfpt.settings]`: any `rtsc.dfpt.EPhSettings` field, for example `do_relax`, `k_spacing`,
`q_spacing`, `el_ph_sigma_Ry`, `el_ph_nsigma`, `use_smearings`, `ecutwfc_Ry`, `timeouts_s`.
`pressure_GPa` defaults to the candidate's target pressure. Each smearing becomes one
alternative spectrum.

### `[certify]` (team CERT; recorded in every verdict)
`allow_unverified_rejection` / `allow_unverified_certification` (false), `T_target_K` (300; must be finite and > 0), `mu_hi` (0.15; finite, ≥ 0), `omega_c_factor`
(10; ≥ 1: omega_c = factor × max ω of each spectrum), `omega_c_K` (a fixed omega_c in kelvin
that overrides the factor; it must satisfy omega_c ≥ π·T\*), `cutoff_rule` (`"block"`, the
only rule defined), `rejection_model` (`"isotropic"` | `"multiband"`, see §1). Booleans and
strings are rejected. `rtsc certify` applies the same validation. The cutoff rule and omega_c are always written into
`Verdict.assumptions` (F3). The remaining keys are passed unchanged to
`rtsc.certify.CertifyParams`, and an unknown key is an error. They are: `mode`
(`"bracket"`: the tightest certified bracket; `"decide"`: stop as soon as the status is
decided, which is much faster), `row_check` (`"mhat"`, the Lean `CellOK` check), `rel_tol`
(1e-4), `N_cap` (2400), `time_budget_s` (120 s per spectrum and per search; when it runs out
only the refinement stops, never the decision) and `verbose`.

### `[validity]` (flags only, S6)
`enabled`, `nonadiabatic_ratio` (0.1), `anharmonic_h_fraction` (0.25), `anharmonic_pressure_GPa`
(50), `soft_mode_THz` (1.0), `magnetic_moment_warn_muB` (0.005), `stoner_dos` (1.5),
`anisotropy_ratio` (2.5), `mu_star_reference` (0.2), `spectrum_spread_rel` (0.1),
`lambda_large` (3.0).

### `[report]`
`outdir` (= workdir), `top_n` (25).

## 6. Running on a cluster

**One node, several materials.** Set `[run] workers = W` and `[dfpt.launcher] kind = "mpirun"`,
`nprocs = P`. W materials run at the same time and each one starts its own `mpirun -np P`.
Choose W·P ≤ the number of cores.

**Slurm job arrays.** Stages 1–5 are cheap and deterministic, so every array task repeats them
and then computes its own share:

```sh
rtsc run -c config.toml --plan-only            # optional; required once if [ml] is enabled
# array.sbatch
#SBATCH --array=0-31 --nodes=1 --ntasks=64
rtsc run -c config.toml --shard ${SLURM_ARRAY_TASK_ID}/32    # uses [dfpt.launcher] kind="srun"
# afterwards
rtsc collect -c config.toml
```

Shard K computes positions K, K+N, ... of the global priority order and writes
`inputs_shard_K_of_N.json`. `collect` merges the union of all these input lists.

- If any input list was written under a different configuration hash, `collect` **refuses**
  (ConfigError): the verdicts would be relabelled with the wrong T\* or μ\*.
- A material whose record is missing, unfinished, or made under another configuration is
  reported as ERROR instead of being dropped, so S4 still holds. You can re-run a shard at any time. Finished materials are
skipped, and QE steps restart from their own checkpoints (team DFT's resumable workflow).
Records are written atomically.

## Discovery loop

`rtsc discover` runs an autonomous generate → ML pre-screen → dedupe → DFT → certify loop
(spec: `DISCOVERY.md`). Every rule of the soundness contract still applies:

- A material **satisfies a goal** only when `certify.classify` returned CERTIFIED_CANDIDATE at
  that goal's T\* for that material at a pressure inside the goal's range. A certificate at 77 K
  never counts for the 300 K goal.
- Generators, the ML pre-screen, the search policy and budgets are HEURISTIC. They only decide what
  is generated and in which order it is computed, or mark a material DEFERRED with a recorded
  reason. The ledger refuses any other status from them (S1/S2/S3 guard in
  `rtsc.discovery.ledger`). ML results are kept in provenance, never in the properties that the
  existence filters read.
- Every generated structure carries the flag HYPOTHETICAL_STRUCTURE.

```sh
rtsc discover run -c examples/discovery/config.toml        # SYNTHETIC demo: no DFT, no ML
rtsc discover continue -w runs/discovery_demo              # resume; fresh stop limits
rtsc discover continue -w runs/discovery_demo --certified-per-goal 3 --max-rounds 20
rtsc discover status -w runs/discovery_demo
rtsc discover list   -w runs/discovery_demo --goal G-RT --status UNDECIDED
rtsc discover export -w runs/discovery_demo [--finetune [--include-synthetic]] [--out DIR]
touch runs/discovery_demo/STOP                             # stop after the current round
```

The demo attaches **SYNTHETIC** α²F spectra to substituted prototypes. They are invented
numbers, not data for the formulas shown. Material ids start with `SYNTH-`, verdicts carry
SYNTHETIC_DATA, and the report prints a banner.

**Round.** Each round runs these phases. Each phase is committed atomically to
`<workdir>/ledger.sqlite`.

1. propose: the bandit policy picks the requests;
2. generate: the structure is validated and fingerprinted, then deduplicated against the whole
   ledger. Duplicates are linked and never re-run;
3. ML pre-screen;
4. select: the top candidates within `dft_per_round` get a DFT job. The rest stay DEFERRED;
5. execute: the existing pipeline runs once per goal in `rounds/<k>/<goal>/`;
6. record: the policy snapshot is saved and the reports are written.

A material runs through the DFT stage once, in the pipeline config of the first enabled goal
whose pressure range contains it. Every further in-scope goal re-runs only `certify.classify`
(at its own T\*) on the same alternative spectra, through the pipeline's own certify-only resume
path. Materials outside a goal's range get DEFERRED (stage `scope`) for that goal.

**Counting, config changes, ledger invariants.** These are specified in `DISCOVERY.md` §3:

- a single predicate `counts_for` decides whether a material satisfies a goal. It checks T\*, the
  material and DFT pressures, and the current `certify_hash` and `dft_hash`;
- non-counting CERTIFIED rows are shown as "does not count";
- at the start of every session the loop reconciles the ledger with the current config. It
  re-certifies on stored spectra, or re-queues the material for DFT;
- `[dfpt.settings]` must not override the pressure;
- append-only and the S1/S2/S3 guard are enforced by SQLite triggers;
- the Slurm array is submitted with `--hold`, its id is recorded, then it is released.

**Crash safety.** A crash loses at most the phase in flight. `continue` resumes the unfinished
round at its phase. The pipeline reuses a record only if its material id, configuration hash and
input digest all match, so finished materials are never recomputed. Jobs with no finished record
become retryable ERRORs and are computed again later.

**Outputs.** These are rewritten every round:

- `discovery_report.md`: assumption banner, goals, per-round progress, certified candidates per
  goal, the open (UNDECIDED) list ranked by certified T_up, and rejected materials with their
  certificates;
- `candidates_<goal>.csv`: every material with its current status for that goal;
- `rounds/<k>/<goal>/`: the standard pipeline reports and certificates.

The ledger also holds:

- materials and duplicate links;
- the ML pre-screen results;
- the DFT jobs;
- the append-only verdict history, with exactly one current status per material and goal;
- rounds, policy snapshots and sessions.

**Configuration.** The TOML is a normal pipeline config plus a `[discovery]` table. `[run]` and
`[[sources]]` are ignored. `[filters]`, `[spectra]`, `[dfpt]`, `[certify]` (except `T_target_K`,
which each goal sets), `[validity]` and `[priority]` apply to every round.

| key | default | meaning |
|---|---|---|
| `discovery.workdir` | `discovery_run` | ledger, reports, `rounds/` |
| `discovery.seed` | 0 | policy RNG seed; with the ledger, it fully determines the choices |
| `discovery.executor` | `local` | `local` (in-process; `workers` processes) or `slurm` |
| `discovery.id_prefix` | `D` | material ids `<prefix><seq>-<formula>` (`SYNTH-` is prepended in demo mode) |
| `discovery.goals.<id>` | G-RT, G-AMB | `T_star_K`, `P_max_GPa`, `pressures_GPa` (pressure bins), optional `chemical_systems`, `space_groups`, `enabled`, `ambient` |
| `discovery.search` | | `chemical_systems` (default seeds for every goal), `space_groups`, `policy` (`thompson` or `ucb`), `requests_per_round` (4), `structures_per_request` (4), `exploration_floor` (0.25), `ucb_c` (1.0) |
| `discovery.generators.<name>` | required | `type` = `prototype_substitution` (`prototypes` = directory, `keep_common`), `directory` (`path`, `format` = local/airss/uspex/calypso, `pressure_tol_GPa`), `mattergen` (MatterGen config keys) or `python` (`class = "module:Name"`, plus constructor kwargs) |
| `discovery.prescreen` | `null` | `type` = `null`, `mattersim` or `python`; `defer_if_ml_unstable` (true), `max_e_above_hull_ml_eV` (0.2), `use_ml_relaxed` (true: the ML-relaxed geometry is the DFT starting point) |
| `discovery.budget.dft_per_round` | 4 | materials sent to the DFT stage per round |
| `discovery.stop` | | `certified_per_goal` (1, counted per session), `max_rounds` (10), `max_wall_hours`, `max_dft_jobs`; a `STOP` file in the workdir also stops the loop |
| `discovery.dedupe` | | `use_structure_matcher` (true when pymatgen is installed), `ltol`, `stol`, `angle_tol`, `confirm_tol_A` (0.02, the absolute pair-distance confirmation when pymatgen is missing). Pressures are never binned |
| `discovery.slurm` | | `partition`, `account`, `time`, `nodes`, `ntasks`, `cpus_per_task`, `gpus`, `extra` (raw `#SBATCH` lines), `setup` (shell lines), `max_array` (64), `poll_seconds` (30), `timeout_hours`, `python`, `sbatch`, `squeue`, `sacct` |
| `discovery.synthetic_spectra` | off | DEMO ONLY: `enabled`, `A_range_K`, `center_range_K`, `width_frac` |

Goals:

- G-RT: T\* = 300 K, P ≤ 200 GPa.
- G-AMB: T\* = 77 K (configurable up to 300 K), P ≤ 1e-4 GPa. A larger P_max is a configuration
  error.

**Policy.** The arms are goal × chemical system × pressure bin × space group × generator. Each
arm has a Beta posterior. The reward comes from DFT-level outcomes:

- CERTIFIED gives 1;
- UNDECIDED gives 0.6 + 0.3·min(1, T_low/T\*);
- REJECTED gives 0.5·min(1, T_up/T\*);
- EXCLUDED gives 0.05;
- ERROR gives 0.

ML-unstable deferrals and requests that produced nothing count as failures. The exploration floor
reserves ceil(floor·n) requests per round for arms never tried. The statistics are recomputed from
the ledger, and the RNG is seeded by (seed, round). So a continued run makes the same choices as an
uninterrupted one.

**Slurm.** The `slurm` executor works in three steps:

1. it writes the plan (`rtsc run --plan-only` semantics);
2. it submits one `sbatch --array` whose task K runs
   `python -m rtsc.discovery.executors shard <job.json> K N`. That is `pipeline.run(shard="K/N")`;
   `[dfpt.launcher]` sets the MPI launch inside each task;
3. it polls `squeue` until the array has left the queue, reads the states from `sacct`, then
   calls `pipeline.collect`.

On `continue`, a recorded job id that is still queued is waited for, not resubmitted.

**Fine-tuning export.** `export --finetune` writes the DFT-level labels (structure → A, λ, ω_log,
certified bracket, status) in MatterGen's CSV dataset layout, via
`rtsc.generate.finetune_export`. SYNTHETIC materials are excluded unless `--include-synthetic` is
given. The export is documented, not run.

## 7. Development notes

* The orchestration tests (`tests/test_pipeline*.py`) need no network, no QE and no pymatgen.
  While `rtsc.certify.classify` does not exist, they patch in a clearly marked, non-rigorous test
  stub (`tests/test_pipeline_stubs.py`). That stub never ships in the package.
* Heavy dependencies (pymatgen, mp-api, ASE, phonopy, MACE/CHGNet) are imported lazily.

## ML adapters (MatterGen generator, MatterSim pre-screen)

Both adapters are optional. The core runs and tests without them. Everything they produce is
**HEURISTIC** (ARCH S2). A generated structure always carries `HYPOTHETICAL_STRUCTURE`. An ML
pre-screen result can only reorder candidates or let the loop mark them DEFERRED with a recorded
reason. The adapters have no API that sets a Status. They never produce EXCLUDED_ASSUMPTION,
REJECTED_PROVEN or CERTIFIED_CANDIDATE.

### Sources checked (2026-10-27)

| What | Source |
|---|---|
| MatterGen | github.com/microsoft/mattergen @ `92423660` (README, MODEL_CARD, `mattergen/scripts/generate.py`, `checkpoints/*/config.yaml`, `common/utils/globals.py`, `common/data/dataset.py`, `conf/`); PyPI `mattergen` 1.0.3; MIT licence |
| MatterSim | github.com/microsoft/mattersim @ `5c84771` (README, MODEL_CARD, `docs/examples/*`, `forcefield/potential.py`, `applications/{relax,phonon}.py`, `utils/download_utils.py`); PyPI `mattersim` 1.2.5; MIT licence |
| phonopy | installed version 4.8.6 (attribute API `.mesh`, `.total_dos`, `.projected_dos`) |

### Two environments

MatterGen and MatterSim cannot share one interpreter, so they run in separate environments.

- **Main environment (Python 3.12+).** This holds rtsc and MatterSim; the deliverable's
  `environment.yml` uses Python 3.12. MatterSim requires Python ≥ 3.12 and numpy ≥ 2. Install it
  with `pip install mattersim`, which pulls in torch, torch_geometric, phonopy, pymatgen and ase.
- **MatterGen environment (Python 3.10, Linux, CUDA GPU).** The official README installs it with
  `uv venv .venv --python 3.10; uv pip install -e .` or
  `uv pip install mattergen --find-links https://data.pyg.org/whl/torch-2.2.0+cu121.html`. Its
  metadata pins `torch==2.2.1+cu118` and `numpy<2`, so it cannot be installed into the main
  environment.

The adapter never imports `mattergen`. It runs the `mattergen-generate` executable as a
subprocess and parses the files that executable writes.

### `rtsc.generate.mattergen.MatterGenGenerator`

**Command line.** The adapter builds this command:

```
<command_prefix...> <executable...> OUTDIR
    --pretrained_name=<name> | --model_path=<dir> --checkpoint_epoch=last
    --batch_size=B --num_batches=K --seed=S --record_trajectories=False
    [--properties_to_condition_on={'chemical_system':'H-La','space_group':225}
     --diffusion_guidance_factor=2.0]
```

- `executable` is set in the config. Examples: `["/opt/mattergen/.venv/bin/mattergen-generate"]`
  or `["conda","run","-n","mattergen","mattergen-generate"]`.
- `command_prefix` is also set in the config, e.g. `["srun","-p","gpu","--gres=gpu:1"]`.
- Set `shell=true` for prefixes such as `ssh host` that pass the command through a remote shell.

**Pretrained checkpoints and their conditioning properties.** Checkpoints come from Hugging Face
`microsoft/mattergen` (by default) or from the repo's Git LFS. This table reproduces
`property_embeddings_adapt` from each checkpoint's `config.yaml`:

| checkpoint | conditions on |
|---|---|
| `mattergen_base`, `mp_20_base` | nothing (unconditional) |
| `chemical_system` | `chemical_system` (`"El1-El2"`) |
| `space_group` | `space_group` (1..230) |
| `dft_mag_density` | `dft_mag_density` |
| `dft_band_gap` | `dft_band_gap` |
| `ml_bulk_modulus` | `ml_bulk_modulus` |
| `dft_mag_density_hhi_score` | `dft_mag_density`, `hhi_score` |
| `chemical_system_energy_above_hull` | `chemical_system`, `energy_above_hull` |

**Mapping a request to conditioning.**

- A request's chemical system and space group, plus any `extra["mattergen_conditioning"]`
  values, are passed to the model only when the checkpoint supports them.
- Anything else raises `UnsupportedConditioningError`. The exception is a field listed in
  `allow_unconditioned`: it is dropped, the drop is recorded, and the outputs are filtered to the
  requested elements.
- Elements outside MatterGen's 76 `SELECTED_ATOMIC_NUMBERS` are refused. These include the noble
  gases, Tc, Pm and Z ≥ 84.
- For a fine-tuned `model_path`, the supported properties are read from `<model_path>/config.yaml`.

**Guidance.** `diffusion_guidance_factor` is the classifier-free guidance γ. A value of 0 means
unconditional sampling; larger values follow the conditions more closely but give less diverse
samples. The README examples use 2.0, which is the adapter's default.

**Outputs read.**

- `generated_crystals.extxyz` is preferred.
- The fallback is `generated_crystals_cif.zip`, which holds `gen_<i>.cif`.
- `generated_trajectories.zip` is not used (trajectories are switched off).

**Provenance and reruns.** Each structure records:
- the model and checkpoint
- the conditioning, plus any fields that were dropped
- the guidance factor, seed, batch sizes and command line
- the output file and the structure's index in it

A finished run with the same command line is reused rather than rerun.

**Hardware.** The MODEL_CARD reports that sampling 1,000 structures takes about 2 h on one V100
GPU.

**Pressure.** MatterGen has no pressure conditioning. Its training sets (MP-20 and Alex-MP-20)
are ambient-pressure DFT relaxations. Generated structures are therefore ambient-like proposals,
and the pre-screen relaxes them at the target pressure. Each structure's provenance says this.

### `rtsc.mlff.mattersim.MatterSimPrescreener`

**Model.** The calculator is created with
`from mattersim.forcefield import MatterSimCalculator`, then
`MatterSimCalculator(load_path="MatterSim-v1.0.0-1M"|"MatterSim-v1.0.0-5M"|<path>, device=, dtype=)`.
Both checkpoints use the M3GNet architecture: 1M is faster and 5M is more accurate. They are
cached in `~/.local/mattersim/pretrained_models` and downloaded from raw.githubusercontent.com
when missing. The MODEL_CARD gives the stated range as 0–5000 K and up to 1000 GPa, for bulk
materials only.

**Relaxation.** Cell and positions are relaxed with an ASE `FrechetCellFilter` (or
`ExpCellFilter`) using `scalar_pressure = P·ase.units.GPa`. The optimiser is FIRE by default, or
BFGS. The result records the enthalpy per atom H = (E + PV)/N and the residual pressure computed
from the stress.

**ML hull (`rtsc.mlff.hull`).** Every elemental, binary and higher sub-system is relaxed with the
**same** potential at the **same** pressure.
- Reference phases come from `reference_dir` (any format ASE can read), from Materials Project
  when `use_mp` is set (0 GPa phases re-relaxed at P), or from a recorded elemental fallback.
  The fallback is the `ase.build.bulk` reference state, or a crude dimer crystal for H, N, O
  and other diatomic elements.
- Results are cached as JSON, one file per sub-system, pressure and model.
- Missing high-pressure polymorphs or crude fallback references put the hull too high, which
  makes `e_above_hull_ml` optimistic. For a heuristic that can only defer, this is the harmless
  direction.

**ML phonons (`rtsc.mlff.phonons`).**
- Finite displacements in phonopy (0.01 Å). The supercell is either `"auto"` (shortest length
  ≥ 7 Å, capped at 160 atoms) or given explicitly.
- Force constants are symmetrised. The mesh is Γ-centred and the three acoustic modes at Γ are
  ignored.
- `dynamically_stable_ml` is false when the lowest frequency is below −`imaginary_tol_cm1`
  (default 15 cm⁻¹).
- The result includes the total DOS (in cm⁻¹), the fraction of imaginary modes, and the q-point
  of the minimum.

**Optional heuristic score.** `MatterSimPrescreener.heuristic_score(result)` returns the fraction
of the phonon DOS that sits on light elements (H–F) above 500 cm⁻¹. A phonon-DOS ω_log (not
α²F) is also reported. Both are for ordering only.

**Suggested deferrals.** Notes beginning `ML_SUGGEST_DEFER:` mark cases where the ML hull is above
`suggest_defer_e_hull_eV` (0.2 eV/atom) or where ML imaginary modes were found. The loop owns the
decision and records the reason.

### Fine-tuning export (`rtsc.generate.finetune_export`)

`export_from_ledger(ledger.sqlite, outdir)` writes `outdir/<name>/{train,val}.csv`. This is the
same layout as MatterGen's `datasets/mp_20/*.csv` and is consumed by `csv-to-dataset`.

**Columns.**
- `material_id`, `cif` and `chemical_system`.
- Float labels: `rtsc_A_K`, `rtsc_lambda`, `rtsc_omega_log_K`, `rtsc_Tc_low_K`, `rtsc_Tc_up_K`
  and `rtsc_pressure_GPa`.
- Metadata: `rtsc_status`, `rtsc_goals` and `rtsc_formula`.

**Scope.** Only verdicts from the DFT `pipeline` stage are exported (REJECTED_PROVEN,
CERTIFIED_CANDIDATE and UNDECIDED). When a material has several goals, their brackets are
intersected.

**Templates.** The export also writes `mattergen_conf/property_embeddings/<prop>.yaml` (copied
from `dft_mag_density.yaml`), `mattergen_conf/data_module/<name>.yaml`, and
`FINETUNE_INSTRUCTIONS.txt`.

**Fine-tuning steps.** These follow the README section "Fine-tune on your own property data":
1. Add the new property names to `PROPERTY_SOURCE_IDS`.
2. Run `csv-to-dataset`.
3. Run `mattergen-finetune ... data_module.properties=["rtsc_A_K"]`.

This export is documented but has not been run here.

### What was run in this environment (2026-10-27; 2 CPU cores, no GPU)

**MatterSim: real runs.** These used `MatterSim-v1.0.0-1M` on CPU (mattersim 1.2.5,
torch 2.14.1, phonopy 4.8.6, Python 3.13 venv; the checkpoint was downloaded from
raw.githubusercontent.com).

| Structure | Pressure | Relaxed lattice | Phonons (2×2×2 supercell, 8³ mesh) | ML hull |
|---|---|---|---|---|
| fcc Al | 0 GPa | a = 4.035 Å | stable; DOS ends at 317 cm⁻¹ | 0 vs the fcc reference |
| fcc Al | 10 GPa | a = 3.906 Å; stress-derived P = 9.79 GPa | not computed | not computed |
| MgB2 | 0 GPa | a = 3.077 Å, c = 3.523 Å | stable; DOS ends at 778 cm⁻¹; light-element HF fraction 0.44 | 0 vs fallback Mg (hcp) and B (fcc, crude) |
| H3S Im-3m (from examples) | 150 GPa | relaxed to 149.4 GPa | **unstable**: min −450 cm⁻¹ at the H point, q = (−½,−½,½) in phonopy's primitive reciprocal basis | — |

- Wall time was about 20 s to load the model, then seconds per structure.
- The H3S result shows why ML results are heuristic. The harmonic ML instability may be real
  (harmonic Im-3m H3S is reported unstable below ~180 GPa) or a model artefact. Either way, under
  S2 it can only defer the structure.

**MatterGen: not run.** It is blocked here for these reasons:
1. huggingface.co, download.pytorch.org and data.pyg.org are refused by the egress proxy.
2. `uv pip install mattergen` is unsatisfiable, because the pinned `torch==2.2.1+cu118` exists
   only on download.pytorch.org.
3. There is no GPU.

The checkpoints themselves would be reachable through media.githubusercontent.com (Git LFS);
this was verified with `data-release/mp-20/mp_20.zip`. The adapter was tested against a **fake**
`mattergen-generate` script written by the tests. That script emits the documented file names
with placeholder structures, which are not MatterGen samples. The ASE-written extxyz and CIF-zip
formats that MatterGen uses (`eval_utils.save_structures`) were also parsed in tests.

**Tests.**
- `tests/test_generate_mattergen.py` and `tests/test_generate_finetune.py`.
- `tests/test_mlff_prescreen.py` uses ASE's EMT potential as a stand-in calculator, so its numbers
  are not MatterSim numbers.
- `tests/test_mlff_mattersim_real.py` uses real MatterSim and is skipped unless `mattersim` is
  installed.
