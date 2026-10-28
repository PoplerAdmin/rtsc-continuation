# Discovery loop: autonomous generate → screen → certify → continue

This extends ARCH.md. Every rule in ARCH §3 (the soundness contract) remains binding.

## 1. Goals (config `[discovery.goals]`, any number, each with an id)

| Goal | Default | Target |
|---|---|---|
| G-RT | on | certified Tc_model >= T_rt (default 300 K) at some pressure P <= P_max (default 200 GPa) |
| G-AMB | on | certified Tc_model >= T_amb (default 77 K, configurable up to 300 K) at ambient pressure, i.e. P <= 1e-4 GPa |

- A material **satisfies a goal** only when certify.classify returns CERTIFIED_CANDIDATE at that goal's T* and pressure.
- UNDECIDED materials are listed separately as "open", ranked by certified T_up.
- REJECTED_PROVEN materials are listed with their certificate.

## 2. Loop (`rtsc discover run|continue|status|list|export`)

Each round does the following:
1. **Propose.** The search policy picks generation requests: chemical system, optional space group, pressure, goal, count.
2. **Generate.** A generator adapter produces structures. Supported generators: MatterGen, prototype substitution, CSP-output adapters. Every generated structure carries the HYPOTHETICAL_STRUCTURE flag.
3. **ML pre-screen.** The MatterSim adapter relaxes each structure at the request's pressure. It then computes an ML energy above the hull, against reference phases relaxed with the same potential, and an ML phonon dynamic-stability check. These results are HEURISTIC. Under S2 they may only reorder candidates or mark them DEFERRED, with a recorded reason. They never yield EXCLUDED_ASSUMPTION or a proven status.
4. **Dedupe.** Structures are deduplicated against the whole ledger (all rounds), using a structure fingerprint. A duplicate is linked to its existing entry, never re-run.
5. **Select.** Within the per-round DFT budget, the top-priority candidates go to DFT. The rest stay DEFERRED in the ledger and can be picked up in later rounds.
6. **DFT stage.** This is the existing pipeline: DFT existence filters (EXCLUDED_ASSUMPTION per S3), then DFPT, then certify.classify, then validity flags.
7. **Record.** Results go into the ledger. The policy is updated from the outcomes, using only heuristics.
8. **Stop check.** The loop stops when `stop.certified_per_goal` is reached for every enabled goal, or when `max_rounds`, `max_wall_hours` or `max_dft_jobs` runs out, or when a STOP file is present. The report and candidate list are written after every round.

`continue` resumes from the ledger with fresh stop limits. Nothing already computed is redone, and nothing is lost.

## 3. Ledger (SQLite, `<workdir>/ledger.sqlite`; append-only for verdicts)

- **Tables:** `materials` (fingerprint, structure, source, round, generation request), `duplicates`,
  `ml_screen`, `dft_jobs`, `verdicts` (status, bracket, rule ids, assumptions, flags, certificate
  path, goal, T*, material pressure, DFT pressure used, stage, origin, input digest, pipeline config
  hash, `dft_hash`, `certify_hash`), `rounds`, `policy_state`, `sessions`.
- **Append-only, in SQL.** Triggers reject UPDATE, DELETE, and an INSERT whose id already exists
  (`recursive_triggers` is on, so INSERT OR REPLACE is rejected as well).
- **Guard, in SQL** (trigger `verdicts_guard`; `append_verdict` checks the same matrix first). A
  verdict row is accepted only if:
  - stage = `pipeline`, or status = DEFERRED, or (status = ERROR and stage is `generate` or
    `executor`);
  - REJECTED_PROVEN / CERTIFIED_CANDIDATE ⇒ origin = `rtsc.certify.classify`;
  - EXCLUDED_ASSUMPTION ⇒ origin `rtsc.stability*`;
  - CERTIFIED_CANDIDATE ⇒ in_scope = 1.
- **Current status = the verdict with the largest id per (material, goal).** This is a view, so
  each pair has exactly one current status. History rows are kept. Every material has a row for
  every enabled goal: out-of-range goals get DEFERRED with stage `scope`.
- **Counting (one predicate, `goals.counts_for`).** A material satisfies goal G iff its current
  row for G has:
  - status CERTIFIED_CANDIDATE, stage `pipeline`, origin `rtsc.certify.classify`;
  - T* = G.T*;
  - material pressure and DFT pressure both inside G's current range, and equal;
  - `certify_hash` = the current hash of `[certify]` with T* = G.T*;
  - `dft_hash` = the current hash of every non-certify pipeline section.

  Reports, CSVs, `status` and `list` use this predicate. A CERTIFIED row that fails it is shown
  as "CERTIFIED_CANDIDATE (does not count)" with the reason.
- **DFT pressure.** `[dfpt.settings]` may not set `pressure_GPa` / `press_kbar`, `do_relax = false`
  or a fixed-volume `cell_dofree`: each material is computed at its own pressure. A proven
  verdict whose spectra were computed at another (or unknown) pressure is not recorded. An
  ERROR with flag DFT_PRESSURE_MISMATCH is recorded instead.
- **Config changes (checked at the start of every `run` / `continue`).** The config of session N
  is stored in meta `config_session_N`.
  - Scope changes re-assign `scope` / awaiting statuses.
  - If `dft_hash` changed, the material gets a DEFERRED with stage `config` and is re-queued for
    DFT.
  - If a goal is new, or its T* or `[certify]` changed, `certify.classify` re-runs on the stored
    spectra (certify only; the rows are appended).
- **Reconcile is idempotent.** A (material, goal) pair is left alone when its current row is final
  under the current hashes. That means a pipeline verdict (classify-origin ones also need the
  current `certify_hash`), or a non-retryable executor ERROR (e.g. DFT_PRESSURE_MISMATCH) made
  under the current `dft_hash` and `certify_hash`.
- **Accounting.** Each round records its number of generator outputs (in the same transaction as
  the material and duplicate rows). It must equal the material rows plus the duplicate rows created
  in that round.
- **ACID writes.** Each round phase is committed atomically, so a crash mid-round loses at most
  in-flight work. Unfinished jobs are retried on `continue`.
- **Verdict reuse.** A verdict is reused only when its input digest and config hash match (as in
  pipeline.py resume).
- **Dedupe.** The fingerprint (cell-choice invariant, exact pressure to 1e-6 GPa) only proposes
  candidates. A duplicate needs confirmation, by StructureMatcher or by absolute pair-distance
  spectra out to the reduced-cell size. An unconfirmed hit becomes a new material, noted
  `possible_duplicate_of`.

## 4. Adapters (optional dependencies; the core runs and tests without them)

- **`rtsc.generate.mattergen`**
  - Wraps Microsoft MatterGen: CLI or Python API, a configured checkpoint, and conditioning on supported properties (e.g. chemical system, space group, energy above hull, and others available in the installed version).
  - Every name used must be verified against the official repository.
  - Supports a local GPU or remote execution through the executor.
- **`rtsc.mlff.mattersim`**
  - Wraps Microsoft MatterSim: an ASE calculator for relaxation under pressure, phonons via phonopy, and ML hull energies against a reference set.
- **Fine-tuning export:** `rtsc discover export --finetune` writes the ledger's DFT-level labels (structure → A, λ, ω_log, certified bracket, status) in a format MatterGen fine-tuning accepts. That enables a property-conditioned generator later. It is documented, not run here.
- **Executors:**
  - `local`: process pool.
  - `slurm`: submits DFPT jobs or job arrays with sbatch, polls squeue/sacct and collects results. The generation step can target a GPU partition.

## 5. Search policy (heuristic; it may only influence what is generated and in which order)

- **Arms:** chemical systems × pressure bins (and optionally space groups), seeded from config.
  - Example seeds: hydrogen-rich binaries and ternaries with alkaline-earth or rare-earth elements; boron/carbon frameworks.
  - Ambient-mode seeds: light-element covalent metals.
- **Bandit** (e.g. Thompson sampling or UCB) over a reward from DFT outcomes:
  - certified or undecided T_up / T*, A_hi, and the multiband κ when available;
  - penalties for ML-unstable and ERROR outcomes.
- **Exploration floor:** a minimum fraction of each round goes to unexplored arms.
- **Determinism:** the policy is deterministic given its RNG seed and the ledger, so a continued run reproduces the same choices. This holds only when `max_dft_jobs` is not the binding limit: that budget is counted per session, so a session that hits it selects fewer materials in a round than an uninterrupted run would.

## 6. Outputs

Rewritten each round:
- `discovery_report.md`: goals, progress, top candidates per goal, open (UNDECIDED) list, and the model-assumptions banner;
- `candidates_<goal>.csv`;
- `ledger.sqlite`;
- `rounds/<k>/` holding that round's pipeline outputs, i.e. the standard reports and certificates.
