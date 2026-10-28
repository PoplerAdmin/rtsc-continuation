# rtsc bundle: a Tc ceiling for phonon-mediated superconductors, its proofs, and a screening program

> **Continuing this work on a new machine?** Read `HANDOFF.md` first. A prompt for a new Claude
> session is in `STARTER_PROMPT.txt`. The autonomous discovery loop (MatterGen → MatterSim →
> DFT → certify, run with `rtsc discover run|continue|status|list`) is described in
> `rtsc/DISCOVERY.md`. Its ledger invariants are proved in the Lean library `DiscoveryLogic`.

This bundle contains:

* **a theorem.** In isotropic Migdal–Eliashberg theory,
  `k_B Tc ≤ A / 4.33148145 ≈ 0.23086789 · A`, where `A = ∫ α²F(ω) dω`. The bound holds for every
  phonon spectrum, every μ* ≥ 0 and every Coulomb cutoff, and the constant is sharp to a
  relative 1e-8 (`proof/PROOF.md`);
* **its checks.** Rigorous ball-arithmetic certificates for the numerical step
  (`certificates/`) and Lean 4 formalisations (`lean/`) of the reduction, of the screening
  program's decision rules, and of a multiband version of the bound;
* **`rtsc`**, a Python program that screens candidate compounds against a target Tc (default
  300 K). It gives each compound a verdict whose "proved" outcomes rest on the Lean-checked rules
  and on Arb interval checks (`rtsc/`).

Contents:

1. [What is claimed and what is not](#1-what-is-claimed-and-what-is-not)
2. [Bundle layout](#2-bundle-layout)
3. [Verifying the bundle](#3-verifying-the-bundle-verify_allsh)
4. [Quickstart on a workstation](#4-quickstart-on-a-workstation)
5. [Running on a cluster](#5-running-on-a-cluster)
6. [Hardware sizing](#6-hardware-sizing-rough)
7. [Rebuilding and checking the Lean proofs](#7-rebuilding-and-checking-the-lean-proofs)
8. [Re-running the certificates for the constant](#8-re-running-the-certificates-for-the-constant)
9. [What was tested while preparing this bundle](#9-what-was-tested-while-preparing-this-bundle)
10. [Known gaps](#10-known-gaps)

---

## 1. What is claimed and what is not

### Claimed (proved)

* **Theorem** (`proof/PROOF.md`). Take isotropic Migdal–Eliashberg theory with the linearised gap
  equation, any α²F ≥ 0 with total weight A, any μ* ≥ 0 and any Coulomb cutoff. Then
  `k_B Tc ≤ A/4.33148145`, so Tc ≥ 300 K requires A ≥ 1299.44 K (about 112 meV). The best
  universal constant c* satisfies `1/4.33148149 ≤ c* ≤ 1/4.33148145`. The proof has three parts:
  * an analytic reduction to single Einstein modes;
  * a computer-assisted certificate for the Einstein family (python-flint/Arb, 128-bit balls;
    `certificates/wE1`, independently `certificates/wB`; both audited, see `proof/audits/`);
  * a sharpness witness.
* **Lean 4, library `TcR1`.** Formalises the reduction and the main theorem for *finite atomic*
  spectra. The numerical certificate enters as an explicit hypothesis `H(a0, η)`, and `2π` is a
  free real parameter κ. There is no `sorry`, and the only axioms are `propext`,
  `Classical.choice` and `Quot.sound`.
* **Lean 4, library `PipelineLogic`.** Proves that the screening program's decision rules are
  sound: the universal ceiling (R-UNIV), μ* monotonicity, the Schur/Collatz–Wielandt cell test
  (R-CELLS), the Rayleigh lower bound (C-RAYLEIGH), and the `classify` function that combines
  them.
* **Lean 4, library `MultiBand`.** Proves the multiband bound `k_B Tc ≤ κ/a0`, with κ the Perron
  root of the symmetrised band-coupling matrix.
* **`rtsc` verdicts.** `REJECTED_PROVEN` means Tc_model < T\*. `CERTIFIED_CANDIDATE` means
  Tc_model ≥ T\*. Each is proved for every spectrum in the stated weight envelope, every
  alternative spectrum (for example every DFPT smearing) and every μ\* in the stated range. The
  certified direction also depends on the configured cutoff rule. Only `rtsc.certify.classify`
  can produce these two verdicts. Its booleans are Arb checks, and its control flow mirrors the
  Lean `classify`. `rtsc recheck` re-verifies any verdict from its stored certificate data.

### Not claimed

* **Not a prediction of a measured Tc.** Every verdict concerns the model: isotropic,
  adiabatic (Migdal), harmonic, the given α²F, and constant μ\*. The model leaves out
  anharmonicity and quantum nuclear motion (both large for hydrogen), non-adiabatic and vertex
  corrections, anisotropic/multiband gaps, spin fluctuations, the accuracy of DFT/DFPT
  (functional, pseudopotentials, convergence) and synthesizability. The program marks these
  risks with *flags* (`ANHARMONIC_H`, `NONADIABATIC`, `MULTIBAND_RISK`, ...). A flag never
  changes a verdict.
* **Physical filters are not proofs.** `EXCLUDED_ASSUMPTION` comes from convex hull,
  metallicity, magnetism or imaginary phonon modes, and is not proved. `DEFERRED` means that the
  material was not computed. `UNDECIDED` means that the certified bracket contains T\*.
* **Multiband materials.** A single-band rejection uses the DOS-averaged A. It carries the
  assumption `M-ISO-SINGLEBAND` and the flag `MULTIBAND_RISK`, because the averaged A can
  underestimate Tc (a counterexample in `certificates/wMB` exceeds it by 9 %). The multiband
  rejection rule R-MULTI is formalised in Lean, but `rtsc/ARCH.md` §4c still marks it as
  *pending audit*.
* **Scope of the Lean proofs.** Lean covers finite atomic spectra only; general measures are
  proved on paper. The interval-arithmetic certificates (Prop. 6/7) and the Python code that
  implements the rules are **not** verified in Lean. Their correspondence to the Lean
  statements is documented in `rtsc/docs/CORRESPONDENCE_LEAN.md`.
* **The trusted base** is the Lean 4 kernel and Mathlib, python-flint/Arb 0.9 (ball
  arithmetic), the Python certificate scripts as written, and the Python implementation of the
  rules.
* **Synthetic data.** `rtsc demo` uses invented spectra labelled `SYNTHETIC-*`. They test the
  software and say nothing about any real material. The bundled QE outputs under
  `rtsc/tests/data_dfpt` are real program output from deliberately unconverged smoke tests and
  have no physical meaning.

---

## 2. Bundle layout

```
rtsc_bundle/
  README.md               this file
  verify_all.sh           end-to-end check (section 3)
  VERSIONS.txt            when/where the bundle was assembled, toolchain and Mathlib revisions
  MANIFEST.sha256         sha256 of every file
  EXCLUDED.txt            what was deliberately left out (.lake, caches, install logs)
  KNOWN_GAPS.txt          gaps detected automatically at assembly time
  .dockerignore           for `docker build -f env/Dockerfile .`
  env/
    environment.yml       conda env: Python stack + QE 7.x with EPW (conda-forge)
    Dockerfile            Ubuntu 24.04 + micromamba env + rtsc   (not built; see section 9)
    slurm/
      rtsc_array.sbatch   job-array template: rtsc run --shard K/N
      rtsc_collect.sbatch rtsc collect after the array
      submit.sh           submits both, collect depends on the array
      config_cluster.toml cluster config template (srun launcher, pseudos, settings)
  lean_checks/
    AllAxioms.lean        #print-axioms audit over EVERY declaration of the three libraries
  rtsc/                   the Python package (pyproject.toml, README.md, ARCH.md, docs/, tests/, examples/)
  lean/                   Lean 4 project: TcR1, PipelineLogic, MultiBand; Check*.lean; *_axioms_output.txt
  proof/
    PROOF.md              the theorem and its proof (final write-up)
    SPEC.md               the model definitions
    audits/               referee reports and break-test logs
  certificates/
    wE1/                  primary certificate for a0 = 4.33148145 (run_all.sh, CERT.md, logs)
    wB/                   independent route-B proof and certificate (PROOF.md, cert_*.py, logs)
    wMB/                  multiband bound: proof, Arb check of Proposition N, counterexamples
```

**Path map.** The component documents were written in the development workspace and use
absolute paths. In the bundle they correspond to:

| path in component docs | in the bundle |
|---|---|
| `/home/claude/rtsc` | `rtsc/` |
| `/home/claude/tc/lean` | `lean/` |
| `/home/claude/tc/final/PROOF.md` | `proof/PROOF.md` |
| `/home/claude/tc/SPEC.md` | `proof/SPEC.md` |
| `/home/claude/tc/wE1`, `wB`, `wMB` | `certificates/wE1`, `wB`, `wMB` |
| `/home/claude/tc/audit_*` | `proof/audits/audit_*` (reports and logs only) |
| `wA`, `wC`, `wD`, `wE2`, `wG`, `wH` (other routes, cited in PROOF.md §3) | not included |

---

## 3. Verifying the bundle (`verify_all.sh`)

```sh
cd rtsc_bundle
./verify_all.sh                    # python + pytest + demo + recheck + Lean   (~20-30 min on 2 cores)
./verify_all.sh --cert             # + the wE1 certificate for the constant      (+6-10 min)
./verify_all.sh --cert-all         # + route-B certificate and multiband Prop. N (+5 min)
./verify_all.sh --skip-lean        # machines without a Lean toolchain (e.g. the Docker image)
./verify_all.sh --help
```

| step | what it checks |
|---|---|
| python | `pip install --no-deps` of `rtsc/` into a throw-away venv works (offline build first; if the local setuptools cannot build it, an isolated build that fetches setuptools from the index); core imports; `rtsc` entry point exists |
| pytest | the full `rtsc` test suite (`tests/`, including the oracle cross-checks). The real-QE test runs only with `RTSC_RUN_QE=1` |
| demo | `rtsc demo` (SYNTHETIC). The status counts add up to the number of inputs, there is one verdict per input, and the SYNTHETIC banner is present |
| recheck | `rtsc recheck` re-verifies every classified demo verdict, and at least one REJECTED_PROVEN and one CERTIFIED_CANDIDATE must be among them. `rtsc certify` on a bundled QE `a2F.dos5` must give REJECTED_PROVEN |
| lean | `lake build`. The `.lean` sources contain no `sorry`, `admit`, `axiom`, `native_decide`, `implemented_by` or `extern`. `Check.lean`, `CheckPipeline.lean` and `CheckMultiBand.lean` print axioms for every listed theorem. `lean_checks/AllAxioms.lean` collects the axioms of **every** declaration in the three libraries. Only `propext`, `Classical.choice` and `Quot.sound` may appear |
| cert (opt.) | `certificates/wE1/run_all.sh` must print `FINAL VERDICT: CERTIFIED`. It runs on a copy, so the shipped logs are untouched |
| cert-B (opt.) | wB `cert_ends.py`, `cert_upper.py` (a = 4.3314814 on w ∈ [1.3, 40]) and `cert_lower.py` (sharpness at a = 4.33148149); wMB `cert_neg.py` |

The report is in `rtsc_verify/verify_report.txt` and the logs in `rtsc_verify/logs/`. The exit
status is non-zero if any step fails.

**The Lean step and Mathlib.** The bundle does not ship `lean/.lake` (1.3 GB). The first
`lake build` must obtain Mathlib (section 7). `--lean-seed DIR` copies `DIR/.lake` from an
existing build with the same `lake-manifest.json`. `--lean-dir DIR` builds and checks in `DIR`
after confirming that its sources are byte-identical to `lean/`.

---

## 4. Quickstart on a workstation

```sh
# Python >= 3.11. Option A: pip
python3 -m venv ~/venvs/rtsc && . ~/venvs/rtsc/bin/activate
pip install './rtsc[dft,mp]'          # core (numpy, scipy, python-flint) + pymatgen, ASE, phonopy, spglib, mp-api
pip install './rtsc[test]' && (cd rtsc && pytest -q)

# Option B: conda (includes Quantum ESPRESSO 7.x with EPW)
conda env create -f env/environment.yml && conda activate rtsc
pip install --no-deps ./rtsc

rtsc demo --workdir rtsc_demo         # SYNTHETIC spectra end to end, ~10 s; read rtsc_demo/summary.md
rtsc recheck rtsc_demo/SYNTHETIC-strong-coupling/verdict.json

# certify one alpha^2F file (two columns: omega, a2F; give the omega unit)
rtsc certify my_a2F.dat --units meV --target 300 --mu-hi 0.15 --envelope 0.10 --out verdict.json
# a QE matdyn a2F.dos file (units read from the file header), e.g. the bundled smoke-test output:
rtsc certify rtsc/tests/data_dfpt/qe67_al_smoke/a2F.dos5 --format qe_a2f_dos --mode decide
# EPW output: --format epw
```

`rtsc run -c config.toml` runs the full pipeline: sources, dedupe, physical filters, ordering,
DFPT, certify, flags and reports. Example configs are in `rtsc/examples/`. All configuration
keys are documented in `rtsc/README.md` §5. On one machine, set `[run] workers = W`,
`[dfpt.launcher] kind = "mpirun"` and `nprocs = P` with W·P ≤ the number of cores.

---

## 5. Running on a cluster

### 5.1 Software

Use either the conda environment (`env/environment.yml`: QE 7.x incl. `epw.x`, OpenMPI 5,
generic OpenBLAS/FFTW build) or the site's QE module plus a Python venv with
`pip install './rtsc[dft,mp]'`. A site QE build tuned to the hardware (MKL/AOCL, vendor MPI) is
typically 1.5–3× faster for DFPT and is preferable for production. Point
`[dfpt.launcher] bin_dir` at it, or put it on `PATH`.

**Ubuntu 24.04 `apt install quantum-espresso` (QE 6.7) aborts at start-up** with
`buffer overflow detected`. The package is built with `_FORTIFY_SOURCE`, and the md5 helper
calls `snprintf(buf+2i, sizeof(char*)=8, ...)` on a 33-byte buffer. `rtsc` ships a workaround,
`rtsc/rtsc/dfpt/contrib/snprintf_chk_shim.c`, which you compile and preload for the QE
processes only:

```sh
gcc -O2 -shared -fPIC -o ~/lib/libsnprintf_chk_shim.so rtsc/rtsc/dfpt/contrib/snprintf_chk_shim.c
```
```toml
[dfpt.launcher]
env = { LD_PRELOAD = "/home/me/lib/libsnprintf_chk_shim.so" }
mpi_args = ["-x", "LD_PRELOAD"]        # OpenMPI: export to all ranks
```

The same source line is unchanged up to QE 7.6. A from-source build compiled with
`-D_FORTIFY_SOURCE` (the gcc default on Ubuntu 24.04) can therefore abort in the same way. Build
with `CFLAGS="-O2 -U_FORTIFY_SOURCE"`, or use the shim. The conda-forge `qe` 7.6 binaries do not
reference `__snprintf_chk`, so the problem does not arise with them. This was checked with
`nm -D` but the binaries were not run, see section 9. Do not export `LD_PRELOAD` globally,
because it weakens FORTIFY for every program.

### 5.2 MPI launcher for QE / EPW

`rtsc` starts one MPI launch per QE step (`relax → scf_fit → scf → ph → q2r → matdyn`) with the
`[dfpt.launcher]` settings:

| key | meaning |
|---|---|
| `kind` | `"srun"` inside a Slurm allocation, `"mpirun"` on a single node, `"plain"` serial |
| `nprocs`, `npool` | MPI ranks and k-point pools (`-nk`) for `pw.x`/`ph.x`. `nprocs` must be divisible by `npool`. Typical: 32–128 ranks, `npool` = 4–16 |
| `mpi_args` | e.g. `["--exact", "--cpu-bind=cores"]` for srun (needed for concurrent steps when `workers > 1`, Slurm ≥ 21.08); `["--bind-to","core"]` for OpenMPI |
| `serial_codes` | `q2r.x`, `matdyn.x` run serially (default) |
| `omp_threads` | default 1 (pure MPI) |
| `bin_dir`, `env` | QE location; extra environment for QE only (e.g. the shim) |

**EPW.** The automated workflow uses `ph.x` with `electron_phonon='interpolated'`. EPW is
supported through input writers (`rtsc.dfpt.epw_inputs`) and the parser (`--format epw`, or
`spectra_files = [{format = "epw", ...}]` in a source manifest). Run EPW yourself and attach the
`*.a2f` output. `epw.x` requires `-nk` equal to the number of ranks.

A full template is in `env/slurm/config_cluster.toml`. It was checked here: it loads, the plan
stage runs, and the launcher command is `srun -n 32 --exact --cpu-bind=cores pw.x -nk 8 -in ...`.

### 5.3 Slurm job arrays

```sh
cp env/slurm/config_cluster.toml mycfg.toml     # edit workdir, sources, pseudos, launcher
# optional, required once if [ml] is enabled:
rtsc run -c mycfg.toml --plan-only
env/slurm/submit.sh mycfg.toml 32 8             # 32 shards, <= 8 running; collect runs afterwards
```

`rtsc_array.sbatch` runs `rtsc run -c $RTSC_CONFIG --shard $SLURM_ARRAY_TASK_ID/$RTSC_NSHARDS`.
Shard K computes positions K, K+N, ... of the global priority order. `rtsc_collect.sbatch`
(dependency `afterany`) runs `rtsc collect`. It merges the records and writes `candidates.csv`,
`certificates/` and `summary.md`. A material without a finished record is reported as ERROR,
never dropped. Shards can be re-run safely, because finished materials are skipped and QE
restarts from its checkpoints. Set `--ntasks` = `[run] workers` × `nprocs`. Put `[run] workdir`
on a parallel file system: QE scratch (`*.wfc*`, `_ph0/`, dvscf) can reach 10–100 GB per
material in flight.

### 5.4 Pseudopotentials

`[dfpt.pseudos]` needs `dir`, plus either `sssp_json` (file names and recommended cutoffs per
element) or `files = { H = "...", S = "..." }` (optionally with `cutoffs_Ry`). Recommended: SSSP
(Standard Solid-State Pseudopotentials, PBE efficiency or precision) from Materials Cloud
(https://www.materialscloud.org/discover/sssp). Unpack the tarball, then point `dir` at the
directory and `sssp_json` at the JSON next to it. `rtsc` reads the UPF headers, records each
file's md5 and functional in the provenance, refuses missing files or element mismatches, and
warns about mixed functionals. At megabar
pressures, check that the pseudopotential core radii do not overlap. The SSSP tests were done
near ambient pressure.

### 5.5 Materials Project

```sh
export MP_API_KEY=...        # https://next-gen.materialsproject.org/api ; never put the key in config.toml
pip install './rtsc[mp]'
```

The `materials_project` source reads the key from the variable named in `api_key_env` (default
`MP_API_KEY`). MP structures and energies are for 0 GPa. At a target pressure the DFT stage
vc-relaxes, and the hull filter does not use MP's 0 GPa `energy_above_hull`.

### 5.6 Containers

`env/Dockerfile` (build context = bundle root): `docker build -f env/Dockerfile -t rtsc:0.1.0 .`.
It produces Ubuntu 24.04 with a micromamba env from `env/environment.yml` and `rtsc` installed,
running as an unprivileged user. Mount pseudopotentials at `/pseudo` and work at `/work`. For
multi-node MPI on a cluster, prefer Apptainer with an MPI that matches the host, or native
modules. **The image has not been built or run** (section 9).

---

## 6. Hardware sizing (rough)

Almost all of the cost is DFPT electron–phonon work (`ph.x` over the irreducible q-points with
a dense k grid). Certification, the demo and the Lean checks are negligible next to it. The
numbers below are **order-of-magnitude estimates from general QE experience. They are not
measurements of this pipeline.** The only DFPT runs made while building it were tiny smoke
tests (fcc Al at 15 Ry and q 2×2×2, minutes on 2 cores). Expect factors of 3–10 either way,
depending on cutoffs (hydrides at megabar pressure often need 60–100 Ry), k/q grids, symmetry,
pseudopotentials, the number of empty bands, and how well the hardware and QE build match.

| material class (primitive cell) | typical grids | core-hours per material (DFPT + e-ph, all smearings) | memory |
|---|---|---|---|
| elemental metal, 1–2 atoms | q 6³–8³, k 16³ / fit 32³ | 10 – 100 | 1–2 GB/rank |
| small hydride, 4–8 atoms, high symmetry (e.g. Im-3m H3S) | q 4³–6³, k 24³ / fit 48³ | 500 – 5,000 | 2–4 GB/rank |
| clathrate hydride, 10–20 atoms (e.g. Fm-3m LaH10, 11 atoms) | q 4³–6³ | 5,000 – 50,000 | 2–4 GB/rank |
| 20–50 atoms or low symmetry | q 3³–4³ | 10⁵ – 10⁶ (usually impractical for screening) | 4+ GB/rank |

* Rule of thumb: cost scales as N_atoms³⁻⁴ × (irreducible q-points). Lower symmetry raises both
  factors. The vc-relax adds 5–20 %. Extra smearings are post-processing and almost free.
* Disk: 10–100 GB scratch per material in flight. Records and reports are small.
* Screening budget example: 200 small hydrides × ~2,000 core-hours ≈ 4×10⁵ core-hours
  (≈ 3,000 node-hours on 128-core nodes). Use `[budget]` and the priority ordering to spend it
  on the most promising candidates. Everything beyond the budget is DEFERRED and listed.
* Certify: seconds to minutes per spectrum on one core with `mode = "decide"`. The default
  `time_budget_s = 900` per spectrum caps the worst case (10 smearings → ≤ 2.5 core-hours).
* Lean from scratch, without the Mathlib cache: about 1 hour on 2 cores for the Mathlib subset
  used. With `lake exe cache get`: minutes plus a download of a few GB. The wE1 certificate takes
  about 6 core-minutes.

---

## 7. Rebuilding and checking the Lean proofs

Lean `leanprover/lean4:v4.35.0` (`lean/lean-toolchain`), Mathlib tag `v4.35.0`, commit
`56a00585...` (`lean/lake-manifest.json`).

```sh
# 1. elan (Lean version manager); it picks the toolchain from lean/lean-toolchain automatically
curl -sSf https://raw.githubusercontent.com/leanprover/elan/master/elan-init.sh | sh -s -- -y --default-toolchain none
export PATH="$HOME/.elan/bin:$PATH"

# 2. dependencies + Mathlib binaries (needs GitHub and the Mathlib cache server)
cd lean
lake exe cache get          # downloads prebuilt Mathlib .olean files for this exact commit
#   If the cache is unreachable (it was blocked where this bundle was prepared), skip it:
#   `lake build` then compiles the needed Mathlib modules from source (~1 h on 2 cores).

# 3. build the three libraries (TcR1, PipelineLogic, MultiBand)
lake build

# 4. axiom checks
lake env lean Check.lean            # TcR1: 34 theorems
lake env lean CheckPipeline.lean    # PipelineLogic: 56 theorems
lake env lean CheckMultiBand.lean   # MultiBand: 67 theorems
lake env lean ../lean_checks/AllAxioms.lean   # every declaration (464 at assembly time)
grep -rnwE 'sorry|admit' --include='*.lean' . --exclude-dir=.lake   # must print nothing
```

Every line must read `depends on axioms: [propext, Classical.choice, Quot.sound]`, a subset of
it, or `does not depend on any axioms`. `sorryAx` or any other axiom means failure.
`verify_all.sh` automates all of this. Reference outputs from the development build are in
`lean/axioms_output.txt`, `lean/pipeline_axioms_output.txt` and
`lean/multiband_axioms_output.txt`. The statement-by-statement mapping from Lean to the paper
and to the Python code is in `lean/README.md` and `rtsc/docs/CORRESPONDENCE_LEAN.md`.

---

## 8. Re-running the certificates for the constant

Needs Python ≥ 3.11 with `python-flint` (Arb) and `numpy`. NumPy only proposes the Cholesky
factors and test vectors, which are then re-checked in Arb.

```sh
cd certificates/wE1
./run_all.sh                          # a0 = 4.33148145, eta = 1e-9, sharpness at a1 = 4.3314815
#   ~6 min on one core; prints "FINAL VERDICT: CERTIFIED ..." and writes log_4.33148145.txt,
#   cells_4.33148145.{log,json}.  Coarser and faster: ./run_all.sh 4.33148 1e-7 4.331482 150 600
cd ../wB                              # independent route B (PROOF.md there, "Lemma 8'/9")
python3 cert_ends.py 43314814 10000000 13/10 40
python3 cert_upper.py 43314814 10000000 13/10 40 80 1000      # ~2 min
python3 cert_lower.py 433148149 100000000 75748661 10000000 300  # sharpness used in PROOF.md Prop. 7
cd ../wMB && python3 cert_neg.py      # multiband Proposition N scalar inequalities -> PASS
```

`run_all.sh` writes its logs into its own directory and overwrites the shipped
`log_4.33148145.txt`. `verify_all.sh --cert` runs it on a copy instead. See
`certificates/wE1/CERT.md` for the proof structure: ends, 126 cells, the odd sector, and
sharpness.

---

## 9. What was tested while preparing this bundle

Tests were run on Ubuntu 24.04, 2 shared cores, Python 3.13.17, python-flint 0.9.0, numpy 2.5.3,
scipy 1.18.1 and Lean 4.35.0. pymatgen, ASE, spglib and mp-api were not installed.

| item | status |
|---|---|
| `assemble.sh` → `bundle/`, `rtsc_bundle.tar.gz` | run |
| `verify_all.sh` (python, pytest, demo, recheck, lean, `--cert-all`) on the assembled bundle | run; results in the hand-over report |
| `pip install --no-deps ./rtsc` and the `rtsc` entry point | tested (venv inside verify_all). The offline `--no-build-isolation` build failed here because this machine's only setuptools is Debian's patched copy for its system Python 3.12 (`install_layout` error) and was picked up by Python 3.13. The normal isolated build (pip fetches setuptools) worked |
| `env/environment.yml` | **solved** with `micromamba create --dry-run` (conda-forge, linux-64): python 3.12, qe 7.6 `mpi_openmpi`, openmpi 5.0, python-flint 0.8, pymatgen, mp-api, phonopy. **Not installed or run** |
| conda-forge `qe` 7.6 package contents | `epw.x`, `pw.x`, `ph.x`, `q2r.x` and `matdyn.x` are present. Their dynamic symbols contain no `__snprintf_chk` (Ubuntu's QE 6.7 binaries do). **Not executed** |
| `env/Dockerfile` | linted with hadolint 2.15 (with shellcheck), clean. **Not built** |
| Slurm scripts | shellcheck clean. `config_cluster.toml` loads and plans in rtsc, and the launcher command was inspected. **Not run under Slurm** (no Slurm here) |
| real QE runs | only those already in the components (QE 6.7 + shim, Al and Pb smoke tests). Reproduce with `RTSC_RUN_QE=1 RTSC_QE_SHIM=1 pytest rtsc/tests/test_dfpt_qe_real.py` |
| Lean from a clean checkout with `lake exe cache get` | **not tested** here (cache host blocked). The bundle's sources were built against the development `.lake` |
| macOS / ARM | not tested |

---

## 10. Known gaps

The list below was written when the bundle was assembled. `KNOWN_GAPS.txt` is regenerated by
`assemble.sh` for the automatically detectable ones.

* `rtsc/docs/CORRESPONDENCE.md` is referenced by `rtsc/README.md` ("to be written") and
  `rtsc/ARCH.md`, but it does not exist. Only `docs/CORRESPONDENCE_LEAN.md` is present.
* The multiband rejection rule R-MULTI is marked "pending audit" in `rtsc/ARCH.md` §4c.
* `PROOF.md` §3 cites independent certificates from routes A, E2 and H and an analytic backstop
  C. Their directories are not part of this bundle; only E1, B and MB are included.
* The pipeline's DFPT stage has no `[dfpt.shim]` switch. With Ubuntu's QE 6.7 package, set the
  shim through `[dfpt.launcher] env` as shown in section 5.1. The helper
  `rtsc.dfpt.workflow_from_config` knows `[dfpt.shim]`, but `rtsc run` does not use it.
* `rtsc run --plan-only` prints `inputs: 0` and all-zero status counts, although `inputs.json`
  lists the planned inputs. This is cosmetic.
