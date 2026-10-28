# HANDOFF: rigorous room-temperature-superconductor screening pipeline

Read this file first. It is written for a fresh Claude session, or a person, continuing the work
on a machine with an NVIDIA GPU. Everything referenced here is inside this bundle.

Prepared: 2026-10-27 by the coordinating session (cloud workspace, 2 CPU cores, no GPU, restricted
network). The goal of the handoff is to run the parts that need a GPU and an unrestricted network:
MatterGen, MatterSim at scale, and real Quantum ESPRESSO/EPW.

---

## 0. What this project is, in one paragraph

The project screens candidate materials for phonon-mediated superconductivity at a target Tc
(default 300 K) with proof-carrying verdicts. The proofs rest on two theorems:

- **Single band.** In isotropic Migdal–Eliashberg theory, k_B·Tc ≤ A/4.33148145 ≈ 0.2308679·A,
  where A = ∫α²F(ω)dω. This holds for every spectrum, every μ* ≥ 0 and every Coulomb cutoff, and
  is sharp to a relative 1e-8.
- **Multiband.** The multiband generalization is Tc ≤ κ/4.33148145, where κ is the Perron root of
  the band coupling-weight matrix 𝒜_ij = ∫α²F_ij dω.

A Python program (`rtsc`) computes or ingests α²F spectra and classifies each material. It uses
rigorous ball arithmetic (python-flint/Arb) and decision logic that mirrors Lean 4 definitions;
Lean proves that logic sound. An autonomous discovery loop proposes candidates (MatterGen),
pre-screens them (MatterSim), sends the best to DFT, and certifies them. It keeps a persistent
ledger and can be continued indefinitely.

## 1. Bundle map

| Path | What it is |
|---|---|
| `HANDOFF.md` | this file |
| `README.md` | install, verification, cluster use, path map |
| `verify_all.sh` | one-shot verification: tests, demo, recheck, Lean build and axiom audit, optional certificates |
| `rtsc/` | the Python package; see `rtsc/README.md` for the full config reference |
| `rtsc/ARCH.md` | architecture and **soundness contract** (S1–S6), errata §4b, multiband rule §4c |
| `rtsc/DISCOVERY.md` | discovery-loop spec and ledger semantics |
| `rtsc/docs/CORRESPONDENCE.md` | every Python decision → Lean theorem, plus ranked list of gaps |
| `rtsc/docs/CORRESPONDENCE_LEAN.md` | Lean-side statements, file:line maps (TcR1, PipelineLogic, MultiBand, DiscoveryLogic) |
| `rtsc/rtsc/certify/PROOFS.md` | hand proofs of the engine's analytic lemmas (not in Lean) |
| `lean/` | Lean 4.35 + Mathlib v4.35.0 project: `TcR1`, `PipelineLogic`, `MultiBand`, `DiscoveryLogic` |
| `proof/PROOF.md`, `proof/SPEC.md` | the single-band theorem and model definition |
| `proof/audits/` | referee reports and logs |
| `certificates/wE1` | primary certificate for the constant (`run_all.sh` → FINAL VERDICT: CERTIFIED) |
| `certificates/wB` | independent second certificate, including sharpness |
| `certificates/wMB` | multiband theorem (`PROOF.md`) and Prop. N certificate (`cert_neg.py`) |
| `env/` | `environment.yml` (main env), Dockerfile, Slurm templates |

## 2. Status at handoff (all verified)

**Theory**
- **Single-band theorem.** Proved and certified. Five independent certificates were produced by
  three numerical methods. Three referees found it sound: two audited certificates and one audited
  the global chain.
- **Multiband theorem.** Audited sound. Two independent numerical attacks (thousands of
  configurations, up to 8 bands) found no violation.

**Lean** (no `sorry`; only the axioms `propext`, `Classical.choice`, `Quot.sound`)

| Library | Theorems | Contents |
|---|---|---|
| `TcR1` | 34 | reduction lemmas and main theorem; numerical certificate H(a0, η) as a hypothesis |
| `PipelineLogic` | 82 | Schur row test, envelopes, cells/cover, Rayleigh, classify and multiband combination soundness |
| `MultiBand` | 67 | multiband identity and theorem; Prop. N as a hypothesis |
| `DiscoveryLogic` | 69 | ledger state machine: heuristics never prove, append-only, accounting, goal soundness, crash/continue |

**Python**
- 442 tests pass (8 skipped for optional packages).
- Four adversarial audit suites are included: `tests/audit`, `tests/audit_certify`,
  `tests/audit_discovery`, and the forgery tests.
- The independent float oracle (`tests/oracle`) agrees on 132/132 cases.
- A clean-room install was tested on aarch64 Linux (Mac mini VM). The demo is identical and the
  constant certificate re-certifies in 44 s.

**Real external codes, as run in the cloud workspace**
- Quantum ESPRESSO 6.7 (Ubuntu package) ran an fcc-Al electron–phonon smoke test end to end.
  - It needs an LD_PRELOAD shim for a FORTIFY abort; see `rtsc/rtsc/dfpt/contrib/`.
  - The conda-forge QE 7.x binaries do not have this issue.
- EPW 5.3 produced per-state λ for fcc Pb. Its Eliashberg solver then crashed, which is not needed.
- MatterSim ran on CPU: Al, MgB2, and H3S at 150 GPa.
- **MatterGen has NOT been run.** Network policy blocked it and there was no GPU. Its adapter is
  tested only against a fake CLI.

## 3. Rules a continuing session must not break

1. **Only `rtsc.certify.classify` may yield REJECTED_PROVEN or CERTIFIED_CANDIDATE (S1).** Its
   control flow is a verbatim mirror of Lean `classifyAlt` / `classify` / `combineMultiband` /
   `classifyFull`, and `tests/test_correspondence.py` diff-checks the copied Lean text. If you
   change classify logic, update the Lean mirror and re-prove it; never change only one side.
2. **Heuristics only reorder or defer (S2).** This covers MatterGen, MatterSim, descriptors, the
   search policy and budgets. The SQL triggers in the discovery ledger enforce this; do not
   weaken them.
3. **Every decision-relevant comparison is exact or Arb.** Floats only propose. Every infinite sum
   needs a proved tail bound, which goes in `PROOFS.md` and gets audited.
4. **Uncertainty is conservative (S5).**
   - Alternatives combine as a conjunction.
   - An empty evidence list never proves anything.
   - Inputs that could make A too small (truncated files, dropped weight) block rejection unless
     the user explicitly opts in.
5. **Verdicts are about the model (S6).** Out-of-model risks are flags and never upgrade a verdict.
   Single-band rejections carry M-ISO-SINGLEBAND and MULTIBAND_RISK.
6. **Label synthetic data.** Never present synthetic or mock data as real.
7. **Proof-relevant changes get an independent audit.** Every new proof-relevant change gets an
   adversarial audit and, where feasible, a Lean mirror. Both kinds of review found real bugs
   every time in this project, including after "done" reports.

## 4. Environment setup on the GPU machine

Keep everything in virtual environments; do not touch the system Python.

**4.1 Main environment (Python ≥ 3.12)**
```bash
conda env create -f env/environment.yml        # or: uv venv -p 3.12 .venv && source .venv/bin/activate
pip install -e "rtsc[dft,mp,mattersim,test]"    # python-flint is pinned to 0.9.0 (verified trust base)
./verify_all.sh                                 # add --cert to re-run the constant certificates (minutes)
```

**4.2 MatterGen (separate environment: Python 3.10, Linux, CUDA)**

MatterGen pins `torch==2.2.1+cu118`, so it cannot share the main environment.
```bash
git clone https://github.com/microsoft/mattergen && cd mattergen
uv venv .venv --python 3.10 && source .venv/bin/activate && uv pip install -e .
# checkpoints download from Hugging Face (microsoft/mattergen) on first use, or: git lfs pull
mattergen-generate /tmp/mg_test --pretrained-name=chemical_system_energy_above_hull \
  --batch_size=16 --num_batches=1 \
  --properties_to_condition_on="{'chemical_system':'Mg-B','energy_above_hull':0.05}" \
  --diffusion_guidance_factor=2.0
```
Then point the rtsc adapter at that executable. The config keys `[discovery.generators.<name>]`
with `type="mattergen"` and `executable`/`command_prefix` are documented in `rtsc/README.md`,
"ML adapters".

Notes:
- The checkpoints support conditioning on: `chemical_system`, `space_group`, `dft_mag_density`,
  `dft_band_gap`, `ml_bulk_modulus`, and the joint `chemical_system_energy_above_hull` /
  `dft_mag_density_hhi_score`.
- There is **no pressure conditioning**. Structures are relaxed at the target pressure by
  MatterSim.
- **First real-run task:** check that the adapter parses real `generated_crystals.extxyz` /
  `generated_crystals_cif.zip`. It was only tested against the documented format.

**4.3 MatterSim** is installed in the main environment by the `mattersim` extra. Use the GPU if
available (`device="cuda"`). The `MatterSim-v1.0.0-5M` model is more accurate than the default 1M.

**4.4 Quantum ESPRESSO + EPW**
- Prefer conda-forge `qe` 7.x (it includes `epw.x`).
- Get pseudopotentials from SSSP (efficiency or precision) and map them in
  `[dfpt.pseudos]`.
- Set the launcher with `[dfpt.launcher]` (mpirun/srun, np, -nk).

**4.5 Lean**
```bash
cd lean && lake exe cache get && lake build   # the cache makes this minutes, not ~1 hour
lake env lean CheckPipeline.lean; lake env lean CheckMultiBand.lean; lake env lean CheckDiscovery.lean; lake env lean Check.lean
```
`verify_all.sh` also runs the whole-library axiom audit, `lean_checks/AllAxioms.lean`.

## 5. Next steps, in priority order

1. **Verify the bundle on the new machine.** Run `./verify_all.sh --cert`; everything should
   pass. Re-run `rtsc demo`.
2. **Benchmark against known superconductors before any discovery run.** This "retrodiction"
   step is the most important missing validation. The certify engine has only seen synthetic
   spectra and one coarse Al smoke test.
   - Compute α²F with DFPT for: Al (Tc 1.2 K), Pb (7.2 K), Nb (9.3 K), MgB2 (39 K; the
     multiband case), and H3S at ~150–200 GPa (~200 K).
   - Use converged grids and several smearings as alternatives.
   - Check that the certified brackets of Tc_model are consistent with the literature, given μ*
     in [0.10, 0.15].
   - Check that MgB2 shows the expected single-band underestimate and is correctly flagged
     MULTIBAND_RISK.
   - Record the results in a benchmark report. This is also the end-to-end test of
     `rtsc run` with real QE: so far the workflow was smoke-tested only on its own.
3. **Run MatterGen for real** (§4.2) and run a structure set through the MatterSim prescreener.
   Inspect the outputs by hand.
4. **Change the policy for ML phonon instability.** MatterSim (1M) calls H3S at 150 GPa
   dynamically unstable (minimum −450 cm⁻¹). DFT and anharmonic studies find it stable, and it
   is the known ~200 K superconductor. So "ML-unstable ⇒ DEFERRED" can bury real candidates,
   especially hydrides.
   - Change it to a priority penalty with a guaranteed exploration quota.
   - It must still never exclude a candidate (S2).
   - This is a heuristic change only: no proofs are affected. Keep `tests/test_discovery_soundness.py`
     passing.
5. **Run a first small discovery** with `examples/discovery/config.toml` adapted:
   - goals G-RT (300 K, P ≤ 200 GPa) and G-AMB (e.g. 77 K, ambient);
   - `dft_per_round` 2–4, `max_rounds` 3;
   - `rtsc discover run`, then `continue`, then `status` / `list`.
6. **Multiband data path.** EPW only provides per-state λ_nk, so the multiband κ from EPW is a
   loose row-sum upper bound (about 3× the isotropic A for Pb). Investigate band-resolved α²F_ij
   from EPW internals or other codes for MgB2-like materials.
7. **Optional: fine-tune MatterGen** on the ledger's DFT labels (`rtsc discover export --finetune`)
   once enough real data exists, after hundreds of DFT-labelled structures, not before.

## 6. Known gaps and open items (none is a known soundness bug in the shipped decision paths)

Details are in `rtsc/docs/CORRESPONDENCE.md` (gaps G1–G16) and the audit reports.

**Proof coverage**
- **Lean scope.** Lean covers finite atomic spectra and the decision logic. Not in Lean: general
  measures, λ(0) = ∞, the numerical certificates H(a0, η) and Prop. N (they are Arb certificates,
  audited), and the engine's hand-proved tail lemmas (S1, T0, R1, R2, T1, S2 in `PROOFS.md`,
  audited). The `cform` row check is hand-proved only and is off by default.
- **Discretization.** Verdicts are about the parsed atomic spectrum. Discretizing a continuous
  α²F is outside the proofs and is labelled DISCRETISATION in every verdict.

**Multiband**
- Detailed balance (per frequency) and positive semidefinite N_i μ*_ij are assumed, not checked.
  The cw route does not check the Kolmogorov cycle condition. All of this is recorded as
  assumptions.
- Never supply band data in the "contributions normalized to the isotropic α²F" convention; it
  is refused when tagged.

**Discovery loop**
- A continued run reproduces the same choices only when `max_dft_jobs` is not the binding limit.
- Schema-v1 ledgers are not migrated.
- Two last fixes landed after the `DiscoveryLogic` Lean library was written: reconcile is now
  idempotent (D5), and accounting checks per-round generator counts (D3). Both only remove writes
  or add checks, so the proved invariants still hold. However, the file:line references in
  `CORRESPONDENCE_LEAN.md` § DiscoveryLogic for `_reconcile` and `accounting()` are slightly
  stale. Refresh them next time the Lean model is touched.
- Structure de-duplication uses pymatgen StructureMatcher when installed. Install it: the fallback
  is conservative but slower.

**Practical limits**
- Hardware sizing in `README.md` is a rough order of magnitude. Electron–phonon DFPT dominates
  cost: expect tens to thousands of core-hours per material depending on size and convergence.
- Mac mini drive (`Claude_Code_Workspaces/rtsc_test`): two half-created venv folders (`.venv`,
  `venv312`) were left by a failed attempt on the FUSE mount. They are harmless; delete them when
  convenient.

## 7. Suggested working style (it worked well here)

- **Coordinator plus specialists.** One coordinating session; specialist agents for engine, DFT,
  orchestration and Lean; and independent adversarial auditors who never see the author's
  reasoning.
- **Write the contract first.** Write the spec and soundness contract before code. Formalize the
  decision logic in Lean as it is built: it found real holes (an empty-evidence vacuous truth and
  a mis-specified matrix diagonal).
- **Keep the machine-checked chain short.** Python decisions are verbatim Lean mirrors, plus a
  correspondence test that fails on drift.
