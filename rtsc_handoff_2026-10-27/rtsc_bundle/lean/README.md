# TcR1 — Lean 4 formalization of Route R1 (L1–L6, MAIN)

Source: `TcR1/Basic.lean` (one file, ~530 lines, no `sorry`, no `axiom`).
The numerical hypothesis `H(a0, η)` is a hypothesis (`TcR1.H`) of the main theorems.

## Versions
* Lean `leanprover/lean4:v4.35.0` (from `lean-toolchain`; installed from the GitHub release
  `lean-4.35.0-linux.tar.zst` into `~/.elan/toolchains/leanprover--lean4---v4.35.0`, because
  elan's default download host is blocked).
* Mathlib tag `v4.35.0` (commit `56a005855ea17e0b0958f73ec801ae140b0a1679`), whose
  `lean-toolchain` is `leanprover/lean4:v4.35.0`.
* Mathlib modules built **from source** (binary cache is blocked), only:
  `Mathlib.Basic.Real.Basic`, `Mathlib.Algebra.Order.BigOperators.Group.Finset`,
  `Mathlib.Algebra.BigOperators.Ring.Finset`, `Mathlib.Tactic.{Linarith,Positivity,FieldSimp,Ring}`
  and their dependencies (~1040 jobs, ~50 min at `nice -n 19` on the shared 2-core machine).

## Build / check
```sh
export PATH=$HOME/.elan/bin:$PATH
cd /home/claude/tc/lean
MATHLIB_NO_CACHE_ON_UPDATE=1 lake update          # only first time (clones deps)
nice -n 19 lake build Mathlib.Basic.Real.Basic Mathlib.Algebra.Order.BigOperators.Group.Finset \
  Mathlib.Algebra.BigOperators.Ring.Finset Mathlib.Tactic.Linarith Mathlib.Tactic.Positivity \
  Mathlib.Tactic.FieldSimp Mathlib.Tactic.Ring    # from source
nice lake build                                   # builds TcR1
lake env lean Check.lean                          # #print axioms for every theorem
```
Logs: `build_mathlib.log`, `build_mathlib2.log`; axioms output: `axioms_output.txt`.

## `#print axioms` output (all 34 theorems; excerpt)
```
'TcR1.main' depends on axioms: [propext, Classical.choice, Quot.sound]
'TcR1.main_even' depends on axioms: [propext, Classical.choice, Quot.sound]
'TcR1.main_odd' depends on axioms: [propext, Classical.choice, Quot.sound]
'TcR1.unstable_imp_lt' depends on axioms: [propext, Classical.choice, Quot.sound]
'TcR1.stable_of_ge' depends on axioms: [propext, Classical.choice, Quot.sound]
```
Every theorem in the file gives exactly `[propext, Classical.choice, Quot.sound]` (see
`axioms_output.txt`).

## Modelling conventions
* Finitely supported `x`: a pair `N : ℕ`, `x : ℕ → ℝ`; all sums over `Finset.range N`.
  `nsq N x = ∑_{n<N} x n ^ 2` is `|x|²`.
* `κ : ℝ` stands for `2π` in `(2π k T)²`. Every theorem holds for every real `κ`, in
  particular `κ = 2 * Real.pi`. (Importing `Real.pi` would pull in ~1600 analysis modules.)
* Cutoff `Nc : WithTop ℕ` (`⊤` = no cutoff); μ*-term `2 μ (∑_{n<N, n<Nc} x n)²`.
* Finite atomic spectrum `Spectrum ι`: Finset `S`, weights `p j ≥ 0`, frequencies `Ω j > 0`;
  `A = ∑_{j∈S} p j`.

## Definitions (namespace `TcR1`)
| Math | Lean |
|---|---|
| `|n-m|` | `adist n m` |
| `C(x) = -∑(2n+1)x_n²` | `C N x` |
| `∑_{k=1}^n L(k)` | `Ssum L n` |
| even / odd kernel with λ(0) cancelled | `Kp L n m`, `Km L n m` |
| λ-part of Q⁺ / Q⁻ for coupling sequence `L` | `Fp L N x`, `Fm L N x` |
| `g_T(k,w) = 2w/(w²+(κkT)²)` | `g κ T k w` |
| `λ_T(k) = ∑_j p_j g_T(k,Ω_j)` | `lam κ sp T k` |
| `2μ*(∑_{n<N_c} x_n)²` | `muTerm μ Nc N x` |
| `Q⁺_T(x)` | `Qp κ sp μ Nc T N x` |
| `Q⁻_T(x)` | `Qm κ sp T N x` |
| `2v/(v²+(κk)²)` | `ge κ v k` |
| `Qe(a,v;x)` (= Q⁺ at T=1, μ*=0, single atom a at v) | `Qe κ a v N x` |
| `H(a0,η)` | `H κ a0 η` := `∀ v>0, ∀ N x, Qe κ a0 v N x ≤ -η * nsq N x` |
| unstable (ρ ≥ 0) | `Unstable Q` := `∀ ε>0, ∃ N x, (∃ n<N, x n ≠ 0) ∧ -ε*nsq N x < Q N x` |

## Lemmas
| Route item | Lean name | Statement (informal) |
|---|---|---|
| form = SPEC form (any finite value `L0` for λ(0)) | `Fp_spec_form`, `Fm_spec_form`, `Qp_spec_form`, `Qm_spec_form` | the cancelled forms equal the SPEC forms with `λ(0) := L0` |
| form = ROUTE display (off-diag + diag) | `Fp_route_form`, `Fm_route_form` | |
| linearity in the spectrum | `Kp_sum`, `Kp_smul`, `Fp_sum`, `Fp_smul` | |
| **L1** | `L1` | `Qp = C - muTerm + ∑_j p_j Fp(g_T(·,Ω_j))` |
| **L2** | `L2_abstract`, `L2` | pointwise `C + A·F(w) ≤ -η|x|²` ∀w>0 ⇒ `C + ∑ p_j F(Ω_j) ≤ -η|x|²` (A>0) |
| **L3** | `L3` (+ `muTerm_zero`) | `-muTerm μ Nc N x ≤ 0` for `μ ≥ 0` |
| **L4** | `L4_pointwise`, `L4` | `A g_T(k,w) = (A/T) ge(w/T,k)`; `C + A·Fp(g_T(·,w)) = Qe(A/T, w/T)` (T>0, w>0) |
| **L5** | `C_le`, `Qe_linear`, `L5` | `C ≤ -|x|²`; `Qe(a)=C+a·G`; `Qe(a0)≤-η|x|²`, `0≤a≤a0`, `a0>0` ⇒ `Qe(a)≤-min(η,1)|x|²` |
| **L6** | `K_abs`, `Fm_le_Fp_abs`, `C_abs`, `nsq_abs`, `lam_nonneg`, `L6` | `Qm(x) ≤ Qp(μ=0)(|x|)`, `| |x| | = |x|` |
| **MAIN** | `main_even`, `main_odd`, `main` | `H κ a0 η`, `A>0`, `μ≥0`, any `Nc`, `T>0`, `A/T ≤ a0` ⇒ ∀N x, `Qp ≤ -min(η,1)|x|²` and `Qm ≤ -min(η,1)|x|²` |
| corollary | `not_unstable_of_bound`, `unstable_imp_lt`, `stable_of_ge` | with `a0>0`, `η>0`: (Unstable Qp ∨ Unstable Qm) at `T>0` ⇒ `T < A/a0`; `T ≥ A/a0` ⇒ stable in both sectors |

## Not done
Stretch goal (general-measure L2 via `MeasureTheory` integrals): not attempted. Its
import closure (`Mathlib.MeasureTheory.Integral.Bochner.Basic`) is ~2300 Mathlib modules,
~1600 more than built here, i.e. several more hours of source build on the shared 2-core
machine, far over the 20-minute job limit in SPEC.md.

---

# PipelineLogic: validity of the RTSC pipeline decision logic

Files: `PipelineLogic.lean` (root) and `PipelineLogic/{Common,Schur,Monotone,Cover,Rayleigh,Assembly,Classify}.lean`,
about 960 lines in total. The library is registered in `lakefile.toml` as a second `lean_lib`
and imports `TcR1`. `TcR1/Basic.lean` is unchanged, and no new Mathlib modules are needed.
There is no `sorry` and no `axiom`. It covers ARCH §3–§4 of `/home/claude/rtsc/ARCH.md`.

```sh
nice lake build                        # builds TcR1 and PipelineLogic (9 new jobs)
lake env lean CheckPipeline.lean       # #print axioms for all 56 theorems
```

The output is in `pipeline_axioms_output.txt`. Every theorem uses a subset of
`[propext, Classical.choice, Quot.sound]`. The pure-Bool facts about `classify` use fewer axioms
or none.

| Item | Lean names |
|---|---|
| P1 weighted Schur / Collatz–Wielandt row test | `schur_row_test` (finite N), `schur_row_test_inf` (truncation form), `rowSum_mono`, `rows_of_eventually` |
| P2 interval envelope | `Mq_le_Mhat`, `rowSum_le_envelope`, `envelope_form_bound`, `stable_of_even_bound`, `R_MU` |
| P3 monotonicity and envelope | `g_antitone`, `lam_antitone`, `lam_mem_cell`, `lam_env`, `lam_box`, `A_env` |
| P4 cover | `Tiles` (def), `cover`, `cover_points` |
| P5 Rayleigh, box minimum, μ* | `rayleigh_unstable`, `Qp_anti_mu`, `boxMinQ` (def), `boxMin_le`, `boxMin_attained`, `sel_mul_eq_min`, `rayleigh_cert` |
| P6 assembly | `stable_of_A_le`, `R_UNIV`, `CellOK` (def), `cell_stable`, `cells_stable`, `R_CELLS` |
| P7 classify | `Status`, `AltEvidence`, `classifyAlt`, `classify`, `EvidenceSound`; `classify_cases`, `status_exclusive`, `classify_eq_rejected_iff`, `classify_eq_certified_iff`, `classify_rejected_sound(')`, `classify_certified_sound`, `evidence_consistent` |

Plain-math statements, the ARCH rule each theorem backs, the obligations on the Python side,
deviations and flaws found in ARCH are listed in `/home/claude/rtsc/docs/CORRESPONDENCE_LEAN.md`.

---

# MultiBand: multiband bound `k_B Tc ≤ κ/a0` (wMB/PROOF.md)

Files: `MultiBand.lean` (root) and `MultiBand/{Defs,Algebra,Identities,Odd,Main,Pipeline,Faithful}.lean`
(about 1170 lines). The library is registered in `lakefile.toml` as a third `lean_lib` and imports
`TcR1` only. `TcR1` and `PipelineLogic` are unchanged, and no new Mathlib modules are needed:
a clean rebuild compiles exactly the 8 `MultiBand` modules (`multiband_build.log`).
There is no `sorry` and no `axiom`.

```sh
nice lake build                        # builds TcR1, PipelineLogic, MultiBand
lake env lean CheckMultiBand.lean      # #print axioms for all 67 MultiBand theorems
```
Output: `multiband_axioms_output.txt`. 66 theorems depend on `[propext, Classical.choice, Quot.sound]`,
and `adist_comm'` depends on `[propext, Quot.sound]`.

## Model (namespace `MultiBand`)
* Bands `Fin B`; gap `φ : Fin B → ℕ → ℝ` with truncation `N` (finitely supported); DOS `dos i = N_i`.
  `|x|² = nrm2 dos N φ = Σ_i N_i Σ_{n<N} φ_{n,i}²` (`nrm2_eq_x`: equals `Σ_i |√N_i φ_i|²`).
* `MSpec B α`: atoms `s ∈ S`, frequencies `Ω s > 0`, weights `W s i j ≥ 0` with `W s i j = W s j i`
  (`W_ij = N_i α²F_ij`). `Amat sp dos i j = (Σ_s W s i j)/N_i = 𝒜_ij`.
* `tp` is `2π` (any real). `κ` is the Perron root; it enters only through a positive vector `v`
  with `𝒜v = κv` (Lemma 4.1) or `𝒜v ≤ κv` (Theorem 1, Collatz–Wielandt form).
* `Cw` = wMB's `C(y) = Σ(2n+1)y_n²` (`= -TcR1.C`); `Bp`, `Bm`, `Dn`, `Dq`, `Gp`, `Gm` = `B^±_w`, `D_w`, `G^±_w`;
  `Qplus`, `Qminus` = `Q^±` of wMB §1 (Z with row sums, pairing with full matrix, Coulomb
  `2Σ_ij N_iμ*_ij S_iS_j`, `S_i = Σ_{n<min(N,N_c)} φ_{n,i}`); `Pterm` = `P`; `theta/psi/chi/mw` = `θ, ψ, χ, m`.
* Hypotheses: `TcR1.H tp a0 η` (Prop. 6), `PropN tp a0 := ∀ w>0, ∀ N y, Qe tp (-a0) w N y ≤ 0`,
  `CoulPSD dos μ := ∀ S, 0 ≤ Σ_ij N_iμ_ij S_iS_j`.
* `MUnstable`, `MStable` = multiband analogues of `TcR1.Unstable` / `PipelineLogic.Stable`.

## Theorems
| wMB item | Lean | Statement |
|---|---|---|
| §1 faithfulness | `Qplus_lambda_form`, `nrm2_eq_x` | `Q⁺ = -Σ N_iC(φ_i) + Σ_i N_i Σ_j[Σ_nm(λ_ij(|n-m|)+λ_ij(n+m+1))φ_{n,i}φ_{m,j} - Σ_n(λ_ij(0)+2Σ_{k≤n}λ_ij(k))φ_{n,i}²] - Coulomb`, `λ_ij(k)=Σ_s (W_s,ij/N_i) g_T(k,Ω_s)` |
| single-band link | `Gp_self`, `Qe_eq` | `G⁺_w(y,y) = TcR1.Fp`, `Q^E_{a,w} = TcR1.Qe tp a w` (any real `a`) |
| Lemma 1.1 | `lemma11_even`, `lemma11_odd`, `Pterm_nonneg` | `Q^± = -ΣN_iC + (1/T)Σ_s Σ_ij W G^±(φ_i,φ_j) - P (- Coulomb)`; `P ≥ 0` for `T>0` |
| Lemma 4.1 | `lemma41` (+ `gauge_G`, `gauge_C_eig`) | `𝒜v=κv`, `κ,T,v_i,N_i ≠ 0`: `Q⁺ = (1/κ)Σ_s Σ_ij m_ij[Q^E_{κ/T,w}(ψ_ij)+Q^E_{-κ/T,w}(χ_ij)] - P - Coulomb` |
| §4.3 step | `QE_pm` | `Q^E_{±a0}≤0`, `0≤a≤a0` ⇒ `Q^E_{±a} ≤ -(1-a/a0)C` |
| Theorem 1 | `thm1_even` | H, `η≥0`, `a0>0`, PropN, `N_i>0`, `v>0`, `𝒜v≤κv`, CoulPSD, `T>0`, `κ≤a0T` ⇒ `Q⁺ ≤ -(1-κ/(a0T))|x|²` |
| Theorem 1 (literal route) | `thm1_even_via_identity` | same conclusion via `lemma41`, for `𝒜v=κv`, `κ>0` |
| Theorem 3 | `thm3_odd` | `T>0`, `N_i≥0` ⇒ `Q⁻ ≤ -Σ N_iC(φ_i) ≤ -|x|²` (no H, no Prop N) |
| Tc ≤ κ/a0 | `mstable_of_gt`, `unstable_imp_le` | `T>0`, `κ<a0T` ⇒ `MStable`; unstable `T>0` ⇒ `T ≤ κ/a0` |
| §7 rule | `R_UNIV_MB`, `R_UNIV_MB_eig`, `R_UNIV_MB_rowsum`, `sub_env`, `sub_of_K` | `κ_hi < a0T*` ⇒ `κ_hi/a0 < T*` and stable for all `T ≥ T*`, every PSD Coulomb matrix and cutoff; envelope `𝒜' ≤ (1+e)𝒜`; symmetric-gauge certificate `Ku ≤ κ_hi u`; max-row-sum fallback |

Plain statements, deviations and pipeline obligations: `/home/claude/rtsc/docs/CORRESPONDENCE_LEAN.md`
(section "Library `MultiBand`").

---

# PipelineLogic.Multiband: mirror of `combine_multiband` (ARCH §4c)

File: `PipelineLogic/Multiband.lean`. It imports `PipelineLogic.Classify` and `MultiBand`, and
is included in the `PipelineLogic` root. There is no `sorry` and no `axiom`, and no new
Mathlib modules are needed. `CheckPipeline.lean` now prints the axioms of all 82
`PipelineLogic` theorems (56 + 26) to `pipeline_axioms_output.txt`. Each theorem depends on a
subset of `[propext, Classical.choice, Quot.sound]`.

| Item | Lean names |
|---|---|
| Mirror | `RejectionModel`, `combineMultiband` (verbatim `combine_multiband`, including `n_spectra` as `nSpectra : Option ℕ`), `classifyFull evs f m = combineMultiband (classify evs) f m (some evs.length)` |
| Combinatorics | `combine_mb_rejected_iff`, `combine_mb_certified_iff`, `combine_mb_cases`, `classifyFull_cases`, `classifyFull_exclusive`, `combine_mb_empty_not_rejected`, `classifyFull_empty_not_rejected`, `combine_mb_rejected_needs_mb`, `combine_certified_imp`, `combine_mb_zero_spectra`, `classifyFull_zero_spectra(_not_rejected)` |
| R-MULTI semantics | `MBAlt`, `MBEnv`, `CWCore`, `RMultiCert`; routes `route_cw`, `route_rowsum`, `route_rowsum_declared`, `route_eig` |
| Soundness | `classifyFull_mb_rejected_sound(')` (multiband model; every alternative, envelope member, PSD Coulomb matrix, cutoff, `T ≥ T*`); `classifyFull_iso_rejected_sound`; `classifyFull_certified_sound` (single-band model only, i.e. M-ISO-SINGLEBAND) |
| Conditional lift | `Qplus_const`, `nrm2_const`, `munstable_of_iso`, `certified_transfers_to_multiband` |

Plain statements: `/home/claude/rtsc/docs/CORRESPONDENCE_LEAN.md`, section
"`PipelineLogic.Multiband`".

---

# DiscoveryLogic: the discovery-loop ledger as a state machine

Files: `DiscoveryLogic.lean` (root) and `DiscoveryLogic/{Ledger,Loop,Continue,Soundness}.lean`
(about 1400 lines). The library is registered in `lakefile.toml` as a fourth `lean_lib`. It
imports `PipelineLogic.Classify`, `PipelineLogic.Multiband` and `Mathlib.Data.List.Nodup`, all
of them already built, so a rebuild compiles only the 5 `DiscoveryLogic` modules. There is no
`sorry` and no `axiom`. The library models the post-audit Python code of `rtsc/discovery/`
(`ledger.py` schema v2, with the SQL guard trigger, the `current_status` view, `counts_for` and
`_reconcile`).

```sh
nice lake build                        # builds TcR1, PipelineLogic, MultiBand, DiscoveryLogic
lake env lean CheckDiscovery.lean      # #print axioms for all 69 DiscoveryLogic theorems
```
The output is in `discovery_axioms_output.txt`. 68 theorems depend on
`[propext, Classical.choice, Quot.sound]`; `le_foldr_max` depends on `[propext, Quot.sound]`.

## Model (namespace `DiscoveryLogic`)
* `Event`: one `verdicts` row (material, goal id, round, `Stage`, `Origin`, `Status` from
  `PipelineLogic`, T*, material pressure, `in_scope`, DFT pressure, `dft_hash`, `certify_hash`).
  `Origin` abstracts the origin string: `classify` / `stability*` / other.
* `appendAllowed`: the guard of `append_verdict`. `sqlGuard`: the SQL trigger
  `verdicts_guard`. `commit L es`: one transaction (all of `es`, or nothing).
* `currentEvent` / `currentStatus`: the latest event per (material, goal), i.e. the
  `current_status` view.
* `State`: ledger tables (`Data`), the session `Config` (goals, `certify_hash`, `dft_hash`), the
  round and its phase. `DOp` / `Op` / `step`: one operation per Python write path (generate one
  structure, pre-screen one material, add a job, select-DEFER, record one pipeline outcome,
  reconcile writes, start a round, set a phase, start a session with a possibly new config).
  `Reachable`: every state that any sequence of operations reaches from an empty ledger. This
  includes every crash point, because the operations are at least as fine-grained as the
  Python transactions.
* `Env`: the outside world (generators, dedupe, ML pre-screen, policy, pipeline).
  `runGen/runPre/runSel/runExec/runDone`, `finish` and `resume` model the round phases of
  `_run_round` and `run`/`continue`. `EnvAgree k e1 e2` holds the explicit determinism
  hypotheses.
* `countsFor` / `Satisfies`: `goals.counts_for` and `certified_ids`. `PipelineSound` is the
  hypothesis that links the ledger rows to `classifyFull` and `EvidenceSound`.

## Main theorems
| Item | Lean names |
|---|---|
| guard = SQL trigger | `appendAllowed_eq_sqlGuard` |
| S2-ext / S1 / S3 per row | `appendAllowed_heuristic`, `appendAllowed_proven`, `appendAllowed_excluded`, `appendAllowed_certified_inScope` |
| S2-ext / S1 / S3 on reachable ledgers | `S2_ext`, `no_proven_from_heuristic`, `proven_only_from_classify`, `excluded_only_from_stability`, `current_certified_from_classify`, `reachable_proven_pressure` |
| append-only | `step_prefix`, `step_ext`, `dataStep_ext`, `resume_prefix`, `snap_prefix` |
| accounting | `accounting` (with `reachable_inv`, `Inv`) |
| goal soundness | `satisfies_sound` (`counts_for` + `PipelineSound` + `classifyFull_certified_sound`) |
| continuation | `continue_after_crash`, `continue_currentStatus`, `continue_with_session`, `newSession_same`, `continue_prescreen_partial` |

Plain statements, the Python file:line mapping, the hypotheses and the discrepancies found are
in `/home/claude/rtsc/docs/CORRESPONDENCE_LEAN.md`, section "Library `DiscoveryLogic`".
