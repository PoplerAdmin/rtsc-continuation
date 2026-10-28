# Lean ↔ ARCH correspondence: library `PipelineLogic`

Lean project: `/home/claude/tc/lean` (Lean 4.35.0, Mathlib v4.35.0). Library `PipelineLogic`
(files `PipelineLogic/{Common,Schur,Monotone,Cover,Rayleigh,Assembly,Classify}.lean`) imports
`TcR1` (`TcR1/Basic.lean`, unchanged). It has no `sorry` and no `axiom`. Every theorem below
depends only on `propext`, `Classical.choice` and `Quot.sound`, or on fewer of them; see
`/home/claude/tc/lean/pipeline_axioms_output.txt`. To re-check:

```sh
export PATH=$HOME/.elan/bin:$PATH; cd /home/claude/tc/lean
nice lake build && lake env lean CheckPipeline.lean
```

## Conventions (same as TcR1)

* Finitely supported `x`: a pair `N : ℕ`, `x : ℕ → ℝ`. All sums run over `n < N`, and
  `|x|² = nsq N x = Σ_{n<N} x_n²`.
* `κ` stands for `2π`. Every theorem holds for every real `κ`, so in particular for `κ = 2π`.
  The theorems that need Proposition 6 take `H(a0, η)` (`TcR1.H κ a0 η`) as a hypothesis, with
  `a0 > 0` and `η > 0`.
* Spectrum: `TcR1.Spectrum ι` is a finite atomic measure with atoms `j ∈ S`, weights
  `p_j ≥ 0` and frequencies `Ω_j > 0`. `A = Σ_j p_j`, and `λ_T(k) = lam κ sp T k = Σ_j p_j g_T(k, Ω_j)`.
* `Unstable Q` (TcR1) means: for every `ε > 0` there is an `x ≠ 0` with `Q(x) > −ε|x|²`.
  This is exactly `ρ ≥ 0`.
* **Stable** (`PipelineLogic.Stable κ sp μ Nc T`) means
  `¬Unstable(Q⁺_{T,μ,Nc}) ∧ ¬Unstable(Q⁻_T)`, i.e. ρ⁺ < 0 and ρ⁻ < 0 (ARCH §4 notation).
* **Envelope** (`InEnv e sp sp'`): `sp'` has the same atoms and frequencies as `sp`, and
  `(1−e)p_j ≤ p'_j ≤ (1+e)p_j`. This is `Spectrum.rel_envelope = e`.
* Cutoff: `Nc : WithTop ℕ` (`⊤` = none). The μ* term is `2μ*(Σ_{n<min(N,Nc)} x_n)²`. In the
  classifier the model's cutoff rule is an arbitrary function `Ncf : ℝ → WithTop ℕ` (T ↦ N_c(T)).
  Python's `N_c(T) = #{n ≥ 0 : (2n+1)πT ≤ ω_c}` is one instance.
* Matrices: `Mq L n m` is the matrix of the μ* = 0 even form, with `C(x) + F⁺_L(x) = Σ_{n,m} Mq_nm x_n x_m`:
  - off-diagonal: `Mq_nm = L(|n−m|) + L(n+m+1)`;
  - diagonal: `Mq_nn = L(2n+1) − 2Σ_{k=1}^{n} L(k) − (2n+1)`.

  `Mhat hi lo` is the interval envelope:
  - off-diagonal: `hi(|n−m|) + hi(n+m+1)`;
  - diagonal: `hi(2n+1) − 2Σ_{k=1}^{n} lo(k) − (2n+1)`.

  `rowSum M v N n = Σ_{m<N} M_nm v_m`.

## Theorems

### P1: weighted Schur / Collatz–Wielandt row test (backs R-CELLS)

| Lean | Statement |
|---|---|
| `schur_row_test` | Finite N. Suppose that, on `n, m < N`, `M` is symmetric, `M_nm ≥ 0` for `n ≠ m`, `v_n > 0`, and `Σ_{m<N} M_nm v_m ≤ −η v_n` for every `n < N`. Then for every `x`: `Σ_{n,m<N} M_nm x_n x_m ≤ −η Σ_{n<N} x_n²`. |
| `schur_row_test_inf` | Infinite matrix, truncation form. Suppose `M` is symmetric with `M_nm ≥ 0` for `n ≠ m`, `v_n > 0` for every `n`, and for every `N` and every `n < N`: `Σ_{m<N} M_nm v_m ≤ −η v_n`. Then for every `N` and `x`: `Σ_{n,m<N} M_nm x_n x_m ≤ −η|x|²`. |
| `rowSum_mono` | If `M_nm ≥ 0` for `n ≠ m` and `v ≥ 0`, then for `n < N ≤ N'`: `Σ_{m<N} M_nm v_m ≤ Σ_{m<N'} M_nm v_m`. Partial row sums increase. |
| `rows_of_eventually` | Same sign hypotheses. Suppose that for every `n` there is an `N0` with `Σ_{m<N} M_nm v_m ≤ −η v_n` for all `N ≥ N0`. Then the bound holds for every `N` and every `n < N`. |
| `term_le` | AM–GM step: if `M ≥ 0` and `v_n, v_m > 0`, then `M x_n x_m ≤ (M/2)(x_n² v_m/v_n) + (M/2)(x_m² v_n/v_m)`. |

**What Python must establish.** For every row `n` (head rows explicitly, tail rows `n ≥ N` by
the analytic bound), certify the full infinite row sum `Σ_{m≥0} Mhat_nm v_m ≤ −η v_n` in Arb.
`rows_of_eventually` and `rowSum_mono` turn this into the truncation hypothesis. The test vector
must satisfy `v_n > 0` for every `n`, including the analytic tail.

### P2: interval envelope (backs R-CELLS)

| Lean | Statement |
|---|---|
| `Mq_le_Mhat` | If `lo(k) ≤ L(k) ≤ hi(k)` for all `k`, then `Mq L n m ≤ Mhat hi lo n m` entrywise. |
| `rowSum_le_envelope` | If in addition `v ≥ 0`, then `Σ_{m<N} Mq_nm v_m ≤ Σ_{m<N} Mhat_nm v_m`. |
| `envelope_form_bound` | Suppose `L(k) ≥ 0` for `k ≥ 1`, `lo ≤ L ≤ hi`, `v > 0`, and for every `N` and `n < N`: `Σ_{m<N} Mhat_nm v_m ≤ −η v_n`. Then for every `N` and `x`: `C(x) + F⁺_L(x) ≤ −η|x|²`. This holds for every L in the box. No `lo ≥ 0` is needed. |
| `Mq_symm`, `Mq_offdiag_nonneg`, `form_eq_Mq` | `Mq` is symmetric. Its off-diagonal entries are `≥ 0` when `L ≥ 0`. Its quadratic form is `C + F⁺_L`. |
| `stable_of_even_bound` | Let `η > 0`. If `C(x) + F⁺_{λ_T}(x) ≤ −η|x|²` for every `N` and `x` (the even sector at μ* = 0), then Stable at T for every μ* ≥ 0 and every cutoff. The odd sector follows via TcR1.L6. |
| `R_MU` | For a fixed cutoff, if `¬Unstable(Q⁺ at μ* = 0)` and μ* ≥ 0, then `¬Unstable(Q⁺ at μ*)`. |

### P3: monotonicity in T, cell box, weight envelope (backs R-CELLS)

| Lean | Statement |
|---|---|
| `g_antitone` | If `0 ≤ T1 ≤ T2` and `ω > 0`, then `g_{T2}(k, ω) ≤ g_{T1}(k, ω)`. This holds for any `k` and any `κ`. |
| `lam_antitone` | If `0 ≤ T1 ≤ T2`, then `λ(k; T2) ≤ λ(k; T1)`. |
| `lam_mem_cell` | If `0 ≤ T1 ≤ T ≤ T2`, then `λ(k; T2) ≤ λ(k; T) ≤ λ(k; T1)`. |
| `lam_env` | If `sp'` is in the envelope `e` of `sp`, then `(1−e)λ_sp(k; T) ≤ λ_sp'(k; T) ≤ (1+e)λ_sp(k; T)`. |
| `lam_box` | If `sp'` is in the envelope and `0 ≤ T1 ≤ T ≤ T2`, then `(1−e)λ_sp(k; T2) ≤ λ_sp'(k; T) ≤ (1+e)λ_sp(k; T1)`. This is exactly the R-CELLS box. |
| `A_env` | `A(sp') ≤ (1+e)A(sp)` for `sp'` in the envelope. |
| `InEnv_refl` | If `e ≥ 0`, then `sp` is in its own envelope. |

### P4: cover lemma (backs R-CELLS contiguity check)

| Lean | Statement |
|---|---|
| `Tiles lo hi cells` (def) | `cells = [(a_0, b_0), …, (a_K, b_K)]` is nonempty, `a_0 = lo`, `b_i = a_{i+1}` exactly, and `b_K = hi`. |
| `cover` | If `Tiles lo hi cells` and `lo ≤ T ≤ hi`, then some cell `(a, b)` in the list has `a ≤ T ≤ b`. |
| `cover_points` | Breakpoint form. If `K ≥ 1` and `t_0 ≤ T ≤ t_K`, then there is an `i < K` with `t_i ≤ T ≤ t_{i+1}`. |

**Python obligation.** Check `cells[0][0] == T_up`, `cells[i][1] == cells[i+1][0]` and
`cells[-1][1] == t_end` with exact (rational / dyadic) arithmetic. Then check `A_hi ≤ a0 · t_end`
in Arb. See deviation D1.

### P5: Rayleigh certificate, box minimum, μ* monotonicity (backs C-RAYLEIGH)

| Lean | Statement |
|---|---|
| `rayleigh_unstable` | For any form `Q`: if `x ≠ 0` (some `n < N` with `x_n ≠ 0`) and `Q(x) ≥ 0`, then `Unstable Q`. |
| `not_stable_of_unstable` | If `Unstable(Q⁺_{T,μ,Nc})`, then `¬Stable(T)`. |
| `Qp_anti_mu` | For the same cutoff, `μ ≤ μ_hi` implies `Q⁺_{μ_hi}(x) ≤ Q⁺_{μ}(x)`. |
| `boxMinQ` (def) | `boxMinQ = C(x) − 2μ(Σ_{n<Nc} x_n)² + Σ_j w_j F_T(Ω_j; x)`, where `F_T(ω; x) = F⁺(g_T(·, ω); x)` (`atomF`). The weight is `w_j = (1−e)p_j` if `F_T(Ω_j; x) ≥ 0`, and `(1+e)p_j` otherwise. |
| `Qp_env_eq` | `Q⁺` of an envelope spectrum equals `C − μ-term + Σ_j p'_j F_T(Ω_j; x)`. It is affine in the weights (TcR1.L1). |
| `affine_box_min`, `affine_box_min_mem` | On the box `Π[lo_j, hi_j]`, `Σ_j c_j F_j ≥ Σ_j sel_j F_j`. The sign-chosen point `sel` lies in the box. |
| `sel_mul_eq_min` | If `lo ≤ hi`, then `sel(F) · F = min(lo · F, hi · F)`. |
| `boxMin_le` | For every `sp'` in the envelope: `boxMinQ ≤ Q⁺_T[sp'](x)`. |
| `boxMin_attained` | If `0 ≤ e ≤ 1`, some `sp'` in the envelope has `Q⁺_T[sp'](x) = boxMinQ`. So `boxMinQ` is the exact minimum over the envelope. |
| `rayleigh_cert` | If `x ≠ 0` and `boxMinQ(μ_hi, Nc, T, x) ≥ 0`, then `Unstable(Q⁺_{T, μ, Nc}[sp'])` for every `sp'` in the envelope and every `μ ≤ μ_hi`. |

**Python obligation.** Compute an Arb lower bound of
`C(x) − 2μ_hi(Σ_{n<N_c(T)} x_n)² + Σ_j min((1−e)p_j F_j, (1+e)p_j F_j)`, and check that it is `≥ 0`.
`sel_mul_eq_min` makes this valid even when a ball `F_j` contains 0. The same `N_c(T)` must be
used as in the model's cutoff rule. Check `x ≠ 0` exactly.

### P6: assembly (backs R-UNIV and R-CELLS)

| Lean | Statement |
|---|---|
| `stable_of_A_zero` | If `A = 0`, the spectrum is stable for every `T`, every `μ* ≥ 0` and every cutoff. |
| `stable_of_A_le` | Under `H(a0, η)`: if `A ≤ a0 T`, the spectrum is stable at `T` for every `μ* ≥ 0` and every cutoff. This is TcR1.stable_of_ge plus the `A = 0` case. |
| `univ_stable` | Under `H`: if `(1+e)A ≤ A_hi` and `A_hi ≤ a0 T`, every spectrum in the envelope is stable at `T` for every `μ* ≥ 0` and every cutoff. |
| `R_UNIV` | Under `H`: suppose `(1+e)A ≤ A_hi` and `A_hi < a0 T*`. Then `A_hi/a0 < T*`, and for every spectrum in the envelope, every `μ* ≥ 0`, every cutoff and every `T ≥ A_hi/a0` (in particular every `T ≥ T*`), the spectrum is stable. |
| `CellOK κ sp e T1 T2` (def) | `T1 ≥ 0`, and there exist `η > 0` and `v` with `v_n > 0` for all `n` such that for every `N` and `n < N`: `Σ_{m<N} Mhat(hi, lo)_nm v_m ≤ −η v_n`. Here `hi(k) = (1+e)λ_sp(k; T1)` and `lo(k) = (1−e)λ_sp(k; T2)`. |
| `cell_stable` | If `CellOK(T1, T2)` holds, every spectrum in the envelope is stable at every `T ∈ [T1, T2]`, for every `μ* ≥ 0` and every cutoff. |
| `cells_stable` | Under `H`: suppose `Tiles T_up t_end cells`, every cell satisfies `CellOK`, `(1+e)A ≤ A_hi` and `A_hi ≤ a0 t_end`. Then every spectrum in the envelope is stable at every `T ≥ T_up`, for every `μ* ≥ 0` and every cutoff. |
| `R_CELLS` | The same hypotheses plus `T_up < T*` give stability at every `T ≥ T*`. |

R-ODD: no separate check is needed. Every stability conclusion above is derived from an
even-sector μ* = 0 bound valid for all `x`, and TcR1.L6 (`Q⁻(x) ≤ Q⁺_{μ*=0}(|x|)`) transfers it
to the odd sector (`stable_of_even_bound`). The standalone row-sum identity (PROOF Lemma 2) is
not formalised because no rule needs it.

### P7: `classify` (backs S1, S4, S5 and the §4 `classify`)

Lean data types:

```lean
inductive Status | rejectedProven | certifiedCandidate | undecided
                 | excludedAssumption | deferred | error
structure AltEvidence where          -- one per alternative spectrum
  univFires : Bool                   -- R-UNIV:   A_hi < a0·T*
  cellsBelowTarget : Bool            -- R-CELLS:  certified cover and T_up < T*
  rayleighAtOrAboveTarget : Bool     -- C-RAYLEIGH at some T_low ≥ T*

def classifyAlt (ev) := if ev.univFires then rejectedProven
  else if ev.cellsBelowTarget then rejectedProven
  else if ev.rayleighAtOrAboveTarget then certifiedCandidate else undecided

def classify (evs : List AltEvidence) :=
  if evs.isEmpty then undecided
  else if evs.all (classifyAlt · == rejectedProven) then rejectedProven
  else if evs.all (classifyAlt · == certifiedCandidate) then certifiedCandidate
  else undecided
```

Python `certify.classify` must have exactly this control flow. Note that **an empty list of
alternatives gives UNDECIDED**.

Each boolean has the following meaning (`EvidenceSound`). Each `True` must be backed by an Arb
proof of:

* `univFires` → `UnivCert`: there is an `A_hi` with `(1+e)A ≤ A_hi` and `A_hi < a0 T*`.
* `cellsBelowTarget` → `CellsCert`: there exist `T_up`, `t_end`, `A_hi` and `cells` with
  `T_up < T*`, `Tiles T_up t_end cells`, `CellOK` for every cell, `(1+e)A ≤ A_hi` and
  `A_hi ≤ a0 t_end`.
* `rayleighAtOrAboveTarget` → `RayleighCert`: there exist `T ≥ T*`, `N` and `x ≠ 0` with
  `boxMinQ(μ_hi, Ncf(T), T, x) ≥ 0`.

| Lean | Statement |
|---|---|
| `classify_cases`, `classifyAlt_cases` | The result is always REJECTED_PROVEN, CERTIFIED_CANDIDATE or UNDECIDED. It is never EXCLUDED, DEFERRED or ERROR. |
| `status_exclusive` | The three statuses are pairwise distinct, so exactly one holds. |
| `classify_eq_rejected_iff` | `classify evs = REJECTED` if and only if `evs ≠ []` and every `classifyAlt ev = REJECTED`. |
| `classify_eq_certified_iff` | `classify evs = CERTIFIED` if and only if `evs ≠ []` and every `classifyAlt ev = CERTIFIED`. |
| `classify_singleton` | `classify [ev] = classifyAlt ev`. |
| `classifyAlt_rejected_sound` | Under `H`, with sound evidence and `classifyAlt = REJECTED`: there is `T_b < T*` such that every spectrum in the envelope is stable for every `μ* ≥ 0`, every cutoff and every `T ≥ T_b`. |
| `classifyAlt_certified_sound` | With sound evidence and `classifyAlt = CERTIFIED`: there is `T_low ≥ T*` such that, for every spectrum in the envelope and every `μ* ≤ μ_hi`, `Q⁺` at `T_low` (cutoff `Ncf(T_low)`) is Unstable, and therefore `T_low` is not stable. |
| `evidence_consistent` | Under `H`, with `e ≥ 0` and `μ_hi ≥ 0`, sound evidence never has (univ ∨ cells) ∧ rayleigh. So the rule order does not matter. |
| `classify_rejected_sound` | Under `H` (`a0, η > 0`), with sound evidence for every alternative and `classify = REJECTED_PROVEN`: the list is nonempty, and for every alternative there is `T_b < T*` such that every spectrum in its envelope is stable at every `T ≥ T_b`, for every `μ* ≥ 0` and every cutoff. Hence `Tc_model < T*`. |
| `classify_rejected_sound'` | The same, stated for every `T ≥ T*` (ARCH S5 wording). |
| `classify_certified_sound` | With sound evidence (no `H` needed) and `classify = CERTIFIED_CANDIDATE`: the list is nonempty, and for every alternative there is `T_low ≥ T*` such that every spectrum in its envelope and every `μ* ≤ μ_hi` is unstable at `T_low`. Hence `Tc_model ≥ T_low ≥ T*`. |

Rule IDs:
- R-UNIV: `R_UNIV`, `univ_stable`, `stable_of_A_le`.
- R-CELLS: `cells_stable`, `R_CELLS`, `cell_stable`, `envelope_form_bound`,
  `schur_row_test_inf`, `rows_of_eventually`, `lam_box`, `cover`.
- R-MU: `R_MU`, `stable_of_even_bound`.
- R-ODD: TcR1.L6 via `stable_of_even_bound`.
- C-RAYLEIGH: `rayleigh_cert`, `boxMin_le`, `boxMin_attained`, `sel_mul_eq_min`, `Qp_anti_mu`.
- classify / S5: `classify_*`.

## Deviations from ARCH (all conservative or clarifying)

* **D1. Cell end point.** ARCH says the cells tile `[T_up, A_hi/a0]`. But `A_hi/a0` is an Arb
  ball, not an exact number, so exact contiguity up to it cannot be checked. Lean instead uses
  cells tiling `[T_up, t_end]` with exact endpoints, plus the Arb inequality `A_hi ≤ a0·t_end`.
  R-UNIV then covers `T > t_end`.
* **D2. Cutoff.** The cutoff rule is an arbitrary function `Ncf : T ↦ N_c(T)`. REJECTED holds for
  every cutoff and every `μ* ≥ 0`, not only for `[0, μ_hi]`. CERTIFIED holds for the run's cutoff
  rule only (see F3).
* **D3. Empty list.** `classify [] = UNDECIDED` (see F1).
* **D4. Stronger quantifier order for CERTIFIED.** A single `T_low` (per alternative) works for
  every envelope spectrum and every `μ* ≤ μ_hi`.
* **D5. A = 0.** Handled separately: a spectrum with `A = 0` is always stable. TcR1.stable_of_ge
  needs `A > 0`, and the envelope with `e = 1` can reach `A = 0`.

## Flaws and pitfalls found in ARCH's logic

* **F1. Empty conjunction.** S5's "for every alternative spectrum" is vacuously true when the run
  produced no spectrum, for example when all DFPT smearings failed. A literal implementation
  would then return REJECTED_PROVEN or CERTIFIED_CANDIDATE with no evidence. `classify` must
  return UNDECIDED, or ERROR, for an empty list.
* **F2. Diagonal of Mhat.** ARCH writes `lower(2·sum_{k<=n} λ(k))`. The sum must start at
  `k = 1`, because `λ(0)` has cancelled. Including `k = 0` subtracts `2λ(0) = 4∫α²F/ω > 0` from
  the diagonal, which is unsound: the envelope would no longer dominate the true matrix.
* **F3. Certified verdicts depend on the cutoff convention.** C-RAYLEIGH is monotone in μ* only
  at a fixed cutoff, and there is no monotonicity in `N_c`. A CERTIFIED_CANDIDATE is therefore a
  statement about the model with the block cutoff `N_c(T, ω_c)` for the configured `ω_c`. It is
  not a statement about other cutoffs or about the column-only convention. The verdict's
  `assumptions` should record `ω_c` and the cutoff rule. REJECTED verdicts hold for every cutoff.
* **F4. "R-CELLS bound T_up < T\*".** This is correct as written, with a strict inequality,
  which gives `Tc_model ≤ T_up < T*`. If `T_up = T*` were allowed, the conclusion would weaken
  to "stable for all `T ≥ T*`", which gives only `Tc_model ≤ T*`, not `< T*`.
* **F5. Envelope scope.** The `e`-envelope perturbs atom weights only. It does not cover
  uncertainty in frequency positions. Frequency shifts must come from the alternative spectra.
* **F6. Box minimum with undetermined signs.** If an Arb ball `F_j` straddles 0, the
  "sign-chosen endpoint" is ambiguous. Use `min((1−e)p_j F_j, (1+e)p_j F_j)` (`sel_mul_eq_min`).
* **F7. Reported T_low across alternatives.** Each alternative has its own `T_low`. The combined
  certified lower bound is the minimum over alternatives, and the combined upper bound `T_up` is
  the maximum over alternatives. For R-UNIV the bound is `A_hi/a0`.

---

# Library `MultiBand`: Lean ↔ wMB/PROOF.md and ARCH §4c (R-MULTI)

Lean project `/home/claude/tc/lean`, files `MultiBand/{Defs,Algebra,Identities,Odd,Main,Pipeline,Faithful}.lean`.
The library imports `TcR1` only, and `TcR1` and `PipelineLogic` are unchanged. It has no `sorry`
and no `axiom`. All 67 theorems depend on `[propext, Classical.choice, Quot.sound]` or fewer; see
`/home/claude/tc/lean/multiband_axioms_output.txt`. To re-check:
`nice lake build && lake env lean CheckMultiBand.lean`.

## Conventions
* Bands `i : Fin B`, DOS `N_i = dos i > 0`. The gap is `φ_{n,i} = φ i n`, truncated to `n < N`.
  `|x|² = Σ_i N_i Σ_n φ_{n,i}²`, i.e. `x = √N φ` (`nrm2_eq_x`).
* Spectrum `MSpec`: finitely many atoms `s` at `Ω_s > 0` with symmetric nonnegative matrices
  `W_s` (`W_ij = N_i α²F_ij`). `𝒜_ij = Σ_s W_s,ij / N_i` (`Amat`).
* `tp` = 2π (any real). The forms `Qplus`, `Qminus`, `Pterm` are those of wMB §1 and §1.1
  verbatim, with `l_w(0)` kept in both `B` and `D`. `Qplus_lambda_form` rewrites `Qplus` with the
  physical `λ_ij(k)`: the Z term uses the row sum `Σ_j λ_ij`, the pairing uses the full matrix.
* Hypotheses:
  - `TcR1.H tp a0 η`: Prop. 6 of final/PROOF.md, with `η ≥ 0` and `a0 > 0`.
  - `PropN tp a0`: `Q^E_{-a0,w}(y) = -C(y) - a0 yᵀG⁺_w y ≤ 0` for every `w > 0` and every `y`.
    This is wMB Prop. N verbatim.
  - `CoulPSD`: `Σ_ij N_i μ*_ij S_i S_j ≥ 0` for all `S`.
  - A positive vector `v` with `𝒜v = κv` (Perron vector) or `𝒜v ≤ κv` (Collatz–Wielandt).
    Perron–Frobenius is not formalized.
* `MStable` means neither sector is `MUnstable`, where `MUnstable` means `ρ ≥ 0` in the norm `|x|`.

## Theorems (plain statements)
| Lean | Statement |
|---|---|
| `lemma11_even` / `lemma11_odd` | Lemma 1.1: `Q^± = -Σ_i N_iC(φ_i) + (1/T)Σ_s Σ_ij W_s,ij G^±_{Ω_s/T}(φ_i,φ_j) - P(φ)` (minus the Coulomb term for `Q⁺`). This holds for every `T`. |
| `Pterm_nonneg` | `P ≥ 0` for `T > 0`. |
| `lemma41` | Lemma 4.1 (exact): if `𝒜v = κv` and `κ, T, v_i, N_i ≠ 0`, then `Q⁺ = (1/κ)Σ_s Σ_ij v_iv_jW_s,ij [Q^E_{κ/T,Ω_s/T}(ψ_ij) + Q^E_{-κ/T,Ω_s/T}(χ_ij)] - P - Coulomb`. |
| `QE_pm` | If `Q^E_{±a0,w}(y) ≤ 0` and `0 ≤ a ≤ a0`, then `Q^E_{±a,w}(y) ≤ -(1 - a/a0)C(y)`. |
| `thm1_even` | Theorem 1, even sector. Assume H(a0,η) with `η ≥ 0` and `a0 > 0`, Prop N, `N_i > 0`, `v > 0` with `𝒜v ≤ κv`, `N_iμ*_ij` PSD, any cutoff `N_c`, `T > 0` and `κ ≤ a0T`. Then `Q⁺(φ) ≤ -(1 - κ/(a0T))|x|²`. |
| `thm1_even_via_identity` | The same conclusion, proved along wMB §4.3 through `lemma41`, for `𝒜v = κv` with `κ > 0`. |
| `thm3_odd` | Theorem 3: for `T > 0` and `N_i ≥ 0`, `Q⁻(φ) ≤ -Σ_i N_iC(φ_i) ≤ -|x|²`. There are no other hypotheses. |
| `mstable_of_gt` | If `T > 0` and `κ < a0T`, both sectors are stable, with `c = 1 - κ/(a0T) > 0` (even) and `c = 1` (odd). |
| `unstable_imp_le` | An unstable `T > 0` satisfies `T ≤ κ/a0`, so `Tc ≤ κ/a0`. |
| `sub_env` | If `𝒜v ≤ κv`, `𝒜' ≤ (1+e)𝒜` entrywise (any atoms and frequencies), `1+e ≥ 0` and `(1+e)κ ≤ κ_hi`, then `𝒜'v ≤ κ_hi v`. |
| `sub_of_K` | Let `r_i > 0` with `r_i² = N_i`, `K_ij = W^tot_ij/(r_ir_j)`, `u > 0` and `Ku ≤ κu`. Then `v = u/r` satisfies `𝒜v ≤ κv`. |
| `R_UNIV_MB` | Multiband R-UNIV. Under H and Prop N, take `v > 0` with `𝒜v ≤ κv`, `(1+e)κ ≤ κ_hi`, `κ_hi < a0T*` and `T* > 0`. Then `κ_hi/a0 < T*`, and every spectrum with `𝒜' ≤ (1+e)𝒜` is stable at every `T ≥ T*`, for every PSD `N_iμ*_ij` and every cutoff. |
| `R_UNIV_MB_eig` | The same with the Perron eigenvector, `κ ≤ κ_hi` and no envelope. |
| `R_UNIV_MB_rowsum` | Fallback: `max_i Σ_j 𝒜_ij ≤ κ_hi < a0T*` gives stability for all `T ≥ T*` (`v = 1`). |

## What Python must establish for R-MULTI
1. **W must be symmetric.** The theorem is about a model with `N_iα²F_ij = N_jα²F_ji`. If the
   data are symmetrized (`W̄`), the verdict is a statement about the symmetrized model.
2. **A certificate for `κ_hi`.** Lean does not formalize Perron–Frobenius. Either
   - produce `u > 0` and check `K u ≤ κ_hi u` componentwise in Arb (`sub_of_K`), with
     `K = N^{-1/2} W̄ N^{-1/2}`. For any `κ_hi > ρ(K)` such a `u` exists: for example
     `u = (κ_hi − K)^{-1}1`, which also works for reducible `K`. Or
   - use the max-row-sum fallback (`R_UNIV_MB_rowsum`), which needs no vector.

   An Arb enclosure of `λ_max(K)` alone is sound mathematically, through Perron–Frobenius, but
   that step is outside Lean.
3. **Margin.** Use a strict inequality `κ_hi < a0·T*`. If `𝒜` is uncertain by a relative `e`
   entrywise, use `(1+e)κ ≤ κ_hi` (`sub_env`).
4. **Coulomb.** `N_iμ*_ij` must be PSD (CE 4 of wMB shows the rule fails otherwise).

## Deviations from wMB/PROOF.md
* **MB-D1. Collatz–Wielandt instead of irreducibility.** Theorem 1 is proved for any `v > 0`
  with `𝒜v ≤ κv`, not only for the Perron vector of an irreducible `𝒜`. This also covers the
  reducible case of §4.4, which wMB treats by a block argument, and the envelope `sub_env`.
  `lemma41` keeps the exact equality version.
* **MB-D2. Boundary `T = κ/a0`.** With Prop N stated without margin (`≤ 0`), the hypotheses
  give `Q⁺ ≤ -(1-κ/(a0T))|x|²`, and the constant is 0 at `T = κ/a0`. Stability is therefore
  proved for `T > κ/a0`, as in wMB §4.3, which assumes `T > κ/a0`. This is enough for
  `Tc ≤ κ/a0` and for the strict pipeline rule.
* **MB-D3. Theorem 3 proved differently.** wMB uses positivity of the Toeplitz symbol and
  Lemma 2 of final/PROOF.md. Lean uses an elementary Schur test instead. The odd kernel
  `l(|n−m|) − l(n+m+1)` is `≥ 0` because `l_w` decreases. Its truncated row sums equal
  `D_n − Σ_{k=N−n}^{N+n} l(k) ≤ D_n` (`odd_row_le`, with the `N ≤ n` case handled separately).
  The conclusion is the same, `Q⁻ ≤ -Σ N_iC ≤ -|x|²`.
* **MB-D4. Finite atomic spectra, `x` in φ-variables, `Fin B` bands.** General measures,
  `λ_ij(0) = ∞`, the operator domain and anisotropic continuum limits are not formalized.
  Theorem 2 (`Λ_A`), sharpness (§4.6) and the ordering of functionals (§4.7) are not formalized.
* **MB-D5. PSD hypothesis.** Only the quadratic form `Σ U_ij S_iS_j ≥ 0` is used. Symmetry of
  `U` is not needed.

## Flaws found in wMB/PROOF.md
None in the parts formalized: Lemma 1.1, Lemma 4.1, the assembly in §4.3, Theorem 3 and the
§7 rule. All of them hold as stated. Minor notes:
* In §4.5, "Σ_j W_ij φ_iᵀDφ_i = Σ_ij W_ij ½(…)" should have `Σ_ij` on the left. This is a typo.
* §4.3 needs only `η ≥ 0` from Prop. 6.
* At `T = κ/a0` the argument gives only `ρ⁺ ≤ 0`, not `ρ⁺ < 0` (MB-D2). This is harmless.

---

# `PipelineLogic.Multiband`: Lean mirror of `combine_multiband` (ARCH §4c)

File `/home/claude/tc/lean/PipelineLogic/Multiband.lean`. It imports `PipelineLogic.Classify`
and `MultiBand`. It has no `sorry` and no `axiom`. All 26 theorems depend on
`[propext, Classical.choice, Quot.sound]` or fewer; see `pipeline_axioms_output.txt`, which
now has 82 lines.

## The Lean definition, verbatim from `rtsc/certify/classify.py`

```lean
inductive RejectionModel | isotropic | multiband
def combineMultiband (sb : Status) (mbFires : List Bool) (model : RejectionModel)
    (nSpectra : Option ℕ := none) : Status :=
  match model with
  | .isotropic => sb
  | .multiband =>
    if nSpectra == some 0 then .undecided
    else
    let mbRej : Bool := !mbFires.isEmpty && mbFires.all id
    if mbRej && sb != .certifiedCandidate then .rejectedProven
    else if sb == .certifiedCandidate && !mbRej then .certifiedCandidate
    else .undecided
def classifyFull (evs : List AltEvidence) (mbFires : List Bool) (model : RejectionModel) :
    Status :=
  combineMultiband (classify evs) mbFires model (some evs.length)
```

`n_spectra : int | None` maps to `nSpectra : Option ℕ`, with `None` ↦ `none` (the default).
`classifyFull` passes `some evs.length`, because Python passes `n_spectra=len(spectra)` and
there is one `AltEvidence` per spectrum. The assumption `BAND-DATA-NOT-PER-ALTERNATIVE` is an
annotation only. The Python `rejection_model` string maps as follows: `"isotropic"` ↦ `.isotropic`, and any other
string ↦ `.multiband`. `classify` accepts only `"multiband"`. `recheck` does not validate the
string, but an unknown value falls into the multiband branch, which is the conservative one.

## Meaning of an R-MULTI boolean (`RMultiCert a0 T* a`)

For one band-coupling alternative `a = (sp, N, e, fires)`, `fires = True` must mean:

* `T* > 0`;
* `e ≥ 0`;
* `N_i > 0`;
* there is a `κ_hi < a0·T*` with `CWCore`, i.e. some `v > 0` and `κ` with `𝒜v ≤ κv` and
  `(1+e)κ ≤ κ_hi`.

Here `sp` is the band model the check is about. For the eig route with non-symmetric data, its
total weights are `W̄`. The four `check_route` routes each produce `CWCore`:

| Python route | Lean | Exact check |
|---|---|---|
| `"cw"` | `route_cw` | `v > 0` and `(1+e)(𝒜v)_i ≤ κ_hi v_i` |
| `"rowsum"`, source `A_K` | `route_rowsum` | `(1+e)Σ_j 𝒜_ij ≤ κ_hi` |
| `"rowsum"`, source `row_sum_upper_K` | `route_rowsum_declared` | `Σ_j 𝒜_ij ≤ R_i` (data assumption) and `(1+e)R_i ≤ κ_hi` |
| `"eig"` | `route_eig` | `r_i > 0`, `r_i² = N_i`, `u > 0` and `(1+e)Σ_j W^tot_ij/(r_ir_j) u_j ≤ κ_hi u_i` (via `MultiBand.sub_of_K`) |

The envelope is `MBEnv a sp'`: `𝒜'_ij ≤ (1+e)𝒜_ij` entrywise, with the same `N`. Atoms and
frequencies may differ.

## Theorems

| Lean | Statement |
|---|---|
| `combineMultiband_isotropic` | In isotropic mode the result is the single-band status, unchanged. |
| `combine_mb_rejected_iff` | Multiband mode: REJECTED if and only if `nSpectra ≠ some 0`, (`mbFires ≠ []` and every entry is True), and `sb ≠ CERTIFIED`. |
| `combine_mb_zero_spectra` | Multiband mode with `nSpectra = some 0`: UNDECIDED. |
| `classifyFull_zero_spectra`, `classifyFull_zero_spectra_not_rejected` | **Zero spectra.** `classifyFull [] mbFires model = UNDECIDED` in both modes, whatever the band evidence, so it is never REJECTED. |
| `combine_mb_certified_iff` | Multiband mode: CERTIFIED if and only if `nSpectra ≠ some 0`, `sb = CERTIFIED`, and not (`mbFires ≠ []` and all True). |
| `combine_mb_cases`, `classifyFull_cases` | **Exhaustive.** The result is REJECTED, CERTIFIED or UNDECIDED. In multiband mode this holds for every input status. |
| `classifyFull_exclusive` | **Exclusive.** No two of the three statuses hold at once. |
| `combine_mb_empty_not_rejected`, `classifyFull_empty_not_rejected` | **An empty band-coupling list never gives REJECTED in multiband mode**, whatever the single-band status. |
| `combine_mb_rejected_needs_mb` | A multiband-mode REJECTED requires R-MULTI on every alternative. A single-band rejection alone is never enough. |
| `combine_certified_imp` | In either mode, CERTIFIED requires the single-band status to be CERTIFIED. |
| `classifyFull_mb_rejected_sound` | Assume H(a0,η) with `η ≥ 0` and `a0 > 0`, Prop N, and every True backed by `RMultiCert`. If the multiband-mode result is REJECTED, then both the spectrum list and the band-coupling list are nonempty. Moreover, for every alternative there is a `κ_hi` with `κ_hi/a0 < T*` such that every envelope member, every PSD `N_iμ*_ij`, every cutoff and every `T > 0` with `T > κ_hi/a0` is `MStable`. This is the reported `T_up = κ_hi/a0`, rounded up. |
| `classifyFull_mb_rejected_sound'` | The same, stated for every `T ≥ T*`: stable in the multiband model for every band-coupling alternative, every envelope member, every PSD Coulomb matrix and every cutoff. |
| `classifyFull_iso_rejected_sound` | Isotropic-mode REJECTED is exactly the single-band REJECTED (`classify_rejected_sound`). This is a statement about the single-band model only (M-ISO-SINGLEBAND, MULTIBAND_RISK). |
| `classifyFull_certified_sound` | **What CERTIFIED means, in either mode.** The single-band alternative list is nonempty, and for every alternative spectrum there is a `T_low ≥ T*` such that every spectrum of its envelope and every `μ* ≤ μ_hi`, at cutoff `Ncf(T_low)`, is unstable **in the single-band (DOS-averaged, isotropic) model**. Nothing is proved about the multiband model; this is the content of the `M-ISO-SINGLEBAND` tag. |
| `Qplus_const`, `nrm2_const`, `munstable_of_iso` | Constant-gap trial vector. For `T > 0` and `N_tot = Σ N_i > 0`: `Q⁺_multi(φ_i ≡ y) = N_tot · Q⁺_iso(y)`, where `α²F_iso = Σ_ij W_ij/N_tot` (`isoAvg`) and `μ*_iso = Σ_ij N_iμ*_ij/N_tot` (`muAvg`). Single-band instability of the average therefore implies multiband instability. |
| `certified_transfers_to_multiband` | **Conditional lift of M-ISO-SINGLEBAND.** Suppose the single-band C-RAYLEIGH certificate holds for `sp0` at `T_low > 0`, the true band model's DOS average `isoAvg(sp, N)` lies in the e-envelope of `sp0` (same atoms and frequencies), and `μ*_iso ≤ μ_hi`. Then the multiband model is not stable at `T_low`, with the same cutoff. **The pipeline checks neither hypothesis.** |

## Discrepancies and notes on the Python logic

The decision logic is sound as written. I found no case where the Python returns a status that
the theorems above do not support. Notes:

* **MB-N1. Bracket labelling.** In multiband mode with status CERTIFIED or UNDECIDED, the
  reported `[T_low_K, T_up_K]` is the single-band bracket. `T_up` from R-CELLS is not an upper
  bound in the multiband model, because the multiband `Tc` can exceed the DOS-averaged one (the
  wMB counterexample gives 9 %). It is valid only under `M-ISO-SINGLEBAND`. The verdict carries
  that tag for CERTIFIED, but not for UNDECIDED. In multiband REJECTED the bracket is reset to
  `[None, κ_hi/a0]`, which is correct.
* **MB-N2. Detailed balance assumption.** The assumption string for the cw/rowsum routes reads
  "detailed balance with some N > 0 (N_i A_ij = N_j A_ji)", which is the integrated form. The
  Lean model (`MSpec.W_symm`) needs `N_iα²F_ij(ω) = N_jα²F_ji(ω)` at every frequency (per
  atom). This is the physical detailed balance, but it is stronger than the integrated
  equality. The wording should say "for every ω".
* **MB-N3. Multiband REJECTED and the single-band spectra (fixed).** Before the `n_spectra`
  guard, `spectra = []` with R-MULTI firing gave REJECTED. That was logically sound for the
  multiband model, but it contradicted F1. It now gives UNDECIDED (`classifyFull_zero_spectra`). A single-band
  CERTIFIED together with a firing R-MULTI gives UNDECIDED. By
  `certified_transfers_to_multiband`, that conflict is impossible when the single-band spectrum
  is the DOS average of the band data and `μ*_iso ≤ μ_hi`, so it signals inconsistent inputs.
* **MB-N4. CERTIFIED could be lifted to the multiband model.** This would need two exact
  checks: that each single-band spectrum (with its envelope) contains the DOS average of the
  band-coupling model, and that `Σ N_iμ*_ij/N_tot ≤ μ_hi`. The pipeline does neither today, so
  the `M-ISO-SINGLEBAND` tag is required and correct.

---

# Library `DiscoveryLogic`: Lean ↔ `rtsc/discovery/` (ledger state machine)

Lean files: `/home/claude/tc/lean/DiscoveryLogic/{Ledger,Loop,Continue,Soundness}.lean`. There is
no `sorry` and no `axiom`. `#print axioms` output: `/home/claude/tc/lean/discovery_axioms_output.txt`.
68 of the 69 theorems depend on `[propext, Classical.choice, Quot.sound]`; `le_foldr_max` depends
on `[propext, Quot.sound]`. The model is the post-audit code: `ledger.py` schema
`rtsc.discovery.ledger/2`, `loop.py` with `_reconcile`, and `goals.counts_for`. The line numbers
are for the files as of 17:45 UTC on 2026-10-27.

## Definitions
| Lean | Python | Notes |
|---|---|---|
| `HStage`, `Stage` | `ledger.py:51-58` (`HEURISTIC_STAGES`, `ALL_STAGES`) | closed types; an unknown stage raises `ValueError` at `ledger.py:308-309` and is rejected by the SQL trigger at `ledger.py:117` |
| `Origin` (`classify` / `stability` / `other`) | origin strings; `CLASSIFY_ORIGIN` at `pipeline.py:55`, `startswith("rtsc.stability")` | the two predicates are disjoint, so this is a function of the string |
| `heuristicAllowed` | `ledger.py:51-57` | |
| `isProven` | `status.py:27` `PROVEN_STATUSES` | `Status` is `PipelineLogic.Status` |
| `Goal`, `inScopeB` | `goals.py:67-79` (`Goal`, `in_scope`: `P ≤ P_max`) | goal ids are `ℕ` |
| `Event` | `verdicts` columns, `ledger.py:99-106` | bracket, rule ids, flags, paths, digests omitted |
| `appendAllowed` | `append_verdict` guard, `ledger.py:307-318` | exact clause-by-clause mirror |
| `sqlGuard` | trigger `verdicts_guard`, `ledger.py:115-126` | |
| `commit` | insert at `ledger.py:319-331` inside `Ledger.tx`, `ledger.py:160-175` | all or nothing (rollback) |
| `currentEvent`, `currentStatus` | view `current_status`, `ledger.py:127-128`; `Ledger.current`, `ledger.py:333-345` | latest event per (material, goal) |
| `Phase` | `ROUND_PHASES`, `ledger.py:59` | |
| `hEv`, `initialEvent` | `_initial_status_for`, `loop.py:366-384` | `in_scope` is `g.in_scope(P)` in every branch |
| `prescreenEvent` | `loop.py:550-557` | |
| `selectEvent` | `loop.py:597-606` | |
| `configEvent` | `_reconcile`, `loop.py:205-218` | |
| `execErrorEvent` | `_record`, `loop.py:762-772` | |
| `mismatchEvent`, `pOk`, `pTol` | `_record`, `loop.py:779-789` (`abs(dp - P) > 1e-9`) | |
| `pipelineEvent` | `_record`, `loop.py:811-818` | T* is `g.T*` (`goals.py:170`, check at `loop.py:775-778`); hashes are the session's |
| `alreadyRecorded` | `loop.py:791-795` | record path, digest and config hash are abstracted to the round |
| `Config` | `DiscoveryConfig.enabled_goals`, `certify_hash`, `dft_hash`, `goals.py:174-186` | hashes are `ℕ` |
| `MatRow`, `Data.fresh` | `materials` at `ledger.py:76-82`; `next_seq` at `ledger.py:201-203` | id = `seq` |
| `GenItem`, `DOp.genItem` (with dedupe decision `dd`), `newMat` | `_generate`, `loop.py:411-465`; `_match`, `loop.py:347-360` | a hit is a duplicate (`loop.py:437-442`), otherwise a new material with its initial statuses (`loop.py:443-459`) |
| `DOp.prescreen` | `_record_prescreen`, `loop.py:521-558` | `defer` = `_ml_decision`, or a duplicate after ML relaxation |
| `DOp.addJob`, `DOp.selectDefer` | `_select`, `loop.py:589-607` | |
| `DOp.exec` | `_record` (all branches), `loop.py:758-818` | |
| `DOp.reinit`, `DOp.configDefer` | `_reconcile`, `loop.py:192-218` | |
| `addGoalEvents`, `Op.newSession` | `__init__` at `loop.py:113-124`, `_add_goal_statuses` at `loop.py:386-391` | a missing old goal id is a `ConfigError` and changes nothing |
| `Reachable` | any sequence of the above from an empty ledger | covers every crash point |
| `Env`, `runGen` … `runDone`, `finish`, `resume` | `_run_round` at `loop.py:318-336`; `run` at `loop.py:144-157` | |
| `toScreen`, `roundMats` | `loop.py:494-495`, `loop.py:595` | |
| `EnvAgree` | (hypothesis) | determinism assumptions, stated explicitly |
| `countsFor`, `Satisfies` | `counts_for` at `goals.py:211-233`; `certified_ids` / `counts` at `loop.py:169-175`; `report.py` `goal_rows` | |
| `CertParams`, `PipelineSound` | (hypothesis) | what a `certify_hash` stands for, and the PipelineLogic contract on each classify row |

## Theorems (plain statements)
* `appendAllowed_eq_sqlGuard`: for every row, the Python guard and the SQL trigger agree. A raw
  SQL insert therefore cannot store anything that `append_verdict` would refuse.
* `appendAllowed_heuristic`, `appendAllowed_proven`, `appendAllowed_excluded`,
  `appendAllowed_certified_inScope`: an accepted row from a heuristic stage is DEFERRED or
  ERROR. A proven status needs `pipeline` and `rtsc.certify.classify`. EXCLUDED_ASSUMPTION needs
  `pipeline` and `rtsc.stability*`. CERTIFIED needs `in_scope = 1`.
* **S2-ext**, `S2_ext` and `no_proven_from_heuristic`: in every reachable ledger, every row from
  generate, scope, config, prescreen, select, policy or executor is DEFERRED or ERROR. So no
  reachable ledger holds REJECTED_PROVEN, CERTIFIED_CANDIDATE or EXCLUDED_ASSUMPTION from a
  heuristic stage. `proven_only_from_classify` and `excluded_only_from_stability` state S1 and
  S3 on the ledger.
* `current_certified_from_classify`: a current CERTIFIED row comes from classify and has
  `in_scope = 1`. `reachable_proven_pressure`: every proven row has a DFT pressure within
  `1e-9` GPa of the material pressure (audit F3).
* **Append-only**: `step_prefix` / `step_ext` say that every operation, including
  `newSession`, extends the verdict list by rows the guard accepted. `resume_prefix` says the
  same for a whole round. `snap_prefix` says every crash snapshot is a prefix of the final
  ledger.
* **Accounting**, `accounting`: in every reachable state, every registered material has exactly
  one current status (`∃!`) for every goal of the current config. That status is the latest
  row. `generated = materials + duplicates` holds, every duplicate links to an existing
  material, and material ids are unique.
* **Goal soundness**, `satisfies_sound`: suppose `counts_for(g, current row of m)` holds with
  the current hashes `ch`, `dh`, and suppose `PipelineSound` holds. Then:
  * the material pressure and the DFT pressure are both `≤ P_max(g)` and differ by at most
    `1e-9`;
  * there is a nonempty list of alternative spectra of `m`, computed under DFT config `dh` at
    that pressure;
  * for each alternative there is `T_low ≥ T*(g)` at which every spectrum of its envelope is
    unstable in the even sector, for every `μ* ≤ μ_hi` of the `[certify]` config `ch`. This
    means `Tc_model ≥ T*(g)` in the isotropic single-band model; it is the conclusion of
    `PipelineLogic.classifyFull_certified_sound`.
  
  `counts_for` is self-contained: the theorem does not need the guard, and it holds whatever
  config changes happened between sessions.
* **Continuation**:
  * `continue_after_crash`: take a crash after `start_round` or after any committed phase of a
    round, followed by `resume` under `EnvAgree`. The result is exactly the uninterrupted state.
    The same holds for a crash before `start_round`.
  * `continue_currentStatus`: the corollary for current statuses.
  * `newSession_same`: `continue` with the same config adds no rows.
  * `continue_with_session`: the same statement including the session start, assuming
    `_reconcile` (`R`) writes nothing at the crash point.
  * `continue_prescreen_partial`: a crash after any number of the per-material pre-screen
    transactions, then a re-run of the phase (which skips screened materials), gives the
    uninterrupted phase result.

## Hypotheses (explicit, not proved)
* `EnvAgree k e1 e2`: the original and the continued process see the same generator outputs
  for the stored requests of round `k` and the same dedupe decisions. They also see the same
  per-material pre-screen decisions (independent of batching), the same selection (policy
  deterministic given seed, round and ledger, **and the same DFT budget**), and the same
  pipeline outcomes (pipeline resume reuses finished records).
* `continue_with_session`: `_reconcile` writes nothing at the crash point.
* `PipelineSound`: every classify row with hashes `ch`, `dh` and DFT pressure `pd` was produced
  by `classifyFull` at the row's T*, under the `[certify]` parameters of `ch`, on spectra
  computed for that material under `dh` at `pd`, with `EvidenceSound` evidence. On the Python
  side this rests on several facts:
  * `goal_pipeline_config` sets `T_target_K = g.T*` (`goals.py:170`);
  * `record_reusable` compares the full pipeline config hash (`pipeline.py:259`), so a reused
    record was made under the current `[certify]` and DFT config;
  * `_record` writes the session's `certify_hash(g)` and `dft_hash` (`loop.py:811-818`);
  * `_reconcile` re-queues materials whose rows carry another `dft_hash` (`loop.py:205-218`).
* Goal ids are unique, and a `newSession` keeps every old goal id (`loop.py:118-120`).

## Modelling simplifications
* Strings (ids, origins, hashes) are abstracted to `ℕ` or to finite types. Floats are `ℝ`.
* The dedupe decision (confirmed fingerprint / StructureMatcher, `loop.py:347-360`) is an
  arbitrary function of the ledger. A hit that is not a registered material creates a new
  material; in Python every hit is a registered material.
* The select-DEFERRED set is "valid materials of round `k` not chosen". Python also excludes
  ML-relaxed duplicates (`dup_of`, via `_pool`). The difference adds only heuristic DEFERRED
  rows in the model.
* `_record`'s already-recorded test compares record path, digest and config hash. The model
  abstracts these to the round.
* Job states (PENDING/SUBMITTED/DONE/RETRY/STALE), `ml_screen` contents, reports and policy
  state carry no verdicts and are not modelled, except as `jobs` and `screened`.

## Discrepancies found (Python vs. the invariants)
**Pre-fix code (17:07 version).** These were reported before the audit fixes and are resolved by
the current code:
* D1, the same as audit F5. On `continue -c` with a changed T* or P_max for an existing goal
  id, `__init__` compared only goal ids (old `loop.py:117-124`). The report and the CSV listed
  every current CERTIFIED row under "Tc_model ≥ <new T*> proved" without checking the row's
  T* (old `report.py:118`). Only `certified_ids` re-checked T* and the scope. The fix is
  `counts_for` plus `_reconcile`.
* D2. The `append_verdict` guard trusted the caller-supplied `in_scope` and `T_star_K`.
  Goal soundness held only through the call sites (`_record`, old `loop.py:651-654, 675-681`).
  The fix makes `counts_for` recompute both.
* D3. `accounting()` computes `generated` as `len(materials) + n_duplicates`
  (`ledger.py:432-434`), so "generated = materials + duplicates" holds by construction there.
  The real invariant, that each generator output yields exactly one material row or one
  duplicate row, is the `count` field of `Inv` (theorem `accounting`). It is not checked
  against `counts["generated"]` in Python.

**Current code.**
* D4, continuation and the DFT budget. `continue` starts fresh stop limits, so `budget_left`
  (`loop.py:152-153`) can differ from the uninterrupted run when `max_dft_jobs` is binding.
  `_select` then chooses a different set, and the continued run is not "the same choices".
  This breaks `EnvAgree.select`. It is by design (fresh limits), but the claim that
  `continue` reproduces the same choices holds only when `max_dft_jobs` is not binding.
* D5, `_reconcile` is not idempotent for DFT-pressure mismatches. A DONE material can have
  one goal with a pipeline row (for example UNDECIDED at 300 K) while another in-scope goal
  has the ERROR row `DFT_PRESSURE_MISMATCH` (certified at 77 K on spectra with an unknown or
  different pressure). That second row is non-retryable, so the job stays DONE. At every
  session start, `_reconcile` puts the second goal in `todo` (`loop.py:224-232`), because its
  current row is not a pipeline row. It then re-runs `certify.classify` on the stored spectra
  in a new `reconcile/s<sid>/` directory and appends another identical mismatch ERROR row
  (`loop.py:779-789`). The current status is unchanged. However, computation is redone on
  every `continue` and the history grows without bound, against DISCOVERY.md §2 ("Nothing
  already computed is redone"). The theorem `continue_with_session` therefore holds with
  `hR` only up to these duplicate rows. This finding comes from reading the code; it was not
  executed.
* D6, a minor edge case that errs on the safe side. `_record` accepts a proven row when
  `|dp − P| ≤ 1e-9` (`loop.py:779`). `counts_for` additionally requires `in_scope(dp)`
  (`goals.py:227`). A row with `P ≤ P_max < dp ≤ P + 1e-9` is stored as CERTIFIED but never
  counts. The report shows it as "certified, does not count". This is consistent and safe.
* No write path bypasses the guard. The only `INSERT INTO verdicts` is at `ledger.py:320-321`;
  there is no other write to `verdicts`, and `current_status` is a view.
