# CORRESPONDENCE: Python decisions ↔ Lean theorems (audited)

This file maps every decision the Python program makes to the result that justifies it. That
result is either a Lean theorem or a proof outside Lean, and proofs outside Lean are labelled
**[non-Lean]**. A ranked list of gaps follows the map. The Lean statements were read from the
`.lean` sources in `/home/claude/tc/lean`. Python line numbers refer to the package as audited
on 2026-10-27. `docs/CORRESPONDENCE_LEAN.md` holds the Lean-side theorem list and ERRATA
F1–F7 / D1–D5.

Re-check:

```sh
cd /home/claude/rtsc && nice python -m pytest -q tests/test_correspondence.py      # 50 passed (~4 s)
cd /home/claude/tc/lean && nice lake build && lake env lean CheckPipeline.lean && lake env lean CheckMultiBand.lean
```

---

## 0. Trusted base

| Item | What is trusted | Where it enters |
|---|---|---|
| Lean 4 kernel v4.35.0 and Mathlib v4.35.0 | Checking of `TcR1`, `PipelineLogic` and `MultiBand`. No `sorry` and no `axiom`. The only axioms are `propext`, `Classical.choice` and `Quot.sound` (`pipeline_axioms_output.txt`, `multiband_axioms_output.txt`). | Every Lean row below |
| **H(a0, η)**, with a0 = 433148145/10⁸ and η = 10⁻⁹ | Lean hypothesis `TcR1.H κ a0 η` with κ = 2π. It is discharged **[non-Lean]** by the wE1 certificate `/home/claude/tc/wE1` (126 dyadic w-cells on [1,400] plus analytic ends, Arb at 128 bits) = PROOF.md Prop. 6, audited SOUND (`tc/audit_E1/REPORT.md`). wE1's matrix `M+(a,w)` is exactly Lean's `Qe κ a w` with κ = 2π (diagonal `l(2n+1) − (2n+1) − 2S_n`, off-diagonal `l(|n−m|) + l(n+m+1)`, `l(k) = 2aw/(w²+(2πk)²)`). | R-UNIV, R-CELLS (through `cells_stable`), R-MULTI |
| **Prop N**, for the same a0 | Lean hypothesis `MultiBand.PropN tp a0`, i.e. `Qe tp (−a0) w N y ≤ 0`. It is discharged **[non-Lean]** by wMB §4.2 (an analytic argument plus the 1-D Arb check `tc/wMB/cert_neg.py`, verdict PASS), audited SOUND (`tc/audit_MB/REPORT_SUMMARY.md`). | R-MULTI only |
| python-flint **0.9.0** (FLINT **3.6.0**) | `arb` ball arithmetic, `arb.pi()`, `arb_poly` products, `arb.min`/`union`, certain comparisons (`<`, `<=`, `>=` return True only if they hold for every point of the balls), and exact `fmpq` arithmetic. `pyproject.toml` only pins `>=0.6`. | Every Arb/`fmpq` check |
| CPython floats | Exact conversion double → `arb`/`fmpq` (`q_of_float` uses `as_integer_ratio`) and exact `==`/`<` on doubles. | Tiling, T comparisons |
| Hand proofs in `rtsc/certify/PROOFS.md` **[non-Lean]** | Lemmas T0, M, S1, S2, Rw, R1, R2, T1(p), Theorem C ("check_cell ok ⇒ CellOK"), Theorem R, Theorem C2 (cform). | R-CELLS, C-RAYLEIGH |
| Paper proofs **[non-Lean]** | `tc/final/PROOF.md` (general measures, λ(0) = ∞, Tc := sup{unstable T}, Foundations wG Thm 2), `tc/wMB/PROOF.md`. | Scope statements only (§5) |

---

## 1. Dictionary of definitions (Lean ↔ Python)

| Concept | Lean (file:line) | Python | Match? |
|---|---|---|---|
| Spectrum | `TcR1.Spectrum` (Basic.lean:60): Finset of atoms, `p_j ≥ 0`, `Ω_j > 0` | `Spectrum.__post_init__` (spectrum.py:56–74) enforces `omega > 0` finite, `weight ≥ 0` finite, `0 ≤ e < 1`. `ArbSpectrum` (arbcore.py:82–107) drops atoms with `w = 0` | Yes. Dropping zero-weight atoms changes no λ, A or Q value, and their envelope is {0}. `A = 0` raises `ValueError` in Python, so the result is ERROR (conservative; Lean `stable_of_A_zero` would allow a rejection). |
| κ | free real `κ`. All theorems hold for every κ, and H is certified for κ = 2π | `TWO_PI = 2*arb.pi()` (arbcore.py:21–22) | Yes. The same κ = 2π is used in H (wE1), in λ, in g and in the cutoff. |
| g, λ | `g κ T k w = 2w/(w²+(κkT)²)`, `lam = Σ p_j g` (Basic.lean:71–75) | `a_j/(o2_j + (2πkT)²)` with `a_j = 2p_jΩ_j`, `o2_j = Ω_j²` (arbcore.py:147–153). For 2πkT ≥ 2Ω_max the moment series with a rigorous remainder is used (arbcore.py:154–177, **Lemma S1 [non-Lean]**) | Yes |
| λ(0) | Never read: `Kp` uses `L(2n+1)`, `L(adist n m)` with n ≠ m, and `Ssum L n = Σ_{k<n} L(k+1)` | `out[0] = 0` (arbcore.py:180), `hi[0] = lo[0] = 0` (cells.py:91,94). The diagonal sum starts at k = 1 (cells.py:149–151) | Yes (ERRATA F2). Tested in §7. |
| Even form Q⁺ | `Qp = C + Fp(lam) − muTerm` (Basic.lean:82–84). `Kp` diagonal is `L(2n+1) − 2Ssum L n`, off-diagonal is `L(|n−m|)+L(n+m+1)` | C-RAYLEIGH: `Q = Σ_j min(...)·F_j − Σ(2n+1)x_n² − 2μ(Σ_{n<min(N_c,N_x)}x_n)²` with `F_j = Σ_k g(k,Ω_j)γ_k` (rayleigh.py:44–59, 108–130). `γ_k` is the exact regrouping of `Σ_{n,m} Kp_nm x_n x_m` by k | Yes. Literal comparison with Lean `boxMinQ` in §7 (T5). |
| Matrix of Q⁺ at μ* = 0 | `Mq` / `Mhat hi lo` (Common.lean:32–39). Off-diagonal `hi(adist)+hi(n+m+1)`, diagonal `hi(2n+1) − 2Ssum lo n − (2n+1)` | `check_cell` mhat row: `TH + rem − 2v_n S_lo − (2n+1)v_n` (cells.py:152–159), where `TH = P1[n]+P2[M_v−1−n]+P2[M_v+n]` is the exact Σ over m < M_v of the off-diagonal and Hankel-diagonal terms | Yes for m < M_v. The columns m ≥ M_v are covered by Lemmas R1/R2 **[non-Lean]**. Literal comparison in §7 (T6). |
| Envelope | `InEnv e sp sp'` (Common.lean:26): same atoms and Ω, `(1−e)p ≤ p' ≤ (1+e)p` | `rel_envelope` per Spectrum. `hi = (1+e)λ(T1)`, `lo = (1−e)λ(T2)` (cells.py:90–94). `A_hi_q = (1+e)A` exact (arbcore.py:103–107) | Yes. e ≥ 0 is required for the nominal spectrum to lie in its own envelope (`InEnv_refl`), and Python enforces 0 ≤ e < 1. |
| Cutoff | `Nc : WithTop ℕ`. `muTerm` sums `n < N` with `n < Nc`. In `classify` it is an arbitrary function `Ncf : ℝ → WithTop ℕ` | `N_c(T) = #{n ≥ 0 : (2n+1)πT ≤ ω_c}`, certified in Arb (rayleigh.py:30–41). The μ sum runs over `range(min(N_c, N_x))` (rayleigh.py:127–129) | Yes. Python's rule is one instance of `Ncf`. N_c = 0 is allowed in both (the paper's model statement lists N_c ≥ 1; the pipeline flags `COULOMB_INACTIVE`). Tested in §7. |
| Stable / Unstable | `Stable = ¬Unstable(Qp) ∧ ¬Unstable(Qm)`. `Unstable Q` ⇔ ρ ≥ 0 | not computed directly | — |
| Status | `inductive Status` with 6 constructors (Classify.lean:18–25) | `rtsc.status.Status` with 6 members | Bijection (tested). |
| AltEvidence | `structure AltEvidence` with 3 Bools (Classify.lean:32–36) | `@dataclass AltEvidence` (classify.py:58–62), same field names and order | Yes (tested). |
| a0, η | free reals with `0 < a0`, `0 < η` | `A0_FRACTION = 433148145/10⁸`, `A0_ETA_FRACTION = 10⁻⁹` (constants.py:17–19). `A0_Q` exact `fmpq` (search.py:18) | Yes. η > 0 is required by `R_UNIV`/`cells_stable` (they use `stable_of_ge`, which needs `min η 1 > 0`). The certified η is 10⁻⁹ > 0. MultiBand needs only η ≥ 0. |

---

## 2. Decision map, single-band (default `rejection_model = "isotropic"`)

### 2.1 `classifyAlt` and `classify_evidence`: control flow

| Python | Lean | Correspondence |
|---|---|---|
| `classifyAlt` classify.py:65–73 | `PipelineLogic.classifyAlt` Classify.lean:40–44 | Literal: same `if/elif` order, same results. |
| `classify_evidence` classify.py:76–84 | `PipelineLogic.classify` Classify.lean:48–52 | Literal. `len(evs) == 0` ↔ `evs.isEmpty`. `all(...)` ↔ `List.all` (pure, so short-circuiting is irrelevant). `==` on the Enum ↔ `DecidableEq` `==`. |

Verified by `tests/test_correspondence.py`:

- The 5 transcribed Lean snippets still occur verbatim in Classify.lean (drift test).
- The two functions agree with the transcription on all 8 `AltEvidence` values, on all
  37 449 lists of length ≤ 5, and on 3 000 random lists of length 6–60.
- The 5 Lean `example`s hold in Python.
- The Lean theorems `classify_cases`, `classifyAlt_cases`, `status_exclusive`,
  `classify_eq_rejected_iff`, `classify_eq_certified_iff` and `classify_singleton` hold
  for the Python functions on all lists of length ≤ 4.

### 2.2 Empty-list rule (ERRATA F1)

| Python | Lean | Note |
|---|---|---|
| classify.py:77–78 `if len(evs) == 0: return UNDECIDED` | first branch of `classify`. `classify_eq_rejected_iff` and `classify_eq_certified_iff` both carry `evs ≠ []` | Literal. |
| `classify([])` (public entry) classify.py:283–338 | — | Returns UNDECIDED with an empty bracket for both rejection models (tested). |
| pipeline.py:479–480 `certify_spectra` raises on zero spectra, giving ERROR via `compute_material` | ARCH F1 ("UNDECIDED, or ERROR") | Conservative. |

### 2.3 Boolean `univFires` (rule R-UNIV)

**Python:** `evidence_for_alternative` classify.py:136:
`univ = A_hi_q < A0_Q * q_of_float(T*)`, all in exact `fmpq` arithmetic. `A_hi_q = (1+e)·Σ_j q(w_j)`
(arbcore.py:103–107). The same expression appears in `rtsc.certify.check_univ` (`__init__.py`:29–34)
and in the search at classify.py:294, but only the line-136 value reaches the decision.

**Lean semantics:** `UnivCert e a0 T* sp := ∃ Ahi, (1+e)·sp.A ≤ Ahi ∧ Ahi < a0·T*`
(Classify.lean:132). It is consumed by `R_UNIV` (Assembly.lean:56):

| Hypothesis of `R_UNIV` | Discharged by | Rigour |
|---|---|---|
| `hH : H κ a0 η` | wE1 certificate (§0) | Arb, **[non-Lean]**, audited |
| `ha0 : 0 < a0` | constant 433148145/10⁸ | exact |
| `hη : 0 < η` | η = 10⁻⁹ | exact |
| `hAhi : (1+e)·A ≤ Ahi` | `Ahi := A_hi_q = (1+e)A`, equality | exact `fmpq` |
| `hU : Ahi < a0·T*` | classify.py:136, **strict** | exact `fmpq` (tested at exact equality: does not fire) |
| `sp` is a `TcR1.Spectrum` | `Spectrum.__post_init__` validation | exact |

Conclusion used: `Ahi/a0 < T*`, and Stable at every `T ≥ Ahi/a0` for every envelope spectrum,
every μ* ≥ 0 and every cutoff (`R_UNIV` → `univ_stable` → `stable_of_A_le` → `TcR1.stable_of_ge`
or `stable_of_A_zero`).

### 2.4 Boolean `cellsBelowTarget` (rule R-CELLS, with R-MU and R-ODD inside)

**Python:** classify.py:139–187. `cells_below = cells_ok ∧ T_up_cert < T*` (line 185), where
`cells_ok` requires all of D1, `Tiles` and every `check_cell` ok (lines 180–182).

**Lean semantics:** `CellsCert κ e a0 T* sp := ∃ Tup tend Ahi cells, Tup < T* ∧ Tiles Tup tend cells ∧
(∀ c ∈ cells, CellOK κ sp e c.1 c.2) ∧ (1+e)A ≤ Ahi ∧ Ahi ≤ a0·tend` (Classify.lean:137–139). It is
consumed by `cells_stable` / `R_CELLS` (Assembly.lean:94, 109).

| Hypothesis | Discharged by (Python) | Rigour |
|---|---|---|
| `Tup < T*` | classify.py:185 `T_up_cert < T_star` (doubles, **strict**, ERRATA F4) | exact (tested at equality) |
| `Tiles Tup tend cells` (Cover.lean:10). Nonempty, `a₀ = Tup`, `bᵢ = aᵢ₊₁`, `b_K = tend` | classify.py:153–164. Empty list ⇒ no evidence, matching `Tiles _ _ [] = False`. `ends[0][0]==T_up`, contiguity and `ends[-1][1]==t_end` with exact double equality | exact (gap case tested) |
| `(1+e)A ≤ Ahi` | `Ahi := A_hi_q`, equality | exact |
| `Ahi ≤ a0·tend` (D1, **non-strict**) | classify.py:149 `A0_Q*q(t_end) >= A_hi_q and t_end > 0` | exact (tested at equality) |
| `H`, `0<a0`, `0<η` | as in 2.3 | — |
| `CellOK κ sp e T1 T2` (Assembly.lean:73), per cell: | `check_cell` (cells.py:54–200), mhat variant | see the next block |
| … `0 ≤ T1` | cells.py:73 `T1 > 0 and T1 <= T2` (stronger) | exact |
| … `∃ η > 0` | cells.py:79 `eta > 0` (double) | exact |
| … `∃ v, ∀ n, 0 < v n` | head `u_m > 0` finite (cells.py:77). Tail `v_m = C(2m+1)^{−p}` with `C > 0`, `p ∈ {2,3}` (cells.py:79–82, 99) | exact. Holds for **every** n. |
| … `hi = (1+e)·lam κ sp T1`, `lo = (1−e)·lam κ sp T2` | cells.py:89–94. λ from direct atom sums, or the series of **Lemma S1 [non-Lean]** with remainder | Arb |
| … `∀ N, ∀ n < N, rowSum (Mhat hi lo) v N n ≤ −η v n` | Head rows n < N_e (cells.py:148–170): `TH` (exact Σ over m < M_v) + `rem` (min of **Lemma R1** and **Lemma R2 [non-Lean]**, using **Lemma T0 [non-Lean]** `hi(k) ≤ b/k²`) − `2v_n S_lo(n)` − `(2n+1)v_n` + `η v_n` < 0 in Arb. Tail rows n ≥ N_e (cells.py:171–197): **Lemma T1(p) [non-Lean]**, i.e. `D(N_e) > 0` and `G(N_e) < 1` in Arb, plus monotonicity of G in n (hand proof). Passing from the certified infinite-row bound to every truncation N uses `rowSum_mono` / `rows_of_eventually` (**Lean**, Schur.lean:88, 100), plus the elementary fact that a nondecreasing partial sum is ≤ any bound on its tail-completed value (**[non-Lean]**, PROOFS.md Thm C) | Arb + hand proofs. Literal Lean `Mhat`/`rowSum` cross-check for all N ≤ 400 in §7 (T6). |

Lean chain used: `cell_stable` (Assembly.lean:81) ← `envelope_form_bound` (Schur.lean:135) ←
`schur_row_test_inf`, `rowSum_le_envelope`, `Mq_le_Mhat`, `form_eq_Mq`; `lam_box` (Monotone.lean:73,
which is g antitone in T plus the weight envelope); `stable_of_even_bound` (Schur.lean:146), which is
**R-MU** (Lemma 3 via `L3`) and **R-ODD** (`TcR1.L6`). `cover` (Cover.lean:17) handles T in a cell.
`univ_stable` handles T > t_end.

Bracket value: `info["T_up"] = T_up` when cells_ok. Otherwise it is `t_end` if D1 holds (valid
upper bound by `univ_stable`). Otherwise it is None (classify.py:151–152, 180–181).

The **cform** variant (`row_check = "cform"`) does not establish Lean `CellOK`. It rests on
**PROOFS.md Thm C2 [non-Lean]**. See gap G9.

### 2.5 R-MU and R-ODD

No Python boolean. They are used *inside* the R-CELLS proof (`stable_of_even_bound`, `R_MU`,
`TcR1.L6`) and inside R-UNIV (`TcR1.main_odd`). `rule_ids` lists them whenever a valid cell
cover exists (classify.py:356–357), which is a label only.

### 2.6 Boolean `rayleighAtOrAboveTarget` (rule C-RAYLEIGH)

**Python:** classify.py:188–209. `ray_above = ray_ok ∧ T_low ≥ T*` (line 209). `ray_ok` requires the
certificate's `mu_hi` and `omega_c` to equal the run's values exactly (line 194) and
`check_rayleigh` ok (rayleigh.py:62–136).

**Lean semantics:** `RayleighCert κ e μhi T* Ncf sp := ∃ T, T* ≤ T ∧ ∃ N x, (∃ n<N, x n ≠ 0) ∧
0 ≤ boxMinQ κ sp e μhi (Ncf T) T N x` (Classify.lean:143–146). It is consumed by `rayleigh_cert`
(Rayleigh.lean:137).

| Hypothesis | Discharged by | Rigour |
|---|---|---|
| `T* ≤ T` (**non-strict**) | classify.py:209 `T_low >= T_star` | exact (tested at equality) |
| `x ≠ 0` | rayleigh.py:74 `any(v != 0)` on doubles | exact |
| `Ncf T` = the run's cutoff rule | rayleigh.py:76–78 `certified_n_c(T, ω_c) == N_c claim`, with both defining inequalities certified in Arb. classify.py:194 checks ω_c equals the run's `omega_c` (one ω_c for all alternatives: `resolved_omega_c` = factor·max over **all** spectra, classify.py:102–107) | Arb (tested against the literal definition) |
| `0 ≤ boxMinQ` | rayleigh.py:124–131. The per-atom `min((1−e)p_jF_j, (1+e)p_jF_j)` is Lean `sel_mul_eq_min` (ERRATA F6). `F_j = Σ_k g γ_k` is an exact regrouping of `atomF = Fp(g(·,Ω_j))`. For k ≥ k_s, g is evaluated by **Lemma S2 [non-Lean]** (alternating series with a certified remainder) | Arb (literal Lean `boxMinQ` cross-check in §7 (T5)) |
| μ ≤ μ_hi for every μ in the claim | `Qp_anti_mu` (Lean) | — |

Conclusion used: for every envelope spectrum and every μ ≤ μ_hi, `Unstable(Q⁺ at T, cutoff N_c(T))`,
hence `¬Stable(T)` (`not_stable_of_unstable`). This holds **only** for the configured block-cutoff
rule (ERRATA F3, D2). No H is needed.

### 2.7 Combination over alternatives (S5)

| Python | Lean |
|---|---|
| classify.py:330–338: one `AltEvidence` per spectrum, then `classify_evidence` | `classify_rejected_sound` / `classify_rejected_sound'` (Classify.lean:214, 228) and `classify_certified_sound` (Classify.lean:244), with `EvidenceSound` (Classify.lean:149) for each alternative |

REJECTED_PROVEN ⇒ (with H, a0 > 0, η > 0) the list is nonempty, and every alternative has some
`T_b < T*` such that every envelope spectrum is Stable for all `T ≥ T_b`, all μ* ≥ 0 and all cutoffs.
CERTIFIED_CANDIDATE ⇒ the list is nonempty, and every alternative has some `T_low ≥ T*` such that
every envelope spectrum is unstable for all μ ≤ μ_hi at cutoff `Ncf(T_low)`.
`evidence_consistent` (Classify.lean:193) shows that, with e ≥ 0 and μ_hi ≥ 0, sound evidence
never has both a rejection rule and C-RAYLEIGH true, so the rule order cannot matter.

Mismatch: the Lean list theorems use **one** `e` for all alternatives. Python allows a different
`rel_envelope` per spectrum (`pipeline._with_envelope` takes max(own, config), pipeline.py:382–384).
The per-alternative theorems `classifyAlt_*_sound` apply unchanged, but the list statement needs
a one-line generalisation (G11).

Evidence provenance: in `classify` the booleans of search-produced certificates use the
`verified` set (classify.py:324–329). This is a content hash of a certificate that already passed
`check_cell`/`check_rayleigh` on the same `ArbSpectrum`. `recheck` (classify.py:419–475) recomputes
every check with `verified=None`.

### 2.8 Bracket combination (ERRATA F7)

`_combine_bracket` classify.py:218–224 computes `[min T_low_i, max T_up_i]`, and None if any
alternative lacks the bound. **No Lean theorem.** The intended statement is: if every alternative
satisfies `T_low_i` unstable (cutoff `Ncf(T_low_i)`, μ ≤ μ_hi) and Stable for all `T ≥ T_up_i`
(every μ* ≥ 0, every cutoff), then every alternative has an unstable point ≥ min T_low and is
stable above max T_up. This is a trivial min/max lemma, but Lean does not define `Tc_model := sup{T :
¬Stable T}` at all (G10). The two ends carry different quantifiers: the lower end is for
μ* ∈ [0, μ_hi] at the block cutoff, the upper end is for every μ* ≥ 0 and every cutoff. The
min/max/None semantics are tested exhaustively for up to 3 alternatives.

### 2.9 Status assignment outside `rtsc.certify`

| Python | Guarantee |
|---|---|
| `pipeline._pipeline_verdict` pipeline.py:262–272 | Only DEFERRED or ERROR (assertion). EXCLUDED_ASSUMPTION comes from `rtsc.stability` (S3), never from Lean-backed rules. |
| `pipeline._finish` pipeline.py:275–277 and the run-summary check pipeline.py:735–736 | REJECTED_PROVEN or CERTIFIED_CANDIDATE only with origin `rtsc.certify.classify` (S1). |
| `certify_spectra` pipeline.py:471–492, `annotate_model`, `ensure_cutoff_assumption`, `validity.apply_flags` | Add assumptions and flags only. Each asserts that status and bracket are unchanged. |

These are runtime assertions with no Lean counterpart. Lean `classify_cases` shows that
`classify` itself never emits EXCLUDED, DEFERRED or ERROR.

---

## 3. Decision map, multiband (`rejection_model = "multiband"`, opt-in)

### 3.1 R-MULTI, one band-coupling alternative

**Python:** `rmulti_evidence` multiband.py:275–318 re-derives `fires` from a stored block. It
checks:

- `T* > 0` finite (278–279);
- the digest (281–282);
- the normalisation refusal (283–286);
- `check_route` (170–247);
- `κ_hi < a0·T*` in exact `fmpq`, **strict** (295).

`T_bound` = κ_hi/a0, rounded up to a double (297–301).

| Route | Python check | Lean theorem | Hypotheses not checked by Python (recorded as assumptions) |
|---|---|---|---|
| `rowsum`, source A_K | `(1+e)·Σ_j A_ij ≤ κ_hi` for every row, exact `fmpq` (187–204) | `MultiBand.R_UNIV_MB_rowsum` (Pipeline.lean:89), or `R_UNIV_MB` with v = 1 and `sub_env` | an `MSpec` exists, i.e. symmetric `W_s` = detailed balance for **some** N > 0; `CoulPSD`; the true 𝒜' ≤ (1+e)𝒜 entrywise |
| `rowsum`, source `row_sum_upper_K` | `(1+e)·R_i ≤ κ_hi` exact | the same | additionally: the declared R_i bound the true row sums (data assumption "DATA BOUND/CAVEAT") |
| `cw` | `v > 0` doubles and `(1+e)(A v)_i ≤ κ_hi v_i` exact `fmpq` (205–221) | `R_UNIV_MB` (Pipeline.lean:55) with `sub_env` (Pipeline.lean:20) | detailed balance for some N; `CoulPSD`; envelope |
| `eig` | N_i > 0, u > 0. `W = N_i A_ij` exact, symmetrised to W̄ if not exactly symmetric, `r_i = √N_i` in Arb, `(1+e)(K u)_i ≤ κ_hi u_i` in Arb (222–246) | `sub_of_K` (Pipeline.lean:36) ⇒ `Ā v ≤ (κ_hi/(1+e)) v` with v = u/r, then `sub_env` and `R_UNIV_MB` | W symmetric, or else the verdict is about W̄; the envelope covers 𝒜, not N; `CoulPSD` |
| all | `0 ≤ e < 1` (173–174), `κ_hi > 0` | `he : 0 ≤ 1+e`, `hT0 : 0 < T*`, `H` (η ≥ 0 suffices), `PropN`, `0 < a0` | H and Prop N from the certificates in §0 |

Bracket: `T_up = max_i T_bound_i`, `T_low = None` (classify.py:390–392). This is justified by
`mstable_of_gt` (Main.lean:176) / `unstable_imp_le` (Main.lean:193): stable for every
T > κ_hi/a0. `R_UNIV_MB` itself only states T ≥ T*.

### 3.2 `combine_multiband` ↔ Lean `PipelineLogic.combineMultiband` / `classifyFull`

| Python | Lean (`tc/lean/PipelineLogic/Multiband.lean`) | Correspondence |
|---|---|---|
| `combine_multiband(status_sb, mb_fires, rejection_model, n_spectra=None)` classify.py:262–277 | `combineMultiband sb mbFires model (nSpectra : Option ℕ := none)` | Literal. `n_spectra is not None and n_spectra == 0` ↔ `nSpectra == some 0` (→ UNDECIDED, checked before `mb_rej`; `None` ↦ `none`). `mb_rej = len>0 and all(...)` ↔ `!mbFires.isEmpty && mbFires.all id`. The two `if` branches are in the same order. The string `"isotropic"` ↦ `.isotropic`, and every other string ↦ `.multiband`, as stated in the Lean header. |
| `status = combine_multiband(status_sb, mb_fires, model, n_spectra=len(spectra))` classify.py:413 and :489 (recheck) | `classifyFull evs mbFires model := combineMultiband (classify evs) mbFires model (some evs.length)` | Literal: one `AltEvidence` per spectrum, so `len(spectra) = evs.length`. The assumption `BAND-DATA-NOT-PER-ALTERNATIVE` (band-coupling count ≠ spectrum count) is an annotation only. The decision still needs every coupling to fire, which Lean requires anyway. |
| `rmulti_evidence(...)[0]`, one boolean per band alternative (classify.py:227–239) | `MBAlt.fires`, with semantics `RMultiCert a0 T* a` = `0<T* ∧ 0≤e ∧ ∀i, 0<N_i ∧ ∃κhi<a0T*, CWCore` | The route lemmas `route_cw`, `route_rowsum`, `route_rowsum_declared` and `route_eig` turn each Python route check (§3.1) into `CWCore`. The rowsum `R_i` bound is a data assumption, and detailed balance is part of the model (`MSpec.W_symm`). |

Truth table (tested exhaustively against the Lean transcription). `sb` is the single-band
status and `mb_rej := mb_fires ≠ [] ∧ all fire`:

| model | sb | mb_rej | result |
|---|---|---|---|
| isotropic | any | any | sb (any n_spectra) |
| multiband, n_spectra = 0 | any | any | UNDECIDED |
| multiband | ≠ CERTIFIED | true | REJECTED |
| multiband | REJECTED / UNDECIDED / other | false | UNDECIDED |
| multiband | CERTIFIED | true | UNDECIDED (the two models conflict) |
| multiband | CERTIFIED | false | CERTIFIED (single-band model, `M-ISO-SINGLEBAND`) |

(Rows after the second assume n_spectra ≠ 0, i.e. `None` or ≥ 1.)

Lean theorems that back it:

- `combineMultiband_isotropic`, `combine_mb_rejected_iff`, `combine_mb_certified_iff`. Both iff
  theorems now carry the conjunct `n ≠ some 0`.
- `combine_mb_zero_spectra` (multiband, `n_spectra = 0` ⇒ UNDECIDED),
  `classifyFull_zero_spectra` and `classifyFull_zero_spectra_not_rejected`. With zero spectra
  the full decision is UNDECIDED in both modes, whatever the band evidence, so it is never
  REJECTED.
- `combine_mb_cases` and `classifyFull_cases` (only the three statuses are possible), and
  `classifyFull_exclusive`.
- `combine_mb_empty_not_rejected` and `classifyFull_empty_not_rejected` (an empty band list
  never rejects).
- `combine_mb_rejected_needs_mb`.
- `combine_certified_imp`.
- `classifyFull_mb_rejected_sound` and `classifyFull_mb_rejected_sound'`. Their hypotheses are
  `H(a0,η)` with η ≥ 0, `PropN`, and `RMultiCert` for every True. Their conclusion: the spectrum
  list and the band list are both nonempty, and for every band alternative there is `κ_hi/a0 < T*` such that every
  envelope member, every PSD `N_iμ*_ij`, every cutoff and every `T > κ_hi/a0` give `MStable`.
  This backs the reported `T_up = κ_hi/a0`, rounded up.
- `classifyFull_iso_rejected_sound` (the isotropic-mode REJECTED is the single-band REJECTED).
- `classifyFull_certified_sound` (CERTIFIED in either mode is a statement about the single-band
  model only).
- `certified_transfers_to_multiband`. This is a conditional lift whose hypotheses (the DOS
  average lies in the envelope, and `μ*_iso ≤ μ_hi`) the pipeline does **not** check.

Per CORRESPONDENCE_LEAN.md, there is no `sorry` or `axiom`. The axiom output
(`pipeline_axioms_output.txt`) lists these theorems with only `[propext, Classical.choice, Quot.sound]`.
After the n_spectra change, Lean was rebuilt (`nice lake build`, 1057 jobs, OK) and
`CheckPipeline.lean` was re-run. It covers 82 theorems, each depending on a subset of
`[propext, Classical.choice, Quot.sound]`.

What remains outside Lean on this path:

- **Prop N** and **H** (certificates, §0).
- Detailed balance is not checked for the `rowsum`/`cw` routes (G3). Lean's `MSpec.W_symm`
  needs it per frequency, which is stronger than the integrated `N_iA_ij = N_jA_ji`
  (CORRESPONDENCE_LEAN MB-N2).
- B = 0 makes the Lean conclusion vacuous (G12).
- With `spectra = []` and firing band data, the public `classify` returns REJECTED. This is
  sound, because the claim is about the band data only (MB-N3). The pipeline refuses empty
  spectra anyway.

---

## 4. Mismatch audit (Lean definition vs Python computation)

| Topic | Finding |
|---|---|
| **Q form vs matrix entries** | Lean's `Qp`/`Kp` is the cancelled (reduced) form. `TcR1.Qp_spec_form` shows it equals SPEC's form for any finite λ(0). R-CELLS uses Lean `Mhat` entry by entry (F2: Σ from k = 1). C-RAYLEIGH uses `γ_k`, the exact regrouping of `Σ Kp x x`. Both are checked against literal Lean definitions in §7 (T5, T6). **Match.** |
| **κ = 2π** | The Lean theorems hold for every κ. H is certified with 2π. Python uses `2*arb.pi()` everywhere (λ, g, N_c). **Match**, provided the reader instantiates κ := 2π (the Lean statements never fix it). |
| **Cutoff semantics** | Lean filters `n < Nc`. Python sums `n < min(N_c, N_x)` with N_c certified from its definition. REJECTED holds for every cutoff (`Nc` universally quantified). CERTIFIED holds for `Ncf(T_low)` only, and Python records ω_c and the rule (classify.py:366–368). One ω_c serves all alternatives, matching Lean's single `Ncf`. Pipeline wording mismatch: `cutoff_description` (pipeline.py:528–538) says "per spectrum", but it is only used when no lower bound is claimed (G15). |
| **Finite atomic spectra vs the parsed spectrum** | Lean covers exactly the finite atomic measure that Python passes. `Spectrum.from_grid` (spectrum.py:127–173) turns the DFPT density into midpoint atoms. It drops mass at ω ≤ 0, clips negative α²F, and converts units in floating point. That map is an input convention outside every proof (recorded as `DISCRETISATION(...)` in the assumptions). See G2. |
| **η > 0** | Required by `R_UNIV`, `cells_stable` and `classify_rejected_sound` (`hη : 0 < η`). Certified η = 10⁻⁹. Each cell carries its own η > 0 (a separate quantity, `CellOK`'s ∃η). MultiBand needs η ≥ 0. **Match.** |
| **Strictness** | R-UNIV `<` ✓; D1 `≤` ✓; R-CELLS `T_up < T*` ✓; C-RAYLEIGH `T ≥ T*` ✓ and `boxMinQ ≥ 0` ✓; the row test certifies `row + ηv_n < 0`, which is stronger than Lean's `≤` ✓; R-MULTI `<` ✓. All boundaries are tested at exact equality (§7 T4), except R-MULTI's, which is read off multiband.py:295. |
| **Envelope** | `0 ≤ e < 1` is enforced (spectrum.py:72–74). e ≥ 0 matters: for e < 0, Lean's envelope does not contain the nominal spectrum and the soundness theorems would be vacuous. Per-alternative e: see G11. |
| **T1 > 0, T1 ≤ T2** | Stronger than Lean's `0 ≤ T1` (no ordering needed). ✓ |
| **Zero-weight atoms / A = 0** | Dropped or ERROR. Conservative. ✓ |

---

## 5. Scope that the paper covers but Lean does not

- **General (non-atomic) measures.** PROOF.md §1 takes α²F to be any finite positive Borel measure.
  The reduction to Einstein spectra (Lemma 2) is proved there for arbitrary measures (audited in
  `audit_B`). Lean `TcR1.Spectrum` and `MultiBand.MSpec` are **finite atomic** only (MB-D4). The
  pipeline only ever passes atomic measures, so the Lean coverage is complete for what Python
  decides. The step from the physical (continuous) α²F to the parsed atoms is G2.
- **λ(0) = ∞** (∫α²F/ω = ∞). PROOF.md defines the model by the reduced form in that case (wG G5.1).
  Lean never meets it, because atomic spectra with Ω > 0 have finite λ(0), and `Qp_spec_form`
  covers finite λ(0) only.
- **Tc := sup{unstable T}** and its relation to the linearised gap equation (Foundations, wG
  Thm 2), the operator domain, and N_c = ∞ via a limit are in the paper only. Lean states
  pointwise stability and instability in T.
- **Multiband:** Theorem 2 (Λ_A), sharpness §4.6, the ordering §4.7, the continuum anisotropic
  limit, Perron–Frobenius, and Prop N itself (a hypothesis in Lean) are not in Lean (MB-D4).

---

## 6. Gaps, ranked by severity

Severity means how much of a REJECTED or CERTIFIED verdict rests on something that is not a
Lean theorem or a certified Arb check, and whether the gap can be reached in production.

| # | Severity | Gap | Effect / remedy |
|---|---|---|---|
| **G1** | ~~High~~ → **Closed (logic)** | `combine_multiband` now has a Lean mirror (`PipelineLogic.combineMultiband` / `classifyFull`, Multiband.lean), and the Python function is literally equal to it (tested exhaustively, §7 T3). Still open on this path: Prop N and H are certificates (G5), detailed balance is unchecked (G3), B = 0 is vacuous (G12), and the bracket is unlabelled (G7). | — |
| **G2** | **High** (scope) | The parsed atomic spectrum ≠ the physical α²F. Discretisation, ω ≤ 0 drop, negative clipping and unit rounding are outside every proof. | Every verdict is about the parsed measure (labelled `DISCRETISATION`). For R-UNIV the paper's general-measure theorem would transfer the rejection to the continuous α²F if `∫α²F ≤ A_hi` were certified. Python does not certify that, because midpoint mass ≠ integral and dropped ω ≤ 0 mass is not counted. Remedy: a certified quadrature bound into `rel_envelope`, or state the rejection for the atomic model only (current behaviour). |
| **G3** | **Medium** | R-MULTI `rowsum`/`cw` routes: the hypothesis "detailed balance for some N > 0" (Lean `MSpec.W_symm`) is assumed, not checked, although it can be checked from A_K (the Kolmogorov cycle condition). `CoulPSD` is assumed (no μ*_ij data). | Recorded in the assumptions. Remedy: check the cycle condition exactly when A_K is given. |
| **G4** | **Medium** | R-CELLS rests on hand proofs **[non-Lean]**: Lemma S1 (λ series), T0 (`hi ≤ b/k²`), R1/R2 (columns m ≥ M_v), T1(p) with G monotone (rows n ≥ N_e), and Theorem C's passage from the certified infinite row sum to Lean's ∀N truncation hypothesis. | Every R-CELLS rejection depends on them. They are cross-checked here on one certificate (T6, every N ≤ 400 against literal Lean `Mhat`). Remedy: formalise T0/R1/R2/T1 (pure real inequalities), or stop using R-CELLS for decisions beyond R-UNIV. |
| **G5** | **Medium** (trusted base) | H(a0, η) and Prop N are Lean hypotheses discharged by Python/Arb certificates (wE1, wMB `cert_neg.py`), audited but not machine-checked. python-flint 0.9.0 / FLINT 3.6.0 is trusted, and `pyproject` pins only `>=0.6`. | Remedy: pin `python-flint==0.9.0` and record the FLINT version in the certificate. |
| **G6** | **Medium** | C-RAYLEIGH: the `γ_k` regrouping (exact algebra) and Lemma S2 (per-atom series for g with a remainder) are hand proofs. | Cross-checked against literal Lean `boxMinQ` on 14 cases (both branches, ok = True and ok = False, e ∈ {0, 0.05, 0.1, 0.3}). |
| **G7** | **Medium-low** (labelling) | Multiband mode, verdict not REJECTED: `T_low`/`T_up` are the **single-band** bracket with no `M-ISO-SINGLEBAND` label and no `MULTIBAND_RISK` flag. Reproduced: a weak spectrum gives UNDECIDED with `T_up = 12.12 K < T* = 300 K` and no single-band label. | The status is conservative, but the reported bound is a single-band statement presented unlabelled. Remedy: label the bracket, or clear `T_up` in that case. |
| **G8** | **Low-medium** (labelling) | `evidence_for_alternative` accepts `cform` cells whatever `params.row_check` is (classify.py:168), but the `R-CELLS[cform]` assumption depends on `params.row_check` (classify.py:371–372). `recheck` does not compare assumptions in the single-band case. | Only reachable with an edited or mixed certificate. Remedy: derive the assumption from the cell variants actually present, or reject cform cells when `row_check="mhat"`. |
| **G9** | **Low-medium** | The `cform` variant (opt-in) is not formalised (PROOFS.md Thm C2). | It is labelled when chosen through the params (see G8). |
| **G10** | **Low** | There is no Lean definition of `Tc_model` and no lemma for the bracket combination `[min T_low, max T_up]`. | It is a trivial min/max argument, tested here. |
| **G11** | **Low** | The Lean list theorems use a single `e` (and `Ncf`, `μhi`) for all alternatives, while Python allows a different `rel_envelope` per spectrum. | The per-alternative theorems apply unchanged. Remedy: restate `classify_*_sound` with `e : alts → ℝ`. |
| **G12** | **Low** | `rmulti_evidence` / `recheck` accept an empty band matrix (B = 0): every route check passes vacuously and R-MULTI "fires" (reproduced). In Lean the conclusion over `Fin 0` is vacuous. `rmulti_for` raises on B = 0, so production cannot reach this, but `recheck` with supplied empty band data can. | Remedy: require B ≥ 1. |
| **G13** | **Low** | R-MULTI `T_up = κ_hi/a0` rests on `mstable_of_gt` / `unstable_imp_le` (strict T > κ_hi/a0, MB-D2), not on `R_UNIV_MB`. | Sound. Documented here. |
| **G14** | **Info** | `recheck` does not validate `rejection_model`: any string other than `"isotropic"` gets multiband semantics. | This matches Lean's string ↦ `.multiband` map, which is the conservative branch. `classify` validates it. |
| **G15** | **Info** | `pipeline.cutoff_description` says "ω_c per spectrum", but `classify` uses one ω_c = factor·max over all spectra. It is shown only when no lower bound is claimed. | Wording only. |
| **G16** | **Info** | The paper's model lists N_c ≥ 1, but Python can use N_c = 0 (πT > ω_c). | Lean covers N_c = 0, and the pipeline flags `COULOMB_INACTIVE`. |

No mismatch was found between the Lean `classify`/`classifyAlt`/`combineMultiband`/`classifyFull`
and the Python control flow.
None was found in the strictness of any rule, and none in the Q / Mhat / boxMinQ / cutoff
formulas.

---

## 7. Tests (`tests/test_correspondence.py`)

| Id | Test | Checks |
|---|---|---|
| T1 | `test_lean_source_verbatim` (5 cases), `test_status_constructors_bijective`, `test_altevidence_fields_match_lean` | The transcribed Lean text still occurs verbatim in Classify.lean. Status and AltEvidence correspond. |
| T2 | `test_classifyAlt_equals_lean` (8), `test_classify_equals_lean_exhaustive_len_le_5` (37 449 lists), `..._random_long_lists` (3 000), `test_lean_examples_hold_in_python`, `test_lean_classify_theorems_hold_in_python`, `test_top_level_classify_empty_list_is_undecided` | The control flow is literally Lean `classifyAlt`/`classify`. The Lean theorems about `classify` hold. F1 holds at the public entry point. |
| T3 | `test_lean_multiband_source_verbatim` (3), `test_lean_multiband_examples_hold_in_python`, `test_combine_multiband_equals_lean_exhaustive` (6 statuses × all fire lists of length ≤ 6 × 3 model strings, including an unknown one, × n_spectra ∈ {None, 0, 1, 2, 5}), `test_classifyFull_equals_lean` (all evidence lists of length ≤ 3 × fire lists of length ≤ 3 × 2 models), `test_combine_multiband_isotropic_is_identity`, `test_zero_spectra_never_rejected`, `..._characterisation` (with n_spectra), `..._soundness_shape` | `combine_multiband` and the composed decision are literally Lean `combineMultiband`/`classifyFull`. The Lean characterisations `combine_mb_rejected_iff`, `combine_mb_certified_iff`, `combine_mb_cases` and `combine_certified_imp` hold in Python. |
| T4 | `test_combine_bracket_F7`, `test_R_UNIV_strict`, `test_D1_nonstrict_and_no_cells`, `test_R_CELLS_strict`, `test_C_RAYLEIGH_nonstrict` | F7 semantics, and every boolean's strictness at exact boundaries: `A_hi = a0·T*` exactly via `a0·5⁸ = 433148145/2⁸`; `T_up = T*`; `T_low = T*`; a tiling gap; a μ_hi or ω_c mismatch. |
| T5 | `test_certified_n_c_matches_definition` (300), `test_check_rayleigh_Q_equals_lean_boxMinQ` (12 seeds), `test_check_rayleigh_accepting_case_equals_lean_boxMinQ` (2) | `check_rayleigh`'s Q ball overlaps Lean `boxMinQ` evaluated literally from `C`, `muTerm`, `Kp`/`Fp`, `g`, `sel`, and the decision agrees. The cases cover the direct and series branches, zero-weight atoms, e > 0, and N_c from 0 to 517. |
| T6 | `test_cell_certificate_satisfies_lean_CellOK_rows` | For an accepted mhat certificate: `rowSum (Mhat hi lo) v N n ≤ −η v_n` for N ∈ {N_e, 400} and every n < N (so for every N ≤ 400, by `rowSum_mono`), with λ from direct atom sums. Rows n ≥ N_e exercise Lemma T1. |

Result (python 3.13, python-flint 0.9.0, `nice`, 2 threads): **50 passed in 3.5 s.**

Mutation check (in memory, package files untouched). The tests catch each of these:

- an empty list mapped to REJECTED;
- "any alternative rejects" instead of "all alternatives reject";
- a corrupted `γ_k` (it differs from Lean's `boxMinQ` in every seed with N_x ≥ 2).
