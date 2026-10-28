# Adversarial review, round 2: α²F convergence protocol v2 (EPH_CONVERGENCE_PROTOCOL.md)

Reviewer: independent referee, 2026-10-27. Nothing was changed.

**Read:** HANDOFF §0 and §3; ARCH §2–4c; DISCOVERY §3; the v1 review; protocol v2; `epw76_interface_spec.md`.

**EPW source checked:**
- `EPW_src_selfen.f90`: `selfen_phon_q`, lines 889–1405.
- `EPW_src_io_io_selfen.f90`: `selfen_ph_write`.
- `EPW_src_use_wannier.f90`: E_F, restart, qwindow, allocation, write order.
- `EPW_src_utilities_bzgrid.f90`: q-grid generation and weights.

**rtsc source checked:**
- `pipeline.py`: `blocking_tags`, `apply_input_policy`, `_with_envelope`, omega_c handling.
- `certify/classify.py`: `CertifyParams.resolved_omega_c`.
- `discovery/loop.py`: `_reconcile`, `_final_under_current_config`, `_record`.
- Lean `DiscoveryLogic/Ledger.lean`: the guard.

Severity uses the v1 scale:
- **CRITICAL:** can yield a wrong proved verdict, meaning the recorded E-ENV(e) and alternative set are false for the converged DFT α²F in the direction that matters.
- **MAJOR:** the result is unreliable or vacuous, or a contract or plumbing violation.
- **MINOR:** everything else.

---

## 0. What is verified correct in the source (so it is not re-litigated)

**Normalisation of P1 is right.** In `selfen_phon_q` with `delta_approx=.true.`:
- weight = π ω_q w_k δ̃(ε_k) δ̃(ε_k+q). Here δ̃ = `w0gauss(x,0)/σ` is always a Gaussian, so the weight is ≥ 0.
- g2 = |epf17|²/(2ω).
- γ_qν = Σ weight·g2.
- λ_qν = γ_qν / (π ω_qν² N_F).
- N_F = `dos_ef(ngaussw, σ, ef0, all bands, full fine k grid)/2` (per spin). `wkf` carries the spin factor.

So λ_qν is EPW's standard double-delta λ, and λ = Σ_q w_q Σ_ν λ_qν. It is the same quantity `a2f_main` sums ("Summed el-ph coupling"). With w_q = 1/N_q:
- α²F = Σ (1/N_q) λ_qν ω_qν/2 · δ(ω − ω_qν) gives 2∫α²F/ω = λ exactly.
- A = Σ λ_qν ω_qν/(2N_q).

There is no factor-2 error, provided ω comes from the linewidth file. Its last column is Γ = 2γ (FWHM), and it must never be used to rebuild λ.

**Degenerate modes.** λ is averaged over modes with |Δω| < 1e-6 Ry. The average preserves Σλ, and Σλω to about 1e-6 relative, within each degenerate set. It is harmless.

**Temperature.** With `delta_approx`, λ does not depend on `temps`. The Fermi–Dirac factors are computed but unused.

**Fermi window.** EPW drops a (k, k+q) pair only if no band at k (or at k+q) lies within fsthick of E_F. With fsthick ≥ 8σ, every dropped term carries a factor ≤ exp(−64). The C4 fsthick hole is closed.

**One smearing per file.** The file holds **smearing index 1 only**: `lambda_all(imode, iqq, 1, itemp)`. With `nsmear=1`, ef0 = efnew, the fine-grid E_F at degaussw with `ngaussw`. So each σ needs its own EPW run.

---

## CRITICAL

### C1. The q direction has no convergence control (v1 C2 is fixed for k only)
**Scenario.**
- Q_1 = Q_c is the coarse DFPT grid itself. Interpolation is the identity there, so it is the least informative grid.
- Q_2 = 2·Q_c. These are the only two q grids.
- The only q term in e is 3·|A(Q_2) − A(Q_1)|/A.

Example: λ(Q_1) = 1.20 and λ(Q_2) = 1.00 with A equal to 2 %. This is plausible: near-Γ acoustic and Kohn-anomaly weight is under-sampled on a coarse grid, and λ ∝ Σw/ω is the quantity that moves while A barely does.
- The protocol gives e = 0.10.
- If λ converges in q with p = 1, then λ_∞ ≈ λ(Q_2) − (λ(Q_1) − λ(Q_2)) = 0.80. That is 20 % below the weakest alternative, outside its 10 % envelope.
- C-RAYLEIGH at the envelope minimum certifies a coupling 12 % above the truth. The lower bound is the tight one, so the result is a wrong CERTIFIED_CANDIDATE.
- For rejection the mirror case applies: λ and A rise with q, and R-CELLS/R-UNIV use an A_hi that is too low.

**Why.**
- v2 applies Richardson rates, windowed weights and the three-grid monotonicity test to k only.
- In q there are two points, so no rate can be estimated. λ, ω_log and W_b are not measured, and nothing is admitted or refused.
- The factor 3 on |ΔA| covers p ≥ 0.42 for refinement ratio 2, but only for A.
- F5 says frequency-position uncertainty must enter as alternatives. Alternatives at Q_1 and Q_2 cover the Q_1→Q_2 change, not Q_2→∞.

**Remedy.**
- Treat q exactly like k: three nested q grids, with Q_1 > Q_c so that at least one interpolated grid is in the chain. For example 2Q_c, 4Q_c, 8Q_c at fixed K.
- Apply Richardson to λ, A, ω_log and W_b, with the same refusal rules (see M1).
- Put R_q into e.
- If cost forbids this, make q non-convergence an explicit verdict-blocking assumption rather than a 3·|ΔA| term.

### C2. G6 can never fail (circular), and the coarse-grid / g-interpolation bias stays uncovered
**Scenario.**
- P6: e = max(0.10, 3·max(…, G6 gap)) ≥ 3·gap.
- P5: G6 passes iff gap ≤ e/2.
- Since gap ≤ e/3 < e/2, G6 always passes unless e > 0.5, which is the separate envelope error.

G6 is therefore not a gate, only an envelope term. As an envelope term its factor is unjustified:
- Q_c′ = Q_c − 2 is **coarser**, a refinement ratio of about 1.33 (for example 6 → 8).
- The remainder beyond Q_c is gap/((4/3)^p − 1): 1.3·gap for p = 2, 3·gap for p = 1, 6.5·gap for p = 0.5.
- The term measures λ and A only.
- G5 checks interpolated ω only, at 3 q points. The g / γ interpolation is not checked at all, although v1 C3 asked for |Δγ| or |Δλ_q| at off-grid q.

So a bias in the dvscf / Wannier-g interpolation that is common to Q_c and Q_c′ goes undetected. One that differs is under-covered for p < 1. Either case reproduces v1 C3: a wrong verdict in either direction.

**Remedy.**
- Make the check grid **denser** than Q_c, not coarser. Alternatively, use a direct DFPT el-ph calculation at ≥ 2–3 off-grid q, compared on γ_qν/λ_qν as well as on ω.
- Write G6 as a real gate with a threshold independent of e. For example: gap ≤ 0.05, else ERROR.
- Separately, add a rate-justified term to e.

### C3. Silent weight-loss channels in the phself route are open, and the A-DROP plumbing does not block
**(a) The tag does not block as specified.** `pipeline.blocking_tags` blocks rejection only for:
- `A-FILE-INCOMPLETE…`
- `A-FILE-COMPLETE(unverified)`
- `A-DROP-NONPOS(<m>K)` **with m > 0**.

For an excluded mode, EPW has already written λ = 0. `g2_tmp = 0` for ω ≤ eps_acoustic, including every imaginary mode.

So the dropped mass computed from the file is 0. The natural implementation emits `A-DROP-NONPOS(0K)`, or a new tag name such as `A-DROP-EXCLUDED`. Neither blocks. The excluded coupling disappears and R-UNIV/R-CELLS can reject with an A that is too small: a wrong REJECTED_PROVEN, contrary to HANDOFF §3.4.

The near-Γ imaginary exemption (radius 0.05|G_min|) also lets a genuinely soft optical branch near Γ drop out. This is the H3S-type case, where the dropped coupling is large.

**(b) Restart writes zeros.** `use_wannier.f90:897–928`: with `restart=.true.`, only the electron self-energy, spectral, transport and absorption data are re-read. `lambda_all` is re-zeroed (line 817) and **not restored**. A resumed phonselfen run therefore writes λ = 0 for every q before `iq_restart`. These are valid-looking non-negative zeros with no tag.
- ARCH §5 says the launcher "resumes".
- `epw_inputs.py` passes `restart` and `restart_step` through.

**(c) q weights are not checked.** `lambda.phself` has no q coordinates and no weights.
- With `mp_mesh_q=.true.`, the q points are the IBZ with non-uniform `wqf`.
- With `rand_q`/`filqf`, the points are arbitrary.

Applying 1/N_q silently mis-weights A and λ. Nothing in P1 pins these settings or checks row count = nqf1·nqf2·nqf3.

**(d) The qwindow may drop q points.** `qwindow_wrap` reduces `totq` to the q points that have (k, k+q) pairs inside fsthick. `lambda_all` is allocated with that `totq`, but `selfen_ph_write` loops over `nqtotf`. If totq < nqtotf, which can happen for small pockets (the risky materials), either:
- the rows are legitimately zero, or
- indexing goes out of the allocated range.

**Not verifiable from the excerpts** (`image_division` is not included). Also, `selecqread=.true.` reuses a `selecq.fmt` made with a different fsthick/σ.

**(e) nimage > 1 may double ω.** With `a2f=.true.`, `write_phdos` recomputes `wf` for **all** q on every image. Then `mp_sum(wf, inter_image_comm)` runs before `selfen_ph_write`, so with nimage > 1, ω would be multiplied by nimage. That scales A and leaves λ unchanged. **Plausible from the source, not verified.**

**Remedy.**
- Specify the exact tag string and require a blocking semantics that does not depend on a known mass. For example, a new `A-DROP-EXCLUDED(n modes)` added to `blocking_tags`, with a unit test.
- Pin these settings and refuse anything else: `restart=.false.`, `mp_mesh_q=.false.`, no `filqf`/`rand_q`, `selecqread=.false.`, `shortrange=.false.`, `ngaussw=0`, nimage = 1, nstemp = 1, `lsda`/noncollinear unsupported (or separately validated).
- Parser assertions:
  - row count = nqf1·nqf2·nqf3;
  - every λ ≥ 0;
  - per-q Σ_ν λ_qν equals the stdout `lambda___( tot )` for that q;
  - the global Σ matches `Summed el-ph coupling`;
  - at the Q_c points, ω equals the coarse dyn ω, which is a free check because interpolation is exact there;
  - the stdout line "We only need to compute N q-points" shows N = nqtotf, or else zero rows are proven legitimate.
- Each failed assertion ⇒ ERROR.

---

## MAJOR

### M1. The admission rule still admits false plateaus and under-covers slow convergence
- **"Monotone or converged" escape.** A non-monotone X(K_1), X(K_2), X(K_3) is admitted if |X3 − X2| ≤ 0.2 %, and then R = 4.4·0.2 % = 0.9 %.
  - Example: 1.00, 0.90, 0.901 passes. The sequence has shown a 10 % oscillation, which is the size of error to expect, and agreement of the last pair is the aliasing / false-plateau case of v1 C3. This is still open.
- **Fallback factor.** For "p < 0.5 or undefined" the factor is 4.4. That value was derived for ratio 1.5 at p = 0.5. For ratio 2:
  - the true factor is 1/(2^p − 1), which is 2.4 at p = 0.5, 6.2 at p = 0.2 and 14 at p = 0.1, and unbounded as p → 0;
  - "undefined" includes **diverging** sequences (|Δ32| > |Δ21|, so p < 0, pre-asymptotic). These get R = 4.4|Δ32| and pass when |Δ32| ≤ 0.45 %.
- **Dead term again for λ, A and ω_log.** Admission requires R/X ≤ 2 %, so 3R ≤ 6 % < 0.10 floor, and the λ/A/ω_log remainders never move e (the v1 C2a pattern). Only W_b, q and G6 can bind. This is acceptable only if the 0.10 floor is argued to cover the admitted worst case, and with the fallback factor above it is not for p < 0.15.
- **ω_log not in e.** R_ωlog is admitted but omitted from the e formula.

**Remedy.**
- Non-monotone ⇒ R ≥ max(|Δ21|, |Δ32|), or refuse.
- p < 0.5 or diverging ⇒ refuse the σ. Do not substitute a fixed factor.
- State the worst-case coverage argument for the floor.

### M2. The σ → 0 closure rests on an untested assumption; G-ANCHOR is nearly tautological and G4 is loose
- **λ_0 comes from the same data as the alternatives.**
  - λ_0 is extrapolated from admitted σ ≥ σ_min (≤ 20 meV). A quadratic extrapolation moves λ by at most about ⅓ of the λ(2σ_min) − λ(σ_min) difference.
  - Admitted λ values sit in [min λ_alt, max λ_alt] by construction, and the interval is widened by e ≥ 0.10.
  - So the gate fails only for enormous curvature. In practice the protection is just "σ_min,adm ≤ 20 meV", together with the unproven belief that λ(20 meV) is within e of λ(0).
- **The λ_0 + aσ² model is unjustified when a band edge or van Hove feature lies within about 1–8σ of E_F.** Examples are the MgB2 σ-band top and the H3S van Hove. There the smearing correction is not a low-order series over 20–160 meV. Also, E_F moves with σ (efnew at degaussw).
- **Unspecified cases.** The rule for fewer than 3 admitted σ is not given: skip, or ERROR? λ_0 is fitted at Q_1 but compared against λ_alt at Q_2. Only λ is anchored, not A.
- **The gate's form is necessary but not sufficient.** λ_0 ∈ [min(1−e), max(1+e)] can lie in a gap between alternative envelopes when the alternatives spread more than 2e. Then the truth is in no single alternative's envelope. The protocol relies on monotonicity in α²F, which the theorems do not state.
- **G4 (tetra) is weak.**
  - Its window is the whole alternative range ± e.
  - It compares λ only.
  - The tetra λ comes from a coarse **shifted** q grid with no q-convergence check: q2r cannot be used with `lshift_q`, and only two k grids are run.
  - It runs on a separate DFPT with `tetrahedra_lin` occupations.
  - `alpha2f.x` aborts on any imaginary frequency, so the result is an ERROR, not a check.
  - So G4 catches only gross (> 20–30 %) errors.

**Remedy.**
- Test the anchor assumption directly. Require σ_min,adm ≤ 10 meV, or require that λ at the two smallest admitted σ agree to e/3.
- Fail G-ANCHOR if the fit residual or curvature is large.
- Make G4 directional and tight: |λ_tetra − λ_0| ≤ max(e/2, q-error of the tetra grid). Give the tetra route its own q check, or state its q grid explicitly.
- Fewer than 3 admitted σ ⇒ ERROR.

### M3. The Γ point makes A-DROP fire on every material, or creates a spurious dominant atom
Every unshifted uniform q grid contains Γ, and both Q_1 and Q_2 do.
- **Acoustic modes at or below eps_acoustic.** With ASR enforced, the three acoustic modes at Γ have ω ≈ 0 ≤ eps_acoustic, so they are excluded and A-DROP is set. **Every** protocol material then has its rejections blocked unless the user sets the global `allow_unverified_rejection`. That flag also waives truncated-file and dropped-weight protections. So the protocol produces either no rejections, or rejections with all input protections off.
- **Acoustic modes just above eps_acoustic.** If ASR is imperfect, an acoustic mode at Γ sits at, say, 0.5 cm⁻¹ > eps_acoustic. It becomes an atom with λ_qν = γ/(πN_F ω²), where γ comes from the ASR-violating g at Γ. Its λ contribution is about 10⁵ times that of a 20 meV mode with the same γ, so it dominates λ, ω_log, the λ admission, G-ANCHOR and G4. G4 is computed on a shifted grid, so it disagrees.
  - Result: ERROR floods, or meaningless λ-based gates.
  - Its contribution to λ(k ≥ 1; T) is about γ/(πN_F N_q(2πkT)²), the same as one ordinary mode. That is a small upward bias in the CERTIFIED direction.

**Physics.** The q = Γ acoustic contribution is exactly 0 in the continuum: ω = 0, and g = 0 by translational invariance. A sound treatment exists, and it is not A-DROP:
- either use shifted q grids that exclude Γ (`filqf` with uniform weights, verified);
- or record exact-Γ acoustic exclusion as a documented, non-blocking identity, *after* checking that the Γ acoustic modes have ω < eps and |g| ≈ 0.

Off-Γ excluded modes (imaginary near Γ) have no sound bound on their dropped weight. For those, ERROR is more honest than a tag that blocks rejection silently or, under opt-in, globally.

### M4. Per-material omega_c conflicts with the config / hash model
P0 fixes omega_c_K per material (10 × coarse max ω). In rtsc, `omega_c_K` is a **global** `[certify]` key. It feeds `certify_hash`, `counts_for` and the reconcile path. Reconcile runs `certify.classify` on stored spectra with the current global `[certify]`.

A per-material value has no carrier:
- If it is written into `[certify]` per material, the hashes diverge, and counted CERTIFIED rows show "does not count".
- If it is not, reconcile and the other-goal re-certification revert to `omega_c_factor × max ω over alternatives`. That brings back v1 M4 silently.

This is not a soundness break (F3 records the value used), but certified statements become inconsistent between paths.

**Remedy.** Add an explicit per-spectrum or per-record omega_c field that classify consumes and the certificate records. Alternatively, keep one global fixed `omega_c_K` for protocol runs.

### M5. D2 plumbing: "ERROR stage executor" does not match where the gate fails, or the re-queue semantics
- The gates fail inside the DFPT stage of `pipeline.compute_material`. Those verdicts are recorded by `_record` with stage `pipeline` and a non-classify origin. That is guard-legal, and reconcile reuses non-classify outcomes "as is", so stored spectra are never re-certified.
- If the protocol instead writes a stage-`executor` row, that is also guard-legal (Lean `Ledger.lean:60`, and the SQL guard). But `_reconcile` checks `dft_hash` **only on stage-`pipeline` rows**. A non-retryable executor ERROR is therefore never re-queued when the config changes, for example a new protocol version or K_4. The "until the config changes" claim is false, and the hard cases are buried permanently, which is the M5 bias concern.
- The protocol thresholds live in a document. Unless the thresholds and the protocol version are part of the hashed config, a code change to a threshold leaves old verdicts counting.

**Remedy.**
- Record gate failures as stage `pipeline`, status ERROR, origin `rtsc.dfpt.eph_protocol`, `retryable=0`.
- Add ledger tests for: a [certify] change (no re-certification), a dft_hash change (re-queue), and goal addition.
- Hash the protocol version and the thresholds.

### M6. The windowed-weight criterion does not establish the per-atom envelope that certify uses
- E-ENV(e) is a per-atom box. C-RAYLEIGH takes per-atom endpoints, and R-CELLS takes (1+e)·w_j. At fixed q, individual λ_qν converge in k much more slowly than window sums. Per-atom relative changes of 30–100 % are compatible with W_b agreeing to 10 %. So E-ENV(e) is literally false per atom.
- Soundness then needs a stated lemma: a window-sum envelope plus kernel variation inside a window bounds λ(k;T) and Q. With 20 windows of width ω_max/20, that is plausible at T* = 300 K, but not for low-T* goals (77 K, or benchmark-level T*).
- Windows with < 1 % of A are ignored, yet low-ω windows can carry ≥ 10 % of λ. The relevant set should be weighted by λ(k;T*), not by A.
- The k admission is measured at Q_1 only, but applied (through e) to the Q_2 alternatives, whose q points differ (nesting vectors).

**Remedy.**
- Define the measured quantities as the functionals certify consumes. For example λ(k; T) for k ≤ k_max at T ∈ {T*, T_up guess}, and the Rayleigh form at its certificate x. Alternatively, prove the window lemma with explicit kernel-variation terms.
- Repeat the K test at Q_2.

### M7. The validation plan cannot expose wrong verdicts
- For Al at T* = 300 K the correct verdict is REJECTED. A sabotaged run that rejects is not "wrong", so the sabotage tests mostly check that gates fire, not that verdicts are right. Sabotage needs a T* near the converged Tc_model, both above and below it.
- "E_F shifted by +0.3 eV" has no mechanism under `efermi_read=.false.`.
- **Missing code-level sabotage**, the cheapest and most decisive tests:
  - an A-DROP tag that must block;
  - EPW restart zeros;
  - IBZ weights (`mp_mesh_q`);
  - `ngaussw=1`;
  - an oscillating or diverging synthetic sequence fed to the admission code;
  - a q-unconverged case (Q_c fine in k, coarse in q).
- For the benchmarks, "pass every gate or fail informatively" accepts any outcome. There are still no tight comparisons with published converged λ/ω_log (asked for in v1 m5).

---

## MINOR

- **m1.** P2 references "G2", which no longer exists in v2. K_4 = 12·K_c breaks nesting and the ratio 2 assumed in the Richardson formula. Which three grids are used after an extension, and with what ratio formula, is unspecified.
- **m2.** P2 says Q grids run "at the densest admitted k grid". P6 also needs (K_2, Q_2). Inconsistent.
- **m3.** Precision:
  - ω is printed E20.8, which is 8 significant digits, not 9. ω rounding is a position uncertainty (F5) and is formally outside the weight envelope. It is negligible, but it should be stated.
  - λ is printed E15.5. Fortran prints exponents with |exp| > 99 without the "E" (for example `0.12345-100`, possible with Gaussian tails ~exp(−200)). The parser must handle that or raise.
  - If stdout `lambda___(…)` lines (f15.6, absolute) are used for other σ, small λ_qν lose relative precision.
- **m4.** The imaginary-mode radius test needs q coordinates. The phself files carry only the index. Reconstruct from EPW's i, j, k ordering (valid only for `mp_mesh_q=.false.` nqf grids), with minimum-image distance.
- **m5.** G1 tolerance min(20 meV, σ_min/2) over E_F ± 1.28 eV is very strict for SCDM above the frozen window. Expect G1 ERROR floods. In the fallback "scdm μ += σ", which σ is meant is unspecified.
- **m6.** Cost:
  - The file holds smearing 1 only, so 6 σ × (3 K at Q_1 + 2 at Q_2) = 30 phself runs, doubled by G6, plus a shifted-grid DFPT for tetra.
  - Budget pressure invites deviation and causes DEFERRED/ERROR bias. Consider a one-line EPW patch that writes every `ismear` (with audit), or document the run count.
- **m7.** `lpolar` / Born-charge handling (v1 C3) is not mentioned. It is fine for clean metals, but needed for weakly metallic polar candidates.
- **m8.** SOC / noncollinear (the Pb benchmark is often run with SOC): verify EPW's `dosef/2` and `wkf` spin convention before accepting such runs.
- **m9.** For several `temps`, `lambda.phself.<T>K` exists per T (identical under `delta_approx`). Pin nstemp = 1, or define which file is parsed.

---

## v1 → v2 closure check

| v1 | closed? | note |
|---|---|---|
| C1 σ→0 | partly | σ_min ≤ 20 meV helps. G-ANCHOR is nearly tautological and G4 loose (M2). |
| C2 envelope | partly | Windows and Richardson exist for k only. The fallback factor is under-covering and non-monotone sequences are admitted (M1). q has no rate (C1). Window ≠ per-atom (M6). |
| C3 shared bias | partly | E_F per grid: closed. G1 tightened. G5 checks ω only, G6 is circular (C2), and the false plateau remains (M1). |
| C4 dropped weight | partly | fsthick closed (verified). The excluded-mode tag does not block (C3a). New silent channels exist (C3b–e). |
| M1 S5 | OK in principle (D1) | Needs an audited S5 text, and the flag set by protocol runs, not only by `use_smearings`. |
| M2 post-hoc | mostly | The large-σ guard is gone and tetra is mandatory. The extension rule points to the undefined G2 (m1). |
| M3 negative bins | closed | Holds iff `shortrange=.false.` (pin it); assert λ ≥ 0. |
| M4 ω_c | plumbing gap | M4 above. |
| M5 status | plumbing gap | M5 above. |
| m1–m7 | mostly | Except the K_4 nesting (m1) and validation (M7). |

---

## Verdict

**Not ready for implementation-with-sabotage-tests as a proved-verdict protocol.** Three CRITICAL items must be redesigned first:
- C1, q convergence;
- C2, the circular G6 and the missing g-interpolation check;
- C3, tag blocking and EPW pins/assertions.

**What can proceed now** (it is independent of the open design questions):
- the phself → Spectrum parser with all C3 assertions and pins;
- the A-DROP / `blocking_tags` unit tests;
- the ledger D2 tests.

M1 (admission refusal rules) and M3 (Γ handling) should be settled in v3 before any benchmark run. Otherwise the benchmarks will mostly measure ERROR floods or blocked rejections.

D1 and D2 are acceptable in substance once the plumbing in M4 and M5 is specified.
