# Adversarial review: α²F convergence protocol (EPH_CONVERGENCE_PROTOCOL.md, v1)

Reviewer: independent referee, 2026-10-27. Read: HANDOFF §0, §3; ARCH §2–4b; README §1, `[spectra]`,
`[certify]`; the protocol; STATUS (Al case); source `spectrum.py`, `dfpt/parsers.py`,
`dfpt/attached.py`, `pipeline.py` (envelope and tag handling), `certify/classify.py` (omega_c).
Nothing was changed.

Severity: CRITICAL = can produce a wrong proved verdict, meaning a verdict whose recorded
assumptions (E-ENV(e), the alternative set) are false for the converged DFT α²F, in the direction
that matters. MAJOR = model-correct but unreliable, or a status/contract violation. MINOR = the rest.

Scope note. Every proved verdict is technically true *for the parsed spectra with their recorded
envelope*. The protocol's stated purpose, though, is to make the envelope plus alternatives cover
the discretisation error. Where that claim fails, `E-ENV(0.10)` is a false assumption recorded as
if established. The CRITICAL items are failures of that kind.

---

## Background facts used below (checked in the source)

- The envelope is per atom and multiplicative on fixed positions (`Spectrum.rel_envelope`). R-UNIV
  uses only (1+e)·A. R-CELLS uses (1+e)·w_j per atom inside λ(k;T) = Σ w_j g_T(k, ω_j), where
  g_T is largest for small ω. C-RAYLEIGH uses the per-atom endpoint that minimises Q.
- `classify` uses one omega_c = factor × max ω over **all** alternatives
  (`CertifyParams.resolved_omega_c`). README `[certify]` says "per spectrum", which is wrong.
- `pipeline._with_envelope` applies e = max(the spectrum's own e, `[spectra] rel_envelope`).
- `load_attached_spectra(kind="epw")` returns one Spectrum per **phonon**-smearing column of the
  `.a2f` file, with rel_envelope 0. The EPW electron smearing is a single value per run.
- `atoms_to_spectrum` clips negative mass up to 1e-6·A without raising and tags it A-NEG-CLIP,
  which blocks CERTIFIED. Weight at ω ≤ 0 is dropped and tagged A-DROP-NONPOS, which blocks
  REJECTED. The EPW consistency check compares the parsed λ with EPW's own integrated λ from the
  **same** truncated grid.

---

## CRITICAL

### C1. No anchor to the σ → 0 limit. Admission favours the most biased smearings.
**Scenario.** On a small-Fermi-surface metal (MgB2 σ bands, H3S, a hydride candidate) the small σ
(5, 10, 20 meV) fail the 2 % N2/N3 test, because they are not resolved at 72³ or at the extensions.
The large σ (40, 80, 160 meV) converge easily in N, so they pass. G2 is satisfied (3 contiguous)
and G3 is satisfied: λ(σ) varies slowly at large σ, for example 0.95 / 1.00 / 1.04. The true
σ → 0 value is monotone outside that range, say 0.85 or 1.15. G4 could catch this, but it is
optional ("where available"; a skip is not a failure). Its tolerance is also *outside* the range
(±5 %).
**Why it is wrong.** For fixed σ, N-convergence proves nothing about the finite-σ bias. Neither
the alternatives (point samples at finite σ) nor the envelope span the σ → 0 limit when λ(σ) is
monotone. The 10 % envelope only absorbs a bias up to 10 % of every atom's weight.
- If the admitted spectra overestimate coupling, C-RAYLEIGH at the envelope minimum can certify
  T_low ≥ T* for the parsed spectra while the converged α²F gives Tc_model < T*. The lower bound
  is the tight one (Al: T_low 0.58 K vs experiment 1.18 K), so a 10–20 % λ bias moves it across
  T*. Result: wrong CERTIFIED_CANDIDATE.
- If they underestimate (λ and A low), R-CELLS or R-UNIV can reject near threshold. Result: wrong
  REJECTED_PROVEN. This is less likely because the upper bounds are loose.

**Remedy.**
- Require the smallest admitted σ to be at most a fixed small value (for example 20 meV), or
  require an extrapolation σ → 0 (a fit over the admitted block plus the converged part of the
  smaller σ). Widen e by |λ(σ_min) − λ_extrap|/λ, or add the extrapolated spectrum as an
  alternative.
- Make the tetrahedron route mandatory wherever QE supports it, which is all NC cases. When it is
  unavailable, record that the σ-bias is uncovered as an assumption that blocks the verdict, for
  example `A-SIGMA-BIAS(unchecked)`.
- Put G4's ±5 % tolerance *inside* the envelope argument: require
  |λ_tetra − λ_adm| ≤ e·λ_adm − (the known spread).

### C2. The envelope rule does not measure what the envelope means. The 3δ term is dead code.
**(a) Dead term.** Admission requires |ΔA|/A ≤ 2 % for every admitted σ. Therefore δ ≤ 0.02,
3δ ≤ 0.06, and e = max(0.10, 3δ) is always exactly 0.10, the old config default. The envelope
never responds to the measured convergence, and the "envelope derivation" in the provenance is
vacuous.

**(b) Wrong quantity.** δ measures the change of the *total* A. The envelope bounds *each atom's*
weight. Grid refinement redistributes weight between frequencies. Examples: a Kohn anomaly or soft
mode sharpening, nesting weight moving into a low-ω acoustic branch, a pocket appearing. A can
agree to 0.5 % while low-ω bins change by 30–100 %.
- R-CELLS: λ(k;T) weights low-ω atoms most heavily, so a redistribution towards low ω raises the
  true λ(k;T) beyond the (1+e)·w_j bound. The bound T_up is then not an upper bound for the
  converged spectrum. Result: wrong REJECTED, possible even when R-UNIV would not fire.
- C-RAYLEIGH: a redistribution towards high ω lowers λ. Result: wrong CERTIFIED.
- The (N2, σ) alternatives cover only the N2 → N3 shape change, not the N3 → ∞ remainder (F5
  says frequency uncertainty must be covered by alternatives; one step of history is not a cover).

**(c) The extrapolation factor is unjustified.** With error ∝ N^−p and N2 = 48, N3 = 72, the
remainder after N3 is δ/((1.5)^p − 1):
- p = 2: 0.8δ
- p = 1: 2δ
- p = 0.5: 4.4δ, which exceeds 3δ

Gaussian-smeared Fermi-surface integrals are not in the asymptotic regime for the small σ. A
non-monotone sequence has no such bound at all (see C3).

**Remedy.**
- Derive e from per-region (coarse frequency windows), λ- and λ(k;T)-relevant differences, not
  from total A. For example, take the max over windows W of |ΔA_W|/A_W with windows of a few
  meV, weighted toward low ω.
- Use all three grids to estimate p, and refuse (DEFER) when the order is not established or the
  sequence is non-monotone.
- Make the factor depend on the measured p.
- Remove the dead `max` or tighten the admission tolerance so the 3δ term can bind.

### C3. "Two densest grids agree" cannot detect biases shared by N2 and N3.
Every fine grid is interpolated from the *same* coarse DFPT and the *same* Wannier model. Each of
the following biases all of N1–N4 identically, passes the 2 % test and G3, and is not in any
alternative:
- **Coarse q-grid / dvscf adequacy.** No check that interpolated phonons and g match direct DFPT
  at an off-grid q. Kohn anomalies, polar materials without `lpolar`, and long-range g are all
  missed. G1 checks bands only.
- **Wannier band error.** G1 allows 50 meV within E_F ± 1 eV. That is 10× σ_min, and for a
  shallow pocket (MgB2 σ bands, ~0.4–0.6 eV) it is a 10–20 % change in pocket area, hence in
  N(E_F) and λ_σ.
- **Fermi level.** If E_F is read from the coarse nscf (`efermi_read`) or otherwise fixed, every
  fine grid shares the shift. If it is recomputed per grid, the variation is not recorded or
  bounded.
- **False plateau.** Both grids miss the same tiny pocket at small σ, or the sequence oscillates
  and two adjacent members agree by chance. With up to two extensions, there are up to three
  chances to hit an agreeing pair. N1 is computed but never used to test monotonicity.

**Why it is critical.** Each produces a common-mode error in A, λ and the shape, of either sign,
with no alternative or envelope term covering it. The direction analysis is as in C1/C2.

**Remedy.**
- Gate the phonon and g interpolation: direct DFPT at ≥ 2 off-grid q, |Δω| and |Δγ| (or |Δλ_q|)
  tolerances. Repeat the run with a denser coarse q grid on a validation subset.
- Require `lpolar` for non-zero Born charges.
- Tighten G1 for small σ, or scale it to σ_min.
- Record E_F per grid and require |ΔE_F| ≪ σ_min.
- Require N1 → N2 → N3 to be monotone (or within the tolerance at all three) before admission.

### C4. Weight dropped silently in EPW makes A too small (violates HANDOFF §3.4).
None of these produces an A-DROP/A-FILE tag, so none blocks rejection:
- **`fsthick`.** EPW includes only states within fsthick of E_F. With σ up to 160 meV and a
  typical fsthick of 0.3–0.5 eV (≈ 2–3σ), the Gaussian tails of the double delta are truncated.
  The protocol does not tie fsthick to σ_max.
- **`eps_acustic` / imaginary modes.** EPW skips modes with ω below eps_acustic, which includes
  interpolation-induced imaginary modes on the fine grid. Their coupling disappears from α²F.
  Nothing in the protocol applies X-DYNAMIC to the EPW fine grid or records how many modes were
  skipped.
- **Phonon-smearing tails.** The `.a2f` frequency grid starts at dw > 0 and stops at a finite
  wmax, so Gaussian tails (degaussq) below dw or above wmax are cut.

The parser's consistency check compares against EPW's integrated λ computed from the same
truncated data, so it passes. The 2 % admission also passes, because the truncation is the same
on every grid.

**Remedy.**
- Require fsthick ≥ (k_c·σ_max + ω_max) with k_c ≥ 6, or drop σ above fsthick/6 *before* the
  sweep. That choice is a priori, not post hoc.
- Count and record modes excluded by eps_acustic or imaginary frequency on the fine grid, with
  their λ_qν. A non-zero count should set an A-DROP-type tag that blocks REJECTED unless the user
  opts in.
- Check A independently: compare Σ_qν λ_qν ω_qν /(2 N_q) against Σ w_j.
- Preferably (see m6), build the Spectrum directly from the (ω_qν, λ_qν) atoms, with no phonon
  smearing and no binning.

---

## MAJOR

### M1. Dropping produced spectra conflicts with S5 as written and with ARCH §5.
- S5: rejection must hold "for every alternative spectrum the run produced (e.g. every DFPT
  smearing)".
- ARCH §5: "Each smearing gives an alternative spectrum… Never silently pick the best."

The protocol computes 6 σ × 3 grids (more with extensions, and × nqsmear phonon-smearing columns
per EPW file). It forwards only the admitted (N3, σ) and (N2, σ) spectra and the tetra spectra.
- The non-admitted σ and N1 are discarded.
- Which EPW phonon-smearing columns become alternatives is **unspecified**. Picking one column is
  an unrecorded selection.

For Al, the coordinating session itself treated excluding σ₁ as "a choice the user still has to
approve". The protocol automates exactly that choice. This is a contract change, not input
handling.

**Remedy.**
- Amend ARCH S5 with a versioned, audited definition of "alternative produced by the run".
- Set the existing SMEARINGS_SUBSET flag on every protocol verdict.
- Record every discarded spectrum with its reason.
- State that every phonon-smearing column (or a fixed, a-priori one) is forwarded.
- Keep a config opt-in (for example `[dfpt] eph_protocol = "v1"`) so the subset rule is an
  explicit user choice.

### M2. Steps that can act as post-hoc selection.
- **Large-σ guard.** It is ambiguous. Large σ almost always pass the N2/N3 test, so they are only
  ever "above the admitted block" if something else removed them. The only plausible trigger is a
  G3 failure. Read that way, the guard lets an operator drop σ = 160 (then 80, …) until the λ
  spread is ≤ 10 %, turning a DEFERRED into a verdict. **Remedy:** state that admission is
  computed once by the P3 rule. No σ may be removed to satisfy G3, and a G3 failure is final.
- **G4 skippable.** "Where available" with "skip is not a failure" makes the only σ-independent
  check optional. A material whose tetra result would DEFER passes if the route is not run.
  **Remedy:** make availability an a-priori, per-material rule (NC pseudopotential ⇒ required)
  recorded before the run, and make a skip block CERTIFIED, at least (see C1).
- **Tetra check is under-specified.**
  - "Two k grids", but no convergence criterion for the tetra result itself.
  - G4 compares λ only, not A or ω_log.
  - The tetra route (ph.x + matdyn) is the one that produces negative bins via Fourier
    interpolation of γ, so tetra spectra will often fail the negativity rule and be dropped. The
    set of alternatives then changes in a data-dependent way.

  **Remedy:** apply the same N2/N3-style test to the tetra route. Gate on λ, A and ω_log. Record
  tetra spectra that were dropped.
- **Grid extensions.** Up to two extensions raise the chance of a spurious agreeing pair.
  **Remedy:** require agreement across the last three grids, not the last two.

### M3. The negative-bin rule is the wrong signal and is not well defined.
- **Vacuous in the EPW Gaussian route.** With Gaussian electronic and phonon smearing, λ_qν ≥ 0
  and α²F ≥ 0 identically, so the rule never fires and gives no convergence information there. In
  the Al case the negatives came from QE matdyn's Fourier interpolation of γ(q), a q-interpolation
  artefact, not k-under-resolution. Methfessel–Paxton smearing would also produce them.
- **Where it is applied is undefined.** Raw table values or parsed atoms? The parser already
  clips negative mass ≤ 1e-6·A into 0 with A-NEG-CLIP. If the rule is checked after parsing,
  "mass > 0" never sees those values, and the A-NEG-CLIP tag then blocks CERTIFIED with no
  explanation. Print precision (E-format zeros) also hides small negatives.
- **Negatives at ω ≤ 0 are ignored.** They are dropped by the parser and only tagged.

**Remedy.**
- Apply the rule to the raw file values, including ω ≤ 0 rows.
- Treat any non-zero A-DROP-NONPOS as non-admission too.
- Do not use negativity as a convergence signal. It is a sanity check.

### M4. omega_c depends on which alternatives were admitted.
omega_c = 10 × max ω over **all** forwarded alternatives (code). A phonon-smearing tail, a
spurious high-frequency interpolated mode in one N2 spectrum, or a tetra spectrum changes omega_c
for the whole verdict. That shifts the meaning of μ* in the certified statement (F3), which STATUS
already shows matters (Al needs μ* > 0.10 in this convention). Numerically, a larger N_c is
conservative for C-RAYLEIGH at fixed μ*, so this is not a soundness break. But certified
statements become incomparable between runs, and the protocol's selection steps silently change
the cutoff convention.

**Remedy.**
- Require a fixed `omega_c_K` for protocol runs, chosen a priori, for example from the coarse
  DFPT ω_max × 10, rounded.
- Fix README `[certify]` ("per spectrum" → "max over all alternatives").

### M5. DEFERRED + UNCONVERGED_EPH is semantically off and biases discovery.
README defines DEFERRED as "not computed" (budget, ordering, duplicate, DFPT disabled).
UNCONVERGED_EPH is "computed, inputs not trustworthy". This is not an S2 violation (DEFERRED proves
nothing) and S4 holds if exactly one status is written. However:
- In the discovery ledger, DEFERRED means "awaiting selection", so the loop can re-select and
  re-run identical expensive DFT, or re-run it with altered settings until it passes, which is
  cherry-picking by repetition.
- The materials most likely to fail G1–G3 (small Fermi surfaces, strong-coupling hydrides) are
  the interesting ones, which repeats the H3S/ML-instability burying problem.

**Remedy.**
- Either use ERROR with reason UNCONVERGED_EPH(...), or keep DEFERRED with a non-requeue stage
  (for example `stage="eph_protocol"`).
- Amend README §1 and the ledger semantics.
- Forbid automatic re-runs with changed protocol parameters unless a new protocol version is
  recorded.
- Track the fraction of UNCONVERGED_EPH by chemistry as a bias monitor.

---

## MINOR

- **m1. Grids are not nested.** 72 is not a multiple of 48, so the "nested sequence" claim is
  false. N4 = 96 and N5 = 120 need not be multiples of the coarse grid (for example coarse 9³:
  24 is not a multiple). N1 is unused. Use grids that are genuine multiples, and use N1 in the
  monotonicity/order test (C2c, C3).
- **m2. The envelope hand-off mechanism is unspecified.** The EPW loader sets rel_envelope = 0,
  and the pipeline applies max(own, config). If a user lowers `[spectra] rel_envelope`, the
  protocol's e is lost unless it is written into each Spectrum. Specify where e is set, and
  assert it in the pipeline.
- **m3. Admission metrics are not defined precisely.** Compute λ, ω_log and A from the *parsed*
  Spectrum that certify sees, not from EPW's printout. Name the denominator grid. Use an absolute
  floor for small λ and ω_log.
- **m4. Small-ω and ASR artefacts.** Acoustic-sum-rule violations near Γ produce spurious low-ω
  λ_qν, which make the λ and ω_log admission metrics noisy. eps_acustic hides some of them (see
  C4). Enforce ASR in EPW, and record the low-ω weight below a fixed ω.
- **m5. The validation plan cannot detect the failures above.** "Experimental Tc inside the
  bracket" is weak because T_up is ~9× experiment. Add tight λ/ω_log comparisons with published
  converged EPW values, and add *sabotage tests*: a deliberately coarse q grid, a small fsthick,
  a shifted E_F, a σ-biased set. The protocol must DEFER on each of them.
- **m6. Phonon smearing (degaussq) and binning are unspecified.** Both move weight in frequency
  (F5). The cleanest remedy is to build the Spectrum from the atoms (ω_qν, λ_qν ω_qν /(2 N_q)):
  no phonon smearing, no bins, no negative bins, no tail truncation. This is exactly rtsc's
  atomic-measure model.
- **m7. Provenance gaps.** Add fsthick, degaussq, nqsmear, eps_acustic, E_F per grid, the coarse
  k/q grids, the counts of skipped and imaginary modes, the omega_c used, and the list of
  discarded spectra with reasons.

---

## Direction summary (per question asked)

- **REJECTION** (upper bounds; uses A_hi, and (1+e)·w_j in R-CELLS):
  - Safe against a uniform underestimate of A up to 10 %.
  - Unsafe against: weight redistribution toward low ω (C2b); common-mode under-coupling (C1, C3);
    silently dropped weight (C4).
  - R-UNIV alone is robust except for C1, C3 and C4.
- **CERTIFICATION** (Rayleigh at the envelope minimum, conjunction over alternatives):
  - Safe against a uniform overestimate up to 10 %.
  - Unsafe against: finite-σ or common-mode over-coupling (C1, C3); redistribution toward high ω
    (C2b).
  - This is the more exposed direction, because the lower bound is tight.
- **Admission of uniformly biased sets.** Yes. A bias that is flat between N2 and N3 is admitted
  by construction (C1, C3).
- **Cherry-picking.** Possible through the large-σ guard reading, the skippable G4, repeated
  extensions and an unspecified choice of phonon-smearing column (M1, M2). The Wannier retry is a
  fixed rule and acceptable.
- **S2/S4.** No S2 violation. S5 wording is violated by discarding produced spectra (M1).
  DEFERRED semantics need amending (M5).

Bottom line: the protocol is a reasonable procedure, but it is not yet sound for its stated
purpose. As written, the envelope is a constant 10 % unrelated to the measurements. Nothing covers
the σ → 0 limit, the coarse-q/Wannier common mode, or EPW weight truncation. C1–C4 must be
resolved, and S5/ARCH amended for M1, before any protocol output feeds a proved verdict.
