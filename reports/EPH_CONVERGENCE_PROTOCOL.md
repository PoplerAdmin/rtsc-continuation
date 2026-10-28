# α²F convergence protocol (contract), v3. 2026-10-27

History:
- v1 → review 1 (`review_eph_protocol_v1.md`): 4 CRITICAL.
- v2 → review 2 (`review_eph_protocol_v2.md`): 3 CRITICAL, 7 MAJOR.
- v3 answers review 2 (§R).

All versions were written before any protocol result existed. Thresholds change only with a
recorded reason and a new version.

**Status: DRAFT.** Outputs are labelled `EPH_PROTOCOL_DRAFT`. Certify results from them are
benchmark diagnostics, not pipeline verdicts, until all three hold:
- a third review of v3 finds no CRITICAL item;
- the sabotage tests (§V) pass;
- the ARCH S5 amendment (§D1) is written into ARCH.md and reviewed.

## 0. What this protocol can and cannot establish
The theorems are about the spectrum they are given. Discretisation error of α²F is outside every
proof and stays an explicit assumption (DISCRETISATION; E-ENV(e) recorded as *estimated*).
The protocol makes the given spectrum a pre-registered, well-converged, audited representative
and quantifies the residual. It cannot turn convergence evidence into proof, and it does not
claim to.

## P0. Fixed inputs (per material, before any protocol output)
- **Coarse DFPT** (`rtsc run`, QE 7.6, norm-conserving pseudopotentials):
  - dyn + dvscf on Q_c;
  - nscf on the full K_c grid, with the **same symmetry settings as the DFPT** (nosym crashes
    EPW 6.1 `gmap_sym`).
- **EPW pinned settings:**
  - `phonselfen=.true., delta_approx=.true., a2f=.false.`;
  - `restart=.false., selecqread=.false.`;
  - `nstemp=1`, `ngaussw=0`, `efermi_read=.false.`;
  - no `-ni` (nimage = 1), `use_ws=.true.`;
  - `lpolar=.false.` (metals; a polar material ⇒ out of scope, ERROR).
- **ω_c fixed a priori:** `omega_c_K = 10 × max ω` of the coarse-DFPT phonon DOS, rounded up to
  100 K. It is passed explicitly (`rtsc certify --omega-c` / `[certify] omega_c_K`) and is the
  same for every alternative and every run of the material.
  - Pipeline integration needs a per-material carrier that `certify_hash`, `counts_for` and
    reconcile all see (integration task I1).

## P1. Spectra are finite atomic measures built from modes
- **Fine q points** are an explicit **shifted** Monkhorst–Pack list (shift ½ in every direction,
  no Γ) passed with `filqf`. Coordinates and weights 1/N_q are therefore known exactly: no IBZ,
  no random q. **Fine k points** are a uniform full `nkf` grid.
- **Construction:**
  α²F = Σ_{q,ν} (1/N_q) · λ_qν ω_qν / 2 · δ(ω − ω_qν),
  - λ_qν from `lambda.phself.*`;
  - ω_qν from `linewidth.phself.*` column 3 (meV). The Γ = 2γ column is never used to rebuild λ.
  - Review 2 confirmed this reproduces EPW's own λ and A, with no factor-2 error.
- **Parser assertions** (any failure ⇒ ERROR):
  - row count = N_q·n_modes;
  - contiguous q indices;
  - λ ≥ 0 and finite;
  - ω > 0 for every mode;
  - per-q Σλ agrees with the stdout `lambda( sum )` lines within print precision;
  - λ_spectrum = q-average of Σ_ν λ_qν to 1e-6.
  - Exponent fields without an "E" are handled. Precision: ω has 8 significant digits, λ 5;
    both are recorded.
- **Weight loss is impossible by construction, or an ERROR:**
  - every (q,ν) on a shifted grid is an atom;
  - any mode with ω ≤ eps_acoustic or imaginary ⇒ ERROR `UNCONVERGED_EPH(imaginary)`;
  - `fsthick = max(1 eV, 8σ)`, so the Gaussian tail < e^(−64).
  - No A-DROP tag is needed in normal operation. If a future variant must drop modes, it uses
    rtsc's existing blocking tag `A-DROP-NONPOS(<m> K)` with m > 0, covered by a unit test.
- **Atom merging:** atoms with identical printed ω are merged exactly.

## P2. Grids and smearings (ratio 2 everywhere)
- **σ ∈ {5, 10, 20, 40, 80, 160} meV**, fixed for every material.
- **k:** K_1 = 2K_c, K_2 = 4K_c, K_3 = 8K_c. One extension K_4 = 16K_c only if no σ is admitted at
  K_3.
- **q (shifted):** Q_1 = 2Q_c, Q_2 = 4Q_c, Q_3 = 8Q_c.
- **k admission** is run at Q_2. **q admission** is run at the densest admitted k grid with the
  admitted σ. **Alternatives** are taken at the densest k and q grids (P6), so admission is
  measured where it is applied.
- All metrics come from the parsed Spectrum objects that certify will see.

## P3. Admission (k direction for each σ; then q direction), the same rule in both
Applied to the three grids G_1, G_2, G_3 (ratio 2) for each X ∈ {λ, ω_log, A, W_b}, where W_b are
window weights:
- **windows:** 20 equal ω windows over (0, ω_max]. A window is tracked if it holds ≥ 1 % of A
  **or** ≥ 1 % of λ (the 1/ω weighting catches low-ω windows).
- **Monotone with a sane rate:** the differences d_1 = X(G_2) − X(G_1) and d_2 = X(G_3) − X(G_2) have
  the same sign, and the observed rate p = log2(|d_1|/|d_2|) ≥ 0.5.
- **Converged:** if both |d_1| and |d_2| are ≤ 0.2 % of X, the item passes regardless of p. This is
  not an escape for oscillation: both steps must be tiny.
- **Remainder:** R_X = |d_2| / (2^p − 1), with p capped at 4. The factor is at most 2.41 at
  p = 0.5. Sequences with p < 0.5, oscillating or diverging, fail.
- **Admitted iff:**
  - every item passes;
  - R/X ≤ 2 % for λ, ω_log and A, and ≤ 10 % for each tracked window;
  - E_F changes by ≤ σ/4 between G_2 and G_3 (k direction).
- **k:** σ is admitted if it passes at Q_2. Admitted σ must form a contiguous block that includes
  the smallest σ that passes. Otherwise ERROR `UNCONVERGED_EPH(k)`.
- **q:** with the admitted σ at the densest admitted k, the q sequence must pass for every
  admitted σ. Otherwise ERROR `UNCONVERGED_EPH(q)`.

## P4. Independent anchors (σ → 0 and shared bias)
- **Tetrahedron (mandatory for NC inputs):**
  - ph.x `lambda_tetra`, `occupations='tetrahedra_lin'` (non-negative weights), shifted q;
  - converged by P3 rules on two refinements of k and two of q, gated on λ, A and ω_log.
  - **Gate G4:** |λ_tetra − λ(σ_min,adm)| ≤ 3 %, |A_tetra − A(σ_min,adm)| ≤ 3 %, and the same for ω_log.
  - Imaginary modes (alpha2f.x aborts) ⇒ ERROR.
  - The tetrahedron is the σ-free anchor. The λ(σ) extrapolation of v2 is dropped as a gate
    (review 2 M2) and only reported.
- **σ_min,adm ≤ 20 meV is required.** Otherwise ERROR `UNCONVERGED_EPH(sigma)`.
- **Direct DFPT check of the interpolated coupling (G6, replaces v2 G5/G6):**
  - at 3 off-grid q points (fixed rule: the first three points of the Q_1 list that are not on
    Q_c), run ph.x single-q DFPT with `electron_phonon='simple'` on K_3 at σ_min,adm;
  - compare γ_qν with EPW interpolated at the same q, k and σ: |Δγ|/γ ≤ 5 % for every mode with
    γ ≥ 1 % of the maximum, and |Δω|/ω ≤ 2 %;
  - fixed threshold, independent of e.
- **G1 Wannier bands:** max |ΔE| ≤ min(30 meV, σ_min,adm/2) over states within E_F ± 4σ_min,adm,
  on ≥ 100 off-grid k points.
  - One fixed fallback: scdm_μ += scdm_σ (1.5 eV), recorded. Then a final pass or fail.

## P5. Envelope
e = max(0.10, 3·max over admitted σ of { R_λ/λ, R_ωlog/ω_log, R_A/A, max_b R_Wb/W_b } over both
the k and q directions, plus the G4 gaps and the G6 max relative Δγ).

e > 0.5 ⇒ ERROR `UNCONVERGED_EPH(envelope)`. e is written into each Spectrum's
`rel_envelope` (the config floor still applies) and recorded as `E-ENV(e, estimated)`.

## P6. Alternatives (pre-registered set; see §D1)
- **Forwarded:** every admitted σ at (K_densest, Q_3), and at (K_densest, Q_2), and at
  (K_previous, Q_3).
- **Not forwarded:** every other computed spectrum, listed in provenance with its reason.
- **Verdicts** carry the flag `SMEARINGS_SUBSET`.
- **Provenance:** grids, σ, fsthick, eps_acoustic, E_F per grid, ω_c, every gate value, p, R,
  tetra and G6 numbers, EPW and QE versions, input hashes, and the protocol version and
  threshold hash.

## P7. Status on gate failure (§D2, decided: ERROR)
ERROR with reason `UNCONVERGED_EPH(<gate>)`, written as:
- stage `pipeline`;
- origin `eph_protocol` (non-classify);
- `retryable = 0`.

The protocol thresholds are part of `dft_hash`, so a protocol change re-queues the material
through the normal reconcile path. Ledger tests are required (integration task I2). Gate failures
are listed as "hard cases" in reports.

## §D. User decisions
- **D1 (decided 2026-10-27: amend S5).** Rationale: forwarding spectra known to be under-resolved
  adds wrong inputs to the conjunction without making any verdict more correct. A pre-registered
  selection with full disclosure is the scientifically sound practice.
  Proposed ARCH S5 text: "For runs with `[spectra] eph_protocol = "v3"`, the alternatives are
  exactly the pre-registered set defined by that protocol version; every other spectrum computed
  in the run is listed in the verdict provenance with the reason; such verdicts carry
  SMEARINGS_SUBSET."
- **D2 (decided 2026-10-27):** ERROR, as in P7.
- **D3 (decided 2026-10-27): no symmetry reduction of the fine q grid.** Fine q sums stay over
  full uniform grids with equal weights, as in P1/P2. EPW `mp_mesh_q` and IBZ filqf lists are not
  used. Rationale: the protocol must not depend on symmetry that Wannier-interpolated couplings
  only approximately satisfy, and it must stay material-agnostic. The literature check that led
  here is in `lit_ibz_q_precedent.md`. Cost is controlled by grid caps, lfast_kmesh/fsthick and
  scheduling instead.

## §V. Validation and sabotage (review 2 M7)
- **Benchmarks:** Al, Pb, Nb, MgB2, H3S@200 GPa. Report every gate, and the certified bracket at
  T* = experimental Tc for μ*_hi = 0.10 and 0.15.
- **Sabotage, evaluated at T* close to the material's bracket edge**, so a wrong pass could flip a
  verdict. Each case must end in ERROR UNCONVERGED_EPH:
  - Al Q_c = 2³;
  - fsthick = 2σ (hand-edited, outside the pins: must be caught by the parser/pin check);
  - E_F shifted by +0.3 eV;
  - SCDM μ = E_F − 5 eV;
  - K_1 = K_2 = K_3;
  - restart=.true. with a truncated run (λ = 0 rows);
  - `ngaussw=1`;
  - a synthetic oscillating sequence;
  - a q-unconverged synthetic case.
- **Code-level tests:** parser assertions (synthetic fixtures), blocking-tag behaviour, ledger
  status rows.

## §I. Integration tasks (before pipeline use)
- I1: a per-material ω_c carrier.
- I2: ledger/status rows and tests.
- I3: threshold hash in `dft_hash`.
- I4: the S5 amendment in ARCH.md, the README and the correspondence docs.
- I5: an independent audit of the implementation.

## §A. v3.1 amendments (review 3: no CRITICAL; applied 2026-10-27 before any protocol result)
These override the corresponding text above.

- **A1 (review 3 M1): identical or misreported grids.**
  - N_k and N_q are pre-registered per run and checked against EPW stdout (k-mesh size, "Size of
    q point mesh", "We only need to compute N = N_q").
  - Every q in the filqf list must have a stdout block with matching coordinates, and its
    `lambda___( tot )` must match the file.
  - The three grids in any sequence must have distinct input hashes. Otherwise ERROR (pin
    violation).
  - In the converged branch, R_X = max(|d_1|, |d_2|), never 0.
  - Every row must have ω > 0. This is the backstop for EPW's out-of-bounds handling of dropped q
    points.
- **A2 (M7): exact formulas.**
  - e = max(0.10, 3·max(R-terms), G4 gaps, G6 relative Δγ). The G4 and G6 terms enter without the
    factor 3.
  - After a K_4 extension the sequence is (K_2, K_3, K_4).
  - The window edges use ω_max of the coarse DFPT (fixed per material).
  - The alternatives are only (K_densest, Q_3) and (K_densest, Q_2) for the admitted σ. The
    (K_previous, Q_3) set is dropped.
  - The stdout key is `lambda___( tot )`.
- **A3 (M6): EPW input whitelist.**
  - Only the keys the driver writes are allowed.
  - The full epw.in is hashed and the hash enters the threshold hash.
  - Any other key ⇒ pin ERROR.
  - A pin violation is its own ERROR reason (`EPH_PIN_VIOLATION`), distinct from UNCONVERGED_EPH.
- **A4 (M3): assumptions are emitted.**
  - Every protocol Spectrum carries the provenance key `discretisation`, so rtsc emits
    DISCRETISATION, and the assumption string `E-ENV(e, estimated by eph_protocol v3.1)`.
  - A unit test checks both appear in the verdict.
  - The ω position uncertainty is handled by alternatives: for every forwarded spectrum, add two
    copies with all ω scaled by (1 ± δω), where δω = max(G6 |Δω|/ω, 0.5 %). F5 forbids putting
    position uncertainty into e.
- **A5 (M2): G6 redesign.**
  - Test points are the 3 points of Q_1 with the largest Σ_ν λ_qν at σ_min,adm, excluding points
    within 0.1·|G_min| of Γ.
  - The comparison is at 2·σ_min,adm (less noise). Degenerate modes are compared as summed sets.
  - Threshold: |Δγ|/γ ≤ 10 % per set (set γ ≥ 1 % of max), and |Δω|/ω ≤ 2 %.
  - Workflow: reuse the coarse scf settings exactly; ph.x single q with `electron_phonon='simple'`
    and `trans=.true.` on the k grid K_2 (cost accepted).
  - Before any pipeline use, G6 must be shown to FAIL on Al with Q_c = 2³ (in §V).
- **A6 (M4): the tetrahedron is a different DFT model** (tetrahedron occupations vs mv 0.02 Ry).
  - For v3.1, G4 is reported with the tolerance widened to 5 % and is NOT used in e.
  - The scf-model offset is measured on Al/Pb/Nb, and a calibrated G4 enters v4.
  - Until then, verdict use requires v4 (DRAFT status unchanged).
- **A7 (m1, m2):**
  - G1 also requires the same number of states within E_F ± 4σ_min,adm at every validation k
    (DFT vs Wannier).
  - EPW's fine-grid E_F must lie within 50 meV of the DFT E_F.
- **A8 (m3): near-Γ acoustic modes** on shifted grids are kept as atoms, so no weight is dropped.
  If they make the P3 sequence non-monotone, the material ends in ERROR, recorded as a hard
  case. That is accepted as conservative for now and revisited after the Al numbers.
- **A9 (M5): cost** is measured end to end on Al before the other materials are scheduled at
  full protocol size.

## §B. v3.2 amendment (2026-10-27; before any admitted result)
- **A10: Wannier frozen window.** `dis_froz_max = E_F(DFT) + 2 eV`; `dis_froz_min` is left at the
  bottom of the Wannier manifold.
  - **Reason:** the first Al run without a frozen window gave EPW fine-grid E_F = 15.07 eV vs DFT
    11.11 eV, with 3.00 electrons counted correctly. So the SCDM-disentangled 4-band subspace
    did not reproduce the DFT bands. Check A7 (|ΔE_F| ≤ 50 meV) flags exactly this.
  - That run is discarded and listed.
  - A frozen window must hold ≤ nbndsub bands at every k; a W90 error ⇒ ERROR
    `UNCONVERGED_EPH(wannier)`.

## §R. Review 2 → v3
| finding | change |
|---|---|
| C1 no q control | P2/P3: three shifted q grids, same Richardson rules; q remainders in e |
| C2 G6 cannot fail | P4 G6: direct single-q DFPT γ vs EPW γ, fixed 5 % threshold |
| C3 silent weight loss | P0 pins (restart, nimage, selecqread); P1 shifted filqf (known q, weights, no Γ); parser assertions; existing blocking tag |
| M1 false plateaus | P3: p ≥ 0.5 required, both steps tiny for the converged escape, ratio-2 factor, R_ωlog in e |
| M2 tautological anchor | P4: mandatory converged tetra gate on λ, A and ω_log; extrapolation demoted to a report |
| M3 Γ / A-DROP | P1: shifted q grids have no Γ; excluded modes ⇒ ERROR |
| M4 ω_c carrier | P0 + I1 |
| M5 D2 mismatch | P7: stage pipeline, origin eph_protocol, retryable 0, thresholds hashed (I2, I3) |
| M6 windows vs per-atom envelope | §0 states the limit; windows also tracked by λ share; admission measured where applied |
| M7 sabotage at the wrong T* | §V: sabotage at the bracket edge, plus code-level tests |
| m1–m9 | numbering fixed, ratio 2 with K_4 = 16K_c, P2/P6 consistent, precision stated, q coordinates from filqf, G1 tolerance widened and bounded, cost accepted, lpolar pinned, SOC spin convention: check before Pb-SOC |
