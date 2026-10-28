# Pb benchmark ERROR: negative α²F at the smallest smearings

Date: 2026-10-27. Branch: `fix/pb-negative-a2f`. Code: `rtsc/dfpt/workflow.py` (`admit_smearings`), `rtsc/dfpt/parsers.py` (`mode_negative_fractions`), `rtsc/pipeline.py` (`_sigma_exclusion_directions`). Tests: `tests/test_dfpt_smearing_admission.py`.

## What failed

The `rtsc run` verdict for fcc Pb at 0 GPa was **ERROR**:

> negative alpha^2F mass at omega > 0 … refusing to clip (negative_policy='error')

This happened for σ₁ and σ₂. QE 7.6 setup:
- PBE with SG15 ONCV 1.2;
- k grid 16³ for scf and 32³ for the e-ph fit;
- q grid 8³ (29 irreducible points);
- matdyn grid 32³;
- 10 double-delta smearings σ_i = i × 0.01 Ry.

## Where the negative values come from

**Theory.** For each phonon mode ν,

α²F_ν(ω) = (1/N_F) Σ_q γ_qν / (π ω_qν) · δ(ω − ω_qν),  with  γ_qν = 2π ω_qν Σ_k |g|² · (non-negative double-delta weights).

Every factor is non-negative, so **α²F_ν(ω) ≥ 0 for every mode ν, exactly**. Negative values can only be numerical.

**The data** (`tests/data_dfpt/qe76_pb_bench/`, the real matdyn output):

| σ | σ (Ry) | per-mode negative mass / A | negative up to | λ (QE footer) | λ change if clipped |
|---|---|---|---|---|---|
| 1 | 0.01 | 3.1 × 10⁻⁴ (total) | 14.6 K | 1.0452 | +0.22 % |
| 2 | 0.02 | 1.2 × 10⁻⁴ (total) | 13.1 K | 1.0870 | +0.09 % |
| 3–10 | 0.03–0.10 | 0 exactly | – | 1.1108 → 1.1919 | – |

- **All three modes** go negative at low ω for σ₁. Pb has one atom, so all three branches are acoustic.
- **Mode by mode the problem is larger than the total shows.** At σ₁, mode 3 is still negative at 31 K. The total is positive there only because mode 1 is positive and the two cancel.
- **Al shows the same pattern.** In the Al benchmark (σ_i = i × 0.005 Ry), only σ₁ = 0.005 Ry is negative. That was the "under-resolved σ₁" in `STATUS_2026-10-27.md`.

**Mechanism.**
1. For a double-delta smearing that is too narrow for the 32³ k grid, γ_qν on the 8³ q grid is dominated by k-sampling noise.
2. matdyn then Fourier-interpolates γ to the 32³ q grid.
3. Near Γ, the true acoustic γ_qν → 0 and ω_qν → 0, so the interpolation's ringing produces negative γ there. The 1/ω factor puts it into low-ω α²F.

**Why this explanation fits the data:**
- the negatives sit only at the bottom of the spectrum (≤ 15–31 K, below the acoustic phonon range);
- they shrink monotonically with σ and vanish **exactly** from σ₃ on;
- they appear in every acoustic branch;
- the same thing appears in Al at its narrowest smearing.

## Options considered

| option | theory | verdict |
|---|---|---|
| Clip negatives to 0 (`negative_policy="clip"`) | Removes only the visible negative half of the ringing. It misses per-mode negatives hidden by cancellation, and it already fails the QE-footer consistency check on Al. | Rejected |
| Hand-pick `use_smearings = [3..10]` | Correct result, but a choice made after seeing the data, separately for each material | Rejected |
| **Admission rule: keep only smearings with non-negative per-mode α²F** | Tests a property that holds exactly in theory. It is decided by the data, before any certification, and is the same rule for every material. | **Adopted** |
| Fix at the source (denser k grid for the e-ph fit) | The real cure, and part of the EPW protocol (P3) | Future runs. The bench route keeps its fixed grids. |

## The rule (`smearing_admission = "mode_positive"`)

1. **Per-smearing test.** Smearing *i* passes if, for every mode ν, the negative α²F_ν mass at ω > 0 is ≤ 10⁻⁶ of the total positive mass. 10⁻⁶ is rtsc's existing print-noise tolerance.
2. **Admitted set.** The admitted smearings are the **contiguous block** from the smallest passing smearing upward, as in protocol P3.
   - Failures *below* the block are excluded: under-resolution is the only mechanism that produces them on a monotone ladder.
   - A failure *above* a passing smearing is an **ERROR** (`SMEARING_ADMISSION_FAIL`), because nothing explains it.
   - No passing smearing at all is also an ERROR.
3. **Records.** The exclusion is recorded in three places:
   - the run provenance (`smearing_admission`: per-σ fractions, footer λ, admitted and excluded smearings);
   - a warning;
   - the existing verdict flag `SMEARINGS_SUBSET`.
4. **Blocking tag.** Each admitted spectrum carries `A-SIGMA-EXCLUDED(sigma<i>:lam=<λ>:A=<A>K, …)`, holding the excluded spectra's own λ and A.

### Can excluding smearings make a proof unsound?

Excluded spectra break α²F_ν ≥ 0, so they are not admissible α²F functions. They cannot simply be ignored, though: they might point to a region the admitted alternatives don't cover.

rtsc treats each admitted alternative *s* as the box (1 ± e_s) w_j, with e = `rel_envelope` = 0.10 in this project. Over that box, λ and A range exactly over [(1 − e) X_s, (1 + e) X_s]. `blocking_tags` compares each excluded smearing's λ and A with the union of these ranges:

| Excluded λ or A lies … | Effect |
|---|---|
| above the union | blocks **REJECTED**, because it might be harder to reject |
| below the union | blocks **CERTIFIED**, because it might be harder to certify |
| inside the union | recorded, but blocks nothing |

A blocked proof becomes UNDECIDED unless `allow_unverified_*` is set. That is rtsc's existing downgrade-only policy.

**Limitation, stated plainly:** this comparison is at the level of moments (λ and A), not pointwise. A spectrum inside the λ and A ranges can still have a different shape. The envelope is already an assumption recorded in every verdict. This rule adds no new assumption type; it only decides whether exclusion needs a new block.

**On the real data, with no envelope (e = 0):**

| Material | Excluded | Excluded λ | Admitted λ | Without envelope | With e = 0.10 |
|---|---|---|---|---|---|
| Pb | σ₁, σ₂ | 1.045, 1.087 | min 1.111 | would block CERTIFIED | inside the band, nothing blocked |
| Al | σ₁ | 0.428 | max 0.415 | would block REJECTED | inside the band, nothing blocked |

## Physics caveats for Pb that this fix does not address

- **No spin–orbit coupling.**
  - The benchmark uses scalar-relativistic pseudopotentials.
  - Published DFPT studies report that SOC substantially increases λ for Pb, since it softens the phonons. **This is from memory and not yet checked against the paper.**
  - So without SOC the model's Tc is expected to fall *below* experiment, which affects any certification at T* ≈ 7.2 K.
  - A fully relativistic run needs SG15 FR pseudopotentials and `lspinorb`; it is recorded as follow-up work.
- **λ rises with σ across the admitted ladder** (1.111 → 1.192). The σ → 0 limit is not reached on the bench grids.
  - This spread is already part of the alternative set the certificate has to cover.
  - The protocol's tetrahedron anchor (P4) is the σ-free cross-check.

## Results with the rule

See `README.md` § Verification for the validation runs (Al and Pb re-collected from their unchanged QE outputs; QE was not re-run).
