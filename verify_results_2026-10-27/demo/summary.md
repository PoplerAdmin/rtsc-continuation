# rtsc run summary: rtsc-demo (SYNTHETIC)

> **SYNTHETIC DATA**: 7 of 7 inputs use invented (synthetic) spectra. Their verdicts exercise the software and say nothing about any real material.

- config hash: `199c23fdf2fd77ea`; finished 2026-10-27T18:27:59+00:00
- target T* = 300.0 K; mu* in [0, 0.15]; weight envelope e >= 0.05
- Coulomb cutoff: CUTOFF(block: N_c(T)=#{n>=0 : (2n+1)*pi*T <= omega_c}; omega_c=10*max(omega) per spectrum)
- rejection model: isotropic (isotropic: rejections assume M-ISO-SINGLEBAND)
- inputs: **7**, each with exactly one status (checked):

| status | count |
|---|---|
| CERTIFIED_CANDIDATE | 1 |
| UNDECIDED | 1 |
| REJECTED_PROVEN | 2 |
| EXCLUDED_ASSUMPTION | 2 |
| DEFERRED | 1 |
| ERROR | 0 |
| **total** | **7** |

## What the verdicts mean (and do not mean)

**What a verdict is about.** Every REJECTED_PROVEN / CERTIFIED_CANDIDATE / UNDECIDED verdict is a
statement about a MODEL, not directly about the material:

* isotropic Migdal-Eliashberg theory, linearised gap equation (Tc_model := sup of unstable T,
  which is >= the usual highest-solution Tc).  With the default rejection_model = "isotropic"
  every rejection uses the DOS-averaged single-band alpha^2F (assumption M-ISO-SINGLEBAND, flag
  MULTIBAND_RISK): in a multiband material the averaged A can underestimate Tc;
* the given alpha^2F (from DFPT or an attached file), within the stated relative weight
  envelope e (true atom weights in [(1-e)w, (1+e)w]); frequency uncertainty enters only through
  the alternative spectra (e.g. DFPT smearings), all of which must agree (S5);
* Coulomb pseudopotential mu* in [0, mu_hi] with the block-form cutoff
  N_c(T) = #{n >= 0 : (2n+1) pi T <= omega_c} (CERTIFIED depends on this rule and omega_c;
  REJECTED holds for every cutoff and every mu* >= 0).

**What a verdict is not.** It is not a prediction of the measured Tc.  Out-of-model physics -
anharmonicity and quantum nuclear motion (important for hydrogen), non-adiabatic effects,
anisotropic or multiband gaps, spin fluctuations, vertex corrections, the accuracy of DFT itself,
thermodynamic and kinetic stability - is NOT covered.  Such risks appear as flags; flags never
change a status.  EXCLUDED_ASSUMPTION verdicts rest on physical, non-proved criteria (hull,
metallicity, magnetism, imaginary modes) and DEFERRED means "not computed".  REJECTED_PROVEN and
CERTIFIED_CANDIDATE are produced only by rtsc.certify.classify, from rigorous ball-arithmetic
checks of rules proved in Lean (library PipelineLogic; docs/CORRESPONDENCE_LEAN.md) -- except
where a verdict's assumptions say otherwise (e.g. row_check = "cform": hand proof, not in Lean);
such verdicts are listed separately below and marked in the proof_basis column.

| status | meaning |
|---|---|
| CERTIFIED_CANDIDATE | proved: Tc_model >= T* for every spectrum in the envelope, every alternative spectrum and every mu* <= mu_hi (with the stated cutoff) |
| REJECTED_PROVEN | proved: Tc_model < T* for every spectrum in the envelope, every alternative and every mu* >= 0 |
| UNDECIDED | the certified bracket [T_low, T_up] contains T*, or alternatives disagree |
| EXCLUDED_ASSUMPTION | removed by a physical, non-proved criterion (recorded) |
| DEFERRED | not computed (budget, ordering, duplicate, no spectrum) |
| ERROR | a stage failed; see the record |

## Top candidates (certified first, then undecided; by certified T_low)

Brackets are single-band (DOS-averaged) unless the bracket column says otherwise.

| material | formula | status | T_low [K] | T_up [K] | bracket | A_max [K] | lambda | flags |
|---|---|---|---|---|---|---|---|---|
| SYNTHETIC-strong-coupling | (synthetic) | CERTIFIED_CANDIDATE | 413.5845207080677 | 727.2338659097802 | single-band (DOS-averaged) | 3000 | 2.84-3.01 | SYNTHETIC_DATA, MAGNETIC_PROXIMITY, STRONG_COUPLING |
| SYNTHETIC-borderline | (synthetic) | UNDECIDED | - | 412.09919068220876 | single-band (DOS-averaged) | 1700 | 0.85-0.85 | SYNTHETIC_DATA |

## Rejections by rule

- CLASSIFY-S5: 2
- R-CELLS: 1
- R-MU: 1
- R-ODD: 1
- R-UNIV: 1

2 rejection(s) hold in the single-band (DOS-averaged) model only (M-ISO-SINGLEBAND); see rejection_model = "multiband".

## Exclusions by criterion (physical assumptions, not proofs)

- X-HULL: 1
- X-METAL: 1

Thresholds: `{"hull": {}, "metallic": {}, "magnetic": {}, "dynamic": {}}`

## Deferred (not computed)

- budget: 1

## Out-of-model flags (never change a status)

- SYNTHETIC_DATA: 7
- MULTIBAND_RISK: 2
- MAGNETIC_PROXIMITY: 1
- STRONG_COUPLING: 1

Files: `candidates.csv`, `candidates.json`, `certificates/<id>.json`, `<workdir>/<id>/record.json`, `<workdir>/<id>/verdict.json` (input of `rtsc recheck`).
