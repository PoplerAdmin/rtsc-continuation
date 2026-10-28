# Five-material validation protocol (pre-registered)

Written 2026-10-27, 23:05 EDT. At that time:
- the only finished pipeline verdicts were Al and Pb at T* = 300 K;
- the Pb runs at T* = Tc_exp had not finished;
- Nb, MgB2 and H3S had no spectra yet.

The criteria below are fixed for all five materials. They will not be changed after any further results are seen. If a later change is unavoidable, it will be listed below with its date and reason, and the original text will be kept.

## What is being validated

rtsc proves statements about **Tc_model**: the Tc of the isotropic Migdal–Eliashberg model, given the computed α²F alternatives, the envelope e = 0.10, and every μ* in [0, μ*_hi]. A `CERTIFIED_CANDIDATE` at T* proves Tc_model ≥ T* > 0, which is a proof of superconductivity **within the model**. Whether the material superconducts in reality also depends on the model and DFT approximations: isotropic gap, harmonic phonons, no spin–orbit coupling, PBE. Comparing with experiment tests those approximations; it does not test rtsc's logic.

## Materials and reference Tc

| Material | Conditions run | Experimental Tc used |
|---|---|---|
| Al (fcc) | 0 GPa | 1.175 K |
| Pb (fcc) | 0 GPa | 7.19 K |
| Nb (bcc) | 0 GPa | 9.25 K |
| MgB2 (P6/mmm) | 0 GPa | 39 K |
| H3S (Im-3m) | 200 GPa | Measured ≈ 203 K at 155 GPa and lower at higher pressure (Drozdov et al., Nature 2015). Value at 200 GPa is taken as ≈ 190 K, **read from memory of the published Tc(P) curve, not checked**. |

All five use the same pipeline settings:
- QE 7.6, PBE, SG15 norm-conserving, scalar-relativistic;
- the bench grids in each `config.toml`;
- `smearing_admission = "mode_positive"`;
- `rel_envelope = 0.10`.

Spectra come from the unchanged QE outputs. Only `rtsc` collection and certification are rerun.

## Tests per material

| Test | Run | Pass condition | What it shows |
|---|---|---|---|
| V1 Rejection soundness | T* = 300 K, μ*_hi = 0.15 | `REJECTED_PROVEN` for Al, Pb, Nb and MgB2. For H3S it is reported only, because its harmonic Tc_model may truly reach about 300 K. | rtsc does not certify known low-Tc metals at room temperature |
| V2 Superconductivity proof | T* = ½·Tc_exp, μ*_hi = 0.13 | `CERTIFIED_CANDIDATE` | rtsc proves Tc_model > 0 for a known superconductor, with margin |
| V3 Accuracy, reported, not pass/fail for rtsc | T* = Tc_exp, μ*_hi = 0.10 and 0.13 | status and bracket [T_low, T_up] recorded | whether the model with conventional μ* reproduces experiment. Expected failures are listed below. |

**Run order:** V2 and V3 use the same spectra as V1. Each is a separate `rtsc run` that reuses the QE outputs; `bench_tools/jobs/validate_material.sh` does this.

## Outcomes expected for physics reasons

These were listed before the results and are not failures of rtsc:

- **Pb:** no spin–orbit coupling. SOC raises λ substantially in published DFPT work (from memory, unverified), so Tc_model may fall below 7.19 K, giving V3 = UNDECIDED.
- **MgB2:** the isotropic model averages over its two gaps. Isotropic Eliashberg with DFT α²F is known to give Tc well below 39 K, so V3 = UNDECIDED or REJECTED is likely. V2 at 19.5 K is uncertain.
- **H3S:** the calculation is harmonic. Anharmonicity lowers λ substantially in the literature, so a harmonic Tc_model above experiment is expected. V3 should then be CERTIFIED. A V1 rejection at 300 K is **not** expected to be easy.

## Changes after registration

None.
