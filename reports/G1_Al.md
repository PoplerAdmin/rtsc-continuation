# Gate G1 (Wannier band quality): fcc Al

Date 2026-10-27. Protocol v3 P4 G1 plus §A A7 (v3.1) and §B A10 (v3.2). Code:
`bench_tools/ephproto/bandcheck.py` (CLI `python -m ephproto.bandcheck run|compare|kpoints`),
tests `bench_tools/tests/test_bandcheck.py` (21 tests, synthetic fixtures). Run data:
`~/rtsc_work/epw_Al/g1/` (pre-registered k set) and `~/rtsc_work/epw_Al/g1_random1000_s20261027/`
(supplementary k set). Each holds `pw_bands/` and one subdirectory per Wannierisation with
`epw_bands/`, `epw_ef_K16/`, `epw_ef_K24/` and `g1_result.json`.

## Criteria as evaluated
σ_min,adm = 20 meV is assumed for now.
- **G1-a (P4):** max |ΔE| ≤ min(30 meV, σ_min,adm/2) = **10 meV** over the states in
  E_F ± 4σ_min,adm = **E_F ± 80 meV**. E_F(DFT) = 11.109 eV.
- **G1-b (A7):** the same number of DFT and Wannier states in that window at every validation k.
- **G1-c (A7):** |E_F(EPW fine grid) − E_F(DFT)| ≤ 50 meV.
- **Validation k count:** at least 100 k points that are not on the 8³ coarse grid.
- **Diagnostic (not gating):** a wide window E_F ± 1 eV with a 30 meV tolerance.

## Method
- **k points (pre-registered):** the half-shifted 5×5×5 MP grid in crystal coordinates.
  - The rule claimed that no point lies on the 8³ coarse grid. That is false for exactly one point:
    (½, ½, ½) = L is on the coarse grid (8 × ½ = 4).
  - That point is removed by an exact rational test, which leaves **124 points**.
- **DFT bands:** pw.x `bands` with the nscf.in settings (SG15 Al, 60/240 Ry, mv 0.02 Ry,
  nbnd 20). It runs on a copy of the scf charge density (`tmp/` is untouched).
  - Eigenvalues come from data-file-schema.xml (Hartree × 27.211386245988).
- **Wannier bands:** epw.x `band_plot=.true.` with `epwread=.true.`, `filkf` set to the same list,
  `filqf` set to Γ, and `efermi_read=.true.` with `fermi_energy=11.109`. nkf and nqf are not set.
  - Each run uses a private outdir that symlinks `out/pwscf.save` and the epmatwp of that
    Wannierisation, so `out/` is never written.
  - Energies come from `band.eig`.
  - The k lists of the two codes agree to better than 2·10⁻⁶ (2π/alat). This is checked on every run.
- **Matching:**
  - Remove DFT bands 1–4 (semicore). Sort the remaining DFT bands and the 4 Wannier bands at each k.
  - Pair Wannier band j with non-excluded DFT band j. A pair counts if either member is in the
    window, so missing and spurious states both appear as large |ΔE|.
  - Count the states in the window separately for each code (A7).
  - A nearest-neighbour distance is also reported. It does not depend on band order.
- **EPW fine-grid E_F:** the protocol phself key set (σ = 20 meV, Gaussian), nkf = 16³ (the gate
  value) and 24³ (a check), and 1 shifted q point. The value is read from the line
  "Fermi energy is calculated from the fine k-mesh".

## Finding 1: the pre-registered k set cannot test G1-a for Al
- On the 124-point shifted 5³ grid, **no DFT state lies within E_F ± 80 meV**. The closest is
  117 meV from E_F.
- The grid is highly symmetric, so its 124 points fall into only a few star-equivalent classes
  (minimum |E − E_F| values come in groups of 6 identical ones).
- So G1-a is vacuous on this set. The code reports it as `inconclusive_no_states` and does **not**
  count it as a pass.
- A **supplementary, post hoc k set** was therefore also run. It is not pre-registered and is
  labelled as such: 1000 uniform random crystal points from `numpy.random.default_rng(20261027)`,
  rounded to 10⁻⁶, none within 10⁻⁴ of a coarse-grid point. pw.x on it takes 34 s on 4 ranks.
- **Proposed amendment for v3.3:** replace the shifted-5³ rule with this seeded random set, with
  N_k ≥ 1000 and a minimum number of in-window states (for example ≥ 20).

## Results

"Pairs" counts the index-paired states with at least one member in the window. "k mism." counts
validation k points where the DFT and Wannier state counts in the window differ.

### New Wannierisation (frozen window, `dis_froz_max` = 13.109 eV): label `frozen`

| k set | window (tol) | max\|ΔE\| | RMS ΔE | pairs | k mism. | worst pair |
|---|---|---|---|---|---|---|
| shifted 5³ (124, pre-reg.) | ±80 meV (10) | 84.9 meV | 81.6 meV | 4 | 4 / 124 | band 2, DFT 11.262, Wannier 11.177 eV |
| shifted 5³ (124, pre-reg.) | ±1 eV (30) | 379.9 meV | 187.2 meV | 54 | 12 / 124 | band 2, DFT 12.046, Wannier 12.426 eV |
| random 1000 (suppl.) | ±80 meV (10) | 288.8 meV | 132.4 meV | 51 | 39 / 1000 | band 3, DFT 11.346, Wannier 11.057 eV |
| random 1000 (suppl.) | ±1 eV (30) | 421.1 meV | 122.1 meV | 407 | 27 / 1000 | band 2, DFT 11.978, Wannier 12.399 eV |

EPW fine-grid E_F: **11.1709 eV at 16³** (+61.9 meV vs DFT) and 11.1246 eV at 24³ (+15.6 meV).

**Per-band breakdown (random 1000 set).** These are mean / max |ΔE| for Wannier minus DFT, all
inside the frozen window:
- band 1 below E_F − 5 eV: −6 / 102 meV
- band 1 from E_F − 5 to E_F − 1 eV: +10 / 100 meV
- band 2 from E_F − 1 to E_F + 1 eV: −4 / 421 meV
- band 3 from E_F − 1 to E_F + 1 eV: −82 / 303 meV
- bands 2 and 3 from E_F + 1 to E_F + 2 eV: about +300 / 700 meV
- above E_F + 2 eV the errors are 2–4 eV. This is outside the frozen window and is expected.

**Pipeline check.** On 40 coarse-grid k points (every 13th point of the nscf list), band_plot
reproduces the nscf eigenvalues inside the frozen window to **0.05 meV** (72 states). So the parsing,
units and k order are right. The off-grid errors are genuine Wannier-interpolation error of this
4-WF set on the 8³ grid. The frozen window makes the bands exact on the grid but not smooth between
grid points.

**Wannier90 diagnostics** (`wannier/pwscf.wout`):
- Disentanglement stopped at the 200-iteration limit with "convergence criteria not satisfied".
  Ω_I = 5.1792 Å², changing by 1.2·10⁻⁵ per iteration.
- Final Ω = 6.394 Å² (spreads 1.69, 1.56, 1.57, 1.57 Å²).
- The WF centres are not on the symmetric tetrahedral sites.
- The outer window spans all 16 non-excluded bands, up to 63.3 eV.

### Discarded Wannierisation (no frozen window): label `discarded_nofrozen`

| k set | window (tol) | max\|ΔE\| | RMS ΔE | pairs | k mism. |
|---|---|---|---|---|---|
| shifted 5³ (124, pre-reg.) | ±80 meV (10) | no states (inconclusive) | | 0 | 0 / 124 |
| shifted 5³ (124, pre-reg.) | ±1 eV (30) | 6977 meV | 4009 meV | 68 | 52 / 124 |
| random 1000 (suppl.) | ±80 meV (10) | 8803 meV | 4810 meV | 31 | 31 / 1000 |
| random 1000 (suppl.) | ±1 eV (30) | 9403 meV | 4029 meV | 517 | 370 / 1000 |

EPW fine-grid E_F: **15.0696 eV at 16³** (+3960.6 meV) and 15.0586 eV at 24³.
Band 1 is off by about 2.4 eV even 5 eV or more below E_F.

## Pass or fail per criterion

| criterion | frozen: pre-reg. 5³ | frozen: random 1000 | discarded: pre-reg. 5³ | discarded: random 1000 |
|---|---|---|---|---|
| N_k ≥ 100 off-grid | PASS (124) | PASS (1000) | PASS (124) | PASS (1000) |
| G1-a max\|ΔE\| ≤ 10 meV in ±80 meV | **FAIL** (84.9) | **FAIL** (288.8) | INCONCLUSIVE (0 states), counted as fail | **FAIL** (8803) |
| G1-b equal state counts (A7) | **FAIL** (4 k) | **FAIL** (39 k) | vacuous pass (0 states) | **FAIL** (31 k) |
| G1-c \|ΔE_F\| ≤ 50 meV, 16³ (A7) | **FAIL** (+61.9 meV) | same run | **FAIL** (+3961 meV) | same run |
| (info) \|ΔE_F\| at 24³ | 15.6 meV, would pass | | 3950 meV | |
| (diag.) ±1 eV, 30 meV | FAIL (380) | FAIL (421) | FAIL (6977) | FAIL (9403) |
| **G1 overall** | **FAIL** | **FAIL** | **FAIL** | **FAIL** |

- The new Wannierisation is a large improvement: E_F is off by 62 meV instead of 3.96 eV, and errors
  near E_F are about 0.1–0.4 eV instead of 7–9 eV. It still fails G1 by an order of magnitude on the
  band criterion.
- On E_F alone it fails at 16³ and passes at 24³. The 16³ value is also limited by k sampling with
  20 meV Gaussian smearing. The protocol does not yet fix the grid used for the A7 E_F check, and it
  should.
- **Suggestions for the next Wannierisation (not done here):**
  - raise `dis_num_iter` (for example to 2000) so disentanglement converges;
  - limit the outer window (`dis_win_max`, for example E_F + 15 to 20 eV) instead of keeping all 16 bands;
  - try sp3 projections centred on the tetrahedral interstitials instead of SCDM;
  - apply the protocol's one fixed fallback (scdm_μ += scdm_σ) if staying with SCDM;
  - consider a denser coarse k grid.

## Notes
- **Provenance caveat:** `~/rtsc_work/epw_Al/wannier/DONE` did not exist when this was run.
  - The new Wannierisation's epw.x had finished: it wrote epmatwp and the .fmt files at 17:33, then
    stopped with the usual post-epwwrite "Cannot load fine q points" abort. The discarded run ended
    the same way.
  - G1 was run on those files read-only. Their sizes and mtimes were recorded before the run
    (`g1_frozen_inputs_stat.txt`) and were unchanged afterwards.
  - If DONE later refers to a different run, delete `g1*/frozen/` and rerun
    `python -m ephproto.bandcheck run material.json --wannier ~/rtsc_work/epw_Al/wannier --epmatwp
    ~/rtsc_work/epw_Al/out/pwscf.epmatwp --label frozen --np 4 --ef-grids 16 24 [--kset random]`.
- The discarded run fails every criterion by a wide margin, as it should. The worst pairs are band 3
  states near E_F that the Wannier set puts 7–9 eV too high. Its EPW E_F is 15.07 eV at 16³ and
  15.06 eV at 24³, 3.96 eV above DFT, which matches what §B A10 recorded.
- Both E_F runs report the EPW coarse-grid E_F as 11.109035 eV, so the reference is consistent.
- The ±80 meV window holds only about 31 states on 1000 random k. The gate statistic is therefore
  limited by sampling near E_F. If a tighter statement is needed, a Fermi-surface-targeted k set
  would be the next refinement.
