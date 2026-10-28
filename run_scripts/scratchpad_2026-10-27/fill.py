p='/mnt/g/ClaudeWorkspace/RTSC_Project/reports/G1_Al.md'
s=open(p).read()
res='''"Pairs" counts the index-paired states with at least one member in the window. "k mism." counts
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
Band 1 is off by about 2.4 eV even 5 eV or more below E_F.'''
ver='''| criterion | frozen: pre-reg. 5³ | frozen: random 1000 | discarded: pre-reg. 5³ | discarded: random 1000 |
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
  - consider a denser coarse k grid.'''
s=s.replace('RESULTS_TABLE',res).replace('VERDICT_TABLE',ver)
s=s.replace('## Notes\n','''## Notes
- **Provenance caveat:** `~/rtsc_work/epw_Al/wannier/DONE` did not exist when this was run.
  - The new Wannierisation's epw.x had finished: it wrote epmatwp and the .fmt files at 17:33, then
    stopped with the usual post-epwwrite "Cannot load fine q points" abort. The discarded run ended
    the same way.
  - G1 was run on those files read-only. Their sizes and mtimes were recorded before the run
    (`g1_frozen_inputs_stat.txt`) and were unchanged afterwards.
  - If DONE later refers to a different run, delete `g1*/frozen/` and rerun
    `python -m ephproto.bandcheck run material.json --wannier ~/rtsc_work/epw_Al/wannier --epmatwp
    ~/rtsc_work/epw_Al/out/pwscf.epmatwp --label frozen --np 4 --ef-grids 16 24 [--kset random]`.
''')
open(p,'w').write(s)
