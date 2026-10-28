# fcc Al: optimized-tetrahedron α²F route (QE 7.6), cross-check against smearing

Date: 2026-10-27. Status: **PARTIAL**. Only the nk=16 e-ph grid finished. The nk=24 and nk=32 runs were killed by
two WSL failures; the host C: drive was full, so the WSL vhdx could not grow. Those runs are pending a relaunch.

Labels used here:
- **[PROVED]** means only a `rtsc certify` verdict that `rtsc recheck` re-verified.
- **[FLOAT]** means a numerical diagnostic, never a proof.
- **[SRC]** means verified in the QE 7.6 source (raw.githubusercontent.com/QEF/q-e/qe-7.6/...).

All data are real DFT. None are synthetic.

## 1. Workflow (verified against the QE 7.6 source)

1. **pw.x scf** with `occupations='tetrahedra_opt'` or `'tetrahedra_lin'`, and `K_POINTS automatic`.
   - Blöchl `'tetrahedra'` is rejected.
   - `lambda_tetra` without tetrahedron occupations is rejected.
   - Source: [SRC phq_readin.f90:982-988].
2. **ph.x run 1 (DFPT)**: `ldisp=.true., lshift_q=.true., nq1-3, fildyn, fildvscf`.
   - The shifted q grid is required: at q=Γ, `elph_tetra_lambda` zeroes the weights and prints
     "λ_qν is singular at q = Γ, you should use shifted q grid" [SRC elph_tetra_mod.f90:122-131].
   - `alpha2f.x` also aborts on any ω ≤ 0 [SRC alpha2f.f90 read_lam].
   - With `lshift_q`, q2r/matdyn cannot be used, so there is **no Fourier interpolation**. The q integration is done by
     tetrahedra directly on the DFPT q grid [SRC INPUT_PH.def lshift_q; alpha2f.f90 compute_a2F].
3. **ph.x run 2 (e-ph)**: the same input plus `electron_phonon='lambda_tetra'` and `nk1-3`.
   - `lambda_tetra` forces `trans=.false.`: the dyn files and dvscf are read from run 1 [SRC phq_readin.f90:464-467].
   - ph.x runs its own nscf for k and k+q on the dense `nk1-3` grid, using symmetry-reduced tetrahedra. No separate pw
     nscf is needed.
   - The bands used are only those crossing E_F (`elph_nbnd_min..max` from phq_setup.f90:412-424). The `nbnd` default
     (or the scf `nbnd`) is enough. I used nbnd=8.
   - `lambda_tetra` uses double-delta weights δ(ε_k)δ(ε_k+q) with ω neglected, the same approximation as the smearing
     'interpolated' route. `gamma_tetra` keeps ω in the step functions.
4. **alpha2f.x** (serial only) takes the same input as run 2, plus `&INPUTa2F nfreq=500 /`.
   - It writes `<prefix>.a2F.dat` with columns ω [Ry], α²F, phonon DOS, PDOS.
   - The grid is ω_i = i·ΔE for i = 1..nfreq, with ΔE = ω_max/nfreq.
   - α²F = ½ Σ λ_qν ω_qν δ(ω−ω_qν), computed by tetrahedron DOS on the q grid.
   - It prints λ and ω_ln two ways: from α²F as Σ 2·α²F_i·ΔE/ω_i, and directly from the tetrahedron q-weights.
   - λ_qν is read back with `(a19,f8.4)` [SRC], so 4 decimals.

### 1a. `tetrahedra_lin` is broken in QE 7.6 (evidence)

The coordinator asked for `tetrahedra_lin` as the primary route, because its weights are non-negative. It is unusable:

- **ph.x crashes with SIGSEGV** on the first q, right after the nscf. I ran it at the production settings
  (`q6s_lin/ph.err`).
- I reproduced the crash on a reduced test: 30 Ry, k 8³, q 2³ shifted, run both serial and with 4 ranks.
- The `tetrahedra_opt` control run with otherwise identical inputs finishes (`~/rtsc_work/tetra_Al/linbug/`).
- Symbolized backtrace (ASLR off, `addr2line`): `dfpt_tetra_mod::dfpt_tetra_main` ← `setup_alpha_pv` ←
  `phq_setup` ← `initialize_ph` ← `do_phonon`.
- **Cause [SRC]:**
  - `LR_Modules/dfpt_tetra_mod.f90` hard-codes `DO ii = 1, 20` at lines 283, 345, 429, 626, 746, 945, 1004 and 1201.
  - For `tetra_type=1`, `opt_tetra_init` allocates `tetra(4,ntetra)` and `wlsm(4,4)` (`nntetra=4`)
    [PW/src/tetra.f90:380-392], so these loops read out of bounds.
  - `alpha2f.f90:320` (`compute_lambda`) has the same hard-coded `20`.
  - `elph_tetra_mod.f90` correctly uses `nntetra`.
- So in 7.6 the tetrahedron e-ph route only works with **`tetrahedra_opt`**. Its smoothing matrix has negative entries
  [SRC tetra.f90:394-427], so negative α²F values are possible in principle. That is why I checked every bin (§3).

## 2. Inputs (as run)

Settings:
- PBE, SG15 `Al_ONCV_PBE-1.2.upf` (norm-conserving), 60/240 Ry.
- Cell from the smearing run's scf.in (vc-relaxed, a = 4.022 Å).
- Work dir: `~/rtsc_work/tetra_Al/q6s_opt/`.

**Choice of q grid: 6³ shifted.** A shifted 8³ grid has 60 irreducible q (spglib count; unshifted 8³ has 29). DFPT on
the shifted 6³ grid (28 q) already took 3 h 52 min on 4 ranks. The machine was heavily oversubscribed by other jobs
(load average 25-100). The shifted 8³ grid would have taken more than 8 h.

`scf.in`:
```
 &control
    calculation = 'scf'
    prefix = 'al'
    outdir = './tmp/'
    pseudo_dir = '/home/user/pseudo/sg15/'
 /
 &system
    ibrav = 0
    nat = 1
    ntyp = 1
    ecutwfc = 60.0
    ecutrho = 240.0
    occupations = 'tetrahedra_opt'
    nbnd = 8
 /
 &electrons
    conv_thr = 1.0d-12
    mixing_beta = 0.5
 /
ATOMIC_SPECIES
 Al  26.98153860 Al_ONCV_PBE-1.2.upf
CELL_PARAMETERS angstrom
     -0.00000000000000     2.01105392100000     2.01105392100000
      2.01105392100000    -0.00000000000000     2.01105392100000
      2.01105392100000     2.01105392100000    -0.00000000000000
ATOMIC_POSITIONS crystal
 Al      0.00000000000000     0.00000000000000     0.00000000000000
K_POINTS automatic
 16 16 16 0 0 0
```
`ph.in` (DFPT; after the forced restart, the run was resumed with `recover=.true.`):
```
Al phonon tetra_opt shifted q
 &inputph
    prefix = 'al'
    outdir = './tmp/'
    tr2_ph = 1.0d-16
    fildyn = 'al.dyn'
    fildvscf = 'dv'
    ldisp = .true.
    lshift_q = .true.
    nq1 = 6
    nq2 = 6
    nq3 = 6
    amass(1) = 26.9815386
 /
```
`elph.in`, which is also the `alpha2f.x` input (NK = 16; 24 and 32 pending):
```
Al electron-phonon lambda_tetra
 &inputph
    prefix = 'al'
    outdir = './tmp/'
    tr2_ph = 1.0d-16
    fildyn = 'al.dyn'
    fildvscf = 'dv'
    ldisp = .true.
    lshift_q = .true.
    nq1 = 6
    nq2 = 6
    nq3 = 6
    amass(1) = 26.9815386
    electron_phonon = 'lambda_tetra'
    nk1 = NK
    nk2 = NK
    nk3 = NK
 /
 &inputa2f
    nfreq = 500
 /
```
Commands:
- `mpirun -np 4 --bind-to none -x OMP_NUM_THREADS=1 pw.x|ph.x -nk 4`
- `alpha2f.x -in elph.in` on 1 rank.

## 3. Results per grid [FLOAT]

The q grid is 6³ shifted (28 irreducible q) and the DFPT k grid is 16³ for every row.

| e-ph k grid | λ (α²F, bins) | λ (direct q-tetra) | ω_log (α²F) | ω_log (direct) | A = ∫α²F dω | N(E_F) states/spin/Ry | negative bins | ph.x e-ph wall |
|---|---|---|---|---|---|---|---|---|
| 16³ | **0.37942** | 0.38032 | 305.29 K | 303.63 K | **61.47 K** | 2.7439 | **0 of 500** (min α²F = 0.0) | 47 min (2817 s) |
| 24³ | pending (killed by the WSL crash) | | | | | | | |
| 32³ | pending | | | | | | | |

Notes on the nk=16 row:
- The trapezoid rule on the same points gives λ = 0.37935 and A = 61.46 K.
- ω_max = 470.0 K and ΔE = 5.953·10⁻⁶ Ry.
- Wall times: scf 112 s; DFPT 3 h 52 min (28 q, 4 ranks, heavy contention); alpha2f.x 0.4 s.
- Wall time per q (4352 k + k+q points for a general q): about 1.7 min at nk=16 under contention. By k-point scaling
  I estimate about 3× that per q at nk=24 and about 8× at nk=32.

## 4. rtsc parseability

**No converter is needed**, so `bench_tools/tetra_to_two_column.py` was not written. `al.a2F.dat` goes directly through:
```
rtsc certify al.a2F.dat --format two_column --units Ry --rule bins --expected-lambda <alpha2f "lambda :" from a2F>
```
- **Unit:** ω column in Ry; α²F is dimensionless. The header line starts with `#` and is skipped. The DOS/PDOS columns
  are ignored with `--columns 0,1`.
- **Binning rule:** `bins`. The grid is uniform, starting at ΔE, so the rtsc Voronoi cells all have width ΔE, including
  the end cells. The rule then puts an atom at ω_i with mass α²F_i·ΔE. This is exactly alpha2f.x's own λ sum.
  - Check: rtsc float λ = 0.37941609898542084; alpha2f.x printed 0.37941609898542056.
- **Caveat:** α²F_i are point values of a tetrahedron DOS, not bin averages. `bins` is a quadrature of the continuous
  α²F, the same one QE uses. Trapezoid differs by 1.7·10⁻⁴ in λ.
- **Completeness:** a two-column file has no end marker. I supplied `--expected-lambda` from alpha2f.x stdout, which
  removes the `A-FILE-COMPLETE(unverified)` block. The rtsc tolerance is 10⁻³ relative.
- Recheck input: the spectra JSON for `rtsc recheck --spectra` was rebuilt exactly as `cmd_certify` builds it (bins,
  envelope 0.10). The script is `.../scratchpad/dump_spectra.py`. The bundle was not edited.

## 5. Certification (nk=16 spectrum)

All runs use the defaults omega_c = 10·ω_max = 4699.9 K and envelope e = 0.10, in bracket mode.

| run | status | T bracket [K] | recheck |
|---|---|---|---|
| T* = 300 K, mu_hi 0.15 | **REJECTED_PROVEN** [PROVED] | [0.30902, 8.31827] | RECHECK OK (same status and bracket) |
| T* = 1.2 K, mu_hi 0.10 | UNDECIDED | [0.81129, 8.31827] | RECHECK OK |
| T* = 1.2 K, mu_hi 0.15 | UNDECIDED | [0.30902, 8.31827] | RECHECK OK |

- Verdict files: `~/rtsc_work/tetra_Al/q6s_opt/elph_k16/cert/verdict_*.json` and `cert/spectra.json`.
- The T* = 300 K rejection carries rule ids R-UNIV, R-CELLS, R-MU, R-ODD, C-RAYLEIGH and CLASSIFY-S5.
- Every verdict carries the standard assumptions: M-ISO-MIGDAL, M-ISO-SINGLEBAND, U-MU, CUTOFF-BLOCK, E-ENV(0.1) and
  DISCRETISATION(bins).
- At 1.2 K both brackets contain T* (Al experiment: 1.18 K), so UNDECIDED is the expected honest answer.

Smearing comparison: I started certifying the smearing a2F.dos2, a2F.dos4 and a2F.dos10 (copies in
`~/rtsc_work/tetra_Al/smear_ref/`). Six verdict JSONs were written before the second WSL failure, but they were never read
back or rechecked. They are not reported here.

## 6. Comparison with converged smearing [FLOAT]

Smearing values are from a2F.dos1-10 of the existing run (k 16³, fit 32³, q 8³ plus Fourier interpolation to 32³),
recomputed with the bins rule.

| quantity | tetra (nk 16³, q 6³ shifted) | smearing σ = 0.02 Ry | smearing range σ = 0.01-0.05 | tetra − smear(0.02) |
|---|---|---|---|---|
| λ | 0.3794 (direct 0.3803) | 0.4116 | 0.4113-0.4152 | **−7.8 %** |
| ω_log | 305.3 K (direct 303.6) | 308.0 K | 304.3-308.0 K | −0.9 % |
| A | 61.47 K | 67.18 K | 67.07-67.53 K | **−8.5 %** |

**Shape.** I histogrammed both α²F onto 10 K bins, normalised each, and took the L1 distance to the σ = 0.02 histogram:

| spectrum | L1 distance to σ = 0.02 |
|---|---|
| tetra nk=16 | 0.175 |
| smearing σ = 0.015-0.05 | 0.006-0.017 |
| smearing σ = 0.01 | 0.029 |
| smearing σ = 0.005 | 0.094 (the under-resolved one with negative bins) |

**Verdict: no agreement within 2-3 % at the grids completed so far.** ω_log agrees within 1 %, so the phonon frequency
distribution is consistent. λ and A are 8 % lower and the shape differs noticeably. Candidate causes, not yet separated:

1. **e-ph k convergence of the tetrahedron double delta at nk=16.** This is the main reason for the nk=24/32 runs.
2. **q integration.** It uses only 28 irreducible shifted 6³ points, with no Fourier interpolation, against 8³
   interpolated to 32³ for smearing. The shape metric is the most sensitive to this.
3. **DFPT occupations.** The dvscf and phonons come from opt-tetra on a 16³ grid, against MV 0.02 Ry.
4. **λ_qν printed to 4 decimals.** This is a minor effect.

For reference only, with different settings: the official QE tetra_example (US pseudo, 40 Ry, q 4³ shifted, nk 16)
gives λ = 0.399.

## 7. Incidents and resources

- WSL went down twice (about 21:07 SIGBUS, then I/O errors around 21:30). The coordinator traced this to a full host
  C: drive: the WSL vhdx could not grow.
- My runs contributed to that disk use. Each tetra DFPT outdir is about 17 GB (28 q × ~0.9 GB of wavefunctions), and the
  stage-2 script copied it per k grid. Deletable without losing results: `q6s_opt/elph_k16/tmp`, `q6s_opt/elph_k24/`
  (unfinished) and `q6s_lin/tmp`. `q6s_opt/tmp` must be kept for nk=24/32.
- **Relaunch plan, when WSL is back and the coordinator allows it:** run nk=24, then nk=32, **in place** in `q6s_opt`
  with no tmp copy. Estimated about 2.5-3 h and 6-7 h on 4 ranks. Then rerun `analyze.py`, `t_cert.sh` and this table.
  - Prepared script: `scratchpad/t_stage2_inplace.sh`.
  - **Disk caveat:** ph.x keeps a per-q outdir `tmp/_ph0/al.q_N/` holding that q's nscf wavefunctions. At nk=32 each
    could be about 8× the 0.9 GB seen at 16³, which is up to roughly 200 GB over 28 q if nothing is removed.
  - The script runs a guard that deletes the wavefunction files of q points whose `.elph.` file already exists. The
    wfc filename pattern is **unverified**: check it in `al.q_N/` before relaunching.
- Scripts are in `G:\ClaudeWorkspace\RTSC_Project\run_scripts\scratchpad_2026-10-27\`:
  `t_stage1.sh`, `t_resume_opt.sh`, `t_stage2.sh`, `analyze.py`, `dump_spectra.py`, `t_cert.sh`, `t_linbug.sh` and
  `t_linbt.sh`.
- Nothing in `~/rtsc_work/rtsc_bundle` or in the existing Al run directory was modified. Files were only read or copied.
