# Step 2: pseudopotentials, cutoffs and rtsc configs for the Pb, Nb, MgB2 and H3S benchmarks

Date: 2026-10-27. QE 7.6 (conda-forge `qe-7.6-mpi_openmpi_hdb26816_104`, OpenMPI 5.0.11, gfortran 16.2, OpenBLAS 0.3.34, FFTW 3.3.11) in micromamba env `qe`.

**Every number in this report comes from a convergence, timing or cross-check test.** Small k grids, unrelaxed geometries and ≤2 MPI ranks were used. None of them is a physical result. No benchmark run was started, and the Al run was not touched.

---

## 1. ph.x e-ph segfault with PAW **and ultrasoft** pseudopotentials

### 1.1 Result
QE 7.6 `ph.x` crashes with SIGSEGV whenever `electron_phonon` is set and any species is ultrasoft or PAW. Both `'interpolated'` and `'simple'` crash. With the same pseudos and no `electron_phonon`, the phonon run completes. Norm-conserving pseudos work.

All runs below were serial (1 MPI rank, OMP_NUM_THREADS=1), with ldisp and nq=1 (Γ only) unless noted. Run directories are in `~/rtsc_work/crashchar/`.

| case | pseudo (type) | electron_phonon | result |
|---|---|---|---|
| Al fcc | SSSP Al.pbe-n-kjpaw_psl.1.0.0 (PAW) | none | **JOB DONE** |
| Al fcc | same PAW | `simple` | **SIGSEGV** |
| Al fcc | same PAW | `interpolated` | **SIGSEGV** |
| Al fcc, q 2³ | same PAW | `interpolated` | **SIGSEGV** (at first q) |
| Al fcc, q 2³ | same PAW, `start_q=2,last_q=2` (Γ skipped) | `interpolated` | **SIGSEGV** (not Γ-specific) |
| Al fcc | same PAW, two-step: phonons OK, then `trans=.false.` + e-ph | `interpolated` | phonon step OK, e-ph step **SIGSEGV** |
| Al fcc | same PAW, no `fildvscf` | `interpolated` | clean stop "El-ph needs a DeltaVscf file" (expected) |
| hcp Mg | GBRV mg_pbe_v1.4.uspp (**US**) | `interpolated` | **SIGSEGV** |
| hcp Mg | Mg.pbe-n-kjpaw_psl.0.3.0 (PAW) | `interpolated` | **SIGSEGV** |
| MgB2 | GBRV mg_pbe_v1.4 + B_pbe_v1.01 (**US**) | `interpolated` | **SIGSEGV** |
| bcc S (toy) | GBRV s_pbe_v1.4.uspp (**US**) | `interpolated` | **SIGSEGV** |
| fcc Pb | SSSP Pb.pbe-dn-kjpaw_psl.0.2.2 (PAW) | `interpolated` | **SIGSEGV** |
| bcc Nb | SSSP Nb.pbe-spn-kjpaw_psl.0.3.0 (PAW) | `interpolated` | **SIGSEGV** |
| Al fcc, q 2³ | SG15 Al_ONCV_PBE-1.2 (NC) | `interpolated` | **JOB DONE**, 90 λ lines |

The earlier finding also holds: `Al.pz-vbc` (NC) and SG15 run fine serially and with 8 MPI ranks.

### 1.2 Crash site
The `ph.x` binary is not stripped, so the gfortran backtrace can be mapped to symbols. I took the frame offsets relative to `_start` and resolved them with `nm`. Every crash, US or PAW, solve_linter path or `trans=.false.` path, has the same innermost frames:

```
adddvscf_ +0x3d0        (LR_Modules/adddvscf.f90, wrapper: adddvscf(ipert, ik, time_reversed))
elphel_   +0x1d29       (PHonon/PH/elphon.f90)
solve_linter_ / phqscf_ / do_phonon_     (trans=.true.)
elphon_   / do_phonon_                   (trans=.false.)
```

`elphel` calls `adddvscf` only `IF (okvan)`, which is the US/PAW augmentation term `int3 * becp1 * vkb`. This explains why NC pseudos are unaffected. In current develop, `elphel` calls it as `adddvscf(ipert, ik, isolv==2)`, and the wrapper then picks `becp1(ik)` or `becpt(ik)`. Likely causes are that `becp1`/`int3` are not allocated or valid in the elphel context, or a mismatch in the new time-reversal argument. I have not confirmed either, because there are no debug symbols and no source build. The crash happens at the start of elphel, after the DFPT self-consistency for the first representation has converged.

### 1.3 Is it known?
I searched the QE GitLab issues through the API for "elphon", "electron-phonon", "PAW phonon" and "segmentation ph.x", and searched the QE-users mailing list. Nothing matches.
- The closest reports are different failures:
  - 2022: an `a2Fsave` EOF with `electron_phonon='interpolated'` for fcc Co ([QE-users 2022-May](https://lists.quantum-espresso.org/pipermail/users/2022-May/049051.html)).
  - 2025-11: the GPU/OpenACC 7.5 elphon `a2Fsave` read failure; Giannozzi says e-ph is not ported to GPU ([QE-users 2025-Nov](https://lists.quantum-espresso.org/pipermail/users/2025-November/052783.html)).
- GitLab #673 (an elphon MPI_RECV signature problem) is unrelated.
- **It appears to be unreported.** I could not test the develop branch, because there is no source build here.

### 1.4 Minimal reproducer (about 20 s, serial)
```
# scf.in  (fit.in identical except: la2F=.true. in &system and K_POINTS 16 16 16 0 0 0)
&control  calculation='scf', prefix='t', outdir='./tmp/', pseudo_dir='<SSSP_1.3.0_PBE_precision>' /
&system   ibrav=2, celldm(1)=7.65, nat=1, ntyp=1, ecutwfc=30, ecutrho=240,
          occupations='smearing', smearing='mv', degauss=0.02 /
&electrons conv_thr=1e-10 /
ATOMIC_SPECIES
 Al 26.98 Al.pbe-n-kjpaw_psl.1.0.0.UPF
ATOMIC_POSITIONS crystal
 Al 0 0 0
K_POINTS automatic
 8 8 8 0 0 0
# ph.in
ph
&inputph
 prefix='t', outdir='./tmp/', tr2_ph=1e-14, fildyn='dyn', fildvscf='dvscf',
 electron_phonon='interpolated', el_ph_sigma=0.005, el_ph_nsigma=10,
 trans=.true., ldisp=.true., nq1=1, nq2=1, nq3=1
/
# run: pw.x -in fit.in; pw.x -in scf.in; ph.x -in ph.in   -> SIGSEGV in adddvscf <- elphel
# control 1: delete the electron_phonon line                -> JOB DONE
# control 2: Al_ONCV_PBE-1.2.upf (SG15), ecutrho=4*ecutwfc   -> JOB DONE with lambda output
```
Generator script: `/mnt/g/ClaudeWorkspace/RTSC_Project/run_scripts/scratchpad_2026-10-27/crash.sh` (cases A/B).

### 1.5 Consequence and workaround
Use norm-conserving pseudos for all `ph.x` e-ph work in QE 7.6. All benchmark configs below use SG15 ONCV 1.2 (scalar-relativistic), which is also the safest family for EPW.

Possible routes for US/PAW, none of them tested:
- an older QE (for example 7.2) in a separate env;
- EPW, which has its own US/PAW support;
- reporting the bug upstream with the reproducer above.

---

## 2. Pseudopotential and cutoff choice

### 2.1 Choice
| element | file (in `~/pseudo/sg15`) | Z_val | max rc (bohr, l) |
|---|---|---|---|
| Pb | Pb_ONCV_PBE-1.2.upf (SR) | 14 (5d 6s 6p) | 3.60 (p,f) |
| Pb (SOC option) | Pb_ONCV_PBE_FR-1.0.upf | 14 | 3.60 |
| Nb | Nb_ONCV_PBE-1.2.upf | 13 (4s4p semicore) | 1.66 (s) |
| Mg | Mg_ONCV_PBE-1.2.upf | 10 (2s2p semicore) | 2.49 (p) |
| B | B_ONCV_PBE-1.2.upf | 3 | 1.90 (p) |
| S | S_ONCV_PBE-1.2.upf | 6 | 2.13 (s), 1.52 (p) |
| H | H_ONCV_PBE-1.2.upf | 1 | 1.14 (s) |

### 2.2 pw.x ecutwfc convergence tests
Settings: ecutrho = 4·ecutwfc, MV smearing 0.02 Ry, 2 ranks, experimental or approximate geometry (unrelaxed). The reference is the highest cutoff tested. Raw data: `~/rtsc_work/conv/results.txt`.

| system (k grid) | ecut | ΔE (mRy/atom) vs ref | ΔP (kbar) vs ref |
|---|---|---|---|
| **Pb** fcc a=4.95 Å (8³), ref 120 Ry | 30 | 0.50 | 5.7 |
| | 40 | 0.004 | 0.05 |
| | **50** | 0.001 | 0.05 |
| **Pb FR+SOC** (8³), ref 120 Ry | 40 | 0.005 | 0.01 |
| | 60 | <0.001 | 0.00 |
| **Nb** bcc a=3.30 Å (8³), ref 200 Ry | 80 | 2.1 | −3.7 |
| | 100 | 1.6 | −0.9 |
| | **120** | 1.0 | −4.0 |
| | 140 | 0.30 | −2.8 |
| | 160 | 0.04 | −0.6 |
| **MgB2** a=3.086, c=3.524 Å (8×8×6), ref 120 Ry | 70 | 0.35 | −2.6 |
| | 80 | 0.15 | −1.7 |
| | **90** | 0.05 | −0.8 |
| | 100 | 0.01 | −0.2 |
| **H3S** Im-3m a=3.09 Å (8³, P ≈ 138 GPa), ref 140 Ry | 70 | 0.21 | −3.7 |
| | 80 | 0.08 | −1.3 |
| | **90** | 0.03 | −0.3 |
| | 100 | 0.02 | +0.1 |

**Chosen cutoffs (Ry):**
- Pb 50 / 200. A margin above the 40 Ry already converged.
- Nb 120 / 480.
- MgB2 90 / 360.
- H3S 90 / 360.

The SOC variant of Pb converges at the same 50 Ry.

**Nb caveat.** SG15 Nb converges slowly in energy: 1 mRy/atom needs 120 Ry. At k 8³ the pressure scatters non-monotonically by about ±2 kbar between 80 and 160 Ry, so the "~1 kbar" target is only met at ≥160 Ry. A 2–4 kbar pressure error corresponds to Δa/a ≈ 0.05% (B ≈ 170 GPa), which is negligible for the phonons. I judged 120 Ry an adequate trade-off. 160 Ry would cost about 1.5× as much.

### 2.3 H3S high-pressure suitability (core-radius overlap)
Nearest-neighbour distances in Im-3m H3S, compared with the SG15 rc sums:

| a (Å) | P (this test, SG15) | S–H | H–H | S–S |
|---|---|---|---|---|
| 3.09 | ≈138 GPa | 1.545 Å = 2.92 bohr | 2.185 Å = 4.13 bohr | 2.676 Å = 5.06 bohr |
| 2.98 | ≈205 GPa | 1.490 Å = 2.82 bohr | 2.107 Å = 3.98 bohr | 2.581 Å = 4.88 bohr |

- rc(S)+rc(H) = 2.13 + 1.14 = **3.27 bohr > S–H**, so the S–H spheres overlap by 0.35–0.45 bohr (12–16%).
  - The overlap comes from the S s-channel (rc 2.13). With the p-channel rc, the sum is 2.66 bohr, which is below the S–H distance.
- H–H (2·1.14 = 2.27 bohr) and S–S (2·2.13 = 4.26 bohr) do not overlap.

I cross-checked SG15 against a harder set: SSSP-efficiency S GBRV US (rc 1.5) + H pslibrary US (rc 1.0), whose rc sum of 2.5 bohr means no overlap. The US set ran at 80/800 Ry, and its cutoff was not convergence-tested.
- **Pressure** (k 8³): 138.2 vs 136.6 GPa at a=3.09 Å (+1.2%); 205.1 vs 202.9 GPa at a=2.98 Å (+1.1%).
- **Γ phonons** at a=2.98 Å (k 12³, phonon only, no e-ph, so the US run works):

  | | T1u modes (cm⁻¹) |
  |---|---|
  | SG15 | 556, 1200, 1654 |
  | US | 524, 1183, 1655 |
  | difference | +6%, +1.4%, 0.1% |

**Assessment.** SG15 is usable for H3S at 150–200 GPa. The pressure offset is about 1% (≈2 GPa). The H-stretching modes agree to ≤1.4%. The low S–H bending/optical mode is about 6% stiffer with SG15, which means a pseudo uncertainty of a few % in ω_log/λ for H3S. Record this as a systematic in the benchmark report.

No harder NC S pseudo is available offline: PseudoDojo is unreachable and oncvpsp is not installed. The relaxed lattice at 200 GPa will sit about 2 GPa "lower" than with US.

### 2.4 Pb SOC
Scalar-relativistic DFPT underestimates λ for Pb: about 1.2 without SOC versus about 1.5–1.6 with SOC (Heid et al., PRB 81, 174527, 2010).
- The default config is SR (SG15 1.2).
- An SOC variant is documented as a comment in `bench/Pb/config.toml`: `files = { Pb = "Pb_ONCV_PBE_FR-1.0.upf" }` plus `extra_pw_system = { noncolin = true, lspinorb = true }`.
- The FR pseudo was convergence-tested at the same cutoff. In the pw.x tests it costs about 3× SR per scf.
- The SOC e-ph path through rtsc (relax, ph, q2r, matdyn) has **not** been tested.

---

## 3. Configs (written, not run)

| material | config | structure (primitive cell) | pressure | ecut (Ry) | k_scf / k_fit | q (irr) | plan-only |
|---|---|---|---|---|---|---|---|
| Pb | `~/rtsc_work/bench/Pb/config.toml` | Pb_fcc.vasp, a=4.95 Å | 0 | 50/200 | 16³ / 32³ | 8³ (29) | OK: 1 input, 1 to_compute |
| Nb | `~/rtsc_work/bench/Nb/config.toml` | Nb_bcc.vasp, a=3.30 Å | 0 | 120/480 | 16³ / 32³ | 8³ (29) | OK |
| MgB2 | `~/rtsc_work/bench/MgB2/config.toml` | MgB2_P6mmm.vasp, a=3.086, c=3.524 Å | 0 | 90/360 | 18×18×12 / 36×36×24 | 6×6×4 (21) | OK |
| H3S | `~/rtsc_work/bench/H3S/config.toml` | H3S_Im-3m.vasp, **4-atom primitive bcc**, a=3.09 Å | 200 GPa | 90/360 | 18³ / 36³ | 6³ (16) | OK |

Common settings, all as in the Al config:
- 8 ranks with `-nk 8` and the explicit `mpirun`/`bin_dir`.
- `do_relax=true`, MV smearing 0.02 Ry.
- `el_ph_sigma` 0.005 Ry × 10 smearings, `tr2_ph=1e-16`, `k_matdyn = k_fit`.
- Filters: hull off, metallic/magnetic/dynamic on.
- Each `structures/manifest.toml` sets `material_id` and `pressure_GPa`.

The irreducible q counts come from spglib, which reproduces QE's 29 for fcc 8³.

Notes:
- **plan-only** (`rtsc run -c config.toml --plan-only` from each bench dir) exits 0. `run/inputs.json` lists one input to compute per material, with the correct species and pressure. The CLI summary prints "inputs: 0" because it counts only records that have verdicts. Plan-only left `run/inputs.json` and a verdict-less `record.json` in each `run/`. Resume ignores verdict-less records.
- I also checked separately that `EPhSettings.from_config`/`Launcher.from_config` accept every `[dfpt]` block and that every pseudo file exists.
- **H3S**: rtsc writes `ibrav=0` with the cell exactly as given and does not reduce to the primitive cell. The 8-atom conventional prototype would double the atoms and fold the q grid, so I wrote a 4-atom primitive cell. I chose 200 GPa rather than 150 GPa because harmonic PBE Im-3m is dynamically stable only above about 180 GPa (Errea et al., 2016). At 150 GPa the dynamic filter would likely reject it. Expect anharmonic and non-adiabatic validity flags.
- **MgB2**: the dense k_fit targets the small σ-band cylinders. The isotropic result should underestimate Tc (λ ≈ 0.7). The certify step must flag MULTIBAND_RISK. A cheaper fallback is k_scf 12×12×8 / k_fit 24×24×16, at about 1/3 of the cost but with a less converged σ contribution.
- **Nb**: non-magnetic, so the magnetic filter is left on as a no-op check.

---

## 4. Cost and RAM estimates (8 ranks, this machine)

**Method.** I ran ph.x e-ph timing tests at the production cutoffs (k 8³ / fit 16³, q 2³; MgB2 6×6×4 / 12×12×8; H3S 6³ / 12³; 2 ranks). They were calibrated against the live Al production run (8 ranks, k 16³, q 8³): 14 of its 29 q took 21.2 min, which projects to about 44 min for ph.x. For each material,

per-q production time ≈ (test time per q) × (k_scf ratio, production/test) × 1.09.

The factor 1.09 is the Al production/test calibration. pw.x relax, scf_fit, q2r and matdyn are small by comparison: an extrapolated scf_fit is <1 min for every material.

| material | test ph time per q (2 ranks) | k ratio | est. per q (8 ranks) | irr q | **est. ph.x wall** | test peak RSS/rank | est. RAM (8 ranks) | est. disk (outdir) |
|---|---|---|---|---|---|---|---|---|
| Al (calibration) | 10.4 s | 8 | ≈91 s (measured) | 29 | ≈45 min (projected) | 104 MB | 0.9 GB (measured ~110 MB/rank) | ~12 GB |
| Pb | 16.9 s | 8 | ≈2.5 min | 29 | **≈1.2 h** | 112 MB | ~1.5 GB | ~17 GB |
| Nb | 31.3 s | 8 | ≈4.5 min | 29 | **≈2.2 h** | 149 MB | ~2.5 GB | ~35 GB |
| MgB2 | 46.3 s | 27 | ≈23 min | 21 | **≈8 h** | 132 MB | ~3 GB | ~55 GB |
| H3S | 16.1 s | 27 | ≈8 min | 16 | **≈2.1 h** | 114 MB | ~1.5 GB | ~11 GB |

Uncertainty is about ×2. The test q points (2³) are high-symmetry points. The production grids contain more low-symmetry q, which need more k points, especially for hexagonal MgB2. The Al calibration absorbs part of this effect. The SOC variant of Pb would be about 3–4× the SR cost.

**RAM.** ph.x keeps wavefunctions in outdir buffers, so RSS per rank is small: about 110 MB for Al production. All estimates are below 4 GB in total, far under the 15 GB WSL limit.

**Disk.** Disk is the larger footprint: Al uses 5.8 GB at half way. The disk estimates scale Al by npw × nbnd × (k per q) × q count. There are 922 GB free.

**Run order.** Pb and Nb (a few hours) run first, then H3S, then MgB2 (overnight). Run them sequentially, not alongside the Al run: each uses all 8 physical cores.

---

## 5. Files
- Configs: `~/rtsc_work/bench/{Pb,Nb,MgB2,H3S}/config.toml`, `structures/*.vasp`, `structures/manifest.toml`.
- Test data:
  - `~/rtsc_work/crashchar/` (crash cases; `summary.txt`)
  - `~/rtsc_work/conv/results.txt` (cutoffs, H3S transferability)
  - `~/rtsc_work/phtime/results.txt` (timings, H3S Γ cross-check)
- Scripts (scratchpad `/mnt/g/ClaudeWorkspace/RTSC_Project/run_scripts/scratchpad_2026-10-27/`): `crash.sh`, `crash2.sh`, `sym.sh`/`sym2.sh` (backtrace symbolisation), `conv.sh`, `phtime.sh`, `h3sgam.sh`, `irq.py`, `mkcfg.py`.
