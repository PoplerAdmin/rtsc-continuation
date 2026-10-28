# Step 3: MatterSim ML prescreen (H3S check + MatterGen Mg-B set)

Date: 2026-10-27. Machine: WSL2 Ubuntu, RTX 4090 (driver 610.88).

> **S2 statement.** Everything in this report is a **HEURISTIC** (ARCH S2). ML hull energies, ML phonons and
> `ML_SUGGEST_DEFER` notes may only **reorder candidates or DEFER them with a recorded reason**. Nothing
> here rejects (EXCLUDED_*) or certifies any structure. All numbers are real outputs of MatterSim runs
> unless explicitly labelled "side check". No `rtsc.certify` or other proof-relevant code was touched, and
> no file under `~/rtsc_work/rtsc_bundle` was edited.

## 1. The "Failed to find C compiler" fix

**Root cause.** The error did not come from `torch.compile`. MatterSim 1.2.5 only uses `torch.jit.script`.
torch 2.14 has an eager-op override registry (`torch/_native/registry.py`, `eager_router`). It sends some
CUDA ops to Triton kernels. Here the op was `bmm_outer_product`, called inside the autograd backward that
MatterSim uses for forces and stress (`potential.py:847`, `torch.autograd.grad`). Before its first launch,
Triton builds a small C stub (`triton/backends/nvidia/driver.py`, `CudaUtils`), and that build needs a host
C compiler. WSL has none.

**What was tried:**

| Setting | Result |
|---|---|
| `TORCHDYNAMO_DISABLE=1 TORCH_COMPILE_DISABLE=1` | still fails, because dynamo/inductor are not involved |
| **`TORCH_DISABLE_NATIVE_JIT=1`** (torch's own global kill switch, `torch/_native/common_utils.py`) | **works**: the ops fall back to the standard ATen CUDA kernels |

**Fix used:** `export TORCH_DISABLE_NATIVE_JIT=1`, set in the run wrapper
`.../scratchpad/run_prescreen.sh`. Nothing was installed. The micromamba compiler env (option ii) was not
needed.

**Results are unchanged.** I compared rattled H3S on CUDA (ATen fallback) against CPU with the 5M model:
|dE| = 0, max |dF| = 1.1e-6 eV/Å, max |dσ| = 1.2e-7 eV/Å³. The Triton path itself could not be compared,
since it cannot run without a compiler. The fallback is the stock ATen implementation.

## 2. Versions, settings, wall times

Python 3.12 (`~/rtsc_work/.venv`), torch 2.14.1+cu130, triton 3.8.0, mattersim 1.2.5, ase 3.29.0,
phonopy 4.8.6. These are the rtsc adapter defaults: FIRE, FrechetCellFilter, fmax 0.02 eV/Å, ≤400 steps,
no symmetry constraint. Phonons use finite displacements of 0.01 Å, an auto supercell (≥7 Å, ≤160 atoms),
an 8×8×8 mesh, and an imaginary tolerance of 15 cm⁻¹. The `ML_SUGGEST_DEFER` thresholds are E_hull > 0.2 eV/atom or ML-unstable.

| Run | Wall time | Notes |
|---|---|---|
| `ms_prescreen.py h3s` (1M and 5M × 150/200 GPa) | 26 s total (each run 2 to 10 s) | OMP_NUM_THREADS=2, ~1.8 cores, 2.0 GB RSS |
| `ms_prescreen.py mg` (5M, 33 structures × 0/150 GPa) | 7 min 48 s total (229 s at 0 GPa, 237 s at 150 GPa) | ~1.75 cores, 2.1 GB RSS |

All 70 relaxations converged. There were no ML failures.

## 3. H3S Im-3m (KNOWN reference)

| Model | P (GPa) | relaxed a_conv (Å) | V/atom (Å³) | ML min freq (cm⁻¹) | q of min | ML stable | defer |
|---|---|---|---|---|---|---|---|
| 1M | 150 | 3.0566 | 3.570 | **−450.4** (rerun −451.0) | H (½,−½,½) | False | DEFER |
| 1M | 200 | 2.9864 | 3.329 | **−302.7** | H | False | DEFER |
| 5M | 150 | 3.0423 | 3.520 | **+120.7** (acoustic, near Γ) | (⅛,0,0) | True | – |
| 5M | 200 | 2.9666 | 3.264 | **+130.8** (acoustic, near Γ) | (0,0,⅛) | True | – |

- The **5M model does not reproduce the 1M instability.** It finds Im-3m ML-stable at both 150 and 200 GPa.
  Its lowest optical frequencies are at H: +256.6 cm⁻¹ at 150 GPa and +321.1 cm⁻¹ at 200 GPa. The 1M
  instability at H is model-specific. At 200 GPa it also conflicts with the expected harmonic DFT
  picture (Im-3m dynamically stable there), so 1M should not drive deferral of H3S-like candidates.
- **Supercell check (side check, `h3s_supercell_xcheck.json`).** The auto supercell for the 8-atom
  conventional cell was capped to an anisotropic 2×3×3 (6.1 Å). I recomputed with isotropic 2×2×2, 3×3×3
  and 4×4×4 conventional supercells. The 1M minimum stayed at −448 to −450 cm⁻¹ (150 GPa) and about −302
  cm⁻¹ (200 GPa), always at H. At 150 GPa the N point is also imaginary (−280 cm⁻¹). So the 1M result is
  converged with respect to supercell size and is not an artifact. The q coordinates in the notes are in
  phonopy's auto-primitive (bcc) basis, so the handoff's "H point" label is correct.
- The H3S E_hull = 0.000 is **meaningless**. The only references are rtsc's crude fallbacks (H "dimer-sc",
  S "fcc(covalent radius)").

## 4. Experimental MgB2 sanity check (5M)

| | a (Å) | c (Å) | V/atom (Å³) | ML min freq (cm⁻¹) | ML stable |
|---|---|---|---|---|---|
| Experiment | 3.086 | 3.524 | 9.688 | – | – |
| 5M, 0 GPa | 3.0712 (−0.48 %) | 3.5283 (+0.12 %) | 9.607 (−0.8 %) | +63 | **True** |
| 5M, 150 GPa | 2.7298 | 2.8433 | 6.116 | +137 | True |

Good agreement with experiment. MgB2 is ML-stable. It is on the hull by construction, because it is itself a hull reference.

## 5. MatterGen Mg-B set (32 real MatterGen structures plus exp. MgB2), 5M

Labels: `cli#i` = `MgB_raw` (CLI run). `s7#i` = adapter run `G-AMB_B-Mg_sgany_P0_n16_s7_*`. All
MatterGen structures are HYPOTHETICAL. H is the enthalpy E+PV per atom.

### 5a. 0 GPa

| label | N | formula | sg in | sg relaxed | E_hull ML (eV/at) | min ω (cm⁻¹) | ML stable | lightHF | ML_SUGGEST_DEFER | H (eV/at) |
|---|---|---|---|---|---|---|---|---|---|---|
| MgB2-exp(KNOWN) | 3 | MgB2 | P6/mmm | P6/mmm | 0.000 | 63 | True | 0.41 | | −5.0960 |
| cli#0 | 8 | Mg4B4 | P6_3/mmc | P6_3/mmc | 0.000 | 7 | True | 0.27 | | −4.2170 |
| cli#1 | 10 | Mg6B4 | P4/mbm | P4/mbm | 0.280 | 30 | True | 0.07 | hull | −3.3770 |
| cli#2 | 12 | Mg7B5 | Pm | Pm | 0.212 | 8 | True | 0.14 | hull | −3.5350 |
| cli#3 | 4 | Mg2B2 | P-3m1 | P-3m1 | 0.322 | 16 | True | 0.19 | hull | −3.8750 |
| cli#4 | 9 | MgB8 | Immm | C2/m | 0.000* | −138 | False | 0.45 | unstable | −5.7550 |
| cli#5 | 18 | Mg4B14 | P3m1 | P-6m2 | 0.000* | −47 | False | 0.42 | unstable | −5.4770 |
| cli#6 | 9 | Mg5B4 | C2/m | C2/m | 0.115 | 16 | True | 0.19 | | −3.7820 |
| cli#7 | 4 | Mg2B2 | P-3m1 | P-3m1 | 0.322 | 16 | True | 0.19 | hull | −3.8750 |
| cli#8 | 12 | Mg8B4 | P2_1/m | P2_1/m | 0.108 | −33 | False | 0.13 | unstable | −3.1890 |
| cli#9 | 9 | Mg3B6 | Pm | Pm | 0.215 | 29 | True | 0.39 | hull | −4.8800 |
| cli#10 | 4 | Mg2B2 | Cmcm | Cmcm | 0.089 | 48 | True | 0.21 | | −4.1080 |
| cli#11 | 16 | Mg4B12 | P3m1 | C2 | 0.000* | 9 | True | 0.52 | | −5.4960 |
| cli#12 | 8 | Mg3B5 | P3m1 | Cm | 0.218 | 11 | True | 0.33 | hull | −4.6530 |
| cli#13 | 5 | Mg3B2 | Cmmm | Cmmm | 0.082 | −85 | False | 0.18 | unstable | −3.5750 |
| cli#14 | 8 | Mg4B4 | Pnma | Pnma | 0.094 | 20 | True | 0.18 | | −4.1030 |
| cli#15 | 8 | Mg3B5 | Amm2 | Amm2 | 0.226 | −58 | False | 0.28 | hull + unstable | −4.6450 |
| s7#0 | 6 | Mg2B4 | P6/mmm | P6/mmm | 0.000 | 42 | True | 0.42 | | −5.0960 |
| s7#1 | 4 | Mg3B | P-3m1 | P-3m1 | 0.335 | −39 | False | 0.00 | hull + unstable | −2.5130 |
| s7#2 | 6 | Mg2B4 | C2/m | C2/m | 0.245 | −172 | False | 0.35 | hull + unstable | −4.8510 |
| s7#3 | 12 | Mg4B8 | P6/mmm | P6/mmm | 0.000 | 37 | True | 0.43 | | −5.0960 |
| s7#4 | 9 | Mg3B6 | C2/m | C2/m | 0.229 | −44 | False | 0.30 | hull + unstable | −4.8660 |
| s7#5 | 10 | Mg6B4 | P-3m1 | P-3m1 | 0.265 | 6 | True | 0.12 | hull | −3.3930 |
| s7#6 | 8 | Mg4B4 | Cmcm | Cmcm | 0.089 | 39 | True | 0.20 | | −4.1080 |
| s7#7 | 20 | Mg2B18 | P-1 | P-1 | 0.000* | 28 | True | 0.52 | | −5.8660 |
| s7#8 | 9 | Mg4B5 | C2/m | P-1 | 0.240 | 19 | True | 0.25 | hull | −4.2570 |
| s7#9 | 8 | Mg3B5 | P3m1 | P3m1 | 0.420 | −308 | False | 0.37 | hull + unstable | −4.4510 |
| s7#10 | 4 | Mg2B2 | R-3m | R-3m | 0.000 | 16 | True | 0.27 | | −4.2170 |
| s7#11 | 17 | Mg6B11 | P1 | P1 | 0.299 | 23 | True | 0.33 | hull | −4.6910 |
| s7#12 | 12 | Mg2B10 | Imma | P1 | 0.000* | 41 | True | 0.54 | | −5.7910 |
| s7#13 | 12 | Mg6B6 | Cmcm | Cmcm | 0.089 | 19 | True | 0.20 | | −4.1080 |
| s7#14 | 12 | Mg7B5 | P4mm | Pmm2 | 0.325 | −82 | False | 0.08 | hull + unstable | −3.4230 |
| s7#15 | 4 | Mg2B2 | P3m1 | P3m1 | 0.156 | −44 | False | 0.30 | unstable | −4.0410 |

\* B-rich (x_B > 2/3). E_hull is set by the fcc-boron fallback and is not informative. See §6.

### 5b. 150 GPa

| label | formula | sg in | sg relaxed | E_hull ML | min ω | ML stable | lightHF | ML_SUGGEST_DEFER | H (eV/at) |
|---|---|---|---|---|---|---|---|---|---|
| MgB2-exp(KNOWN) | MgB2 | P6/mmm | P6/mmm | 0.000 | 137 | True | 0.62 | | 1.7810 |
| cli#0 | MgB | P6_3/mmc | P6_3/mmc | 0.091 | 31 | True | 0.45 | | 4.1550 |
| cli#1 | Mg3B2 | P4/mbm | P4/mbm | 0.162 | −256 | False | 0.33 | unstable | 5.5953 |
| cli#2 | Mg7B5 | Pm | Pm | 0.344 | 17 | True | 0.33 | hull | 5.5493 |
| cli#3 | MgB | P-3m1 | R-3m | 0.091 | 61 | True | 0.46 | | 4.1548 |
| cli#4 | MgB8 | Immm | P-1 | 0.000* | 65 | True | 0.76 | | −0.2219 |
| cli#5 | Mg2B7 | P3m1 | P-6m2 | 0.000* | 18 | True | 0.64 | | 0.9211 |
| cli#6 | Mg5B4 | C2/m | C2/m | 0.163 | −88 | False | 0.44 | unstable | 4.9879 |
| cli#7 | MgB | P-3m1 | P-3m1 | 0.872 | −618 | False | 0.34 | hull + unstable | 4.9367 |
| cli#8 | Mg2B | P2_1/m | P2_1/m | 0.214 | 41 | True | 0.28 | hull | 6.5616 |
| cli#9 | MgB2 | Pm | **P6/mmm** | 0.000 | 57 | True | 0.60 | | 1.7808 |
| cli#10 | MgB | Cmcm | Cmcm | 0.013 | 80 | True | 0.44 | | 4.0766 |
| cli#11 | MgB3 | P3m1 | C2 | 0.000* | 15 | True | 0.68 | | 0.9508 |
| cli#12 | Mg3B5 | P3m1 | P3m1 | 0.428 | −233 | False | 0.48 | hull + unstable | 2.7794 |
| cli#13 | Mg3B2 | Cmmm | Cmmm | 0.064 | 57 | True | 0.34 | | 5.4981 |
| cli#14 | MgB | Pnma | Pnma | 0.020 | 51 | True | 0.44 | | 4.0850 |
| cli#15 | Mg3B5 | Amm2 | Amm2 | 0.294 | 24 | True | 0.51 | hull | 2.6461 |
| s7#0 | MgB2 | P6/mmm | P6/mmm | 0.000 | 80 | True | 0.60 | | 1.7808 |
| s7#1 | Mg3B | P-3m1 | P-3m1 | 0.984 | −254 | False | 0.21 | hull + unstable | 8.4730 |
| s7#2 | MgB2 | C2/m | **Imma** | 0.005 | 85 | True | 0.60 | | 1.7850 |
| s7#3 | MgB2 | P6/mmm | P6/mmm | 0.000 | 69 | True | 0.61 | | 1.7811 |
| s7#4 | MgB2 | C2/m | C2/m | 0.337 | 64 | True | 0.60 | hull | 2.1177 |
| s7#5 | Mg3B2 | P-3m1 | R-3m | 0.140 | 20 | True | 0.36 | | 5.5749 |
| s7#6 | MgB | Cmcm | Cmcm | 0.013 | 70 | True | 0.44 | | 4.0775 |
| s7#7 | MgB9 | P-1 | P-1 | 0.000* | 46 | True | 0.76 | | −0.2730 |
| s7#8 | Mg4B5 | C2/m | C2/m | 0.508 | 30 | True | 0.42 | hull | 3.8106 |
| s7#9 | Mg3B5 | P3m1 | P3m1 | 0.993 | −349 | False | 0.45 | hull + unstable | 3.3441 |
| s7#10 | MgB | R-3m | R-3m | 0.091 | 62 | True | 0.45 | | 4.1548 |
| s7#11 | Mg6B11 | P1 | P1 | 0.274 | 43 | True | 0.55 | hull | 2.3227 |
| s7#12 | MgB5 | Imma | Imm2 | 0.021* | −150 | False | 0.64 | unstable | 0.6005 |
| s7#13 | MgB | Cmcm | Cmcm | 0.013 | 34 | True | 0.44 | | 4.0768 |
| s7#14 | Mg7B5 | P4mm | P4mm | 0.582 | −647 | False | 0.25 | hull + unstable | 5.7887 |
| s7#15 | MgB | P3m1 | P3m1 | 0.574 | 28 | True | 0.40 | hull | 4.6386 |

### 5c. Counts (the 32 MatterGen structures; exp. MgB2 is ML-stable and not deferred at either pressure)

| | 0 GPa | 150 GPa |
|---|---|---|
| ML-unstable (min ω < −15 cm⁻¹) | 11 | 8 |
| ML_SUGGEST_DEFER (any reason) | 20 | 15 |
| of which by hull (> 0.2 eV/at) / by imaginary modes / both | 15 / 11 / 6 | 12 / 8 / 5 |
| deferred at both pressures | 13 | |
| never deferred | 10: cli#0, cli#10, cli#11, cli#14, s7#0, s7#3, s7#6, s7#7, s7#10, s7#13 | |
| ML failures | 0 | 0 |

Of the 10 never deferred, 2 are MgB2 supercells, 6 are MgB polytypes, and 2 are B-rich phases whose hull
value is fallback-limited (§6).

### 5d. Known phases and duplicates

- **MgB2 P6/mmm.** s7#0 (6-atom, 2× cell) and s7#3 (**12-atom, 4× supercell**) are the known AlB2-type
  MgB2. They match the experimental reference in H to within 0.1 meV/atom. cli#9 (Mg3B6, Pm) relaxes to P6/mmm MgB2 **only at 150 GPa**. At 0 GPa it stays in a Pm local minimum 0.215 eV/atom up.
- s7#2 (MgB2, C2/m) relaxes to **Imma at 150 GPa**, 5 meV/atom above P6/mmm. A KHg2-type Imma high-pressure MgB2
  phase has been proposed in the literature. This is worth a DFT look, but the attribution is not verified here.
- **MgB4 (Pnma) and MgB7 (Imma): none.** No input has those compositions (the B-rich inputs are MgB3,
  Mg2B7, MgB5, MgB8, MgB9 and Mg6B11). s7#12 has MgB7's space group (Imma) but composition MgB5, and
  it lowers symmetry on relaxation (P1 at 0 GPa, Imm2 at 150 GPa).
- **MgB (1:1), 9 structures, is not a known stable phase.** Variants: P6_3/mmc, R-3m, P-3m1, Cmcm, Pnma and P3m1.
  At 0 GPa the 5M model puts cli#0 (P6_3/mmc) and s7#10 (R-3m), which have identical H, on the hull. They
  sit 20 meV/atom below the Mg + MgB2 tie line. These are probably stacking polytypes. Treat this as an ML
  result needing DFT, not as evidence of a stable MgB.
- **Duplicates (same formula, same relaxed space group, ΔH < 3 meV/at, ΔV < 1 %).**
  - 0 GPa: {MgB2-exp, s7#0, s7#3}, {cli#3, cli#7} (MgB P-3m1), {cli#10, s7#6, s7#13} (MgB Cmcm). That leaves 29 distinct of 33.
    cli#0 and s7#10 are degenerate in H (−4.2170) but have different space groups.
  - 150 GPa: {MgB2-exp, cli#9, s7#0, s7#3}, {cli#3, s7#10} (MgB R-3m, with cli#0 P6_3/mmc 0.2 meV/at away),
    {cli#10, s7#6, s7#13} (MgB Cmcm).

## 6. Hull caveat: E_hull is optimistic, strongly so for B-rich compositions

The ML hull uses only **hcp Mg** and the **experimental MgB2** structure (both relaxed with 5M at each P). For
boron it uses rtsc's crude elemental fallback: `ase.build.bulk("B")` fails, so rtsc uses
`fallback:fcc(covalent radius)`. No MgB4, MgB7, MgB12 or α-B12 is in the reference set. Every E_hull above is
therefore a **lower bound** (optimistic), which is the harmless direction for S2. But on the B-rich side
it is far too optimistic to defer anything.

**Side check** (not part of the prescreen JSONs; computed from the same 5M enthalpies): I relaxed α-B12 (R-3m,
12 atoms, literature geometry) with 5M. At 0 GPa its H is 1.51 eV/atom below the fcc-B fallback, and 0.56
eV/atom below at 150 GPa. With α-B as the B reference instead of fcc:

| B-rich structure | E_hull rtsc (fcc-B) 0 GPa | E_hull with α-B 0 GPa | rtsc 150 GPa | with α-B 150 GPa |
|---|---|---|---|---|
| cli#4 MgB8 | 0.000 | 0.394 | 0.000 | −0.027 |
| cli#5 Mg2B7 | 0.000 | 0.146 | 0.000 | 0.128 |
| cli#11 MgB3 | 0.000 | −0.005 | 0.000 | −0.089 |
| s7#7 MgB9 | 0.000 | 0.336 | 0.000 | 0.021 |
| s7#12 MgB5 | 0.000 | 0.095 | 0.021 | 0.302 |

These are still optimistic. MgB4 and MgB7 are missing, and α-B is not the stable boron phase at 150 GPa
(α-Ga-type B is). Recommendation: put α-B12, MgB4 (Pnma) and MgB7 (Imma) into `refs/`. The hull cache
must also be cleared when references change (see §7).

## 7. Suspected adapter issues (described only, not fixed)

1. **Boron fallback is very crude** (`rtsc/mlff/hull.py::fallback_element`). ASE's reference state for B is
   not buildable by `bulk()`, so the code falls back to fcc at covalent radius, which is about 1.5 eV/atom too
   high at 0 GPa. The fallback is recorded, and its direction is the "harmless" one. But it makes the hull
   criterion inert for every B-rich candidate. The same applies to the H (`dimer-sc`) and S fallbacks used for H3S.
2. **Hull cache never invalidates.** The cache key is (sub-system, P, model_id) only. It does not include the
   reference set (`reference_dir` contents, `use_mp`) or the relax settings. Failed reference relaxations are
   written to the cache as `{"failed": ...}` entries and never retried (`MLHull.entries_for_subsystem`).
   This did not bite here: `hull_cache/` did not exist before this run, and all 6 cached Mg-B entries are
   successful and converged. But a transient failure, such as the CC error hitting a reference first, or a
   change of `refs/`, would silently reuse stale or empty references.
3. **Auto supercell can be anisotropic for cubic cells** (`phonons.py::auto_supercell`). For the 8-atom H3S
   conventional cell, 3×3×3 (216 atoms) exceeds the 160-atom cap, so one axis is cut to 2×3×3. Verified
   harmless here: isotropic 2/3/4 supercells agree within 3 cm⁻¹.
4. **q-point basis is unlabelled.** phonopy 4.8.6 now defaults `primitive_matrix='auto'` (it warns that this
   differs from v3). The `min_freq_qpoint` in the notes is in the auto-primitive reciprocal basis, which the
   notes do not state.
5. (Not the adapter.) The scratchpad `ms_prescreen.py` originally wrote H3S 5M@150 GPa and Mg-B 5M@150 GPa
   to the same file name, so the Mg run overwrote the H3S file. I fixed it in the script, which now writes
   `prescreen_{h3s|mgb}_<model>_<P>GPa.json`, and reran H3S (identical results within 0.6 cm⁻¹). The old
   un-tagged files are in `~/rtsc_work/mattergen_test/superseded/`.

## 8. Files

- Data, copied to `reports/data/`: `prescreen_h3s_MatterSim-v1.0.0-{1M,5M}_{150,200}GPa.json`,
  `prescreen_mgb_MatterSim-v1.0.0-5M_{0,150}GPa.json`, `h3s_supercell_xcheck.json` (side check),
  `hull_cache/` (reference enthalpies).
- Source of truth: `~/rtsc_work/mattergen_test/`. Scripts are in the session scratchpad: `ms_prescreen.py`,
  `run_prescreen.sh`, `h3s_xcheck.py`, `alphaB.py`, `analyze_mg.py`.
