# QE 7.6 ph.x segfault: electron_phonon + ultrasoft/PAW (reproduction report)

Date: 2026-10-27. Host: WSL2 Ubuntu on Windows 11. All ph.x runs serial (`nice -n 10`, OMP_NUM_THREADS=1), except one 2-rank `-nk 2` test. Nothing was posted anywhere.

## Verdict

**It is a QE source bug introduced in 7.6, not a conda-forge packaging problem.**

In QE 7.6 the signature of `adddvscf` (LR_Modules/adddvscf.f90) gained a third required argument:

```fortran
SUBROUTINE adddvscf(ipert, ik, time_reversed)
   USE uspp, ONLY : okvan
   USE lrus, ONLY : becp1, becpt
   LOGICAL, INTENT(IN) :: time_reversed
   IF (okvan) THEN
      IF (time_reversed) THEN
         CALL adddvscf_(ipert, ik, becpt(ik))   ! becpt is only allocated for time-reversed (magnetic) cases
      ELSE
         CALL adddvscf_(ipert, ik, becp1(ik))
      ENDIF
   ENDIF
```

The callers in `LR_Modules/response_kernels.f90` were updated (`CALL adddvscf(ipert, ik, time_reversed)`), but `SUBROUTINE elphel` in `PHonon/PH/elphon.f90` (qe-7.6 tag, line 550) still calls it with two arguments:

```fortran
           IF (isolv==1) THEN
              call adddvscf (ipert, ik)          ! <-- missing 3rd argument
```

`adddvscf` is an external procedure with no explicit interface, so the compiler does not catch this. At run time `time_reversed` is read through whatever stale pointer is in the third-argument register (`%rdx`). In our run that register points at `__mytime_MOD_called`, which holds 1, so `time_reversed` is read as `.TRUE.`. The code then passes `becpt(ik)`, but `becpt` is not allocated (its descriptor base address is NULL), and the program crashes inside `adddvscf_` on the first access to `becp1_ik%k`. With norm-conserving pseudopotentials `okvan=.FALSE.`, so that branch is never reached, which is why NC is fine. The same two-argument call is still present on the `develop` branch (fetched from the GitHub mirror today).

**Fix (one line, not tested because no source build was done):** `call adddvscf (ipert, ik, .FALSE.)` in `elphel`. It is saved as `reports/qe_bug_repro/elphon_fix.patch`. Because this is undefined behaviour, a different compiler or build could pick up a zero byte and run without crashing, but that is not guaranteed. Every 7.6 build we tested crashes.

## Outcomes (Γ-only ph.x, serial)

| env | conda-forge build | MPI / gfortran-era deps | Al PAW interp | Al PAW simple | Al PAW no-eph | Al SG15 (NC) interp | hcp Mg GBRV-US interp |
|---|---|---|---|---|---|---|---|
| qe (existing) | qe 7.6 mpi_openmpi_hdb26816_104 | openmpi 5.0.11, elpa 2026.02.002, hdf5 2.2.0 | **SEGV** | **SEGV** | OK | OK | **SEGV** |
| qe76mpich | qe 7.6 mpi_mpich_he3edec4_4 | mpich 5.0.2, elpa 2026.02.002, hdf5 2.2.0 | **SEGV** | **SEGV** | OK | OK | **SEGV** |
| qe76b0 | qe 7.6 h2ec8f75_0 (first 7.6 build) * | openmpi 5.0.11, hdf5 2.2.0 | **SEGV** | **SEGV** | OK | OK | **SEGV** |
| qe75b | qe 7.5 h19104ac_2 | openmpi 5.0.11, elpa 2025.06.001, hdf5 1.14.6 | OK | OK | OK | OK | OK |
| qe75a | qe 7.5 hbc83ec2_0 | openmpi 4.1.6, elpa 2021.11.002, hdf5 1.14.3 | OK | OK | OK | OK | OK |
| qe74 | qe 7.4 hbc83ec2_0 | openmpi 4.1.6, elpa 2021.11.002, hdf5 1.14.3 | OK | OK | OK | OK | OK |
| qe72 | qe 7.2 hac89879_4 | openmpi 4.1.6, elpa 2021.11.002, hdf5 1.14.3 | refused** | refused** | OK | OK | OK |

\* Build `h2ec8f75_0` links `libelpa.so.19` but does not declare `elpa` as a dependency, so pw.x fails to load out of the box (a separate packaging defect, fixed in later builds). I added `elpa 2026.02.002 mpi_openmpi_h33b1e68_2` from cache to make it run.
\** QE 7.2 stops in `phq_readin` with "Electron-phonon calculations with PAW not tested". It still works with US (Mg GBRV).
No conda-forge `qe=7.3` exists (the versions jump from 7.2 to 7.4). No `nompi` variant of 7.6 exists. 7.6 builds are only `mpi_openmpi_*` or `mpi_mpich_*`.
"OK" means `JOB DONE` was reached and the expected `lambda(` lines were written (30 for Al, 60 for Mg).

Other variants on the 7.6 openmpi build (Al PAW, interpolated):

| variant | result |
|---|---|
| `nosym=.true.` (pw) | SEGV |
| `noinv=.true.` (pw) | SEGV |
| no `fildvscf` | not testable: ph.x stops with "El-ph needs a DeltaVscf file" |
| 2 MPI ranks, `-nk 2` | SEGV (7.5 h19104ac_2 with `-nk 2`: OK) |
| q=(1/2)-grid point only (`start_q=2,last_q=2`, from earlier crashchar run) | SEGV |
| two-step (phonon first, then `trans=.false.` e-ph) (earlier crashchar run) | SEGV |

So the crash does not depend on symmetry, pools, the q point, the MPI flavour, the dependency stack or the build number. It depends only on QE 7.6, `electron_phonon` being set, and `okvan` (any US/PAW species). The earlier crashchar runs show it also happens for Pb/Nb PAW, bcc S US and MgB2 US. As a sanity check, NC frequencies agree across versions (Γ acoustic, -1.8250 THz in 7.6 vs -1.8248 THz in 7.5, 7.4 and 7.2), and the US Mg results agree between 7.4 and 7.5 (λ(5)=0.0600, γ=5.72 GHz).

**Regression window:** works in conda-forge 7.4 and 7.5 (all builds tested). Broken in 7.6 (all three build lines tested: openmpi _0 and _104, mpich _4). On the source side, `elphon.f90`'s call did not change, but `adddvscf.f90` changed between the qe-7.5 and qe-7.6 tags (the 7.2/7.4/7.5 adddvscf.f90 files are byte-identical, md5 b0d142d9...). In 7.5, `adddvscf(ipert,ik)` took two arguments and a separate `adddvscf_ph_mag` used `becpt`. In 7.6 these were merged into a 3-argument routine and the `adddvscf_ph_mag` call in elphel was commented out, but the remaining call was not updated.

## Backtrace

gfortran's built-in backtrace has no symbols (`???`). The binaries are not stripped (they have `.symtab`) but carry no DWARF info, so gdb 18.1 (conda-forge, env `dbg`) gives function names and instruction offsets but no line numbers:

```
Thread 1 "ph.x" received signal SIGSEGV, Segmentation fault.
0x00005555557db161 in adddvscf__ ()
#0  adddvscf__ ()      <+977>   imul 0x98(%rax),%rdi   ; rax = 0x128 (NULL base + offset)
#1  elphel_ ()         <+7467>  (return after call adddvscf_)
#2  solve_linter_ ()
#3  phqscf_ ()
#4  do_phonon_ ()
#5  MAIN__ ()
#6  main ()
```
(`adddvscf__` is the gfortran mangling of `adddvscf_`, the inner worker; the faulting load reads `dim[2].stride` of the `%k(:,:)` component of the `bec_type` argument.)

Breakpoint on entry to `adddvscf_` (outer wrapper) when called from `elphel_`:
```
ipert=1 ik=1  3rd-arg ptr (rdx)=0x5555573401a0
0x5555573401a0 <__mytime_MOD_called>:  0x00000001      -> time_reversed reads .TRUE.
__lrus_MOD_becpt descriptor: base_addr = 0x0000000000000000   (unallocated)
__lrus_MOD_becp1 descriptor: base_addr = 0x0000555558183e70   (allocated)
```

## Commands

Inputs are listed below and also stored in `reports/qe_bug_repro/`. Per case:
```
export OMP_NUM_THREADS=1 GFORTRAN_ERROR_BACKTRACE=1
nice -n 10 $QE_BIN/pw.x -in fit.in > fit.out     # 16^3, la2F
nice -n 10 $QE_BIN/pw.x -in scf.in > scf.out     # 8^3
nice -n 10 $QE_BIN/ph.x -in ph.in  > ph.out 2> ph.err
```
Packaged as `QE_BIN=... SSSP_DIR=... SG15_DIR=... reports/qe_bug_repro/run.sh all`. It was verified to give SEGV/SEGV/OK/OK/SEGV on 7.6 and OK on 7.5. gdb: `gdb -batch -ex run -ex bt --args $QE_BIN/ph.x -in ph.in` (scripted breakpoint commands are in the scratchpad `gdb.cmd`).
Envs were created with `micromamba create -n <env> -c conda-forge --override-channels qe=<ver>=<build>`.

### Inputs (fcc Al PAW, minimal case)
fit.in:
```
&control
 calculation='scf', prefix='t', outdir='./tmp/', pseudo_dir='<SSSP_1.3.0_PBE_precision>'
/
&system
 ibrav=2, celldm(1)=7.65, nat=1, ntyp=1, ecutwfc=30, ecutrho=240,
 occupations='smearing', smearing='mv', degauss=0.02, la2F=.true.
/
&electrons
 conv_thr=1e-10
/
ATOMIC_SPECIES
 Al 26.98 Al.pbe-n-kjpaw_psl.1.0.0.UPF
ATOMIC_POSITIONS crystal
 Al 0.0 0.0 0.0
K_POINTS automatic
 16 16 16 0 0 0
```
scf.in: same, without `la2F=.true.`, `K_POINTS 8 8 8 0 0 0`.
ph.in:
```
ph
&inputph
 prefix='t', outdir='./tmp/', tr2_ph=1e-14, fildyn='dyn', fildvscf='dvscf',
 electron_phonon='interpolated', el_ph_sigma=0.005, el_ph_nsigma=10,
 trans=.true., ldisp=.true., nq1=1, nq2=1, nq3=1
/
```
Controls: `electron_phonon='simple'`; no electron_phonon line (and no fit step); SG15 `Al_ONCV_PBE-1.2.upf` with ecutrho=120. US case: hcp Mg, `ibrav=4, celldm(1)=6.05, celldm(3)=1.624, nat=2`, Mg at (1/3,2/3,1/4),(2/3,1/3,3/4), `mg_pbe_v1.4.uspp.F.UPF`, k 6x6x4 scf / 12x12x8 fit.

Pseudopotential md5s:
- `Al.pbe-n-kjpaw_psl.1.0.0.UPF` (SSSP 1.3.0 PBE precision): cfc449ca30b5f3223ec38ddd88ac046d
- `mg_pbe_v1.4.uspp.F.UPF` (SSSP 1.3.0 PBE precision, GBRV): 8ffbd8f729fef71095aac5bd8316fb1f
- `Al_ONCV_PBE-1.2.upf` (SG15): 2d9e7bc02555220d62fad08c3c4858e7

## Package builds (micromamba list)
```
## qe
  qe                 7.6           mpi_openmpi_hdb26816_104  conda-forge
  openmpi            5.0.11        hee33aa4_0                conda-forge
  elpa               2026.02.002   mpi_openmpi_h33b1e68_2    conda-forge
  hdf5               2.2.0         mpi_openmpi_h0b8148c_1    conda-forge
  fftw               3.3.11        nompi_h3b011a4_100        conda-forge
  libopenblas        0.3.34        openmp_h156ec52_2         conda-forge
  scalapack          2.2.0         h606478a_6                conda-forge
  libgfortran5       16.2.0        h6b99dfc_7                conda-forge
## qe76mpich
  qe                 7.6           mpi_mpich_he3edec4_4   conda-forge
  mpich              5.0.2         h76636b3_0             conda-forge
  elpa               2026.02.002   mpi_mpich_h49e544b_2   conda-forge
  hdf5               2.2.0         mpi_mpich_h581622d_1   conda-forge
  scalapack          2.2.0         h13b89aa_6             conda-forge
## qe76b0
  qe                 7.6           h2ec8f75_0              conda-forge
  openmpi            5.0.11        hee33aa4_0              conda-forge
  elpa               2026.02.002   mpi_openmpi_h33b1e68_2  conda-forge   (added manually)
  hdf5               2.2.0         mpi_openmpi_h0b8148c_1  conda-forge
## qe75b
  qe               7.5           h19104ac_2               conda-forge
  openmpi          5.0.11        hee33aa4_0               conda-forge
  elpa             2025.06.001   mpi_openmpi_h335e491_1   conda-forge
  hdf5             1.14.6        mpi_openmpi_h8367ee7_10  conda-forge
  libopenblas      0.3.34        pthreads_hf13c14d_2      conda-forge
## qe75a
  qe               7.5           hbc83ec2_0              conda-forge
  openmpi          4.1.6         hc5af2df_101            conda-forge
  elpa             2021.11.002   mpi_openmpi_h335e491_2  conda-forge
  hdf5             1.14.3        mpi_openmpi_hf5a2fd3_4  conda-forge
## qe74
  qe               7.4           hbc83ec2_0              conda-forge
  openmpi          4.1.6         hc5af2df_101            conda-forge
  elpa             2021.11.002   mpi_openmpi_h335e491_2  conda-forge
  hdf5             1.14.3        mpi_openmpi_hf5a2fd3_4  conda-forge
## qe72
  qe               7.2           hac89879_4              conda-forge
  openmpi          4.1.6         hc5af2df_101            conda-forge
  elpa             2021.11.002   mpi_openmpi_h335e491_2  conda-forge
  hdf5             1.14.3        mpi_openmpi_hf5a2fd3_4  conda-forge
## dbg
  gdb               18.1        py314h996d46b_0       conda-forge
```
All envs: fftw 3.3.11 nompi_h3b011a4_100, libgfortran5 16.2.0 h6b99dfc_7 (runtime).

## Resources used
- New envs: qe76mpich, qe76b0, qe75b, qe75a, qe74, qe72, dbg (about 3.4 GB on disk in total). Env `qe` was not modified.
- Download: about 0.9 GB (47 packages not already in env `qe`, 887 MB compressed). The dry-run worst case was 1.29 GB.
- Run directories: `~/rtsc_work/qebug/<env>/<case>_<variant>/`.

## Follow-ups
- Workarounds for production: use conda-forge **qe 7.5** (h19104ac_2 or hbc83ec2_0) for e-ph with US/PAW, or NC pseudos with 7.6. 7.5 and 7.6 gave identical NC results to about 1e-4 THz. Do not mix ph.x versions on one outdir.
- Upstream report or MR (not done, per instructions): the one-line fix in `elphel`. Other 2-argument `adddvscf` callers in the tree (e.g. EPW, HP, other PHonon routines) should also be grepped; I only checked elphon.f90 and response_kernels.f90.
- No conda-forge debug build exists. A source build of 7.6 ph.x with the patch would confirm the fix, but it was intentionally not done here.
