# DRAFT — not posted. For conda-forge/qe-feedstock (GitHub). Review before posting.

> Before posting: confirm https://gitlab.com/QEF/q-e/-/merge_requests/2963 is merged and still says
> what is summarised below (it could not be re-read from this machine because of GitLab rate limiting).

---

**Title:** qe 7.6: ph.x segfaults for electron-phonon with ultrasoft/PAW pseudopotentials (upstream fix in QEF/q-e!2963) — please backport

**Body:**

### Summary

In every conda-forge `qe` 7.6 build I tested, `ph.x` crashes with SIGSEGV whenever `electron_phonon` is set (`'simple'` or `'interpolated'`) and any species uses an ultrasoft or PAW pseudopotential. Norm-conserving pseudopotentials work, and the same inputs work with conda-forge `qe` 7.4 and 7.5.

This is an upstream source bug in QE 7.6, not a build problem. It is fixed on QE `develop` by https://gitlab.com/QEF/q-e/-/merge_requests/2963 ("segmentation fault in "adddvscf" because of insufficient dummy arguments"). No 7.6.x release contains the fix yet. Could the patch be applied to the feedstock until one does?

### Cause

QE 7.6 changed `LR_Modules/adddvscf.f90` to `SUBROUTINE adddvscf(ipert, ik, time_reversed)`. `PHonon/PH/elphon.f90` (subroutine `elphel`, about line 550) still calls `call adddvscf (ipert, ik)` with two arguments. There is no explicit interface, so this compiles. At run time the missing logical is read from an arbitrary address. In gdb it read as `.TRUE.`, so the code used the unallocated `becpt(ik)`. That branch only runs when `okvan` is true, which is why norm-conserving pseudopotentials are unaffected.

The regression came in with QEF/q-e commit c08649a9 ("LR: Merge adddvscf and adddvscf_ph_mag, remove the latter"). According to the MR description, !2963 also fixes the same two-argument call in `ep_matrix_element_wannier.f90`, `elphel_refolded`, `solve_e_fpol.f90` and `TDDFPT/lr_apply_liouvillian_eels.f90`.

The minimal fix for `elphel` (what !2963 does there, given that the call sits inside `IF (isolv==1)`):

```diff
-              call adddvscf (ipert, ik)
+              call adddvscf (ipert, ik, isolv==2)
```

### Builds tested (linux-64, serial, Γ-only)

| build | Al PAW, e-ph interpolated | Al PAW, e-ph simple | Al PAW, no e-ph | Al NC (SG15), e-ph | hcp Mg US (GBRV), e-ph |
|---|---|---|---|---|---|
| qe 7.6 mpi_openmpi_hdb26816_104 | SEGV | SEGV | OK | OK | SEGV |
| qe 7.6 mpi_mpich_he3edec4_4 | SEGV | SEGV | OK | OK | SEGV |
| qe 7.6 h2ec8f75_0 | SEGV | SEGV | OK | OK | SEGV |
| qe 7.5 h19104ac_2 | OK | OK | OK | OK | OK |
| qe 7.5 hbc83ec2_0 | OK | OK | OK | OK | OK |
| qe 7.4 hbc83ec2_0 | OK | OK | OK | OK | OK |

The 7.6 crash also happens with `nosym`, with `noinv`, at non-Γ q points, with 2 MPI ranks and `-nk 2`, and in a two-step run (phonons first, then e-ph with `trans=.false.`).

### Minimal reproducer (about 20 s, serial)

fcc Al, SSSP 1.3.0 PBE `Al.pbe-n-kjpaw_psl.1.0.0.UPF`, ecutwfc/ecutrho 30/240 Ry, MV smearing 0.02 Ry.

1. `pw.x` scf on a 16³ k grid with `la2F=.true.`.
2. `pw.x` scf on an 8³ k grid.
3. `ph.x` at Γ only, with `electron_phonon='interpolated'`, `fildvscf='dvscf'`, `trans=.true.`.

Inputs and `run.sh` are attached. Backtrace (symbols, no debug info):

```
#0 adddvscf_ <+977>   (NULL base of __lrus_MOD_becpt + offset)
#1 elphel_
#2 solve_linter_  #3 phqscf_  #4 do_phonon_  #5 MAIN__
```

### Separate, minor

`qe 7.6 h2ec8f75_0` links `libelpa.so.19` but does not list `elpa` as a run dependency, so `pw.x` fails to start unless elpa is installed separately. This may already be covered by #88.

---

Attachments to add when posting: the files in `reports/qe_bug_repro/` (`run.sh`, `*.in`). Do not attach UPF files; cite SSSP 1.3.0 by file name.
