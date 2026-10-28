# EPW (QE 7.6) interface spec for the automated α²F convergence protocol

Date: 2026-10-27. Scope: what rtsc must write and parse to drive EPW from QE 7.6, plus the independent
ph.x optimized-tetrahedron route. The rtsc bundle was **not** edited.

**Evidence tags.** **[SRC]** = read in source at tag `qe-7.6` (raw.githubusercontent.com/QEF/q-e/qe-7.6/...).
**[RUN]** = checked against the installed binaries (WSL `~/micromamba/envs/qe`, conda-forge `qe 7.6`).
**[INF]** = my inference or something I recall from the literature; I did not verify it here.
The local copies of the sources are in
`G:\ClaudeWorkspace\RTSC_Project\run_scripts\scratchpad_2026-10-27\wsl\src\{raw,raw2,ts}`.

---

## 1. Version, namelist location, key diff

* **EPW version is v6.1.** [RUN] The `epw.x` banner prints `Program EPW v.6.1` and cites Lee et al., npj Comput Mater 9, 156 (2023).
  Wannier90 is 3.1.0 (`wannier90.x -v`; epw.x links `libqe_wannier90.a`). [RUN]
* **The namelist is now read in `EPW/src/readin.f90`**, at `NAMELIST / inputepw /` (line 177) and
  `READ(qestdin, inputepw, IOSTAT=ios)` (line 885). The file `epw_readin.f90` no longer exists. [SRC]
  Defaults are set at lines 540-900. The consistency checks are at lines 950-1150.
  Units are converted at lines 1310-1352. [SRC]
* **Size.** 331 keys at 7.6, against 168 in rtsc's `EPW53_KEYS`. [SRC] The list is saved as `scratchpad\wsl\src\keys76.txt`.

### Removed or renamed keys from EPW53_KEYS
Any of these is a fatal `Bad line in namelist` error at 7.6. [SRC] I also confirmed this for `eps_acustic` and `kmaps` with epw.x. [RUN]

| 5.3 key | 7.6 status |
|---|---|
| `eps_acustic` | **renamed `eps_acoustic`**. The default changed from 5 cm⁻¹ to **0.1 cm⁻¹** (05/2025 comment, readin.f90:588). |
| `kmaps` | **removed** from the namelist (internal variable only). |
| `cumulant` | now `lcumulant` |
| `longrange` | now `longrange_only` |
| `lunif`, `n_r`, `nswc`, `nswfc`, `pwc`, `wsfc` | removed |

Note: QE's own `EPW/examples/mgb2/epw/epw.in` at tag qe-7.6 still contains `eps_acustic`. It is stale and would fail. [SRC+RUN]

### Keys relevant to us (all present at 7.6) [SRC]
`wannierize, nbndsub, bands_skipped, proj, auto_projections, scdm_proj, scdm_entanglement, scdm_mu,
scdm_sigma, dis_win_min/max, dis_froz_min/max, num_iter, wdata, a2f, phonselfen, delta_approx,
degaussw, ngaussw, nsmear, delta_smear, degaussq, delta_qsmear, nqsmear, nqstep, fsthick, nk1-3, nq1-3,
nkf1-3, nqf1-3, mp_mesh_k, mp_mesh_q, use_ws, epwwrite, epwread, epbwrite, epbread, eliashberg, liso,
laniso, limag, ephwrite, temps, nstemp, band_plot, filkf, filqf, efermi_read, fermi_energy,
assume_metal, eps_acoustic, fila2f, iverbosity, etf_mem, vme, lifc, asr_typ`.

New keys that matter to us:
* `a2f_iso`: isotropic α²F on the fly with image parallelization. It is incompatible with `ephwrite` and `laniso`.
* `dos_tetra`: tetrahedron electronic DOS only. It writes `prefix.dos` and does **not** feed α²F.
* `calc_nelec_wann`: set automatically when `dis_win_min` is given.
* `lfast_kmesh`: set by `etf_mem=3`. It requires `mp_mesh_k`, `efermi_read`, and no `phonselfen`.
* `lopt_w2b`, `interpolate` (default .true.), `epw_memdist` (ignored), `fd`, `gridsamp`, `fbw`, `tc_linear`, `icoulomb`.

### Defaults to note [SRC]

| Key | Default |
|---|---|
| `epwwrite` | .true. |
| `dvscf_dir` | `'.'` (a trailing slash is added by `trimcheck`) |
| `fsthick` | 1e10 eV |
| `degaussw` | 0.025 eV |
| `ngaussw` | 1 (Methfessel-Paxton) |
| `nqsmear` | 10 |
| `degaussq`, `delta_qsmear` | 0.05 meV each |
| `nsmear` | 1 |
| `delta_smear` | 0.01 eV |
| `nqstep` | 500 |
| `num_iter` | 200 |
| `use_ws` | .false. |
| `vme` | `'wannier'` |
| `mp_mesh_k` | .false. |
| `scdm_entanglement` | `'isolated'` |
| `scdm_mu` | 0 eV |
| `scdm_sigma` | 1 eV |
| `assume_metal` | .false. |
| `temps` | 300 K when unset |

### Input constraints we must respect [SRC readin.f90]

**General**
* Either `epwread` or `epwwrite` must be true.
* `epbread` cannot be combined with `epbwrite` or with `epwread`.
* `elph`, `wannierize` and `epwread` together is an error ("must use same w90 rotation matrix").
* `4*degaussw > fsthick` is an error.
* `nbndsub` must be ≤ 200.

**Routes and grids**
* `a2f` requires `phonselfen`.
* `phonselfen` (and `specfun_ph`) is an error with `mp_mesh_k` ("only works with full uniform k-mesh").
* `mp_mesh_k` and `mp_mesh_q` cannot both be true.
* `ephwrite` requires `nkf` to be a multiple of `nqf` and uniform grids. `eliashberg_init` repeats the multiple check.

**Wannierization**
* `auto_projections` requires `scdm_proj`, and the reverse. `auto_projections` together with `proj(1)` is an error.

**band_plot**
* It needs `filkf` and/or `filqf`.
* It is an error to give both `filkf` and `nkf`, or both `filqf` and `nqf`.
* With `filkf`, set `efermi_read`. Otherwise EPW prints a warning.

---

## 2. Output files and formats

There are **three distinct writers** that produce α²F-type files. Their formats differ. Every writer puts its files in the CWD
of epw.x, not in `outdir`, except the files under `<outdir>/<prefix>.ephmat/`.

### 2a. Eliashberg route
This route needs `ephwrite=.true.` and `eliashberg=.true.`, with `liso`, `laniso`, or neither.
The writer is `evaluate_a2f_lambda` in `EPW/src/supercond.f90`, lines 392-950. [SRC]

`liso=laniso=.false.` with `eliashberg=.true.` is a legal "compute α²F only" mode (`supercond_driver.f90`:118-131). [SRC]

**`<prefix>.a2f`**
```
 w[meV] a2f and integrated 2*a2f/w for   NQ smearing values          <- header (leading space)
 <w(meV)> <a2F_1..a2F_NQ> <lam_cum_1..lam_cum_NQ>   format (f12.7, 20f12.7), nqstep rows
 Integrated el-ph coupling                                          <- list-directed (leading space)
  #         l_1 ... l_NQ          (15f12.7)
 Phonon smearing (meV)
  #         s_1 ... s_NQ          (15f12.7), s_i = degaussq + (i-1)*delta_qsmear
Electron smearing (eV)  <degaussw>   <- NO leading space here
Fermi window (eV)  <fsthick>
Summed el-ph coupling  <l_sum>       (Σ w0g w0g g²/ω / N_F, a direct sum, not the integral)
```
* The grid is `wsph(i) = i*dwsph`, with `dwsph = wsphmax/nqstep` and `wsphmax = 1.1*max ω` over the fine q grid.
  Sources: `io_supercond.f90`:1420-1437 and, when the file is read back, 1321-1325. The first point is ω = dw, not 0.
* α²F = ½ Σ_k,q,n,m w_k w_q w0g(ε_nk) w0g(ε_mk+q) |g|² δ_σq(ω−ω_qν) / N_F, built only from modes with ω > `eps_acoustic`.
* `w0g = w0gauss(x,0)/degaussw` is always a **Gaussian**, set at `io_supercond.f90`:1647.
* `N_F = dos_ef(ngaussw, degaussw, …)/2` uses the **`ngaussw`** smearing type (default MP-1) and all `nbndsub` bands.
  It is computed in `write_ephmat` (`io_supercond.f90`:2369) and stored in `ephmat/egnv`.
  → Set `ngaussw=0` so that the numerator and denominator use the same kernel. [SRC; the recommendation is INF]
* The cumulative λ is the rectangle sum Σ 2 a2F_i/w_i · dw. [SRC]

**`<prefix>.a2f_proj`**
Header `w[meV] a2f a2f_modeproj`. Rows are `(f12.7,100f12.7)`: w, a2F (first smearing only), then one column per mode.
The footer is `lambda_int = … lambda_sum = …`. [SRC]

**`<prefix>.phdos`**
* Header ` w[meV] phdos[states/meV] for NQ smearing values`, rows `(f12.7, 15f15.7)`. [SRC]
* `<prefix>.phdos_proj` sits beside it.
* Written by `write_phdos` from Wannier-interpolated frequencies over all fine-grid q.

**`<prefix>.lambda_FS`**
* Header `#  k-point  Band Enk-Ef [eV]  lambda`.
* Rows `(3f12.6, i8, f12.6, f24.15)`: k cartesian (2π/alat), band index inside the window, E−E_F (eV), λ_nk.
* It lists every full-grid k with |E−E_F| < fsthick.

**`<prefix>.lambda_k_pairs`**
* A histogram of ρ(λ_nk).
* Two `#` header lines, then rows `(1x,ES20.10,2(2x,ES20.10))`: λ, distribution scaled to 1, distribution unscaled.

**Written only with `iverbosity=2`**
* `<prefix>.lambda_aniso`: header `# Enk-Ef[eV]  lambda_nk# kpt# band`, rows `(2f12.7,2i7)`.
* `<prefix>.lambda_pairs`.
* `<prefix>.lambda_XX.cube` and `<prefix>.lambda.frmsf`.

**`<outdir>/<prefix>.ephmat/freq`**
* It is written only when neither `a2f` nor `a2f_iso` is set (`io_supercond.f90`:769).
* **It is now `ACCESS='stream'`** (`io_supercond.f90`:911-920), with **no record markers**.
* Layout: `int32×5 (nqtotf,nqf1,nqf2,nqf3,nmodes)`, then for each q `real64×3 xq` followed by `real64×nmodes wf` (Ry).
* `egnv` is also a stream file: `nkftot,nkf1-3,nkfs`, then `nbnd,ef,ef0,dosef,degaussw,fsthick` (Ry), and so on.

**Reuse pitfall (find_a2f).**
If `<prefix>.a2f` (or `fila2f`) already exists in the CWD, EPW **reuses it and does not recompute α²F**. It only prints
"a2f file is found…" (`supercond.f90`:1139-1155). [SRC]
→ The protocol must delete or rename `prefix.a2f` before every run, or run each grid in a fresh directory.

**Electron smearing is frozen into ephmat.**
`read_eigenvalues` overwrites `degaussw`, `fsthick`, `ef0` and `dosef` with the values stored in `ephmat/egnv`
(`io_supercond.f90`:1522-1531). [SRC]
→ Changing `degaussw` in a later `ephwrite=.false.` Eliashberg-only run has **no effect**. Each electron smearing needs its own
`ephwrite` pass. That pass can reuse `epwread=.true.`, so it costs interpolation only.

### 2b. Phonon self-energy route
This route needs `phonselfen=.true.` and `a2f=.true.`. The writer is `a2f_main` in `EPW/src/spectral.f90`, lines 1029-1373. [SRC]

**File names.** For each temperature and each electron smearing `isig = 1..nsmear` it writes:
* `<prefix>.a2f.0<isig>.<T>` (for isig ≥ 10 the name is `.a2f.<isig>.<T>`). `<T>` is `f8.3` adjusted, for example `pb.a2f.01.0.075`.
* `<prefix>.a2f_tr.0<isig>.<T>`, `<prefix>.phdos.0<isig>.<T>` and `<prefix>.res.0<isig>.<T>`.
* There is **no plain `<prefix>.a2f`**.

**Format.**
* **No header line.** Rows are `(f12.7,15f12.7)` = w(meV), then a2F_1..a2F_NQ only. **There are no cumulative-λ columns.**
* The footer is the same five items as in 2a, but every line has a leading space.
* `Electron smearing (eV)` = degaussw + (isig−1)·delta_smear.
* `Summed el-ph coupling` = Σ_q w_q Σ_ν λ_qν.

**Grid and clipping.**
* `ww(i) = i*dw`, with `dw = 1.1*max(wf)/nqstep`.
* λ_qν < 0 is clipped to 0 before it enters α²F.

**Kernel (set by `delta_approx`).**
* Default `delta_approx=.false.`: γ uses −π[f(ε_k+q)−f(ε_k)] δ_σ(ε_k+q−ε_k−ω_q). Here f is Fermi-Dirac at `temps`, so the result is **not**
  the double-delta Allen form and it depends on temperature.
* With `delta_approx=.true.` (as in test-suite `epw_metal/epw2_up.in`): γ uses π ω_q w0g(ε_k) w0g(ε_k+q) with a Gaussian of width
  degaussw0. → **Use `delta_approx=.true.`** for the Allen/McMillan α²F. [SRC]

**Fermi level and N_F.** For each smearing, E_F is recomputed by `efermig(…, degaussw0, ngaussw)` when `nsmear>1`, unless `efermi_read`.
N_F uses `ngaussw` (`selfen.f90`:1135-1160). [SRC]

**Cost and constraints.** This route needs `mp_mesh_k=.false.` (full k grid). Adaptive smearing is forbidden
(`selfen.f90`: `adapt_smearing cannot be used with phonon self-energy`). [SRC]

### 2c. `a2f_iso=.true.`
The writer is `EPW/src/supercond_vertex.f90`, lines 690-728. [SRC]
* The header still says `for NQ smearing values` (NQ = `nqsmear`), but each row holds **one** smearing only: `w, a2F, lam_cum`.
* There is an extra footer line `DOS (eV) …`.

### 2d. Compatibility with rtsc parsers
These conclusions come from reading `rtsc/dfpt/parsers.py` (`parse_epw_a2f`, `epw_a2f_to_spectra`) and `band_coupling.py`.

| File | rtsc expectation | 7.6 reality | Verdict |
|---|---|---|---|
| 2a `<prefix>.a2f`, nqsmear ≤ 10 | header with N, 1+2N columns, `Integrated` block equal to last cumulative row | same | **OK** |
| 2a with nqsmear 11-15 | 1+2N columns per line | `20f12.7` wraps each row onto 2 records | **breaks** (shape error) |
| 2a with nqsmear > 15 | — | footer `15f12.7` wraps as well | **breaks** |
| 2b `.a2f.0i.T` | header required | no header, no cumulative columns, different file name | **breaks** (`missing EPW header`) |
| 2c a2f_iso | header N = 1+2N columns | header says nqsmear, file has 3 columns | **breaks** unless `nqsmear=1` |
| `epw_a2f_to_spectra` dw | dw = w[0] | w_i = i·dw | OK, and matches EPW's rectangle rule |
| `epw_omega_max_meV` | last w = 1.1·max ω | still true for 2a (`io_supercond.f90`:1420) | OK |
| `parse_epw_lambda_fs` | 6 columns, `#` header | same | OK |
| `parse_epw_lambda_aniso` | (E, λ, ik, ibnd) | same `(2f12.7,2i7)` | OK |
| `read_epw_freq` | gfortran sequential records | **stream** | **breaks** (`bad record length`) |
| `parse_epw_summed_lambda` | regex `Summed el-ph coupling` | present in 2a and 2b | OK |

**Precision note for certified bounds.** All α²F and λ values are printed as `f12.7`: absolute rounding of ±5·10⁻⁸ per value, and ω to 10⁻⁷ meV.
rtsc should carry this as an explicit rounding envelope (for example a2F_i ± 5e-8). [SRC format; the envelope is INF]

---

## 3. Smearing control for the double-delta integral

* **Electron smearing.** `degaussw` (eV) sets the width of the Gaussian `w0g` in both routes. `ngaussw` affects only E_F and N_F.
  * Route 2a: one electron smearing per `ephwrite` pass. `nsmear` and `delta_smear` are ignored there.
  * Route 2b: `nsmear` × `delta_smear` gives several electron smearings in one run, one file per `isig`. [SRC]
* **Phonon smearing.** `degaussq` (meV) + (i−1)·`delta_qsmear` for i = 1..`nqsmear`. This broadens δ(ω−ω_qν) when α²F is histogrammed.
  The columns of `.a2f` are the **phonon** smearings. [SRC]
* **Frequency bins.** `nqstep` points from dw to 1.1·ω_max. [SRC]
* **Fermi window.** `fsthick` must be ≥ 4·degaussw (hard check). Only states inside the window enter the numerator; N_F includes all bands. [SRC]
* **Adaptive smearing.** It exists (`adapt_smearing` is triggered by `degaussw=0`, `use_wannier.f90`:991-1002), but it is **forbidden for
  phonselfen and a2f** and is meant for self-energy and transport only. In route 2a, `degaussw=0` would divide by zero in `w0g`.
  [SRC; the route 2a part is INF since there is no explicit guard]
* **Tetrahedron.** EPW 7.6 has **no tetrahedron α²F**. `dos_tetra` only replaces the electronic `prefix.dos` used by fbw and Coulomb
  (`io_supercond.f90` `write_dos_tetra`). [SRC]
  → The independent tetrahedron route must be ph.x `lambda_tetra` plus `alpha2f.x` (§7).

---

## 4. SCDM in EPW 7.6

**Required keys**
* `wannierize=.true.`, `auto_projections=.true.`, `scdm_proj=.true.`, **no `proj(i)`**, `nbndsub=<num_wann>`.
* `scdm_entanglement='erfc'` (metals) or `'gaussian'`. `'isolated'` is for isolated groups of bands.
* `scdm_mu` and `scdm_sigma` in **eV, absolute energy scale (not relative to E_F)**.
* Optional: `bands_skipped='exclude_bands = …'` and `dis_*` windows for the W90 disentanglement step. [SRC readin.f90:1094-1100, 1141-1144]

**What the code does** (`pw2wan.f90` `compute_amn_with_scdm`, line 1215 and following). [SRC]
* f(ε) = ½·erfc((ε − scdm_mu)/scdm_sigma) for erfc, or exp(−(ε−μ)²/σ²) for gaussian, with ε = `et*rytoev` (absolute eV).
* It applies f(ε) at **Γ only**, then runs a QR with column pivoting (ZGEQP3) on ψ_Γ(r)·f.
* It **requires Γ to be the first k point** of the nscf list (`No Gamma point found` otherwise).
* The `.win` file gets `auto_projections = .true.` instead of a projections block (`wannierization.f90`:140).
* `scdm_sigma ≤ 0` is an error.
* Defaults: `scdm_mu = 0 eV`. For a metal with E_F ≈ 8 eV, the default would suppress the states at E_F. **Always set it.**

**Atomic orbitals are not needed.** In EPW, `auto_projections` only works together with SCDM and **does not use PP_PSWFC**.
SCDM builds the projections from the Bloch states themselves. [SRC]

**SG15 pseudopotentials** [RUN]
* Al, Pb, Nb, Mg, B, H and S `_ONCV_PBE-1.2.upf` (and every file in `~/pseudo/sg15`) contain an **empty** `<PP_PSWFC></PP_PSWFC>`
  block, `number_of_wfc="0"`, and no `PP_CHI`.
* → projwfc-based projectability, W90 atomic projections and Vitale-style automatic choices are **impossible with SG15**. SCDM itself works.
* For comparison, SSSP efficiency has PP_CHI for Al (2), Pb (3), Nb (4), Mg (2) and H (1). The GBRV `b_pbe_v1.4.uspp.F.UPF` and
  `s_pbe_v1.4.uspp.F.UPF` files report 0 by this grep. They use an older layout, so check their `PP_PSWFC` text before relying on them.

**Choosing scdm_mu, scdm_sigma and nbndsub automatically** [INF / literature, not verified here]
* Vitale et al., npj Comput. Mater. 6, 66 (2020):
  * Compute the projectability p_nk of the Bloch states onto the pseudo-atomic orbitals (this needs PP_PSWFC).
  * Fit p(ε) with ½·erfc((ε−μ_fit)/σ_fit).
  * Set **σ = σ_fit** and **μ = μ_fit − 3σ_fit**, and num_wann = the number of PAOs.
  * Then judge quality by comparing interpolated bands with DFT bands.
* A protocol that works without PSWFC (SG15):
  1. Pick `nbndsub` from a valence-orbital table per element: for example s+p (4) for sp metals, s+p+d (9) for transition metals,
     s (1) for H. Optionally add semicore states only if they are not excluded.
  2. Set `scdm_mu = E_F(DFT, absolute eV) + Δ` with Δ ≈ 1-3 eV, and `scdm_sigma ≈ 1-2 eV`.
  3. Make sure the nscf `nbnd` exceeds the band index of μ + 3σ at Γ (plus excluded bands).
  4. Accept only if the band check in §6 passes. If it fails, scan Δ ∈ {0,1,2,3} eV × σ ∈ {1,2} eV, then increase `nbndsub`.
* Optional alternative: run one projectability-only calculation with an SSSP or PseudoDojo set that has PSWFC, to get μ_fit, σ_fit and
  the PAO count. Then use those numbers with SG15. Band energies differ slightly between pseudopotential families, so the §6 check
  remains the arbiter.

---

## 5. The pp step: `epw_pp.py` against `collect_ph_for_epw`

* The installed `epw_pp.py` is **byte-identical** to `EPW/bin/pp.py` at qe-7.6. [RUN diff]
* **Interface.** Run it inside the ph.x outdir. It asks once for the prefix through `input()`, so it can be run non-interactively as
  `echo <prefix> | python3 epw_pp.py`. It reads the number of q points from `_ph0/<prefix>.phsave/control_ph.xml`. [SRC]
* **Output in `save/`:**
  * `<prefix>.dyn_q<iq>[.xml]` and `<prefix>.dvscf_q<iq>`.
  * For PAW, `<prefix>.dvscf_paw_q<iq>`.
  * **For DFT+U (new): `<prefix>.dnsscf_q<iq>`, `<prefix>.dnsbare_q<iq>`, `<prefix>.occup`.**
  * `<prefix>.phsave/` (copied once).
  * `ifc.q2r[.xml]` if `<prefix>.fc[.xml]` exists.
  * It also copies `dyn0` → `dyn0.xml` in the XML case, and it **deletes `_ph0/<prefix>.q_*/*wfc*`**.
* **EPW reads** `dvscf_dir//prefix.dyn_q*`, `.dvscf_q*`, `.dvscf_paw_q*`, `.dnsscf_q*`, `.dnsbare_q*` and `ifc.q2r`
  (`io.f90`:932, 1855, 1929, 1995, 2054, 2657). [SRC]
* **rtsc's `collect_ph_for_epw` is functionally equivalent for non-+U runs.** It has the same sequential/parallel detection, XML
  detection, PAW detection, phsave copy and fc handling. Two differences:
  * It lacks the DFT+U files. Add them, or refuse when `lda_plus_u` is set.
  * It does not delete the wfc scratch files (harmless).
  In the sequential XML case, pp.py copies with the glob `prefix.dvscf*`, which could also match `dvscf_paw`. rtsc's exact names are
  safer. Nothing else needs to change.

---

## 6. Wannier quality check

**Files** [SRC `printing.f90` plot_band, lines 1834-1995; `bzgrid.f90` loadkmesh]
* **`band.eig`** is written in the CWD by `band_plot=.true.` together with `filkf` (and/or `filqf`):
  * Line 1: `&plot nbnd=  NW, nks=   NK /`.
  * Then for each k: one line `(10x,3f10.6)` with k in **cartesian 2π/alat**, and one line `(1000f20.12)` with **absolute** energies in eV
    (not shifted by E_F) for all `nbndsub` bands.
  * `phband.freq` is the analogue for phonons (meV, f14.4).
* **`filkf` format:** first line `NK crystal` (or `cartesian`), then `kx ky kz w` per line.
* Wannier90's own interpolation (`wdata(1)='bands_plot = .true.'` plus a `kpoint_path`, as in `EPW/examples/pb`) writes
  `<prefix>_band.dat` and `<prefix>_band.kpt` (gnuplot format) during the wannierize step. [SRC example; that W90 writes these in
  library mode is INF]
* Also parse `<prefix>.wout` for the final spreads and the disentanglement convergence.

**Recommended automatic metric** [INF]
1. Generate `NK ≈ 50-200` validation k points that are **not** on the coarse grid (for example a shifted MP grid, or random points).
   Write them as `filkf` (crystal).
2. Run pw.x `calculation='bands'` on the same list with enough `nbnd`.
3. Run epw.x with `epwread=.true., wannierize=.false., band_plot=.true., filkf=..., filqf=<1 point Γ>, efermi_read=.true.,
   fermi_energy=E_F`. Leave nkf and nqf unset.
4. For each k, take every DFT eigenvalue in W = [E_F − 1 eV, E_F + 1 eV] and match it to the nearest Wannier eigenvalue. Use the
   Hungarian method or sorted matching within the frozen window.
5. Report `max|ΔE|` and the RMS over W.
6. Accept if `max|ΔE| ≤ 20-30 meV` and `RMS ≤ 5-10 meV` with W inside the frozen window. Also require no Wannier band in W that has
   no DFT partner (a spurious band).
7. Tighten the tolerance to about ≤ degaussw/5 when the window is used for α²F.

---

## 7. ph.x optimized-tetrahedron route plus alpha2f.x (independent check)

**Selection** [SRC `phq_readin.f90`:464-471]
* `electron_phonon='lambda_tetra'` (elph_tetra=1): λ_qν by optimized tetrahedron double-delta integration on the nk1..3 grid.
* `electron_phonon='gamma_tetra'` (elph_tetra=2): γ_qν with tetrahedron step-function weights. λ is then derived from γ.
* Both set `trans=.false.`, so ph.x **reads** the dynamical matrices and dvscf from an earlier phonon run.

**pw.x requirements** [SRC]
* scf with `occupations='tetrahedra_opt'` or `'tetrahedra_lin'` (`tetra_type` 2 or 1).
* Blöchl `'tetrahedra'` is rejected ("DFPT with the Blochl correction is not implemented").
* `lambda_tetra` without tetrahedron occupations is an error ("electron_phonon and occupation are inconsistent").
* `K_POINTS automatic`. No `nosym` is needed: ph.x runs its own nscf on `nk1..3` with symmetry-aware tetrahedra (`setup_nscf.f90`:145).

**q grid**
* **Use `lshift_q=.true.`.** At q=Γ, `elph_tetra_lambda` sets the weights to 0 and prints "λ_qν is singular at q=Γ, you should use a
  shifted q grid" (`elph_tetra_mod.f90`). [SRC]
* With `lshift_q`, **q2r.x cannot be used**, so there is no Fourier interpolation of the λ_qν. [SRC INPUT_PH.def, q_points.f90]
  → Converging in q means repeating DFPT on denser shifted q grids. The existing unshifted coarse DFPT cannot be reused for this route.

**Inputs** (taken from `PHonon/examples/tetra_example`) [SRC]
* ph run 1: `ldisp=.true., lshift_q=.true., nq1-3, fildyn, fildvscf`.
* ph run 2: the same plus `electron_phonon='lambda_tetra', nk1=nk2=nk3=<dense>`.
* `alpha2f.x -in <same file as run 2>`, with this namelist appended **after** `&INPUTPH`:
  ```
  &INPUTa2F
   nfreq = 500
  /
  ```
  `nfreq` is the only variable.

**Outputs** [SRC elph_tetra_mod.f90, alpha2f.f90]
* **`<fildyn><iq>.elph.<iq>`**, for each irreducible q:
  * `(3f15.8,i8)`: xq, 3·nat.
  * `(6e14.6)`: ω²(ν) in Ry².
  * `     Tetrahedron method`, then `DOS = … states/spin/Ry/Unit Cell at Ef= … eV`.
  * For each mode `(5x,'lambda(',i5,')=',f8.4,'   gamma=',f8.2,' GHz')`.
  * **λ has only 4 decimals.** λ is set to 0 if ω² < 0.
* **`elph_dir/a2Fq2r.51.<iq>`**: the el-ph matrix on the q star. It is meant for q2r and matdyn, which are unusable with lshift_q.
* **`<prefix>.a2F.dat`** (prefix is the pw prefix):
  * Header ` # Frequency[Ry], a2F, DOS[Ry], Partial DOS(Atom 1), PDOS(Atom 2), ..`
  * Rows `(200e25.15)`: ω (Ry), α²F, phonon DOS, PDOS per atom.
  * The grid is ω_i = i·ΔE for i = 1..nfreq, with **ΔE = max ω/nfreq**. There is no margin above ω_max.
  * α²F comes from a q-space optimized-tetrahedron partial DOS of the phonons, weighted by λ_qν·ω_qν/2.
* **`<prefix>.McMillan.gp`**, plus stdout lines `lambda :` and `omega_ln [Ry] :` computed two ways (from α²F and from the direct sum).

**Constraints**
* `alpha2f.x` is **serial only** (`nproc > 1` gives an error). [SRC]
* It **stops on any imaginary frequency** (`read_lam`: "Imaginary frequency"). [SRC]
* It reads λ with `(a19,f8.4)`, so λ_qν carries ±5·10⁻⁵ quantization. [SRC]
* Ultrasoft pseudopotentials are fine: the official example uses `Al.pbe-n-rrkjus_psl`. [SRC]
* PAW, noncollinear and +U with `lambda_tetra` are not explicitly blocked in `phq_readin`, but they are untested. [INF]
* Known local issue: `reports/qe_eph_uspaw_segfault_repro.md` documents an el-ph US/PAW segfault in this build. Check it before using
  US or PAW here.

**Negative α²F is possible** [SRC `PW/src/tetra.f90`:394-427]
* For `tetra_type=2` (`tetrahedra_opt`) the smoothing matrix `wlsm` has negative entries (for example −38/1260 and −56/1260).
* So both the k-space double-delta weights (λ_qν, which are not clipped) and the q-space α²F histogram can be **locally negative**.
  `alpha2f.x` does not clip either one.
* `tetrahedra_lin` (wlsm = identity) gives standard linear-tetrahedron weights, which should be non-negative. [INF]
  → For a certified **lower or upper bound** use `tetrahedra_lin`, or clip and record the clipping explicitly.

---

## 8. Minimal example inputs

### 8.1 pw.x nscf for EPW
Full coarse grid, Γ first.
```
&control  calculation='nscf', prefix='X', outdir='./', pseudo_dir='~/pseudo/sg15' /
&system   ibrav=..., nat=..., ntyp=..., ecutwfc=..., occupations='smearing', smearing='mv',
          degauss=0.02, nbnd=<>= index(mu+3*sigma at Gamma) + n_excluded + margin>, nosym=.false. /
&electrons diago_full_acc=.true., conv_thr=1e-12 /
ATOMIC_SPECIES ... / ATOMIC_POSITIONS ...
K_POINTS crystal          <- rtsc nscf_kpoints_card(nk): full nk1*nk2*nk3 list, (0,0,0) first, weights 1/N
```

### 8.2 epw.in: SCDM, ephwrite, α²F only (route 2a)
Run one pass per electron smearing. Delete `X.a2f` first.
```
--
&inputepw
  prefix='X', outdir='./', dvscf_dir='./save'
  amass(1)=...
  elph=.true., epwwrite=.true., epwread=.false., epbwrite=.false., etf_mem=1, use_ws=.true.
  wannierize=.true., nbndsub=NW, num_iter=500, iprint=2
  auto_projections=.true., scdm_proj=.true., scdm_entanglement='erfc'
  scdm_mu=<E_F_abs + 2.0>, scdm_sigma=1.5          ! eV, absolute scale
  ! bands_skipped='exclude_bands = 1-5'           ! semicore, if any
  ! dis_froz_max=<E_F_abs + 1.0>                  ! optional inner window
  ephwrite=.true., eliashberg=.true., liso=.false., laniso=.false.
  iverbosity=2                                    ! lambda_aniso / lambda_pairs too
  fsthick=<>= 5*degaussw>, degaussw=0.05, ngaussw=0
  degaussq=0.05, delta_qsmear=0.05, nqsmear=10, nqstep=1000      ! nqsmear<=10 (format)
  eps_acoustic=0.1
  nstemp=1, temps=10                              ! eliashberg needs T>0; unused for a2F
  nk1=..,nk2=..,nk3=.., nq1=..,nq2=..,nq3=..      ! = DFPT coarse grids
  mp_mesh_k=.true.
  nkf1=..,nkf2=..,nkf3=.., nqf1=..,nqf2=..,nqf3=..   ! nkf multiple of nqf
/
```
**Later passes (other degaussw or other fine grids).** Set `wannierize=.false., epwwrite=.false., epwread=.true.` and keep everything
else. For a different `degaussw`, `ephwrite` must stay `.true.` (see §2a).

**Route 2b alternative** (several electron smearings in one run). Replace the Eliashberg keys with:
```
  phonselfen=.true., a2f=.true., delta_approx=.true., mp_mesh_k=.false.,
  nsmear=4, delta_smear=0.025, degaussw=0.025, efermi_read=.true., fermi_energy=<E_F_abs>
```
Files are `X.a2f.0<isig>.<T>`, with no cumulative columns and nqsmear ≤ 15.

### 8.3 Band check pass
```
  epwread=.true., epwwrite=.false., wannierize=.false., elph=.true.,
  band_plot=.true., filkf='kval.txt', filqf='gamma.txt', efermi_read=.true., fermi_energy=<E_F_abs>
  (no nkf*/nqf*)
```

### 8.4 Tetrahedron route
`scf.in` uses `occupations='tetrahedra_lin'` (or `tetrahedra_opt`) and `K_POINTS automatic`.
```
ph1.in                          ph2.in (= alpha2f input)
&inputph                        &inputph
 prefix='X', outdir='./',        prefix='X', outdir='./',
 fildyn='X.dyn', fildvscf='dv',  fildyn='X.dyn', fildvscf='dv',
 ldisp=.true., lshift_q=.true.,  ldisp=.true., lshift_q=.true.,
 nq1=4, nq2=4, nq3=4,            nq1=4, nq2=4, nq3=4,
 tr2_ph=1e-16                    electron_phonon='lambda_tetra',
/                                nk1=24, nk2=24, nk3=24
                                /
                                &INPUTa2F
                                 nfreq=1000
                                /
```
Run order: `pw.x < scf.in`, then `ph.x -in ph1.in`, then `ph.x -in ph2.in`, then **`alpha2f.x -in ph2.in`** (serial).

---

## 9. Required changes to rtsc helpers

These are listed only. Nothing was edited.

**`epw_inputs.py`**
1. Add `EPW61_KEYS` (qe-7.6, readin.f90) and `version="6.1"`, and make it the default when the binary banner says v6.1.
   Sniff the banner with `parse_qe_version` on a short run.
2. In the 6.1 writer, rename `eps_acustic` → `eps_acoustic` and **drop `kmaps`**. Both are fatal at 7.6 [RUN]. The current writer
   always emits both, so **every input it produces fails with EPW 6.1.**
3. Support SCDM: `auto_projections=True, scdm_proj=True, scdm_entanglement, scdm_mu, scdm_sigma`.
   * Relax the "proj required" refusal when SCDM is chosen.
   * Forbid `proj` together with `auto_projections`.
   * Require `scdm_mu` explicitly (no default of 0).
   * Assert that the nscf k list starts at Γ (`nscf_kpoints_card` already does).
4. Add an α²F-only mode (`eliashberg=True, liso=False, laniso=False`) and a phonselfen/a2f mode with `delta_approx=True`.
   * The phonselfen/a2f mode must force `mp_mesh_k=False`.
5. Add validators that mirror EPW's own checks:
   * `4*degaussw ≤ fsthick`; nkf a multiple of nqf.
   * `nqsmear ≤ 10` (route 2a) or ≤ 15 (route 2b), so the output does not wrap.
   * `ngaussw=0` recommended.
   * No mp_mesh_k with phonselfen; `auto_projections` ⇔ `scdm_proj`.
6. Make the run step delete or rename stale `<prefix>.a2f` before each Eliashberg run (find_a2f reuse) and record that in provenance.
   Use one directory per (grid, degaussw).
7. Record in provenance that `degaussw` and `fsthick` come from `ephmat/egnv`, not from the input, when `ephwrite=.false.`.
   Better: parse them from the stdout line `Electron smearing (eV) =`.

**`collect_ph_for_epw`**
8. Add the DFT+U files (`dnsscf_q`, `dnsbare_q`, `<prefix>.occup`), or refuse when `lda_plus_u` is set. It is otherwise equivalent to
   7.6's `epw_pp.py`.

**`parsers.py`**
9. `parse_epw_a2f`:
   * Accept rows that wrap across lines: read all numbers and reshape by 1+2N. Alternatively refuse N > 10 with a clear message.
   * Accept the optional `DOS (eV)` line.
   * Detect the a2f_iso mismatch (header N against the column count) and handle it as N_eff = 1.
10. Add `parse_epw_a2f_selfen(path)` for `<prefix>.a2f.0<i>.<T>`:
    * No header; 1+N columns; same footer.
    * Rebuild the cumulative λ by the rectangle rule and check it against `Integrated el-ph coupling`.
    * Take the electron smearing from the footer and T from the file name.
    * Treat `.a2f_tr.*` as transport and never use it for Tc.
11. Add a rounding envelope (±5e-8 per printed a2F value) to the provenance of `epw_a2f_to_spectra`.
12. Add `parse_alpha2f_dat(path)` for ph.x/alpha2f `<prefix>.a2F.dat`:
    * Ry units, grid i·ΔE, column 2 is α²F.
    * Flag negative values (optimized tetrahedron) and the ±5e-5 quantization of λ_qν.
    * Add a reader for `*.elph.*` that recomputes λ = Σ w_q λ_qν, using star weights from the dyn files.

**`band_coupling.py`**
13. `read_epw_freq`: support the **stream** layout of `ephmat/freq` (20-byte header, then (24 + 8·nmodes) bytes per q), and keep the
    sequential layout for 5.x. Detect the layout by the file size.

**New helper**
14. A Wannier band check:
    * Write `filkf` and a pw `bands` input.
    * Parse `band.eig` (cartesian 2π/alat, absolute eV).
    * Compute max|ΔE| and RMS over [E_F−1, E_F+1] eV, and gate the α²F protocol on the result.
