# Adversarial review, round 3: α²F convergence protocol v3 (EPH_CONVERGENCE_PROTOCOL.md)

Reviewer: independent referee, 2026-10-27. Nothing in the protocol or code was changed.

**Read:** HANDOFF §0 and §3; ARCH §2–4c; DISCOVERY §3; reviews 1 and 2; protocol v3; `epw76_interface_spec.md`.

**EPW 6.1 source checked** (scratchpad `wsl/src/raw/`):
- `readin.f90`: the filqf, phonselfen, mp_mesh_k, shortrange and lfast_kmesh checks; the eps_acoustic default.
- `utilities_bzgrid.f90`: `load_mesh`, `loadqmesh_serial`, `loadqmesh_para`, `qwindow_wrap`, `qwindow`.
- `use_wannier.f90`: E_F (efnew), allocation, restart paths, the main q loop, where wf is set, the image sums.
- `selfen.f90`: `selfen_phon_q`.
- `io_selfen.f90`: `selfen_ph_write`.

**QE 7.6 source fetched** from raw.githubusercontent.com/QEF/q-e/qe-7.6:
- `PHonon/PH/elphon.f90`: `elphsum_simple`.
- `PHonon/PH/phq_readin.f90`: the electron_phonon modes, trans, and the tetra checks.

**rtsc code checked:**
- `pipeline.py`: `blocking_tags`, `apply_input_policy`, `annotate_model`, `_with_envelope`, the omega_c check.
- `certify/classify.py`: `resolved_omega_c` and the E-ENV assumption string.
- `dfpt/workflow.py`: the default grids and smearings.

Severity uses the scale of reviews 1 and 2:
- **CRITICAL:** can yield a wrong proved verdict.
- **MAJOR:** an unreliable or vacuous result, or a contract or plumbing violation.
- **MINOR:** everything else.

---

## 0. Source facts established in this round (so they are not re-litigated)

**F-1. filqf with phonselfen is accepted, and the q weights do not enter λ.**
- `readin` rejects filqf only with lfast_kmesh, lopt_w2b, band_plot (if nqf is also set) and a B field. There is no phonselfen restriction.
- `load_mesh` → `loadqmesh_serial` reads `nqtotf`, the coordinate type, xqf and wqf. If filqf is set it wins silently over nqf1..3.
- wqf is used only in the stdout print (`wt:`) and in a warning if the weights do not sum to 1. It is halved for noncollinear runs. λ_qν = γ_qν/(π ω² N_F) is per q. So the 1/N_q weighting is entirely on the protocol side, which is correct.

**F-2. qwindow with filqf.**
- `qwindow_wrap` sets `homogeneous = .FALSE.` when nqf1..3 = 0, and runs the non-homogeneous path. A q point is kept iff some k with a band within fsthick of ef has a k+q band within fsthick.
- With `selecqread=.false.`, `selecq.fmt` is always recomputed. A stale file is overwritten, not reused.
- **Bug, confirmed (v2 C3d).** `lambda_all` and `gamma_all` are allocated `(nmodes, totq, …)`, but `selfen_phon_q` writes `lambda_all(:, iq, …)` with the global index iq ≤ nqtotf, and `selfen_ph_write` reads iqq = 1..nqtotf. If totq < nqtotf, both are out-of-bounds accesses.
- **But** wf is filled only inside the q loop (`wf(nu, qind)`, use_wannier 1131) and is zero-initialised over nqtotf. A dropped q therefore has ω = 0 in the linewidth file, and the P1 assertion "ω > 0 for every mode" turns it into an ERROR.
- **The same catch-all covers every loop-skipping path** that P0 does not pin, as long as the parser enforces ω > 0 strictly: an `ephwrite` restart from a stale `restart.fmt`, `a2f_iso` (iq_restart = totq + 1) and `scatread`.

**F-3. The Fermi window is centred on the right E_F.** `use_wannier` sets `ef = efnew` (the fine-grid E_F at degaussw and ngaussw) before qwindow and before `selfen_phon_q`. The fsthick test (on ef) and the deltas (on ef0 = efnew, with nsmear = 1) use the same centre. The e^(−64) tail claim holds.

**F-4. Other EPW facts.**
- `mp_mesh_k` with phonselfen is a hard error.
- `shortrange` requires lpolar, which is pinned false. Both are closed.
- The default eps_acoustic is 0.1 cm⁻¹ (0.012 meV), not 5 cm⁻¹.
- With a2f=.false. and nimage = 1, `mp_sum(wf, inter_image_comm)` is the identity. Each image fills its own q, so even nimage > 1 is benign without a2f. v2 C3e is closed.

**F-5. The stdout per-q block exists and is fully checkable.**
- Header: `ismear = … iq = … coord.: x y z wt: w Temp: T`.
- Per-mode lines: `lambda___( ν )`, `gamma___`, `omega`.
- Total: **`lambda___( tot )=`** (f15.6).
- The string `lambda( sum )` named in P1 does not exist.
- The DOS and E_F lines are printed only for iq = 1.
- `We only need to compute N q-points` is printed by `qwindow_wrap`.

**F-6. ph.x `elphsum_simple` (QE 7.6).**
- **Deltas and E_F.** The deltas are always Gaussian (`ngauss1 = 0` is hard-coded). E_F and N_F are recomputed for each σ with `el_ph_ngauss` (default 0) on the k set of the ph run. σ = isig·el_ph_sigma, in Ry.
- **γ.** γ_ν = (π/2) Σ_k w_k Σ_ij δ̃δ̃ |⟨k+q,j| z_ν·dV |k,i⟩|², with the spin factor 2 in w_k. This is the same γ as EPW's `gamma_all`. The EPW linewidth file prints **2γ** in meV, so EPW γ[meV] = col4/2 = ph.x γ[GHz]/241.8.
- **λ.** λ = γ/(π N_F ω²), with N_F per spin. λ is set to 0 for ω < 20 cm⁻¹.
- **Output precision.** γ is printed `f8.2` in GHz and λ `f8.4`.
- **Bands.** All nbnd bands of the run are summed. EPW sums only the Wannier/fsthick bands.
- **trans.** With `electron_phonon='simple'`, `trans` is **not** turned off (`phq_readin`). The run is a full DFPT at the scf k grid. "Single-q DFPT on K_3" therefore means a linear-response calculation at 8·K_c unless dvscf is reused (`trans=.false.`, or `ldvscf_interpolate`).
- **Symmetrisation.** `el_ph_sum` is symmetrised like the dynamical matrix.

Conclusion: γ from ph.x 'simple' and EPW `delta_approx` γ are the same functional at the same σ, k grid and E_F convention. They differ only in band energies (DFT vs Wannier), in the band set, and in g and ω (direct vs interpolated). G6 is well posed in principle, but see M2.

---

## CRITICAL

**None found.** I could not construct a scenario in which v3, as written and implemented with the parser assertions taken literally, records an E-ENV(e) and alternative set that are wrong in the direction that matters, without also tripping a gate or a pin.
- The remaining risks to proved verdicts are estimate-quality risks inside the declared estimated assumption: three-point Richardson, and q under-resolution shared between EPW and tetra.
- Plumbing and spec gaps (M1, M3, M6) could reopen silent channels if implemented loosely. They must be fixed before any verdict is fed to the pipeline.

---

## MAJOR

### M1. Completeness and grid-identity checks are file-internal; the protocol's own sabotage test "K_1 = K_2 = K_3" would be ADMITTED
**Scenario.**
- If the three k runs are identical (the sabotage case, a script bug, a stale `filkf` that silently overrides nkf, or the same filqf list written for every Q level), then d_1 = d_2 = 0.
- The "converged" branch fires: both steps ≤ 0.2 %, pass "regardless of p".
- R is not defined in that branch. The natural implementation gives R = 0 (or p = log2(0/0)), so e falls to the 0.10 floor, with no q or k remainder at all.

**Related holes.**
- If d_2 = 0 exactly and d_1 ≠ 0, then p = ∞, capped at 4, and R = 0. This is plausible for small windows, because λ is printed with 5 significant digits.
- "Row count = N_q·n_modes" is tautological if N_q is read from the same file.

**Remedy.**
- Pre-register the expected N_k and N_q for every run, and verify them against EPW stdout:
  - `Size of q point mesh for interpolation`;
  - the k-mesh size line;
  - `Using q-mesh file: <hash-checked path>`;
  - `We only need to compute N q-points` with N = N_q.
- Require a stdout per-q block for **every** iq = 1..N_q with ismear = 1. Its printed coordinates must match the filqf list to 1e-5, and its `lambda___( tot )` must match the file.
- Require pairwise-distinct input hashes across the grids of a sequence.
- In the converged branch, define R = max(|d_1|, |d_2|).
- If d_2 is below print resolution, set R = max(|d_2|, the print quantum of X) and p undefined ⇒ fail unless the converged branch applies.
- Pin no `filkf`, and nqf1 = nqf2 = nqf3 = 0 (see M6).

### M2. G6 (direct DFPT γ check) is now a real gate, but as specified it tests the wrong q points and cannot be met reliably
1. **The point rule picks near-Γ points.**
   - Q_1 is shifted, so none of its points lies on Q_c. "The first three points not on Q_c" are simply the first three in list order: (1,1,1)/(4N), (1,1,3)/(4N), (1,1,5)/(4N). These are near Γ, on one line, with low γ, and the acoustic branches are the worst conditioned for ω.
   - The check never visits nesting vectors or soft branches, where Wannier g interpolation fails (the v1 C3 / v2 C2 concern).
2. **The 5 % per-mode threshold is below the noise G1 permits.**
   - A single-q γ at σ_min (5–20 meV) on K_3 is a sum over few Fermi-surface pairs.
   - G1 allows Wannier vs DFT band differences of up to min(30 meV, σ/2). A shift of σ/2 changes individual Gaussian weights by factors of roughly e^(±0.5…1).
   - So G6 will either flood ERROR or invite threshold tinkering after the fact. It also conflates band error with g error.
3. **Cost and workflow are unspecified.**
   - ph.x 'simple' keeps trans = .true. (F-6), so as written it is a DFPT at K_3 = 8·K_c.
   - The dvscf must come from a DFPT with **identical** scf settings to the coarse run (smearing `mv`, degauss 0.02 Ry, ecut, pseudopotentials). Otherwise G6 compares two different DFT models.
4. **Mode matching.**
   - Compare per degenerate set, not per mode, because EPW averages γ over modes with |Δω| < 1e-6 Ry.
   - |Δω|/ω ≤ 2 % on near-Γ acoustic modes will fail on ASR noise.
5. **Units and precision traps.**
   - EPW γ = col4/2 in meV; ph.x γ is in GHz, f8.2. Set el_ph_ngauss = 0, el_ph_sigma = σ_min in Ry, isig = 1.
   - Compare γ, not ph.x λ, which is zeroed below 20 cm⁻¹.

**Remedy.**
- Choose the points by a pre-registered data rule that targets coupling:
  - the three Q_1 points with the largest Σ_ν λ_qν at σ_min in the EPW run;
  - plus one generic low-symmetry point.
- Compare at 2·σ_min, or compare γ normalised by the respective double-delta phase space (ph.x prints `double delta at Ef`; EPW would need a one-line print), so that band error is separated from g error.
- Set the threshold from the measured G1 error, or keep 5 % on the normalised quantity.
- Specify the workflow: dvscf at the coarse-run settings, then elph on a K_3 nscf with `trans=.false.` (or `ldvscf_interpolate`).
- Compare by degenerate set, and exclude acoustic modes below a stated ω from the |Δω| test.

### M3. Frequency-position error and the DISCRETISATION tag: the §0 limit is stated but not plumbed
- **Positions are outside E-ENV.** G6 tolerates |Δω|/ω up to 2 % between interpolated and direct ω. That is a position uncertainty. ARCH F5 says it must enter as alternatives; the envelope covers weights only. v3 includes only the Δγ term in e. Note that R-UNIV is position-free, but R-CELLS and C-RAYLEIGH are not.
- **DISCRETISATION is not emitted automatically.** §0 relies on "DISCRETISATION" being an explicit assumption. `pipeline.annotate_model` emits `DISCRETISATION(...)` **only** if `Spectrum.provenance` carries `rule` or `discretisation`. A phself parser that does not set it produces verdicts whose assumption list does not name the discretisation at all.
- **"estimated" is lost.** `classify` writes `E-ENV(e) weights only (ERRATA F5)`, so the "estimated" qualifier promised in P5 is not recorded.

**Remedy.**
- Make the parser set `provenance["discretisation"] = "eph_protocol_v3:<grid>"`, with a unit test that the verdict carries `DISCRETISATION(eph_protocol_v3…)`.
- Add an `E-ENV-ESTIMATED(eph_protocol_v3, e)` assumption.
- For positions, either:
  - forward two extra alternatives with ω scaled by (1 ± δ_ω), where δ_ω = the G6 max |Δω|/ω (weights rebuilt from λ_qν·ω/2 consistently); or
  - state the position error explicitly as part of DISCRETISATION and record δ_ω.

### M4. The tetrahedron anchor (G4) compares two different DFT models, its own remainder is not in e, and its convergence is costly
- **Different scf models.** The EPW route uses dyn and dvscf from scf with `mv` smearing at 0.02 Ry (272 meV, `workflow.py`). The tetra route needs `occupations='tetrahedra_lin'` in the scf (`phq_readin`: elph_tetra requires ltetra). Screening, ω and dvscf therefore differ. This is largest at Kohn anomalies, which are the interesting materials. A 3 % gate then measures the scf-smearing model difference as well as the σ bias, and will ERROR on the materials that matter.
- **Missing remainder and compare point.** The tetra result's own P3 remainder R_tetra is not added to the G4 gap in e. The grid at which λ(σ_min,adm) is compared is unspecified.
- **Cost.** "Two refinements of q" for tetra means direct DFPT on shifted q grids up to about 4·Q_c, with no interpolation. That is the cost EPW is meant to avoid.

**Remedy.**
- Define the compare point as (K_densest, Q_3).
- Add R_tetra to the G4 term.
- Calibrate the scf-occupation model difference on the benchmarks (Al and Nb should agree to < 1–2 %), and record it as a known offset or a widened tolerance, fixed before any candidate runs.
- State the tetra grids explicitly, with a cost estimate.

### M5. Feasibility: run count × grid size invites DEFERRED/ERROR bias or quiet deviation
- With the default spacings (k 0.25 Å⁻¹, q 0.8 Å⁻¹), Al and H3S get K_c ≈ 11–12 and Q_c ≈ 4.
- **One EPW run** at K_3 × Q_3 is about 7–9·10⁵ k × 3.3·10⁴ q, roughly 2–3·10¹⁰ (k, q) pairs.
- **Run count.** About 18 runs for k admission, plus about 8 for q, plus about 4 for (K_prev, Q_3), each separate because the file holds smearing 1 only. K_4 = 16·K_c is 8× more.
- **Plus** the tetra DFPT on refined shifted grids, plus G6 at K_3.
- On the single GPU machine of HANDOFF §4 this is days to weeks per material. That is exactly the budget pressure that produces skipped gates or a smaller K_c.

**Remedy.**
- Measure Al end to end before committing.
- Consider the audited one-line `selfen_ph_write` patch to write every ismear. Note that EPW's extra smearings are **linear**: degaussw0 = degaussw + (ismear−1)·delta_smear. The ratio-2 σ list would need two runs (5, 10, 20 and 40, 80, 160 are not arithmetic) or a changed σ list. Changing the list is a pre-registration decision.

### M6. The P0 pins are a blacklist; several unpinned EPW keys change λ or the q/k set silently
Each of the following was checked in source:
- **`scissor ≠ 0` without `assume_metal`:** shifts every band above floor(nelec/2) in a metal, silently.
- **`eig_read`:** replaces the band energies.
- **`nsmear > 1` / `delta_smear`:** the file still holds smearing 1 only. The stdout has several blocks, so the per-q Σλ check could match the wrong one.
- **`nqf1..3 ≠ 0` together with filqf:** filqf wins, but `qwindow_wrap` takes the *homogeneous* path, whose `kpmq_map` assumes q lies on the k grid.
- **`filkf`:** overrides nkf.
- **`eps_acoustic`:** user-set.
- **`ephwrite`, `a2f_iso`, `scatread`, `plrn`:** loop-skipping. These are caught only by the ω > 0 catch-all of F-2.
- **`lifc`, `asr_typ`, `system_2d`, `lscreen`:** change ω or g.

**Remedy.** Make P0 a **whitelist**: the complete EPW input is hashed, and every key not in the whitelist must be at its EPW default, else ERROR. Pin `eps_acoustic` (with its value hashed), `nsmear = 1`, no `filkf`, and nqf1..3 = 0.

### M7. Pre-registration ambiguities that must be fixed before the first benchmark number exists
- **P5 parenthesisation.** Is it e = max(0.10, 3·max{R…, G4 gaps, G6}), or 3·max{R} + G4 + G6? The two readings differ by up to a factor of 3 on the anchor terms.
- **K_4 extension.** Which three grids, and which ratio, are used after it (K_2, K_3, K_4?)?
- **Window edges.** "(0, ω_max]": the ω_max of which spectrum? Edges must be fixed per material, for example from the coarse-DFPT max ω as for ω_c. Otherwise boundaries move between grids and atoms hop windows, causing spurious non-monotonicity.
- **(K_previous, Q_3).** These alternatives belong to neither the k-admission set (at Q_2) nor the q-admission set (at K_densest). They are extra runs; say so.
- **P1 string.** `lambda( sum )` should be `lambda___( tot )` (F-5).

---

## MINOR

- **m1. G1 does not count states.**
  - "max |ΔE| over states within E_F ± 4σ" does not detect a band or Fermi-surface sheet missing from the Wannier manifold. EPW would then lose that sheet's coupling and N_F, while G4 and G6 (all DFT bands) catch it only indirectly.
  - Require equal DFT and Wannier state counts in the window at every test k.
- **m2. No absolute E_F check.** Add a cheap gate: |efnew − E_F^DFT(σ)| ≤ σ/4 on a comparable grid. It catches nelec / nbndskip / `dis_win_min` errors directly instead of through G4 and G6.
- **m3. Near-Γ acoustic atoms.**
  - On dense shifted grids, the near-Γ acoustic atoms (ω ~ 0.1 meV at Q_3 for Pb) carry λ_qν ∝ |g|²/ω². Wannier-g ASR residue then dominates λ and ω_log, which are admission metrics, while barely affecting the Matsubara kernels certify uses.
  - Expect false ERRORs. Report the near-Γ shell's share of λ, and consider gating on λ(n = 1; T*) in addition to λ.
- **m4. k remainder at Q_3.** The k admission is measured at Q_2 but applied to alternatives at Q_3. The effect is small, but it is unmeasured; say so, or repeat one σ at Q_3.
- **m5. Validation and sabotage.**
  - "E_F +0.3 eV" still has no mechanism except `efermi_read=.true.`, which is a pin violation and so tests the pin check, not physics. Say this explicitly.
  - Expected reasons for pin sabotage (fsthick = 2σ, ngaussw = 1, restart) should be `PIN` errors, not `UNCONVERGED_EPH`.
  - There is still no tolerance against published converged λ and ω_log for Al, Pb and Nb (v1 m5).
- **m6. Print precision.** λ is printed E15.5, and the linewidth Γ E22.10 is unused. Small windows quantise. This links to the d_2 = 0 case in M1.
- **m7. G6 frequency comparison.** ph.x prints ω² (`w2`) in the elph file, and the frequencies themselves in the dyn file. Specify the source and the unit conversion.
- **m8.** HANDOFF §3.4 (input loss blocks rejection) is satisfied *more* strictly by ERROR. That is fine, but the existing `A-DROP-NONPOS` path in `dfpt/parsers.py` must not be reused by the phself parser with m = 0.

---

## Closure of round-2 findings

| v2 | status | note |
|---|---|---|
| C1 q control | **closed** (design) | There are three shifted q grids, the same Richardson rules, and R_q in e. Residual: q admission only at K_densest; q under-resolution shared with tetra (inherent). |
| C2 G6 circular / no g check | **partly** | It is now a real gate with a fixed threshold. It samples near-Γ points, its 5 % threshold is below the band-error noise, the workflow and cost are unspecified, and position error is not in e (M2, M3). G4 tetra now carries most of the interpolation-bias protection. |
| C3a tag not blocking | **closed** | Exclusions are ERROR, and no A-DROP is used. |
| C3b restart zeros | **closed** | `restart` is pinned. Other loop-skip paths are caught by ω > 0 (F-2). Whitelist anyway (M6). |
| C3c q weights | **closed** | The filqf coordinates are known and EPW ignores wqf for λ (F-1). The expected N_q must not be read from the file itself (M1). |
| C3d qwindow | **closed conditionally** | The out-of-bounds access is confirmed in source. A dropped q has ω = 0, so ERROR. Add the `We only need to compute N` = N_q assertion (M1). |
| C3e nimage | **closed** | nimage = 1 and a2f=.false. (F-4). |
| M1 false plateaus | **mostly** | The p ≥ 0.5 rule, sign agreement and the ratio-2 factor are fixed. R is undefined in the converged branch, and identical grids are admitted (M1). |
| M2 tautological anchor | **mostly** | The tetra gate is mandatory and directional. It compares a different scf model, R_tetra is missing from e, and the cost is high (M4). |
| M3 Γ / A-DROP | **closed** | There is no Γ on shifted grids. Near-Γ noise remains (m3). |
| M4 ω_c carrier | **open** (I1) | This is acknowledged as an integration task. |
| M5 D2 plumbing | **closed in spec** | The fields are guard-legal (stage `pipeline`) and dft_hash re-queues. I2 and I3 are still to be implemented. |
| M6 window ≠ per-atom | **partly** | §0 states the limit. DISCRETISATION and "estimated" are not plumbed (M3). |
| M7 validation | **partly** | Sabotage now runs at the bracket edge and code-level tests exist. K_1 = K_2 = K_3 would pass (M1), E_F +0.3 eV has no mechanism, and there are no published-value tolerances (m5). |
| m1–m9 | mostly | Open: the K_4 triple (M7), the `lambda( sum )` string (M7), and SOC (deferred). |

---

## Verdict

**v3 is ready to implement and to run benchmarks and sabotage tests as DRAFT diagnostics** (`EPH_PROTOCOL_DRAFT`). No CRITICAL item was found. Before the first benchmark *number* is produced, which is what keeps pre-registration honest, fix:
- M7: the P5 parenthesisation, the K_4 triple, the window edges, and the stdout string;
- M1: the converged-branch R, the expected-count, stdout-block and coordinate assertions, and distinct grids;
- the M2 point-selection rule and the G6 workflow, so that G6 results are interpretable rather than an ERROR flood;
- M6: the whitelist.

**Before any protocol output may feed pipeline verdicts:**
1. Sabotage passes (§V), including K_1 = K_2 = K_3 ⇒ ERROR (impossible under the current rule text).
2. M2 G6 redesigned (coupling-targeted points, a band-error-robust comparison) and shown on Al, Nb and MgB2 to pass, and to fail on Al with Q_c = 2³.
3. M3: DISCRETISATION and E-ENV-ESTIMATED emitted, with a test; position uncertainty either forwarded as ω-scaled alternatives or recorded explicitly.
4. M4: the tetra scf-model offset calibrated on the benchmarks and fixed in a new protocol version; R_tetra added to e.
5. M5: an end-to-end cost measurement on Al, and a decision on the all-ismear patch and σ list (a new version if the σ list changes).
6. I1–I5: ω_c carrier, ledger rows and tests, threshold hash in dft_hash, the ARCH S5 amendment reviewed, and an independent implementation audit.
