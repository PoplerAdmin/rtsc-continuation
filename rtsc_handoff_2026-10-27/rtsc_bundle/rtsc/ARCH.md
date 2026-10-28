# RTSC pipeline — architecture and soundness contract

Goal: a Python program, run on a large machine, that takes candidate compounds and outputs, for each
one, a verdict about whether its phonon-mediated Tc can reach a target T* (default 300 K). Every
"rejected" or "certified" verdict must be backed by machine-checked logic (Lean 4) plus rigorous
ball-arithmetic checks (python-flint Arb).

Background: the proved theorem in /home/claude/tc/final/PROOF.md says that in isotropic
Migdal–Eliashberg theory k_B Tc <= A/a0 with a0 = 4.33148145 and A = ∫α²F dω, for all spectra,
all mu* >= 0 and all cutoffs. Lean: /home/claude/tc/lean (namespace TcR1, e.g. TcR1.stable_of_ge,
TcR1.unstable_imp_lt). Model definitions: /home/claude/tc/SPEC.md.

## 1. Package layout (root /home/claude/rtsc)

    rtsc/
      pyproject.toml   README.md   ARCH.md   docs/CORRESPONDENCE.md
      rtsc/
        constants.py     certified constants with provenance (A0 = 4.33148145, ...)
        status.py        Status enum + Verdict record (the soundness contract, §3)
        spectrum.py      Spectrum = finite atomic measure (omega_K > 0, weight_K >= 0) + envelope + units
        certify/         rigorous Eliashberg engine (ball arithmetic)                       [team CERT]
        dfpt/            Quantum ESPRESSO / EPW workflow writers, runners, a2F parsers      [team DFT]
        sources/         candidate sources: local CIF/POSCAR dirs, Materials Project (mp-api),
                         CSP-output adapters (AIRSS/USPEX/CALYPSO dirs), seed lists      [team ORCH]
        stability/       hull (pymatgen), metallicity, magnetism, dynamic stability   [team ORCH]
        descriptors/     heuristic priority scores, used for ORDERING ONLY            [team ORCH]
        validity/        model-validity flags (Migdal ratio, anharmonic risk, multiband) [team ORCH]
        pipeline.py      stage orchestration, multiprocessing / MPI launch, resumable per-material
                         JSON records in <workdir>/<material_id>/record.json         [team ORCH]
        report.py        candidates.csv / candidates.json / per-material certificate files
        cli.py           `rtsc run -c config.toml`, `rtsc certify a2F.dat ...`, `rtsc demo`
      tests/   examples/

Heavy dependencies (pymatgen, mp-api, ase, phonopy, mace/chgnet, QE binaries) must be OPTIONAL
imports. `rtsc.certify` and `rtsc demo` need only numpy, scipy and python-flint.

## 2. Shared data types (put in spectrum.py / status.py; all teams use these exact names)

    @dataclass(frozen=True)
    class Spectrum:                      # alpha^2F as a finite atomic measure, energies in kelvin
        omega_K: tuple[float, ...]       # atom positions, each > 0
        weight_K: tuple[float, ...]      # atom masses, each >= 0; A = sum(weight_K)
        rel_envelope: float = 0.0        # assumption: true masses lie in [(1-e)w_j, (1+e)w_j]
        label: str = ""
        provenance: dict = field(default_factory=dict)   # file, smearing index, code version, ...
        # helpers: A(), lambda_(), omega_log(), from_grid(omega, a2f, unit) (midpoint weights), to_json()

    class Status(Enum):
        REJECTED_PROVEN      # Tc_model < T* proved (Lean-backed rule + Arb checks)
        CERTIFIED_CANDIDATE  # Tc_model >= T* proved
        UNDECIDED            # certified bracket [T_low, T_up] straddles T*
        EXCLUDED_ASSUMPTION  # excluded by a physical, non-proved criterion (hull, imaginary modes,
                             # magnetism, insulating). The reason, value and threshold are recorded.
        DEFERRED             # not computed (budget / ordering); never silently dropped
        ERROR

    @dataclass
    class Verdict:
        material_id: str; status: Status
        T_low_K: float | None; T_up_K: float | None      # certified bracket on Tc_model
        reasons: list[str]                    # human-readable
        rule_ids: list[str]                   # e.g. "R-UNIV", "R-CELLS", "C-RAYLEIGH" (§4)
        assumptions: list[str]                # e.g. "M-ISO-MIGDAL", "E-ENV(0.10)", "U-MU[0,0.15]"
        flags: list[str]                      # out-of-model risks: "NONADIABATIC", "ANHARMONIC_H", ...
        certificate: dict                     # data needed to re-check (vectors, cells, constants)

## 3. Soundness contract (non-negotiable)

S1. Only functions in rtsc.certify may produce REJECTED_PROVEN or CERTIFIED_CANDIDATE, and only
    through the single decision function `certify.classify(...)`. Its control flow mirrors the Lean
    definition `classify` in the Lean library PipelineLogic. Every boolean it consumes is a
    rigorous Arb check, so a "True" is a proof that holds whatever the floating-point rounding.
S2. Heuristic stages (descriptors, ML potentials, budgets) may only reorder candidates or mark them
    DEFERRED. They never produce REJECTED_PROVEN or EXCLUDED_ASSUMPTION.
S3. Physical existence filters produce EXCLUDED_ASSUMPTION with the exact criterion, and each can be
    disabled in the config.
S4. Every input material ends in exactly one Status. The run summary must account for all inputs.
S5. Uncertainty handling is conservative:
    - REJECTED_PROVEN must hold for every spectrum in the envelope, for every mu* in [0, mu_hi],
      and for every alternative spectrum the run produced (e.g. every DFPT smearing). This is a
      conjunction.
    - CERTIFIED_CANDIDATE must hold for every spectrum in the envelope, for every mu* <= mu_hi,
      and for every alternative spectrum.
    - Anything else is UNDECIDED.
S6. Verdicts are statements about the MODEL (isotropic Migdal–Eliashberg, linearised gap equation,
    block-form Coulomb cutoff, given α²F). Out-of-model risks are flags; they never upgrade a
    verdict. The README must say this plainly.

## 4. Rules (each has a Lean theorem; see docs/CORRESPONDENCE.md)

Notation: Stable(T) means rho+(T) < 0 and rho-(T) < 0. Tc_model := sup{T : not Stable(T)}.

- R-UNIV (universal ceiling): A_hi < a0 * T*  =>  for all T >= T*, Stable(T).
  Lean: TcR1.stable_of_ge. Here A_hi = (1+e)·A is the envelope maximum.
- R-ODD: the odd sector is always stable (row-sum identity). No check is needed.
- R-MU: Stable at mu* = 0 => Stable for every mu* >= 0 (Lemma 3). The mu* term is <= 0 and the
  instability form is decreasing in mu*.
- R-CELLS (per-material upper bound): for T in [T1, T2], g_T(k, ω) = 2ω/(ω² + (2πkT)²) is
  decreasing in T, so λ(k; T) lies in [λ(k; T2)·(1-e), λ(k; T1)·(1+e)].
  A positive test vector v (explicit head plus an analytic tail formula) whose every row satisfies
      sum_m Mhat_nm v_m <= -η v_n
  proves Stable on the cell. Here Mhat is the interval upper envelope: off-diagonals take their
  upper bound, and the diagonal takes upper(λ(2n+1)) - lower(2·sum_{k<=n} λ(k)) - (2n+1).
  This is the weighted Schur / Collatz–Wielandt test, valid because off-diagonals are >= 0.
  Rows n >= N need rigorous analytic tail bounds. Cells must tile [T_up, A_hi/a0] (contiguity is
  checked exactly), and R-UNIV covers T >= A_hi/a0. Conclusion: Tc_model <= T_up.
- C-RAYLEIGH (per-material lower bound): some x != 0 (finitely supported) has Q+_{T}(x) >= 0,
  evaluated as a minimum over the envelope (Q is affine in the atom masses, so the minimum is
  attained at the endpoint chosen by the sign of each atom's contribution) and at mu* = mu_hi with
  block cutoff N_c = N_c(T, omega_c). Then T is unstable, so Tc_model >= T.
  N_c rule: N_c(T) = number of n >= 0 with (2n+1)·π·T <= omega_c, where omega_c is a fixed
  frequency (config; default 10·max omega).
- classify: if R-UNIV fires => REJECTED_PROVEN. Elif the R-CELLS bound T_up < T* => REJECTED_PROVEN.
  Elif C-RAYLEIGH holds at some T_low >= T* => CERTIFIED_CANDIDATE. Else UNDECIDED. Across
  alternative spectra combine by S5.

## 4b. ERRATA from the Lean formalisation (binding; Lean library PipelineLogic, docs/CORRESPONDENCE_LEAN.md)

- **F1. Empty list.** classify([]) must return UNDECIDED, or ERROR if every DFPT alternative
  failed. Never return REJECTED or CERTIFIED with zero alternatives: the conjunction would be
  vacuously true.
- **F2. Mhat diagonal.** The Mhat diagonal is upper(λ(2n+1)) − lower(2·Σ_{k=1..n} λ(k)) − (2n+1).
  The sum starts at k = 1; λ(0) has cancelled and must never appear.
- **F3. Certified verdicts depend on the cutoff.** CERTIFIED holds for the configured block-cutoff
  rule N_c(T, omega_c). Record omega_c and the cutoff rule in Verdict.assumptions.
  REJECTED holds for every cutoff and every mu* >= 0.
- **F5. Envelope scope.** The envelope e covers atom weights only. Frequency-position
  uncertainty must enter as alternative spectra.
- **F6. Unresolved signs.** If an Arb ball for an atom's contribution F_j contains 0, use
  min((1−e)p_j F_j, (1+e)p_j F_j) (Lean sel_mul_eq_min, boxMin_le).
- **F7. Combined bracket.** The combined bracket over alternatives is
  [min of the T_low values, max of the T_up values].
- **D1. Cell tiling.** Cells tile [T_up, t_end] with exact rational/dyadic endpoints, together with
  an Arb check A_hi <= a0·t_end (R-UNIV then covers T >= t_end). Do not tile up to the Arb ball
  A_hi/a0.
- **D5. A = 0.** A = 0 (possible when e = 1) is stable trivially. Keep e < 1.

## 4c. Multiband rule (audited SOUND by an independent referee and a second numerical attacker; Lean libraries MultiBand + PipelineLogic.Multiband; opt-in; see /home/claude/tc/wMB/PROOF.md)

- **R-MULTI.** Given band-resolved coupling weights 𝒜_ij (K) and DOS weights N_i, compute
  κ = λ_max(N^{-1/2} W̄ N^{-1/2}), where W = N_i 𝒜_ij and W̄ = (W + Wᵀ)/2. If κ_hi < a0·T*, the
  material is REJECTED_PROVEN in the multiband (isotropic-per-band, discretised anisotropic) model.
  The Coulomb matrix must satisfy N_i μ*_ij PSD.
- Rejections that rely only on the single-band (DOS-averaged) A carry the assumption
  "M-ISO-SINGLEBAND" and the flag MULTIBAND_RISK. The DOS-averaged A can underestimate Tc in a
  multiband material: a counterexample in wMB gives a 9 % excess.
- Config `rejection_model = "isotropic" | "multiband"`. With "multiband", a rejection needs
  R-MULTI (or the fallback κ <= max_i Σ_j 𝒜_ij). Without band data the result is UNDECIDED.

## 5. DFPT stage (team DFT)

- Quantum ESPRESSO route: vc-relax at the target pressure, scf, ph.x with
  electron_phonon='interpolated' (or EPW), q2r.x, then matdyn.x with dos=.true., la2F=.true.
- Parse a2F.dos* for each smearing into a Spectrum (verify the units from QE's source/docs; convert
  to kelvin). Also parse λ and ω_log and check them against the Spectrum's own values.
- Each smearing gives an alternative spectrum (S5).
- Write inputs from templates, run via a configurable launcher (mpirun/srun), and resume.
- Convergence policy: report spread across smearings and q-grids. Never silently pick the best.
- EPW parser optional. Also a generic two-column a2F reader (omega, a2F) with explicit units.

## 6. Orchestration (team ORCH)

Stages:
1. source
2. dedupe
3. existence filters (EXCLUDED_ASSUMPTION)
4. heuristic ordering + budget (DEFERRED)
5. optional ML dynamic-stability pre-check (ordering/deferral only, S2)
6. DFPT
7. certify.classify
8. validity flags
9. report

- Config is TOML; every threshold is explicit and recorded in each Verdict.
- Parallelism: process pool over materials; each DFPT job gets its own MPI launch.
- Offline/demo mode: synthetic spectra clearly labelled synthetic. Never present invented numbers
  as data for a real compound.
