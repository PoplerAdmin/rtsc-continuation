# A rigorous upper bound on Tc in Migdal–Eliashberg theory

**Theorem.** In isotropic Migdal–Eliashberg theory, every phonon spectrum α²F ≥ 0, every Coulomb
pseudopotential μ* ≥ 0 and every Coulomb cutoff satisfy

    k_B Tc  ≤  A / 4.33148145  ≈  0.23086789 · A ,        A := ∫ α²F(ω) dω .

The constant is sharp to a relative error of 1e-8. The optimal universal constant c* satisfies

    0.2308678918 ≤ c* ≤ 0.2308678940 ,

By Lemma 4 the supremum is approached by single Einstein modes. The certified near-extremal mode
has Ω/T = 7.5748661. Float exploration (not certified) puts the optimum at Ω ≈ 7.5749·Tc, with
λ = 2A/Ω ≈ 1.1436. Uniqueness of the extremal spectrum is not proved.

**Corollary (screening filter).** Tc ≥ 300 K requires A ≥ 1299.44 K, i.e. A > 111.977 meV.

Prior art: as recalled, not re-checked against the paper, the value ≈ 0.2309 matches the bound
Leavens (1975) obtained by numerical optimisation, using the Allen–Dynes (1975) eigenvalue
formulation. What is new here is a complete proof: the reduction to Einstein spectra is proved
for arbitrary measures, and the infinite matrices and the frequency continuum are handled by
rigorous ball-arithmetic certificates. The result was independently audited. The status of the
Lean 4 formalisation of the reduction is given in §3.

Units: k_B = ħ = 1, so temperatures and frequencies are both in kelvin.

---

## 1. Model

α²F is any finite positive Borel measure on (0, ∞) with mass A ∈ (0, ∞). Densities are a
special case, and ∫α²F/ω = ∞ is allowed. For T > 0 and k ≥ 1 define

    g_T(k, ω) = 2ω / (ω² + (2πkT)²),        λ(k) = ∫ g_T(k, ω) dα²F(ω).

Linearising the isotropic imaginary-axis Eliashberg equations at T and setting
φ_n = Δ̃_n / |ω̃_n| gives an equation that splits into even- and odd-frequency sectors. On
finitely supported real sequences x = (x_n), n ≥ 0, the even sector is (μ* cutoff in block form,
N_c ∈ {1, 2, …, ∞}):

    Q⁺_T(x) = −Σ_n (2n+1) x_n² + Σ_{n≠m} [λ(|n−m|) + λ(n+m+1)] x_n x_m
              + Σ_n [λ(2n+1) − 2 Σ_{k=1}^{n} λ(k)] x_n² − 2μ* (Σ_{n<N_c} x_n)²

The λ(0) terms from the Z-renormalisation and from the kernel cancel exactly. The odd sector
Q⁻_T is the same expression with +λ(n+m+1) and +λ(2n+1) replaced by −λ(n+m+1) and −λ(2n+1),
and with no μ* term.

Two independent re-derivations from the full Eliashberg equations confirmed these forms
(audit D; brute-force comparison on the full Matsubara line agrees to 7e-9). They also reproduce
two known limits: the strong-coupling limit Tc → 0.1827 √λ Ω, and the weak-coupling prefactor
e^{−1/2} · 1.13 Ω e^{−(1+λ)/λ}.

When λ(0) = ∞, ω̃ itself is undefined, and the reduced form above is taken as the definition
(wG G5.1).

Define ρ^±(T) = sup_{x≠0} Q^±_T(x)/|x|². Call T **unstable** if ρ⁺(T) ≥ 0 or ρ⁻(T) ≥ 0, and set
Tc := sup{unstable T}, with Tc := 0 if no T is unstable. Foundations (wG, Theorem 2) show this Tc is at least the highest
temperature at which the linearised gap equation has a nontrivial solution. The only hypothesis
needed is that the series converges absolutely, which holds whenever |Δ̃_n| ≤ C|ω̃_n|. For a
finite cutoff in block form the two temperatures are equal.

## 2. Proof

Throughout, x is finitely supported and |x|² = Σ x_n².

**Lemma 1 (coupling).** g_T(k, ·) satisfies 0 < g_T(k, ω) ≤ 1/(2πkT), and g_T is strictly
decreasing in k. Hence λ(k) ≤ A/(2πkT) < ∞ and λ is decreasing in k. With
L := Σ_{k≥1} λ(k), the identity Σ_{k∈Z} 1/(w² + (2πk)²) = coth(w/2)/(2w) gives
L = ∫ (1/2T)(coth(ω/2T) − 2T/ω) dα²F < A/(2T).

**Lemma 2 (the odd sector is never unstable).** Write S_n = Σ_{k=1}^{n} λ(k) and
T_j = Σ_{k≥j} λ(k). Then for every n:

    Σ_{m≥0, m≠n} λ(|n−m|) = S_n + L,        Σ_{m≥0, m≠n} λ(n+m+1) = T_{n+1} − λ(2n+1).

So the row sum of the odd matrix is

    −λ(2n+1) − (2n+1) − 2S_n + (S_n + L) − (T_{n+1} − λ(2n+1)) = −(2n+1),

using L = S_n + T_{n+1}. The off-diagonal odd entries λ(|n−m|) − λ(n+m+1) are > 0, because
|n−m| < n+m+1 and λ is decreasing. A symmetric matrix with nonnegative off-diagonal entries and
row sums r_n satisfies xᵀMx = Σ r_n x_n² − Σ_{n<m} M_nm (x_n − x_m)² ≤ Σ r_n x_n². Hence

    Q⁻_T(x) ≤ −Σ (2n+1) x_n² ≤ −|x|²    for every α²F and every T.

So ρ⁻ ≤ −1. This was found independently by routes D, E1, G, H and C.

**Lemma 3 (Coulomb repulsion only lowers Q⁺).** For μ* ≥ 0, −2μ*(Σ_{n<N_c} x_n)² ≤ 0, so
Q⁺_T ≤ Q⁺_T|_{μ*=0}.

**Lemma 4 (affine in the spectrum, reduction to Einstein modes).** For a, w > 0 let
Q^E(a, w; x) be the μ* = 0 even form of the Einstein spectrum a·δ(ω − w) at T = 1, that is,
with λ(k) replaced by l(k) = 2aw/(w² + (2πk)²). Since A·g_T(k, ω) = 2(A/T)(ω/T)/((ω/T)² + (2πk)²),

    Q⁺_T|_{μ*=0}(x) = (1/A) ∫ Q^E(A/T, ω/T; x) dα²F(ω).

Proof: Q⁺|_{μ*=0} = C(x) + ∫ F_T(ω, x) dα²F, where C(x) = −Σ(2n+1)x_n² and F_T is linear in
the l(k). This is a finite sum of integrals of bounded functions (Lemma 1), and the identity
C(x) = (1/A) ∫ C(x) dα²F completes it. ∎

Consequently, if Q^E(a, w; x) ≤ −η'|x|² for a = A/T and all w > 0, then
Q⁺_T|_{μ*=0}(x) ≤ −η'|x|². This step is an exact inequality, not a stationarity argument.

**Lemma 5 (convexity in the coupling).** Q^E(a, w; x) = C(x) + a·G(w; x) is affine in a, and
C(x) ≤ −|x|². So if Q^E(a0, w; x) ≤ −η|x|² and 0 < a ≤ a0, then Q^E(a, w; x) ≤ −min(η, 1)|x|².
Proof: if G ≥ 0 the a-term is at most its value at a0; otherwise Q^E ≤ C ≤ −|x|². ∎

**Proposition 6 (certificate; computer-assisted).** Let a0 = 4.33148145 and η = 1e-9. For every
w ∈ (0, ∞) and every finitely supported x,

    Q^E(a0, w; x) ≤ −η |x|².

Proof (certificate wE1/CERT.md; audited in audit_E1):
- **Small w (w ≤ 1).** Unit-weight row sums give Q^E ≤ (a·Lang(w/2) − 1)|x|², where
  Lang(u) = coth u − 1/u is increasing. At a0 this is ≤ −0.2898|x|².
- **Large w (w ≥ 400).** A weighted Schur test with weights (8, 5/2, 1, 1, …) gives every row
  ≤ −0.300, uniformly in w.
- **Compact range [1, 400].** The range is covered by 126 dyadic cells that tile it exactly.
  On each cell the infinite form is bounded as follows:
  - the 150×150 leading block, Taylor-expanded in w with a rigorous second-derivative
    remainder over the whole cell;
  - plus a Schur-complement bound for the coupling to the tail: an exact ball-arithmetic Gram
    term over 1500 further columns and an analytic Cauchy–Schwarz bound beyond them;
  - plus an analytic tail-row bound 2T_{m+1} − (2m+1).
  Negative definiteness of the resulting 150×150 ball matrices is verified by an exact
  power-of-two rescaling, a float Cholesky proposal, and a residual bounded in Arb ball
  arithmetic.

All verdict-relevant arithmetic is done in Arb (python-flint 0.9.0, 128 bits). Every infinite
sum is bounded analytically. ∎

**Main theorem, proof.** Let T ≥ A/a0, so a := A/T ≤ a0.
1. By Proposition 6 and Lemma 5, Q^E(a, w; x) ≤ −η'|x|² for all w, with η' = min(η, 1) = 1e-9.
2. Lemma 4 gives Q⁺_T|_{μ*=0} ≤ −η'|x|², and Lemma 3 gives the same bound for every μ* and N_c.
3. Lemma 2 gives Q⁻_T ≤ −|x|².

Hence ρ^±(T) ≤ −η' < 0, so every T ≥ A/a0 is stable, and Tc ≤ A/a0 (in fact strictly, since
the unstable set is closed by wG G3). ∎

**Other Coulomb conventions.**
- **Column-only ("literal") cutoff.** Here the μ* term multiplies every row, so the operator is
  not symmetric. By wG Proposition 2′, any nontrivial solution forces ρ⁺|_{μ*=0}(T) ≥ 0. Proof
  sketch: above the μ* = 0 Tc, the operator K = D − B is coercive with an entrywise nonnegative
  inverse. Pairing the equation with ψ = K⁻¹e ≥ 0 gives s(1 + 2μ* Σψ) = 0, so s = 0 and then
  x = 0. The step above shows ρ⁺|_{μ*=0}(T) < 0 for T ≥ A/a0. Hence, under the absolute
  convergence hypothesis (H) of wG §G2, the highest temperature at which the literal-convention
  gap equation has a nontrivial solution is < A/a0.
- **N_c = ∞.** The μ* term drops out of the supremum (wG G1), so the μ* = 0 bound applies.

**Proposition 7 (sharpness; computer-assisted).** For the Einstein spectrum A·δ(ω − Ω) with
μ* = 0, A/T = 4.33148149, Ω/T = 75748661/10⁷ and an explicit finitely supported x (support
n < 300, produced by wB/cert_lower.py), Arb evaluation gives Q⁺_T(x) = 4.346e-10 ± 7e-38 > 0. So T is unstable and Tc ≥ A/4.33148149. The audit of wB
confirmed this, and wE1 certified the same at 4.3314815 with Q = 3.17e-9.

Hence 1/4.33148149 ≤ c* ≤ 1/4.33148145. ∎

## 3. Verification record

| Item | Who | Method | Result |
|---|---|---|---|
| Reduction (Lemmas 2–5) | found independently by routes A, B, C, D, H | — | identical |
| Proposition 6 (primary) | E1 | Arb + Cholesky residual, 126 cells | a0 = 4.33148145 |
| Independent certificate | B | Taylor + Cholesky residual, 178 cells | a0 = 4.3314814 |
| Independent certificate | A | positive test vectors (Schur), 161 intervals | c ≤ 0.2308685 |
| Independent certificate | E2 | Collatz–Wielandt vectors, 118 cells | a0 = 4.331481 |
| Independent certificate | H | Schur test vectors, 344 intervals | c ≤ 0.23087 |
| Analytic backstop (no eigen-numerics) | C | explicit tables | c ≤ 0.237 (hand: 1/4) |
| Audit of E1 | referee 2 | re-derivation, reruns, break tests | SOUND |
| Audit of B (+ independent re-certifier, interval Cholesky) | referee 1 | same | SOUND |
| Foundations (operator theory, conventions) | G | proofs | PROVED |
| Global logical chain (Lemma 4 identity at 40 digits, conventions, full rerun, sharpness) | referee 3 | re-derivation + reruns | SOUND, wording fixes applied |
| Lean 4 formalisation of Lemmas 2–5 + main theorem (finite-atomic spectra, 2π as a free parameter κ, Prop. 6 as hypothesis H(a0, η)) | F | Lean 4.35 + Mathlib | builds; 34 theorems, no sorry, axioms = propext, Classical.choice, Quot.sound only (re-checked by coordinator) |

**Break tests.** Both audited certifiers FAIL, as a sound certifier must, in each of these cases:
- at a = 4.33148149, 4.3314815 and 4.3314816, all just above the threshold;
- with deliberately corrupted matrix entries (1e-7 in one entry);
- with an index bug, l(n+m) in place of l(n+m+1);
- with a margin η set between the truncated and the true eigenvalue (this checks that the tail
  bound is actually used).

## 4. Scope

This is a theorem about the model: isotropic, Migdal (adiabatic), a given α²F, constant μ*, and
the linearised gap equation. Anisotropic or multiband pairing, non-adiabatic vertex corrections,
anharmonic renormalisation of α²F, and non-phonon pairing are outside it.

The computer-assisted steps assume that Arb / python-flint 0.9.0 is correct and that the
scripts' assertions are active. The scripts are included and run in about 2–6 minutes each.
