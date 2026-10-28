# Multiband Migdal–Eliashberg: a sharp Tc bound by the Perron root of the coupling-weight matrix

Units: k_B = ħ = 1 (temperatures and frequencies in kelvin). a0 := 4.33148145 is the certified
single-band constant of `final/PROOF.md` (Prop. 6). Everything below is about the *model*:
linearized imaginary-axis Migdal–Eliashberg theory, isotropic within each band (or any finite
k-point discretization of the anisotropic theory, see §6).

## 0. Results at a glance

**Theorem 1 (main).** Let bands i = 1..B have Fermi-level DOS weights N_i > 0 and Eliashberg
functions α²F_ij (finite positive Borel measures on (0,∞)) with N_i α²F_ij = N_j α²F_ji, in the
convention where the band-i mass enhancement is λ_i = Σ_j λ_ij. Let the Coulomb matrix be such
that U_ij := N_i μ*_ij is symmetric positive semidefinite (any block cutoff N_c ≤ ∞; scalar
μ* ≥ 0 with μ*_ij = μ* N_j/N is included, and so is μ* = 0). Put

    𝒜_ij := ∫ α²F_ij(ω) dω ,      κ := Perron root of 𝒜 = λ_max(K),  K_ij := sqrt(N_i/N_j) 𝒜_ij .

Then

    Tc  ≤  κ / a0  =  κ / 4.33148145 .

The constant is sharp to 1e-8 for every B and every DOS vector: equality is approached by any
spectrum with equal row sums Σ_j 𝒜_ij (this includes B = 1) concentrated at Ω ≈ 7.5749 Tc.
**Status: COMPUTER-ASSISTED.** The only computer-assisted input is the single-band certificate
(Prop. 6 of final/PROOF.md). Every new step is analytic. One new step, Proposition N, also uses
a one-dimensional Arb check of three elementary scalar inequalities (`cert_neg.py`).

**Ordering of the candidate functionals.** For every admissible 𝒜,

    𝒜_avg := Σ_ij N_i 𝒜_ij / Σ_i N_i  ≤  κ  ≤  Λ_A := λ_max(½(R + K))  ≤  max_i A_i ,

where A_i := Σ_j 𝒜_ij (row sums) and R := diag(A_i). Hence:
* **Valid** (rigorous) rejection functionals are κ (best), Λ_A, max_i A_i and the maximum column sum.
* The **isotropic/averaged** A_avg is **not** valid. The ratio Tc·a0/A_avg is unbounded
  (Counterexample 1). In an MgB2-like connected two-band model the true Tc beats the averaged
  bound by 9 % (Counterexample 5).

**Theorem 2 (weaker, simpler).** Tc ≤ Λ_A / a0. Its proof needs only Prop. 6 and no Prop. N.
Status: COMPUTER-ASSISTED (through Prop. 6 only).

**Theorem 3.** The odd-frequency multiband sector is never unstable. Status: PROVED.

**Screening rule (§7).** Reject a candidate for T_target when κ < 4.33148145 · T_target. For
300 K this means κ < 1299.44 K = 111.98 meV. κ is the Perron root of the band matrix
𝒜_ij = ∫α²F_ij dω = ½ Σ_qν w_q ω_qν λ_ij(qν).

---

## 1. Equations and symmetric formulation (Question 1)

Standard multiband (Suhl–Matthias–Walker / Eliashberg) imaginary-axis equations, linearized in Δ:

    ω̃_{n,i} = ω_n + πT Σ_{m,j} λ_ij(n−m) sgn ω_m
    Δ̃_{n,i} = πT Σ_{m,j} [ λ_ij(n−m) − μ*_ij θ(N_c−…) ] Δ̃_{m,j} / |ω̃_{m,j}| ,
    λ_ij(k) = ∫ 2ω dα²F_ij(ω) / (ω² + (2πkT)²) ,     N_i λ_ij(k) = N_j λ_ji(k).

Here Z_i ω_n = ω̃_{n,i}: the Z-renormalization of band i uses the **row sum** Σ_j λ_ij, while
the pairing uses the full **matrix**. Set φ_{n,i} = Δ̃_{n,i}/|ω̃_{n,i}|. For n ≥ 0,
|ω̃_{n,i}|/πT = (2n+1) + Σ_j [λ_ij(0) + 2 Σ_{k=1}^n λ_ij(k)]. The problem commutes with
n → −n−1 and splits into even (φ_{−n−1} = φ_n) and odd (φ_{−n−1} = −φ_n) sectors.

**Unit kernels.** For w > 0 let l_w(k) = 2w/(w²+(2πk)²), and on sequences n ≥ 0 define

    (B^±_w)_{nm} = l_w(|n−m|) ± l_w(n+m+1),   D_{w,n} = l_w(0) + 2 Σ_{k=1}^n l_w(k),   G^±_w = B^±_w − D_w ,
    C(y) = Σ_n (2n+1) y_n² ,   G_w(y,z) = yᵀ G⁺_w z  (l_w(0) cancels in G).

The single-band Einstein form of final/PROOF.md is Q^E_{a,w}(y) = −C(y) + a G⁺_w(y) (T = 1, mode
a·δ(ω−w)). λ_ij(k) = T⁻¹ ∫ l_{ω/T}(k) dα²F_ij(ω).

**Symmetric form.** Write W_ij := N_i α²F_ij. This is a symmetric matrix of positive measures.
Use variables x_{n,i} = sqrt(N_i) φ_{n,i}, so |x|² = Σ_i N_i Σ_n φ_{n,i}². The even sector is

    Q⁺(φ) = −Σ_i N_i C(φ_i) + (1/T) Σ_ij ∫ dW_ij(ω) [ φ_iᵀ B⁺_w φ_j − φ_iᵀ D_w φ_i ] − 2 Σ_ij U_ij S_i S_j ,
    w = ω/T,  S_i = Σ_{n<N_c} φ_{n,i}.

In x it is the symmetric matrix −(2n+1)δ + [√(N_i/N_j) λ_ij(|n−m|) + …] − diag(Z), the Allen–Dynes
ρ-operator. The odd sector Q⁻ is the same with B⁻ and no Coulomb term.
ρ^±(T) := sup_{x≠0, fin. supp.} Q^±/|x|². T is unstable if ρ⁺ ≥ 0 or ρ⁻ ≥ 0, and
Tc := sup{unstable T}, exactly as in SPEC.md.

*Check.* `fullline.py` builds the equations on the full Matsubara line, computing |ω̃| by direct
summation, for a 2-band, 2-mode spectrum. The critical coupling it finds agrees with the reduced
even form to 1e-6, the size of the Z-sum truncation (`fullline.log`).

**Lemma 1.1 (Z-split; the Z/pairing mismatch is a pair-breaking term).** Because W is symmetric,
Σ_ij W_ij φ_{n,i}² = Σ_ij W_ij φ_{n,i}φ_{n,j} + ½ Σ_ij W_ij (φ_{n,i}−φ_{n,j})². Hence, with μ* = 0,

    Q⁺₀(φ) = −Σ_i N_i C(φ_i) + (1/T) Σ_ij ∫ dW_ij G_w(φ_i, φ_j) − P(φ),
    P(φ) = (1/2T) Σ_ij ∫ dW_ij Σ_n D_{w,n} (φ_{n,i} − φ_{n,j})²  ≥ 0 .

The same holds for Q⁻ with G⁻. PROVED (algebra). The bilinear G_w contains no l(0), and
|G_w(y,z)| is bounded for finitely supported y, z, so the first integral is finite. If some
λ_ij(0) = ∫2α²F_ij/ω = ∞, the display is taken as the definition, with P ∈ [0,∞] (cf. wG G5.1
for B = 1). The physical reading: interband transfers at zero frequency act like interband
impurity scattering on a band-dependent gap.

## 2. Affinity and extreme points (Question 2)

For fixed T and φ, Q⁺₀ is affine in the matrix measure W, as in Lemma 4 of final/PROOF.md. The
admissible W form the convex cone of symmetric, entrywise-nonnegative matrix-valued measures. The
extreme rays of that cone are

    δ_Ω ⊗ E^(ii)    (one intraband mode)        and       δ_Ω ⊗ (E^(ij)+E^(ji))   (one interband mode, i ≠ j),

that is, **a single frequency with a single "edge" of coupling**, not a single frequency with an
arbitrary matrix. The Lemma-4 argument carries over verbatim to any *linear* normalization ℓ(W):
"ℓ ≤ a T ⇒ stable" holds for all W as soon as it holds on the extreme rays. A linear ℓ that is
valid on both kinds of extreme ray must charge every intraband ray at least 𝒜_ii. The best such
functional is therefore at least trace-like, Σ_i 𝒜_ii + …. Trace-like bounds grow with B even
when the bands are decoupled, where Tc = max_i Tc_i. So the naive extreme-point reduction is
valid but weak. The optimal functional κ is *convex* and not linear, and it needs a different
decomposition.

That decomposition is the identity in Lemma 4.1. In the Perron gauge, Q⁺₀ is a *positive
combination of single-band Einstein forms at the two couplings ±κ/T*, minus P. So the effective
extreme points are single-band Einstein modes again (as in the B = 1 theory), now together with
their sign-reversed partners.

## 3. Numerical test of candidates (Question 3)

`mb.py` computes the critical scale a* for a spectrum *shape* (the coupling at which the
even-sector top eigenvalue crosses 0 at T = 1). For a functional Φ, the claim "Tc ≤ Φ/a0" is
equivalent to the ratio a*·Φ(shape)/a0 being ≥ 1. The table gives minimum ratios over random
spectra with random DOS, random weights and random frequencies 2–40 T (`evidence.py`,
`evidence.log`, `scan_multi.py`):

| B, modes, trials | avg | max_i 𝒜_ii | κ | Λ_A | max row sum |
|---|---|---|---|---|---|
| 2, 1, 80 | 0.005 | 0.0 | 1.007 | 1.007 | 1.007 |
| 2, 3, 60 | 0.011 | 0.050 | 1.075 | 1.086 | 1.223 |
| 3, 1, 60 | 0.0014 | 0.0 | 1.015 | 1.023 | 1.085 |
| 3, 2, 40 | 0.016 | 0.0 | 1.144 | 1.231 | 1.377 |
| 4, 1, 40 | 0.040 | 0.0 | 1.107 | 1.306 | 1.838 |
| 4, 3, 25 | 0.042 | 0.40 | 1.369 | 1.593 | 2.293 |

Nelder–Mead minimization of the κ-ratio was run over 2 bands × (1, 2, 3) modes with free DOS,
weights and frequencies (`scan2.py`, `opt_multi.py`, `pert2.py`). It always converged to ratio
1 + O(1e-8), and only on two families: decoupled bands, and equal row sums (for example
N = (0.962, 0.038), 𝒜 = [[1, 0.0416], [1.041, 0]]). It never went below 1. This evidence pointed
to κ, which is then proved in §4. The table also shows that κ is the best of the candidates.
It also shows that κ is far from tight when the row sums are very unequal: for pure interband
coupling with N = (ε, 1−ε), the minimum over w of the ratio is 1, 1.13, 1.73, 3.37 and 10.1 for
ε = 0.5, 0.2, 0.05, 0.01 and 0.001 (`grid_inter.log`).

## 4. Proofs (Question 4)

### 4.1 The Perron-gauge identity

Let 𝒜 be irreducible (§4.4 handles the reducible case). Let v > 0 with 𝒜 v = κ v. This is the
Perron vector; κ equals λ_max(K), because K = N^{1/2} 𝒜 N^{−1/2} is symmetric and similar to 𝒜.
Set

    θ_i = φ_i / v_i ,   ν_i = N_i v_i² ,   dm_ij = v_i v_j dW_ij ,   ψ_ij = (θ_i+θ_j)/2 ,   χ_ij = (θ_i−θ_j)/2 .

Then Σ_j m_ij((0,∞)) = v_i N_i (𝒜v)_i = κ ν_i (equal "row sums" in the gauge).

**Lemma 4.1 (exact identity). PROVED.** With a := κ/T and w = ω/T,

    Q⁺₀(φ) = (1/κ) Σ_ij ∫ dm_ij(ω) [ Q^E_{+a,w}(ψ_ij) + Q^E_{−a,w}(χ_ij) ] − P(φ).

*Proof.* Start from Lemma 1.1.
1. W_ij G_w(φ_i, φ_j) = m_ij G_w(θ_i, θ_j).
2. Polarization for the symmetric G_w gives G_w(θ_i, θ_j) = G_w(ψ_ij) − G_w(χ_ij).
3. The parallelogram law gives ½(C(θ_i)+C(θ_j)) = C(ψ_ij)+C(χ_ij). Hence
   Σ_i N_i C(φ_i) = Σ_i ν_i C(θ_i) = κ⁻¹ Σ_ij m_ij^{tot} ½(C(θ_i)+C(θ_j)) = κ⁻¹ Σ_ij ∫dm_ij [C(ψ_ij)+C(χ_ij)].
4. Collect terms, using (1/T) = (1/κ)·a. ∎

`identity_check.py` confirms the identity numerically (3 bands, 3 modes, random x): the
difference is ≤ 3e-12 on values of order 6e3 (`identity.log`).

### 4.2 Proposition N (negative coupling is harmless). PROVED (analytic + 1-D Arb check)

For every w > 0 and finitely supported y: Q^E_{−a0,w}(y) = −C(y) − a0 yᵀG⁺_w y ≤ 0.

*Proof.*
(i) Let Y be the even extension of y to Z. Then Σ_{n,m∈Z} l(n−m) Y_n Y_m = 2 yᵀB⁺y. The symbol
σ(θ) = Σ_k l(k) e^{ikθ} = cosh((w/2π)(π−θ)) / sinh(w/2) for θ ∈ [0,2π]. This follows from the
Poisson formula Σ_k e^{ikθ}/(k²+c²) = (π/c) cosh(c(π−θ))/sinh(πc) (`check_lowerB.log` checks it
numerically).
(ii) cosh(cs) − 1 ≥ c²s²/2 ≥ c²(1 − cos s). Hence σ(θ) ≥ σ_π + γ(1 + cos θ), with σ_π = 1/sinh(w/2)
and γ = (w/2π)²/sinh(w/2).
(iii) The even restriction of the symbol 1 + cos θ is the tridiagonal form
yᵀTy = ½Σ_{n≥0}(y_n + y_{n+1})² + y_0² ≥ y_0². Therefore B⁺_w ≥ σ_π I + γ e_0e_0ᵀ.
(iv) Consequently −G⁺_w = D_w − B⁺_w ≤ diag(c + 2S(n)) − γ e_0e_0ᵀ, where c(w) = 2/w − 1/sinh(w/2)
and S(n) = Σ_{k=1}^n l_w(k). So

    Q^E_{−a0,w}(y) ≤ [a0(c−γ) − 1] y_0² + [a0(c + 2l(1)) − 3] y_1² + Σ_{n≥2} [a0(c + 2S(n)) − (2n+1)] y_n² .

(v) Every bracket is ≤ 0.
* For n ≥ 2: c + 2S(n) ≤ c + 2L = tanh(w/4) < 1 < 5/a0. This uses Lemma 1 of final/PROOF.md,
  L = ½(coth(w/2) − 2/w).
* For n = 0 and n = 1 the check is the scalar inequalities a0·g0 ≤ 1, with
  g0 = 2/w − (1+(w/2π)²)/sinh(w/2), and a0(c + 2l(1)) ≤ 3. These are verified for all w > 0 in
  Arb (`cert_neg.py`, verdict PASS). Analytic bounds cover the ends: g0 ≤ c ≤ tanh(w/4) ≤ w/4,
  g0 ≤ 2/w, c + 2l(1) ≤ tanh(w/4), and c + 2l(1) ≤ 6/w. On the compact parts [0.9, 8.7] and
  [3.4, 40], cell enclosures follow from monotonicity of 2/w, of 1/sinh(w/2) and of
  (1+(w/2π)²)/sinh(w/2). The last is decreasing because (w/π²) tanh(w/2) < w/π ≤ 1 + w²/4π².
* The certified maxima are 0.671 for a0·g0 and 2.994 for a0(c + 2l(1)). ∎

Numerically the critical negative coupling is 6.263 (at w ≈ 6.6), well above a0
(`negcoup.log`).

### 4.3 Theorem 1, proof

Let T > κ/a0, so a = κ/T < a0.
* Prop. 6 of final/PROOF.md gives Q^E_{a0,w} ≤ −η|·|² ≤ 0. Prop. N gives Q^E_{−a0,w} ≤ 0. Both
  forms are affine in the coupling and equal −C at coupling 0, so for every w:

      Q^E_{±a,w} = (1 − a/a0)(−C) + (a/a0) Q^E_{±a0,w} ≤ −(1 − a/a0) C .

* Insert this in Lemma 4.1, using P ≥ 0 and step 3 backwards:

      Q⁺₀(φ) ≤ −(1 − a/a0) κ⁻¹ Σ_ij ∫dm_ij [C(ψ_ij)+C(χ_ij)] = −(1 − a/a0) Σ_i N_i C(φ_i) ≤ −(1 − κ/(a0T)) |x|² .

* The Coulomb term −2Σ U_ij S_i S_j is ≤ 0 because U is PSD. So Q⁺ ≤ −(1−κ/(a0T))|x|².
* Theorem 3 gives Q⁻ ≤ −|x|².

Hence ρ^± < 0 for every T > κ/a0, so Tc ≤ κ/a0. ∎

**Solutions, not just forms.** Assume every λ_ij(0) < ∞. Then the even and odd operators are −C
plus a bounded symmetric operator (Toeplitz parts bounded by the symbol, bounded Z diagonal,
finite-rank Coulomb), and they are self-adjoint on dom(C) with finitely supported vectors as a
core. The form bound therefore extends to dom(C). Any ℓ² solution x of the linearized equation
satisfies Cx ∈ ℓ², so x ∈ dom(C), and so x = 0 for T > κ/a0. With some λ_ij(0) = ∞ we rely on
the conventions of wG, as in the single-band proof.

### 4.4 Reducible 𝒜

Because W_ij(ω) ≥ 0, 𝒜_ij = 0 forces W_ij ≡ 0. A reducible symmetric pattern therefore splits
into decoupled blocks. Q^± is then a sum over blocks, each block has its own Perron root
κ_b ≤ κ, and the bound follows blockwise. A band with no coupling at all is a block with κ_b = 0
and Q = −C.

### 4.5 Theorem 3 (odd sector) and Theorem 2

*Odd sector.* Split the Z term symmetrically: Σ_j W_ij φ_iᵀDφ_i = Σ_ij W_ij ½(φ_iᵀDφ_i + φ_jᵀDφ_j).
Then, with ψ = (φ_i+φ_j)/2 and χ = (φ_i−φ_j)/2,

    Q^± = −Σ N_i C(φ_i) + (1/T) Σ_ij ∫dW_ij [ ψᵀ G^±_w ψ − χᵀ(B^±_w + D_w)χ ] − Coulomb.

B^±_w is the even or odd restriction of a Toeplitz form with symbol σ ≥ 0, so it is PSD, and
D_w ≥ 0. Hence the χ term is ≤ 0. Lemma 2 of final/PROOF.md gives a·G⁻_w ≤ 0 for a > 0, so
G⁻_w ≤ 0 and Q⁻ ≤ −Σ N_i C(φ_i) ≤ −|x|². **PROVED.**

*Theorem 2.* For the even sector, Prop. 6 gives G⁺_w(ψ) ≤ C(ψ)/a0. Per Matsubara index n, the
sum Σ_ij W_ij (φ_{n,i}+φ_{n,j})²/4 equals x_nᵀ ½(R+K) x_n (in units of T). Hence
Q⁺₀ ≤ Σ_n (2n+1) x_nᵀ[−1 + ½(R+K)/(a0T)] x_n, and this is < 0 for T > Λ_A/a0. This route
needs no Prop. N, but it is weaker than Theorem 1 because Λ_A ≥ κ.

### 4.6 Sharpness. PROVED (via Prop. 7 of final/PROOF.md)

Take any 𝒜 with equal row sums A_i ≡ r and α²F_ij = 𝒜_ij δ(ω − Ω). Then v = 1 and κ = r. The
band-constant trial vector φ_{n,i} = y_n gives P = 0 and Q⁺₀ = (Σ N_i) Q^E_{r/T, Ω/T}(y). Prop. 7
then yields Tc ≥ r/4.33148149. So 1/a0 is optimal to relative 1e-8 for every B and every N.

Lemma 4.1 also shows that near-equality requires P ≈ 0 and χ ≈ 0, that is φ ∝ v and constant on
coupled bands. That points to equal row sums as the only extremal 𝒜. This characterization is
heuristic and not proved as a theorem.

### 4.7 Ordering of the functionals. PROVED (linear algebra)

* 𝒜_avg = uᵀKu with u = √N/|√N|, so 𝒜_avg ≤ κ.
* R − K ≥ 0, because xᵀ(R−K)x = ½ Σ W_ij(φ_i−φ_j)². Hence ½(R+K) ≥ K and Λ_A ≥ κ.
* Λ_A ≤ ½(max A_i + κ) ≤ max A_i.
* κ = ρ(𝒜) ≤ max row sum and ≤ max column sum of 𝒜 (Perron–Frobenius).

## 5. Counterexamples (Question 5)

All Tc values below are float computations with μ* = 0 (`ce.py`, `ce.log`). Here A = 1000 K and
Ω = 7.574866·A/a0, so a single band would have Tc = A/a0 = 230.87 K.

| # | model | true Tc | failing candidate (bound /a0) | κ/a0 |
|---|---|---|---|---|
| 1 | decoupled, N = (ε, 1−ε), 𝒜 = diag(A, 0); ε = 0.2 / 0.05 / 0.01 | 230.87 K | avg: 46.2 / 11.5 / 2.31 K, so Tc·a0/A_avg = 1/ε → ∞ | 230.87 (sharp) |
| 1b | as 1, connected (𝒜_12 = 50 K), ε = 0.05 | 219.17 K | avg: 12.70 K | 230.90 |
| 2 | pure interband, N equal, 𝒜 = [[0,A],[A,0]] | 230.87 K | max_i 𝒜_ii = 0, so every intraband-only functional fails | 230.87 (sharp) |
| 3 | equal row sums, N = (0.2, 0.8) | 230.87 K | none; κ = Λ_A = max row sum, the bound is attained | 230.87 (sharp) |
| 4 | **Coulomb not PSD**: μ*_12 = μ*_21 = 0.2, μ*_ii = 0, N equal, ω_c = 1000 K, **α²F = 0** | ≈ 7.8 K | every phonon functional (κ = 0) | 0 |
| 5 | MgB2-like (illustrative, not fitted): N = (0.43, 0.57), λ = [[1.0, 0.2],[0.151, 0.45]], Ω_σσ = 700 K, others 800 K | 81.13 K | avg: 74.32 K (Tc exceeds it by 9 %); max 𝒜_ii: 80.80 K | 86.53 |

Notes:
* In #4 the s± vector φ_{n,1} = −φ_{n,2} = 1/(2n+1) on n < N_c turns the interband repulsion
  into attraction. Instability holds iff 2μ Σ_{n<N_c} 1/(2n+1) ≥ 1, so Tc > 0 with A = 0. The
  hypothesis "N_i μ*_ij PSD" therefore cannot be weakened to "μ*_ij ≥ 0 entrywise". With a
  fixed index cutoff N_c, every T is unstable once 2μ Σ 1/(2n+1) ≥ 1.
* In #5 the isotropic-average single-band Tc is 71.6 K, while the multiband Tc is 81.1 K
  (×1.13).
* Other averages over bands fail like #1: min, harmonic or geometric means of row sums, and
  trace/B (decoupled bands with only one band active).
* Valid but weaker functionals: max row sum, Λ_A, and λ_max((𝒜+𝒜ᵀ)/2). The last is ≥ ρ(𝒜) by
  the numerical range, but it is not N-covariant.

## 6. Anisotropic Eliashberg

On a k-mesh, the anisotropic linearized equations used by EPW-type codes (states k on the Fermi
surface with weights N_k, kernels λ(k,k′; n−m) with N_k λ(k,k′) = N_k′ λ(k′,k), and Coulomb
μ*N_kN_k′/N_F, which is rank one and PSD) are exactly the multiband equations with B = number of
Fermi-surface states. Theorem 1 therefore applies to every such discretization (**PROVED**), with
κ = ρ(𝒜_{kk′}) ≤ max_k A_k, where A_k = ∫α²F_k(ω)dω is the k-resolved total coupling.

The continuum limit (κ = sup of the spectrum of the integral operator) is not separately proved
(**OPEN**, presumably routine by approximation).

## 7. Recommendation for a screening pipeline

**Quantity to compute.** Start from band- (or Fermi-sheet-) resolved Eliashberg functions,
using the convention in which band i's mass enhancement is λ_i = Σ_j λ_ij, with
λ_ij = 2∫α²F_ij(ω)/ω dω, and band DOS N_i. Compute

    𝒜_ij = ∫_0^∞ α²F_ij(ω) dω      (= ½ Σ_{qν} w_q ω_{qν} λ_ij(q,ν)  = ½ λ_ij ⟨ω⟩_ij,  ⟨ω⟩ = ∫α²F/∫(α²F/ω)),

in kelvin (1 meV = 11.604518 K).
1. Form W_ij = N_i 𝒜_ij and symmetrize W̄ = (W+Wᵀ)/2. The size of ‖W−Wᵀ‖/‖W‖ is a useful
   data-quality check.
2. Compute κ = λ_max( N^{−1/2} W̄ N^{−1/2} ). This number is always ≥ ρ(𝒜), so the
   symmetrization is conservative.
3. If the DOS weights are unavailable, κ ≤ max_i Σ_j 𝒜_ij (largest per-band total coupling) is a
   safe substitute.

**Rejection rule (safe within the model).**

    reject  ⇔  κ < 4.33148145 · T_target        (e.g. T_target = 300 K: κ < 1299.44 K = 111.977 meV).

The rule is valid for any frequency distribution and any number of bands. It holds for any
Coulomb matrix with N_i μ*_ij PSD, which includes a scalar μ* and μ* = 0, and for any block
cutoff.

**Margin for numerical error.** κ is 1-homogeneous and monotone in 𝒜. So if the computed 𝒜 can
be low by a relative δ entrywise, reject only when κ(1+δ) < a0 T_target.

**Fallback** when only λ_ij and the maximum phonon frequency ω_max are known: 𝒜 ≤ (ω_max/2)λ
entrywise, hence κ ≤ (ω_max/2) ρ(λ). Here ρ(λ) is the Perron root of the λ_ij matrix, or
λ_max of its √(N_i/N_j) symmetrization.

**Do not use** the single-band rule applied to the DOS-averaged A_avg = Σ_ij N_i𝒜_ij/N. It can
reject superconductors whose model Tc exceeds A_avg/a0 by any factor (Counterexample 1), and even
by 9 % for realistic two-band parameters (Counterexample 5). Do not use intraband-only
quantities either. If an effective μ*_ij matrix with N_iμ*_ij not PSD is used, no phonon-only
rule is valid (Counterexample 4).

**Scope.** The rule covers the Migdal approximation, the linearized gap equation, isotropy
within bands or a discretized anisotropic theory, and phonon-mediated pairing with a PSD Coulomb
pseudopotential. Like the B = 1 result, it uses Tc := sup{unstable T}, which is an upper bound on
the physical Tc.

## 8. Status table

| Item | Status | Evidence / files |
|---|---|---|
| Multiband even/odd forms, symmetric x-variables | derived; float cross-check against full Matsubara line, agreement 1e-6 | `fullline.py`, `fullline.log` |
| Lemma 1.1 (Z-split = pair-breaking P ≥ 0) | PROVED | algebra |
| Lemma 4.1 (Perron-gauge identity) | PROVED | algebra; `identity_check.py` (agreement 3e-12) |
| Prop. N (Q^E_{−a0,w} ≤ 0 ∀w) | PROVED (analytic + 1-D Arb check of scalar inequalities) | `cert_neg.py`, `cert_neg.log` (PASS); `check_lowerB.py`; `negcoup.log` (critical value 6.26) |
| Theorem 1: Tc ≤ κ/a0 | COMPUTER-ASSISTED (only through Prop. 6 of final/PROOF.md) | §4.3 |
| Theorem 2: Tc ≤ Λ_A/a0 | COMPUTER-ASSISTED (only through Prop. 6) | §4.5 |
| Theorem 3: odd sector stable | PROVED | §4.5 |
| Sharpness of 1/a0 for κ, all B and N | PROVED via Prop. 7 of final/PROOF.md | §4.6, CE 3 |
| Ordering avg ≤ κ ≤ Λ_A ≤ max row sum | PROVED | §4.7 |
| Counterexamples 1–5 | 1, 2 and 4 analytic; 1b, 3 and 5 float | `ce.py`, `ce.log` |
| Equality only for equal row sums | heuristic | §4.6 |
| Pointwise-optimal functional of (𝒜, N); κ is loose for very unequal row sums | OPEN | `grid_inter.log` |
| Continuum anisotropic limit; literal (column-only) Coulomb convention for B > 1 | OPEN | — |

Exploration-only scripts, which are not part of any proof: `scan2.py`, `scan_multi.py`,
`opt_multi.py`, `pert*.py`, `region.py`, `negan.py`. Every script runs in under a few minutes
on one core.
