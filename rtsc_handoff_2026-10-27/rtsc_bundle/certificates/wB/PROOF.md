# Route B — variational calculus on the space of spectra

Model, notation and the definition of Tc are those of `/home/claude/tc/SPEC.md`.
Everything below is in units k_B = hbar = 1.

## Main theorem

**Theorem.** Let a2F be any finite positive Borel measure on (0,inf) with total weight
A = a2F((0,inf)), let mu* >= 0 and N_c in {1,2,...,inf}. Then

    Tc <= A / 4.3314814  =  0.23086789660... * A      (so Tc <= 0.2308679 A).

The constant is sharp to about 2e-8 (relative). Define c* := sup over all (a2F, mu*, N_c) of Tc/A
(the best universal constant). Then

    0.2308678912 <= 1/4.3314815 <= c* <= 1/4.3314814 <= 0.2308678967 ,

and c* = 1/a*, where a* = min_{w>0} a_c(w) and a_c(w) is the critical coupling of the
Einstein spectrum (definition in Corollary 4). The supremum is attained by an Einstein spectrum
a2F = A delta(omega - Omega) with mu* = 0 and Omega/Tc = w*. Floating-point estimates
(guidance only, not certified): a* = 4.33148148714, w* = 7.57487, lambda = 2a*/w* = 1.14365,
Omega = 1.74879 A.

Corollary (screening filter): Tc >= 300 K requires A >= 300 K / c* >= 1299.4 K
(certified: A >= 300 x 4.3314814 K = 1299.44 K).

### Lemma status

| # | Statement | Status |
|---|-----------|--------|
| 1 | Sector reduction: rho^- <= rho^+_{mu*=0}, rho^+_{mu*} <= rho^+_{mu*=0} | PROVED |
| 2 | Q_T[a2F](x) is affine in a2F; extreme-point (Einstein) bound for rho_T | PROVED |
| 3 | Affine in the coupling a, convexity in a, threshold structure | PROVED |
| 4 | Reduction: if Q_{a0,w} <= 0 for all w > 0, then Tc <= A/a0 for every spectrum | PROVED |
| 5 | Small w: rho(a,w) <= -1 + a(coth(w/2) - 2/w) | PROVED (scalar check `cert_ends.py`) |
| 6 | Large w: weighted Hilbert–Schmidt bound, Q <= 0 for w >= 40, a <= 4.3314815 | PROVED (scalar check `cert_ends.py`) |
| 7 | Tail/Schur reduction of the infinite form to an N x N matrix | PROVED |
| 8 | Taylor-in-w enclosure and Cholesky certificate (method) | PROVED |
| 8' | Q_{a_lo,w} <= 0 for all w in [1.3, 40], a_lo = 4.3314814 | COMPUTER-ASSISTED (`cert_upper.py`) |
| 9 | Sharpness: Einstein a = 4.3314815, w = 7.5748661 is unstable at T = 1 | COMPUTER-ASSISTED (`cert_lower.py`) |
| 10 | Continuity of rho(a,w); c* = 1/a*; the min of a_c is attained | PROVED |
| 11 | Variational structure: subgradient (Bergmann–Rainer derivative), KKT, quasi-convexity of Tc | PROVED (KKT at the optimum: see 11c) |
| 12 | Uniqueness of the maximizer (w* unique; maximizing measure is a single delta) | OPEN |

---

## 0. Notation

Write mu for the measure a2F. At temperature T put
g_k(omega) := 2 omega / (omega^2 + (2 pi k T)^2) for k >= 1, so lam(k) = int g_k dmu.
lam(k) is nonincreasing in k, and 0 <= g_k <= 1/(2 pi k T). After the exact cancellation of the
lam(0) terms (SPEC), for finitely supported real x = (x_n)_{n>=0}:

    Q^±_T(x) = sum_{n != m} [lam(|n-m|) ± lam(n+m+1)] x_n x_m + sum_n (±lam(2n+1) - d_n) x_n^2 - [mu* term],
    d_n = (2n+1) + 2 sum_{k=1}^n lam(k),

where the mu* term, 2 mu* (sum_{n<N_c} x_n)^2 >= 0, appears only in Q^+. Write K^±_{nm} for the
coefficient of x_n x_m without the -d_n part. Then rho^±(T) = sup_x Q^±_T(x)/|x|^2.

**Einstein family in reduced units.** For mu = A delta_Omega one has
lam(k) = 2aw/(w^2 + 4 pi^2 k^2) with a = A/T and w = Omega/T. So Q^+_T, with mu* = 0, depends only
on (a, w). We write it Q_{a,w}(x) and its top value rho(a,w).

---

## 1. Lemma 1 (sector and Coulomb reduction) — PROVED

For every mu, T, mu* >= 0, N_c and finitely supported x:
(i) Q^+_{T,mu*}(x) <= Q^+_{T,0}(x);
(ii) Q^-_T(x) <= Q^-_T(|x|) <= Q^+_{T,0}(|x|).
Hence rho^+(T) <= rho^+_0(T) and rho^-(T) <= rho^+_0(T), where rho^+_0 is the even-sector value
with mu* = 0.

*Proof.* (i) The mu* term is -2mu*(sum x_n)^2 <= 0. (ii) For n != m the odd off-diagonal
coefficient is lam(|n-m|) - lam(n+m+1) >= 0, because |n-m| < n+m+1 and lam is nonincreasing.
So replacing x by |x| can only increase the off-diagonal part, and the diagonal part is unchanged.
Next, Q^+_0(y) - Q^-(y) = 2 sum_{n != m} lam(n+m+1) y_n y_m + 2 sum_n lam(2n+1) y_n^2 >= 0 for y >= 0.
Finally | |x| | = |x|. ∎

So T is unstable only if rho^+_0(T) >= 0. **From here on mu* = 0 and the even sector is understood.**
N_c plays no further role.

## 2. Lemma 2 (affinity in the spectrum; extreme points) — PROVED

For finitely supported x let

    q_x(omega) := sum_{n != m}[g_{|n-m|} + g_{n+m+1}](omega) x_n x_m + sum_n [g_{2n+1} - 2 sum_{k=1}^n g_k](omega) x_n^2 .

This is a finite sum of bounded continuous functions. Then

    Q_T[mu](x) = -sum_n (2n+1) x_n^2 + int q_x(omega) dmu(omega).            (2.1)

Hence:
(a) For fixed x, mu -> Q_T[mu](x) is affine. So rho_T(mu) = sup_x Q_T[mu](x)/|x|^2 is a
    **convex** function of mu (a supremum of affine functions).
(b) (Extreme-point principle.) If mu has mass A and nu = mu/A, then
    Q_T[mu](x) = int Q_T[A delta_omega](x) dnu(omega) = int Q_{A/T, omega/T}(x) dnu(omega),
    so

        rho_T(mu) <= sup_{w>0} rho(A/T, w).                                  (2.2)

*Proof.* (2.1) is the definition with lam(k) = int g_k dmu inserted. Every sum is finite, and every
integrand is bounded by sum |x_n x_m|/(2 pi T), so it is mu-integrable. (b) Write -sum(2n+1)x_n^2 =
int (-sum(2n+1)x_n^2) dnu and A q_x(omega) = Q_T[A delta_omega](x) + sum(2n+1)x_n^2. The scaling to
(a,w) is the identity 2A omega/(omega^2 + 4 pi^2 k^2 T^2) = 2(A/T)(omega/T)/((omega/T)^2 + 4 pi^2 k^2).
Then Q_T[mu](x) <= sup_omega Q_{A/T,omega/T}(x) <= sup_w rho(A/T, w) |x|^2. ∎

This is the global statement that the variational analysis needs. The maximization of the top
eigenvalue at fixed T and A is the maximization of a convex functional over the convex set
{mu >= 0, mu(0,inf) = A}. Its supremum is controlled by the extreme points A delta_Omega
(Einstein spectra). No stationarity analysis is required, and none would be sufficient (see 11d).

## 3. Lemma 3 (dependence on the coupling a) — PROVED

For fixed (w, x), Q_{a,w}(x) = -sum(2n+1)x_n^2 + a * qhat_x(w) is affine in a, and
Q_{0,w}(x) = -sum(2n+1)x_n^2 <= -|x|^2. Consequently, for 0 <= a <= a0:

    Q_{a,w}(x) = (1 - a/a0) Q_{0,w}(x) + (a/a0) Q_{a0,w}(x) <= -(1 - a/a0)|x|^2 + (a/a0) Q_{a0,w}(x).   (3.1)

So rho(.,w) is convex with rho(0,w) = -1. If rho(a1,w) >= 0, then rho(a,w) > 0 for all a > a1.
The unstable set {a >= 0 : rho(a,w) >= 0} is therefore an interval [a_c(w), inf) (closed by
Lemma 10), or it is empty (a_c = inf).

*Counterexample to pointwise monotonicity (sub-claim refuted).* Q_{a,w}(x) is **not**
nondecreasing in a for every fixed x. For x = e_1, Q_{a,w}(e_1) = lam(3) - 3 - 2 lam(1), and its
slope in a is (lam(3) - 2 lam(1))/a < 0. Monotonicity of rho in the coupling comes only from
convexity plus rho(0) = -1, as in (3.1).

## 4. Proposition 4 (reduction to the Einstein family) — PROVED

Suppose that for some a0 > 0, Q_{a0,w}(x) <= 0 for every w > 0 and every finitely supported x.
Then **Tc <= A/a0** for every spectrum of mass A, every mu* >= 0 and every N_c.

*Proof.* Let T > A/a0, so a := A/T < a0. By (3.1), Q_{a,w}(x) <= -(1 - a/a0)|x|^2 for all w and x.
By Lemma 2(b), Q_T[mu](x) <= -(1 - a/a0)|x|^2, so rho^+_0(T) <= -(1 - a/a0) < 0. By Lemma 1,
rho^+(T) < 0 and rho^-(T) < 0. So every T > A/a0 is stable and Tc <= A/a0. ∎

The theorem's upper bound follows from Proposition 4 with a0 = a_lo = 4.3314814, using
Lemma 5 for w in (0, 1.3], Lemma 8' for w in [1.3, 40] and Lemma 6 for w in [40, inf).

## 5. Lemma 5 (small w, Gershgorin) — PROVED

For any spectrum (in particular Einstein), with S := sum_{k>=1} lam(k):

    Q(x) <= (-1 + 2S)|x|^2.   For Einstein: 2S = a (coth(w/2) - 2/w) =: a Lgv(w/2).

*Proof.* Q(x) <= sum_n x_n^2 (K_nn - d_n + sum_{m != n} K_nm), because |x_n x_m| <= (x_n^2 + x_m^2)/2,
K is symmetric and K >= 0. The row sum is computed exactly:
sum_{m>=0} K_nm = [sum_{k<=n} lam(k) + S] + [S - sum_{k<=n} lam(k)] = 2S. The diagonal entry
K_nn = lam(2n+1) is the m = n term of sum_m lam(n+m+1). Hence the row-n coefficient is
-(2n+1) - 2 sum_{k<=n} lam(k) + 2S <= -1 + 2S. For Einstein,
sum_{k>=1} 2aw/(w^2 + 4 pi^2 k^2) = (a/w)((w/2)coth(w/2) - 1) (partial fractions of coth).
Lgv(x) = coth x - 1/x is increasing, since Lgv'(x) = 1/x^2 - 1/sinh^2 x > 0.

*Check* (`cert_ends.py`, arb): a_lo Lgv(0.65) = 0.913074... < 1. So Q_{a,w} <= 0 for all w <= 1.3
and a <= a_lo. (Corollary for any mu: rho <= -1 + 2S, so every spectrum with lam-sum S < 1/2 is
stable.)

## 6. Lemma 6 (large w, weighted Hilbert–Schmidt) — PROVED

Let D = diag(2n+1). For any spectrum, Q(x) <= (HS - 1) sum (2n+1) x_n^2 with
HS^2 := sum_{n,m} K_nm^2/((2n+1)(2m+1)) <= 4 sum_{k>=1} lam(k)^2 H_k / k, where
H_k = sum_{j=0}^{k-1} 1/(2j+1) <= 1 + (1/2) ln(2k).

For Einstein and every w >= W >= 2 pi:

    HS^2 <= (4a^2/pi^2) P(L) e^{-2L},   L = ln(w/(2pi)),
    P(L) = (7/16)L^2 + bL + c0,   b = c1 + (e - 1/2)/2,   c0 = c1(1/2 + e) + 1/32,
    c1 = 1 + (1/2) ln 2,   e = (2pi/W)^2.

*Proof.*
- Since d_n >= 2n+1 and K >= 0: Q(x) <= sum K_nm |x_n||x_m| - sum (2n+1)x_n^2 <= (||D^{-1/2} K D^{-1/2}|| - 1) sum (2n+1)x_n^2, and the operator norm is at most the HS norm.
- Use K_nm <= lam(|n-m|)[n != m] + lam(n+m+1) and (u+v)^2 <= 2u^2 + 2v^2, with the exact sums
  sum_n 1/((2n+1)(2n+2k+1)) = H_k/(2k) and sum_{n+m+1=s} 1/((2n+1)(2m+1)) = H_s/s (partial fractions).
  This gives HS^2 <= 4 sum_k lam(k)^2 H_k/k.
- For Einstein, lam^2 = 4a^2 w^2/(w^2 + 4pi^2 k^2)^2. Put phi(k) = (1 + (1/2)ln 2k)/k * w^2/(w^2 + 4pi^2 k^2)^2. phi is decreasing on [1, inf), so sum_{k>=1} phi(k) <= phi(1) + int_1^inf phi.
- phi(1) <= c1/w^2. Substitute k = wu/(2pi) and put u0 = 2pi/w <= 1. Then
  int_1^inf phi = w^{-2} int_{u0}^inf (c1 + L/2 + (1/2) ln u) / (u(1+u^2)^2) du.
- Exact antiderivative: int du/(u(1+u^2)^2) = ln u - (1/2)ln(1+u^2) + 1/(2(1+u^2)). So the integral from u0 to inf is at most L - 1/2 + u0^2.
- int_{u0}^inf ln u/(u(1+u^2)^2) du <= -L^2/8 + 1/16, using ln u/(u(1+u^2)^2) <= ln u/(4u) on [u0, 1] and <= ln u/u^5 on [1, inf).
- Collecting terms gives HS^2 <= (16a^2/w^2) P(L) with u0^2 <= e. Since 1/w^2 = e^{-2L}/(4pi^2), this is the stated bound.
- The bound is decreasing in L >= 0, because 2P - P' = (7/8)L^2 + (2b - 7/8)L + (2c0 - b) with both coefficients > 0. ∎

*Check* (`cert_ends.py`, arb, W = 40, a = a_lo; also passes for a_hi): 2b - 7/8 = 1.3428 > 0,
2c0 - b = 0.3666 > 0, HS^2 <= 0.80475 < 1. So Q_{a,w} <= 0 for all w >= 40.
(W = 30 fails with this bound: 1.19.)

## 7. Lemma 7 (tail/Schur reduction to a finite matrix) — PROVED

Fix N, L and (a,w), and split x = (u,v) with u = (x_0..x_{N-1}). Let B be the N x N section of M
(entries: B_nm = lam(|n-m|) + lam(n+m+1) for n != m, B_nn = lam(2n+1) - d_n). Let c_m(i) := K_im for
i < N <= m, and beta_m := 2m+1-2S > 0. Then

    Q(x) <= u^T (B + G) u + tau |u|^2,   G = sum_{m=N}^{N+L-1} c_m c_m^T / beta_m,
    tau = 4 q^2 / (6 (L-1)^2 (2N + 2L + 1 - 2S)),   q = a w/(2 pi^2).

*Proof.*
- Write Q = u^T B u + 2 u^T C v + v^T E v.
- Lemma 5's row-sum computation, applied to E, gives v^T E v <= -sum_{m>=N} (d_m - 2S) v_m^2 <= -sum beta_m v_m^2.
- Then 2 u^T C v - sum beta_m v_m^2 = sum_m [2 (c_m^T u) v_m - beta_m v_m^2] <= sum_m (c_m^T u)^2/beta_m.
- For m >= N+L, bound (c_m^T u)^2 <= |c_m|^2 |u|^2. Here |c_m|^2 <= 4 sum_{j >= m-N+1} lam(j)^2 (lam is decreasing) and lam(j) <= q/j^2.
- Next, sum_{j>p} j^{-4} <= 1/(3p^3) and sum_{p>=L} p^{-3}/3 <= 1/(6(L-1)^2).
- Finally, 1/beta_m <= 1/beta_{N+L}. ∎

Remarks.
- (i) A trace bound on all of G would be useless: its corner m ~ N, i ~ N gives about lam(1)^2/(2N), which is about 3e-3 at N = 80. G must be kept as a matrix. The top eigenvector is small in the corner, so B + G has top eigenvalue within about 1e-9 of the true rho at N = 80.
- (ii) The naive weighted Cauchy–Schwarz alternative G <= diag(...) costs about 3.6/N^2 on every row and also fails (float check in `fexp.py`).

## 8. Lemma 8 (Taylor-in-w enclosure and certificate) — PROVED (method)

Let I = [w0-h, w0+h] and t = w - w0. Entrywise Taylor expansion with Lagrange remainder gives
B(w) = B0 + t B1 + t^2 R, where B0 = B(w0), B1 = B'(w0) and R_ij = B''_ij(xi_ij)/2 lies in the ball
(1/2)B''(I) (R is symmetric). The derivatives are
lam' = 2a(c^2-w^2)/(w^2+c^2)^2 and lam'' = 4aw(w^2-3c^2)/(w^2+c^2)^3, with c = 2 pi k. Then for
w in I:

    Q(x) <= max_{s=±1} u^T [B0 + s h B1 + G(w)] u + sum_i u_i^2 (h^2 sum_j |R_ij| + tau),

because u^T(B0 + tB1)u is affine in t. G(w) lies in the ball matrix G(I), obtained by ball
evaluation at the ball w = I. For a ball matrix Mb with midpoints mid and radii rad, and any
symmetric M' in Mb, x^T M' x <= x^T sym(mid) x + sum_i x_i^2 (rowrad_i + colrad_i)/2.

It therefore suffices that P_s := -sym(mid(B0 + s h B1 + G(I))) - diag(radii + h^2 rowabs(R) + tau)
is positive semidefinite for s = ±1. This is certified as follows:
1. Compute a float Cholesky factor F of P_s - eps I, with eps = 1e-11.
2. Bound E := P_s - eps I - F F^T in arb; F is exactly representable, and P_s is symmetric with arb entries.
3. Check max_i sum_j |E_ij| < eps. Then P_s = FF^T + eps I + E ⪰ FF^T ⪰ 0, because E is symmetric and ||E||_2 <= max row sum.

### Lemma 8' — COMPUTER-ASSISTED

**Statement.** For a_lo = 4.3314814 (= 43314814/10^7 exactly), every w in [13/10, 40] and every
finitely supported x: Q_{a_lo,w}(x) <= 0.

**Script.** `python3 cert_upper.py 43314814 10000000 13/10 40 80 1000` (N = 80, L = 1000, arb at
128 bits). Adaptive bisection starts from intervals of width 1/4. The script asserts that the
accepted intervals tile [13/10, 40] exactly (rational endpoints). Output: `log_upper_43314814.txt`
(about 2 minutes on one core).

**What it certifies.** For each accepted interval, both P_+ and P_- are PSD, with all inputs
enclosed rigorously. The near-optimal intervals around w = 7.5749 have width about 2e-4 and
certified lambda_min(P) >= 6.4e-9 (float estimate; the certificate only needs
lambda_min > eps = 1e-11). Runtime: 115 s, 178 intervals, 0 failures.

The run at a = 4.331481 (N = 60, L = 300; `log_upper_4331481.txt`) also passed, but it was made
with an earlier version of `certify_nsd` that did not symmetrize mid(Mb). It is **superseded** and
is not used.

## 9. Lemma 9 (sharpness) — COMPUTER-ASSISTED

**Script.** `python3 cert_lower.py 43314815 10000000 75748661 10000000 300` (output `log_lower.txt`).

**What it certifies.** For the Einstein spectrum A delta(omega - Omega) with A/T = 4.3314815,
Omega/T = 7.5748661 and mu* = 0, the vector x (the float top eigenvector of the 300-section,
exactly representable) gives Q^+_T(x) = 3.17349697911e-9 (+/- 5e-38) > 0, evaluated exactly in
arb. So T is unstable, and Tc >= A/4.3314815 for this spectrum. Scaling Omega proportionally to A
gives the same for every A, so c* >= 1/4.3314815 = 0.23086789127.

## 10. Lemma 10 (continuity; c* = 1/a*; attainment) — PROVED

(a) rho(a,w) is continuous on [0,inf) x (0,inf).

*Proof.* lambda_max(B_N) <= rho <= lambda_max(B_N + G_N^full), by Lemma 7 with L = inf (tau = 0) and N > a. Here
||G_N^full|| <= tr G_N^full <= (4/(2N+1-2S)) sum_j j lam(j)^2, which tends to 0 as N -> inf locally
uniformly in (a,w). So rho is a locally uniform limit of the continuous functions
lambda_max(B_N(a,w)). ∎

(b) a_c(w) := min{a : rho(a,w) >= 0} is lower semicontinuous, because {rho >= 0} is closed.
Lemma 5 gives a_c(w) >= 1/Lgv(w/2) -> inf as w -> 0. Lemma 6, where the bound tends to 0 as
w -> inf for fixed a, gives a_c(w) -> inf as w -> inf. So a* = min_w a_c(w) is attained at some w*.

(c) **c* = 1/a*.**
- Lower bound: the Einstein spectrum with Omega = w* A/a* is unstable at T = A/a*.
- Upper bound: by (2.2), for any mu, rho_T(mu) >= 0 forces sup_w rho(A/T, w) >= 0. By (a), (b) and the end lemmas this sup is attained, so A/T >= a_c(w0) >= a*.
- So the optimal universal constant is the Einstein one. a_lo <= a* (Lemma 8' with Lemmas 3 and 5–6) and a* <= 4.3314815 (Lemma 9) give the enclosure stated in the theorem.

## 11. Variational structure — PROVED unless stated

(a) **Subgradient (rigorous Bergmann–Rainer derivative).** Let x be finitely supported with
|x| = 1 and Q_T[mu](x) >= rho_T(mu) - eps. Then for every mu' >= 0:

    rho_T(mu') >= rho_T(mu) - eps + int q_x d(mu' - mu).

This follows directly from (2.1). Where the top eigenvalue is simple, letting x run to the
eigenvector gives delta rho_T/delta mu(omega) = q_x(omega). The implicit function theorem then gives
the Bergmann–Rainer derivative delta Tc/delta a2F(omega) = -q_x(omega)/(d rho/dT) at T = Tc.

(b) **Second order.** rho_T is convex in mu, so along every direction sigma with sigma(0,inf) = 0
the second variation of rho_T is >= 0. Tc is **quasi-convex** in a2F: {mu : Tc(mu) <= T} is the
intersection over T' > T of the convex sets {rho_{T'}(mu) < 0}. Hence
Tc(sum p_i mu_i) <= max_i Tc(mu_i). Along a mixture segment t -> (1-t)mu_0 + t mu_1, Tc has no strict interior local maximum:
if t1 < t0 < t2 with Tc(t1), Tc(t2) < Tc(t0), quasi-convexity would give
Tc(t0) <= max(Tc(t1), Tc(t2)), a contradiction.
This is the global principle: Lagrange stationarity of Tc at fixed A, including the
Bergmann–Rainer condition "delta Tc/delta a2F = const on supp a2F", identifies candidates. These
candidates are generally not maxima. The maximum is located by Lemma 2(b), which reduces the
problem to the 2-parameter Einstein family.

(c) **KKT at the optimum.** Let x* be the top eigenvector at (a*, w*). For x* in finite section
(float, `fkkt.py`, guidance only), w -> Q_{a*,w}(x*) is maximal at w = w*, with value about 0, and
it is negative elsewhere (e.g. -4.6e-5 at w = 7.5 and 7.65, -0.087 at w = 12). That is the KKT
condition "q_{x*}(omega) <= q_{x*}(Omega*)": adding weight at any other frequency cannot raise
Tc to first order. Rigorously, rho(a*,w) <= 0 for every w. Otherwise continuity in a (Lemma 10a) would give
a_c(w) < a*. Hence Q_{a*,w}(x) <= 0 for all w and all finitely supported x, with equality
approached at w* along eigenvector approximants (rho(a*,w*) = 0 by Lemma 10). The
identification of the maximizing w as the unique 7.5749 is numerical only (OPEN 1).

(d) **Why local analysis is insufficient here.** The admissible set {mu >= 0, mass A} is
non-compact: mass can escape to omega -> 0 or omega -> inf. Lemmas 5 and 6 show that both escapes
are uniformly stable (rho -> -1). Without them, a stationary Einstein point would not exclude a
larger supremum approached at the boundary.

## 12. OPEN items (precise)

1. **Uniqueness of w*** (argmin of a_c). Lemma 8' certifies Q_{a_lo,w} <= 0 for all w. It does not
   certify a_c(w) > a* for w != w*. A certificate of rho(a_hi, w) < 0 for |w - w*| >= delta, of the
   same type, would close this but has not been run.
2. **Uniqueness of the maximizing measure.** If Tc(mu) = c* A, then mu/A should be delta at
   Omega = w* Tc. A proof needs the equality case of Lemma 2(b) for the l^2 eigenvector, which in
   turn needs (2.1) on the form domain (monotone convergence) and item 1. Not written out.
3. Not addressed: the improved constant c(mu*) for fixed mu* > 0. The universal bound is attained
   only at mu* = 0, since Lemma 1(i) is strict for x with sum x_n != 0.

No gap remains in the theorem as stated (upper bound A/4.3314814, sharpness A/4.3314815).

## Trusted base and files

The certificates rely on python-flint (arb 128-bit ball arithmetic, arb_mat products). NumPy is
used only to propose the Cholesky factor F and the test vector x, and their correctness is
re-checked in arb.

Files in /home/claude/tc/wB:
- `rig_common.py` — matrix builders, ball G-block, tail tau, NSD certificate.
- `cert_upper.py` — Lemma 8' driver, with the exact-tiling assertion.
- `cert_lower.py` — Lemma 9.
- `cert_ends.py` — scalar checks for Lemmas 5 and 6.
- `log_upper_43314814.txt`, `log_lower.txt`, `log_ends.txt` — certificate outputs.
- `fexp.py`, `fexp2.py`, `fhs.py`, `fkkt.py` — float exploration (guidance only).
