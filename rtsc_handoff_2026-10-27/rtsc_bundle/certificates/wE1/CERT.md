# wE1 — Rigorous certificate for the Einstein-spectrum reduced statement

## Result

Notation (units of T): `a = A/T`, `w = Omega/T`, `mu* = 0`,
`l(k) = 2aw/(w^2+(2 pi k)^2)` (k >= 1), `S_n = sum_{k=1}^n l(k)`, `T(j) = sum_{k>=j} l(k)`, `L = T(1)`.
M+ and M- are the matrices from the brief. They are the same as SPEC.md's Q+ and Q- with
lam(0) cancelled. `sharp.py` checks this numerically: the arb block and SPEC's literal float
matrix agree to 1e-13.

**Theorem (certified).** Let a0 = 4.33148145 and eta = 1e-9 (exact decimals). Then for every
w in (0, inf) and every finitely supported real x = (x_n), n >= 0:

    x^T M+(a0,w) x <= -eta |x|^2      and      x^T M-(a0,w) x <= -|x|^2 .

So R(a0, eta) holds. It also holds for every a in (0, a0] (see the Remark at the end).

**Sharpness (certified).** For a1 = 4.3314815 and w1 = 7.57486629486084 (a dyadic number), the
dyadic vector x in `sharp_x_4.3314815.txt` (support n < 300) gives

    x^T M+(a1,w1) x = 3.17349635946...e-9 > 0     (arb enclosure radius 3e-38).

Hence rho+(a1, w1) > 0, and R(a1, eta) fails for every eta >= 0.

**Consequence.**  `4.33148145 <= a* < 4.3314815`, where `a* := sup{a : R(a, eta) holds for some eta > 0}`.
The optimal Einstein-spectrum constant `c_E = 1/a*` therefore lies in
`(1/4.3314815, 1/4.33148145] = (0.2308678912.., 0.2308678940..]`. The gap has relative width
1.2e-8. A float extrapolation (not part of the certificate) gives a* ~= 4.33148148714
at w* ~= 7.574866.

A coarser companion certificate was also produced: a0 = 4.33148, eta = 1e-7, N = 150, K = 600,
116 cells (`cells_4.33148.log`). In addition, a1 = 4.331482 and a1 = 4.3314816 were certified
unstable.

## How to re-run

    cd /home/claude/tc/wE1
    ./run_all.sh                         # = ./run_all.sh 4.33148145 1e-9 4.3314815 150 1500
    ./run_all.sh 4.33148 1e-7 4.331482 150 600     # coarser, faster variant

The run takes about 6 min on one core. It prints `FINAL VERDICT: ...` and writes
`log_<a0>.txt`, `cells_<a0>.log` (one line per cell) and `cells_<a0>.json` (exact dyadic cell
endpoints). The individual pieces are `ends.py A0 ETA 1 400`,
`run_cells.py A0 ETA N K 1 400` and `sharp.py A1 300`.

Files: `core.py` (ball-arithmetic primitives, verified PD test), `cell.py` (one w-cell),
`run_cells.py` (adaptive cover), `ends.py` (analytic end ranges), `sharp.py` (sharpness).
`fexplore.py` and `astar.py` are float exploration only and not part of the proof.

All verdict-relevant arithmetic uses python-flint `arb` at 128 bits. Floating point is used
only to *propose* a shift and a Cholesky factor, and these are then checked rigorously
(Lemma 6). Decimal inputs (a0, eta, a1) are arb balls that contain the exact decimal.

## Proof structure

The range of w is split as (0,1] ∪ [1,400] ∪ [400,inf).
- Odd sector, all w: Lemma 9 (analytic).
- Even sector, small w: Lemma 7.
- Even sector, large w: Lemma 8.
- Even sector, [1,400]: 126 cells, each certified by Lemmas 2–6.

### Lemma 1 (sums)
`sum_{k in Z} 1/(w^2+(2 pi k)^2) = coth(w/2)/(2w)`. This is the partial-fraction expansion
`sum_k 1/(k^2+c^2) = (pi/c) coth(pi c)` with `c = w/(2pi)`. Therefore

    L = sum_{k>=1} l(k) = (a/2)(coth(w/2) - 2/w) = (a/2) Lang(w/2),    Lang(u) = coth u - 1/u .

Lang is increasing on (0, inf) because `Lang'(u) = 1/u^2 - 1/sinh^2 u > 0`. It also satisfies
0 < Lang < 1, so L < a/2.

l(k) is positive and decreasing in k >= 1, and `l(k) <= 2a/w`. The identity
`T(j) = L - S_{j-1}` holds exactly. Comparison with an integral gives

    T(j) <= l(j) + int_j^inf l(t) dt = l(j) + (a/pi) arctan(w/(2 pi j)).

`core.T_upper` uses the minimum of the arb upper bounds from these two expressions. Both are
evaluated with a ball w, so the bound is valid for every w in the ball.

### Lemma 2 (weighted Schur test)
Let K be symmetric with off-diagonal entries satisfying `|K_nm| <= Khat_nm`, and let h_n > 0.
For finitely supported x:

    |sum_{n!=m} K_nm x_n x_m| <= sum_n x_n^2 sum_{m!=n} Khat_nm h_m/h_n .

Proof: `2|x_n x_m| <= x_n^2 h_m/h_n + x_m^2 h_n/h_m`, then use symmetry. With h = 1 this
becomes the unit-weight row-sum bound `y^T E y <= sum_n y_n^2 sum_m |E_nm|`, which holds for
any symmetric E.

### Lemma 3 (even-sector row sums)
The off-diagonal entries `K+_nm = l(|n-m|) + l(n+m+1)` are >= 0, and
`D+_n = l(2n+1) - (2n+1) - 2 S_n`. For each n:

    sum_{m!=n} l(|n-m|) = S_n + L,      sum_{m>=0, m!=n} l(n+m+1) = T(n+1) - l(2n+1),
    D+_n + sum_{m!=n} K+_nm = 2 T(n+1) - (2n+1).

### Lemma 4 (tail block)
Fix N and write x = y ⊕ z, with y = (x_n), n < N, and z = (x_m), m >= N. Then
`Q+(x) = y^T A y + 2 y^T B z + z^T C z`. Here A is the leading N×N block,
`B_nm = l(m-n) + l(n+m+1) >= 0` for n < N <= m, and C is the tail block.

For a row m >= N of C, the unit-weight Schur sum is

    D+_m + sum_{m'>=N, m'!=m} C_mm' = -(2m+1) - 2S_m + (S_{m-N} + L) + (T(m+N+1) - l(2m+1)) + l(2m+1)
                                    <= -(2m+1) + T(m+1) + T(m+N+1) <= -(2m+1) + 2T(m+1).

The first inequality uses `S_m - S_{m-N} >= 0` and `L - S_m = T(m+1)`. Hence

    z^T C z <= -sum_{m>=N} c_m z_m^2,    c_m := 2m+1 - 2T(m+1).

### Lemma 5 (Schur complement / coupling)
Choose constants `0 < gamma_m <= c_m(w) - eta` for every w in the cell.
- For N <= m < N+K: `gamma_m = lower(2m+1 - 2 Tbar(m+1) - eta)`.
- For m >= N+K: `gamma_m = gamma_min := lower(2(N+K)+1 - 2 Tbar(N+K+1) - eta)`. This is valid
  because T(m+1) <= T(N+K+1).

Here Tbar is the Lemma 1 upper bound over the cell. By AM–GM,
`2 y^T B z <= sum_m [gamma_m z_m^2 + (b_m . y)^2/gamma_m]`, where b_m is column m of B.
Lemma 4 then gives

    Q+(x) + eta|x|^2 <= y^T (A + eta I + G1 + G2) y,
    G1 = sum_{m=N}^{N+K-1} b_m b_m^T/gamma_m,    G2 = sum_{m>=N+K} b_m b_m^T/gamma_m.

**G1** is computed as an arb matrix product. B1 has ball entries evaluated with w equal to the
whole cell, so `G1(w) ∈ G1_ball` for every w in the cell.

**G2** satisfies `G2 <= diag(d)`. By Cauchy–Schwarz,
`(b_m . y)^2 <= (sum_n B_nm)(sum_n B_nm y_n^2)`. The column sums are bounded by

    sum_{n<N} B_nm <= T(m-N+1) + T(m+1) <= T(K+1) + T(N+K+1)    (m >= N+K),

and the row tails by

    sum_{m>=N+K} B_nm = T(N+K-n) + T(N+K+n+1).

So `d_n = (Tbar(K+1)+Tbar(N+K+1))/gamma_min * (Tbar(N+K-n)+Tbar(N+K+n+1))`, rounded up.

### Lemma 5' (Taylor in w on a cell [wc-h, wc+h])
Each entry of A(w) is an integer combination of l(k) with k <= 2N-1, plus the constant
-(2n+1). Let `mu_k = sup_cell |l''(k)|`, where `l''(k) = -4aw(3b^2-w^2)/(w^2+b^2)^3` and
b = 2 pi k. mu_k is evaluated in arb with a ball w.

Taylor's theorem with Lagrange remainder gives `A(wc+t) = A(wc) + t A'(wc) + E(t)` with
`|E_nm| <= h^2/2 U_nm`. Here `U_nm = mu_|n-m| + mu_{n+m+1}` for n != m, and
`U_nn = mu_{2n+1} + 2 sum_{k<=n} mu_k`. By the unit-weight Schur bound,
`y^T E y <= sum_n rem_n y_n^2` with `rem_n = (h^2/2) sum_m U_nm`.

Define

    P(t) = A(wc) + t A'(wc) + diag(rem_n + d_n + eta) + G1(w).

Then Q+(x) + eta|x|^2 <= y^T P(t) y. For fixed w, P is affine in t, so P(t) is a convex
combination of P(+h) and P(-h). Checking both endpoints with the G1 ball and with A(wc), A'(wc)
as balls (exact centre wc, ball a0) therefore proves P(t) <= 0 for all |t| <= h and all w in
the cell.

### Lemma 6 (verified positive definiteness, `core.rigorous_negdef`)
Let S = -P(±h) be a ball matrix. The steps are:
1. Scale exactly by powers of two, `S' = D S D` with `D_ii = 2^{-round(log2 sqrt S_ii)}`. Then
   S > 0 if and only if S' > 0.
2. In float, take `s = lambda_min(mid S')/2 > 0` and a Cholesky factor `L` with
   `mid S' - sI ~= L L^T`. L has exact float entries.
3. Compute `E = S' - sI - L L^T` in arb. It encloses E for every member of the ball.
4. If `max_i sum_j |E_ij|_upper < s`, then `S' = L L^T + sI + E >= (s - ||E||_2) I > 0`.
   This uses `||E||_2 <= ||E||_inf` for symmetric E.

### Lemma 7 (small w, 0 < w <= 1)
By Lemma 3 and Lemma 2 with h = 1, `Q+(x) <= sum_n (2T(n+1) - (2n+1)) x_n^2`. Since
`T(n+1) <= L`, this gives `Q+(x) <= (2L - 1)|x|^2`. By Lemma 1, `2L - 1 = a Lang(w/2) - 1`,
which is increasing in w. For w <= 1 it is at most `a0 (coth(1/2) - 2) - 1 = -0.28984 < -eta`
(checked in arb by `ends.py`).

### Lemma 8 (large w, w >= W = 400)
Use Lemma 2 with weights h_0 = H0 = 8, h_1 = H1 = 5/2 and h_n = 1 for n >= 2. Bound the rows
with `l(k) <= 2a/W` and `T(j) <= L <= a/2`:
- Row 0: `D_0 = l(1) - 1`, `K_01 = l(1)+l(2)`, and `sum_{m>=2} K_0m = T(2)+T(3) <= a`. So
  `R_0 <= -1 + (2a/W)(1 + 2H1/H0) + a/H0`.
- Row 1: `D_1 <= l(3) - 3`, `K_10 = l(1)+l(2)`, and `sum_{m>=2} K_1m = T(1)+T(4) <= a`. So
  `R_1 <= -3 + (2a/W)(1 + 2H0/H1) + a/H1`.
- Row n >= 2: Lemma 3 plus the corrections `(H0-1)K_n0 + (H1-1)K_n1` give
  `R_n <= 2T(n+1) - (2n+1) + (H0-1)(l(n)+l(n+1)) + (H1-1)(l(n-1)+l(n+2)) <= a - 5 + (4a/W)(H0+H1-2)`.

At a0 the arb values are R_0 <= -0.4234, R_1 <= -1.1071 and R_n <= -0.3003. All are < -eta,
so `Q+(x) <= -0.30|x|^2` for all w >= 400.

### Lemma 9 (odd sector, all a > 0, all w > 0)
The off-diagonal entries `K-_nm = l(|n-m|) - l(n+m+1)` are > 0, because |n-m| < n+m+1 and l is
decreasing. Also `D-_n = -l(2n+1) - (2n+1) - 2S_n`. Using the sums from Lemma 3,

    D-_n + sum_{m!=n} K-_nm = -l(2n+1) - (2n+1) - 2S_n + (S_n + L) - (T(n+1) - l(2n+1)) = -(2n+1).

By Lemma 2 (h = 1), `Q-(x) <= -sum (2n+1) x_n^2 <= -|x|^2`. No numerics are needed.

### Compact range [1, 400] (even sector)
The interval is split into dyadic cells. The initial widths are 1/2 on [1,16], 2 on [16,64]
and 8 on [64,400]. Any cell that fails is bisected. For the tight certificate (N = 150,
K = 1500, eta = 1e-9) this gives 126 certified cells that tile [1,400] exactly (the tiling is
checked in the script). The narrowest cell has width 2^-11, near w ~= 7.575, and none failed.
Each cell checks P(+h) and P(-h) with Lemmas 5, 5' and 6. Per-cell data (scaled margin,
remainder rem_0, coupling d_0, time) is in `cells_4.33148145.log`.

### Sharpness (Lemma 10, `sharp.py`)
For any finitely supported x, `Q+(x)` is a finite sum of closed-form entries. It is evaluated
in arb at dyadic w1 and dyadic x, with a1 as a ball containing the decimal. The enclosure is
> 0, so rho+(a1, w1) >= Q+(x)/|x|^2 > 0.

### Remark (monotonicity in a)
`Q+_a(x) = a q_w(x) - sum (2n+1) x_n^2`, where q_w does not depend on a. Suppose a <= a0. If
q_w(x) >= 0, then `Q+_a <= Q+_{a0} <= -eta|x|^2`. Otherwise `Q+_a <= -|x|^2`. So R(a0, eta)
implies R(a, min(eta,1)) for all a <= a0.

## Limits of the certificate
- Everything stated above is fully certified. No step relies on floats without a rigorous
  check, and every infinite sum is bounded analytically (Lemmas 1, 4, 5).
- Not certified: the exact value of a*, beyond the bracket [4.33148145, 4.3314815), and the
  location w* of the maximum (that comes from float exploration only).
- Not covered: the reduction from general a2F, mu* and N_c to this Einstein, mu* = 0 statement.
  Per the brief, another team handles that.
