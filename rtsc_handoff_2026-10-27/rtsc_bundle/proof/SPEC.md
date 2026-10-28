# Shared specification: the model and the target theorem

Units: k_B = hbar = 1 (temperatures and frequencies both in kelvin).

## Data
* `a2F`: the Eliashberg spectral function, a finite positive Borel measure on (0, inf)
  (densities a2F(w) dw are a special case). Total weight
      A := a2F((0,inf)) = \int a2F(w) dw ,   0 < A < inf.
* `mu*` >= 0 (Coulomb pseudopotential), with an optional Matsubara cutoff N_c in {1,2,...} or N_c = inf.
* Temperature T > 0. Fermionic Matsubara frequencies w_n = (2n+1) pi T, n in Z.

## Coupling function
    lam(k) = lam_T(k) = \int 2 w / (w^2 + (2 pi k T)^2) d a2F(w),   k = 1, 2, ...
(lam(0) = 2\int a2F/w is the usual lambda; it cancels identically below and may be infinite.)

## Linearized isotropic Migdal-Eliashberg equations at temperature T
    wtilde_n = w_n + pi T sum_m lam(n-m) sgn(w_m)
    Dtilde_n = pi T sum_m [lam(n-m) - mu* theta(N_c - |m|...)] Dtilde_m / |wtilde_m|
With phi_n := Dtilde_n/|wtilde_n| this becomes, for n >= 0,
    [ (2n+1) + lam(0) + 2 sum_{k=1}^{n} lam(k) ] phi_n = sum_{m in Z} [lam(n-m) - mu*_m] phi_m .
The problem splits into an even-frequency sector (phi_{-n-1} = phi_n) and an odd-frequency
sector (phi_{-n-1} = -phi_n). On finitely supported real sequences x = (x_n)_{n>=0} define

  Even sector:
    Q^+_T(x) = sum_{n,m>=0} [lam(|n-m|) + lam(n+m+1)] x_n x_m
             - sum_{n>=0} [ (2n+1) + lam(0) + 2 sum_{k=1}^{n} lam(k) ] x_n^2
             - 2 mu* ( sum_{0<=n<N_c} x_n )^2
  Odd sector:
    Q^-_T(x) = sum_{n,m>=0} [lam(|n-m|) - lam(n+m+1)] x_n x_m
             - sum_{n>=0} [ (2n+1) + lam(0) + 2 sum_{k=1}^{n} lam(k) ] x_n^2
(the lam(0) x_n^2 terms cancel exactly, so only lam(k), k >= 1, enter.)

The associated symmetric operators M^±_T (Allen & Dynes 1975 "rho" formulation) are
self-adjoint, bounded above, with compact resolvent. Let
    rho^±(T) := sup_{x != 0 finitely supported} Q^±_T(x) / |x|^2 .

## Definition of Tc used here
T is called *unstable* if rho^+(T) >= 0 or rho^-(T) >= 0. Any nontrivial (l^2) solution of the
linearized gap equation at T makes T unstable. Define
    Tc := sup { T > 0 : T unstable }   (Tc := 0 if no T is unstable).
This Tc is >= the physical (highest-solution) Tc, so an upper bound on it is the stronger statement.

## Target theorem (the hypothesis)
There is a universal constant c (independent of the shape of a2F, of mu* >= 0 and of N_c) with
      Tc <= c * A ,
c is to be made as small as possible, proved rigorously (computer-assisted steps must use
rigorous interval/ball arithmetic and rigorous tail bounds), and shown to be (nearly) sharp.
Corollary (screening filter): Tc >= 300 K requires A >= 300 K / c.

## Rules for all workers
* Return concrete lemmas with proofs, equations, code, certificates or explicit counterexamples.
  No status reports, no "routine" claims for unproved global statements.
* Never assume the statement being proved (or an equivalent of it) inside its own proof.
* The machine has only 2 CPU cores shared by several workers: keep numerics light
  (matrices <= ~600, no long brute-force sweeps; use `nice`), and do not run jobs > 20 min.
* Work only inside your own subdirectory of /home/claude/tc/ (named in your brief).
