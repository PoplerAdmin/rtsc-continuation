# Audit of Route B (transcribed by coordinator from the referee's final message)

Overall verdict: SOUND. Tc <= A/4.3314814 rigorously established for all spectra, mu* >= 0, all N_c,
both sectors, conditional on python-flint 0.9.0 / Arb ball arithmetic. Sharpness c* >= 1/4.3314815
(and in fact >= 1/4.33148149) rigorous.

1. Reduction lemmas: SOUND. mu* term <= 0 for every N_c; odd off-diagonals >= 0 (g_k decreasing in k);
   Q+_0(y) - Q-(y) = 2 sum lam(n+m+1) y_n y_m >= 0 for y >= 0. Affine-in-measure lemma valid for general
   finite measures (0 <= g_k <= 1/(2 pi k T)). Non-strict certificate suffices: for T > A/a0,
   Q <= -(1 - a/a0)|x|^2 strictly; T = A/a0 may be unstable but Tc is a sup, so Tc <= A/a0.
2. End lemmas (w <= 1.3, w >= 40): SOUND; re-derived every inequality and closed form.
3. Tail lemma (N = 80, L = 1000): SOUND uniformly over each cell.
4. Cell certificate: SOUND. All inputs exact or enclosed; Taylor remainder evaluated over whole cell;
   matrix symmetrized before residual; residual argument correct; 178 cells tile [13/10, 40] exactly.
   Perturbing the factor by 1e-12 or using a shifted factor makes it reject.
5. Reruns reproduce. Break tests: certifier FAILS at a = 4.3314816, 4.3314815, 4.33148149, and with a
   deliberate 1e-7 error in B_00. Independent re-certifier audit_cell.py (no shared code, N = 90,
   different tail bound, interval Cholesky directly in arb) passes all 178 cells at a = 4.3314814 and
   fails at 4.3314816 and 4.33148149.
6. Sharpness: SOUND; fully in arb; Q(x) = 3.17e-9 +- 5e-38 > 0.

Minor (non-affecting): cell-width statement in PROOF.md inaccurate (min width 1/1024); the "6.4e-9
certified margin" is a float estimate, rigorous margin ~9.9e-12; quasi-convexity remark should use
both sectors; Bergmann–Rainer derivative section assumes unproved differentiability (unused);
superseded log correctly marked; docstring lemma names mismatch.
