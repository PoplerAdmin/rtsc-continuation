# Audit of certifier E1 (transcribed by coordinator from the referee's final message)
Overall verdict: SOUND. No mathematical or implementation error affecting the verdict.
1. Odd sector row-sum identity (rows of M- sum to -(2n+1)) and deduction Q- <= -|x|^2: SOUND (mpmath exact checks).
2. Small-w (w <= 1) and large-w (w >= 400, weights 8,5/2,1,1,...) lemmas: SOUND, uniform in w.
3. Truncation (150 block + analytic tail + 150x1500 coupling + analytic beyond 1650): SOUND; tail handling is
   necessary (eta = 1.2e-8, between truncated and true eigenvalue, correctly FAILS).
4. Cell certificate: SOUND; Taylor remainder over whole cell; endpoint check valid (affine in offset, NSD cone convex);
   power-of-two scaling exact; residual bounded by inf-norm in arb; 126 dyadic cells tile [1,400] exactly.
5. Reruns identical. Break tests FAIL as required: a = 4.3314816, 4.3314815, 4.33148149; eta = 1.2e-8;
   corrupted entries ((0,0)+1e-7, (0,1)+1e-6, index bug l(n+m), (5,5)+1e-4). Independent scan finds no w with rho+(a0,w) >= 0.
6. Code matrix equals SPEC M+- exactly (difference 1.14e-13 vs float SPEC-literal matrix).
Minor: w1 is the double 7942823/2^20 (not a short decimal); sharp_x parsed as doubles; Lemma 6 wording.
Not audited here: reduction from general a2F (audited separately in audit_B).
