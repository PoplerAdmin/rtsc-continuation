# Multiband audit summary (transcribed by coordinator)
Overall: SOUND. Equations re-derived; Lemma 1.1, Lemma 4.1 (40-digit check, worst rel. diff 2.5e-32),
Prop. N (cert_neg.py rerun PASS, tiles exact, ends analytic; FAILs at a0 = 4.35 as required),
Theorem 1/2/3 assembly, Coulomb PSD condition, sharpness all SOUND.
Break attempts: 1620 random configs (B up to 8, DOS spread 1e4) + optimizers: min ratio a0*Tc/kappa
satisfies r >= 1.00000001 (truncation only raises r). PSD Coulomb tests stable.
Recipe notes: (a) kappa = rho(A) without DOS possible; (b) PITFALL: band "contributions" normalised to
the isotropic a2F give too-small kappa (non-conservative) -> must be documented/refused; (c) envelope
1/(1-delta) vs (1+delta) semantics; (d) data violating detailed balance: theorem not formally applicable.
Minor wording: sec 4.2(v) "certified maximum" is a crude bound; fullline agreement 1.3e-7; sec 4.6 heuristic.
