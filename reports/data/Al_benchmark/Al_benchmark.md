## Benchmark: Al (Al-fcc@0GPa[sigma2-10])

> Data: real DFPT data (convergence as stated by the producer of the input).

Input: `/home/user/rtsc_work/bench_out/Al_inputs_s2-10` (qe_dir); smearings: 9; envelope e = [0.1]; literature T* = 1.2 K; generated 2026-10-27T19:44:24Z (python-flint 0.9.0).

### 1. Spectra per smearing (parsed Spectrum; QE lambda file for comparison)

| # | broadening (Ry) | λ | ω_log (K) | A (K) | ω_max (K) | QE λ (lambda file) | Δλ | QE ω_ln (K) | Δω_log |
|---|---|---|---|---|---|---|---|---|---|
| 2 | 0.01 | 0.4152 | 305.9 | 67.53 | 471.5 | 0.4153 | -0.03% | 305.6 | +0.09% |
| 3 | 0.015 | 0.4127 | 307.8 | 67.32 | 471.5 | 0.4128 | -0.01% | 307.5 | +0.09% |
| 4 | 0.02 | 0.4116 | 308.0 | 67.18 | 471.5 | 0.4117 | -0.03% | 307.6 | +0.12% |
| 5 | 0.025 | 0.4113 | 307.5 | 67.07 | 471.5 | 0.4115 | -0.05% | 307 | +0.14% |
| 6 | 0.03 | 0.4118 | 306.8 | 67.08 | 471.5 | 0.412 | -0.04% | 306.3 | +0.16% |
| 7 | 0.035 | 0.4127 | 306.0 | 67.14 | 471.5 | 0.4129 | -0.05% | 305.5 | +0.18% |
| 8 | 0.04 | 0.4136 | 305.4 | 67.21 | 471.5 | 0.4139 | -0.07% | 304.8 | +0.19% |
| 9 | 0.045 | 0.4144 | 304.8 | 67.28 | 471.5 | 0.4147 | -0.07% | 304.2 | +0.21% |
| 10 | 0.05 | 0.4151 | 304.3 | 67.33 | 471.5 | 0.4154 | -0.08% | 303.7 | +0.22% |

Spread over smearings: lambda 0.4113–0.4152 (0.9%); omega_log_K 304.3–308 (1.2%); A_K 67.07–67.53 (0.7%); omega_max_K 471.5–471.5 (0.0%)

### 2. PROVED results (rtsc.pipeline.certify_spectra → rtsc.certify.classify; re-verified with rtsc.certify.recheck)

Statements about Tc_model (isotropic Migdal–Eliashberg, linearised gap equation, block Coulomb cutoff, given α²F, weight envelope e). T_up holds for every μ* ≥ 0; T_low holds for every μ* ≤ mu_hi under the stated cutoff. Combined over smearings per ERRATA F7.

| run | T* (K) | mu_hi | mode | status | combined T_low (K) | combined T_up (K) | recheck | certify wall (s) | recheck wall (s) |
|---|---|---|---|---|---|---|---|---|---|
| T300 | 300.0 | 0.15 | decide | **REJECTED_PROVEN** | - | 17.1497 | OK | 0.0 | 0.0 |
| Tlit_mu0.10 | 1.2 | 0.1 | bracket | **CERTIFIED_CANDIDATE** | 1.29688 | 10.3021 | OK | 13.3 | 1.0 |
| Tlit_mu0.15 | 1.2 | 0.15 | bracket | **UNDECIDED** | 0.577899 | 10.3021 | OK | 20.7 | 1.1 |

Per-smearing proved brackets, Tlit_mu0.10 (from recheck's re-derived evidence):

| # | T_low (K) | T_up (K) | classifyAlt | certify s |
|---|---|---|---|---|
| 2 | 1.35411 | 10.3021 | CERTIFIED_CANDIDATE | 1.6 |
| 3 | 1.32365 | 10.233 | CERTIFIED_CANDIDATE | 1.48 |
| 4 | 1.30489 | 10.1767 | CERTIFIED_CANDIDATE | 1.68 |
| 5 | 1.29688 | 10.1426 | CERTIFIED_CANDIDATE | 1.75 |
| 6 | 1.30165 | 10.149 | CERTIFIED_CANDIDATE | 0.401 |
| 7 | 1.31303 | 10.1747 | CERTIFIED_CANDIDATE | 1.71 |
| 8 | 1.32463 | 10.2041 | CERTIFIED_CANDIDATE | 1.49 |
| 9 | 1.33498 | 10.2289 | CERTIFIED_CANDIDATE | 1.66 |
| 10 | 1.34362 | 10.2504 | CERTIFIED_CANDIDATE | 1.54 |

Per-smearing proved brackets, Tlit_mu0.15 (from recheck's re-derived evidence):

| # | T_low (K) | T_up (K) | classifyAlt | certify s |
|---|---|---|---|---|
| 2 | 0.613224 | 10.3021 | UNDECIDED | 2.37 |
| 3 | 0.59356 | 10.233 | UNDECIDED | 2.37 |
| 4 | 0.582385 | 10.1767 | UNDECIDED | 2.59 |
| 5 | 0.577899 | 10.1426 | UNDECIDED | 1.17 |
| 6 | 0.581359 | 10.149 | UNDECIDED | 2.75 |
| 7 | 0.58832 | 10.1747 | UNDECIDED | 2.54 |
| 8 | 0.595794 | 10.2041 | UNDECIDED | 2.24 |
| 9 | 0.602471 | 10.2289 | UNDECIDED | 2.38 |
| 10 | 0.607915 | 10.2504 | UNDECIDED | 2.24 |

T* = 300 K: **REJECTED_PROVEN**, flags ['MULTIBAND_RISK'], rules ['R-UNIV', 'CLASSIFY-S5'].

omega_c used: 4715.39 K. Verdict files and re-check commands:

- `rtsc recheck /home/user/rtsc_work/bench_out/Al_s2-10/Al_verdict_T300.json --spectra /home/user/rtsc_work/bench_out/Al_s2-10/Al_spectra.json`
- `rtsc recheck /home/user/rtsc_work/bench_out/Al_s2-10/Al_verdict_Tlit_mu0.10.json --spectra /home/user/rtsc_work/bench_out/Al_s2-10/Al_spectra.json`
- `rtsc recheck /home/user/rtsc_work/bench_out/Al_s2-10/Al_verdict_Tlit_mu0.15.json --spectra /home/user/rtsc_work/bench_out/Al_s2-10/Al_spectra.json`

### 3. Float reference, not a verdict

Allen–Dynes (with f1·f2) and an independent float solution of the isotropic linearised Eliashberg equations on the same spectrum (nominal weights, same block cutoff ω_c). **These numbers are NOT proved and are NOT verdicts.**

| # | AD μ*=0.10 | AD μ*=0.15 | Eliashberg μ*=0.00 | Eliashberg μ*=0.10 | Eliashberg μ*=0.15 |
|---|---|---|---|---|---|
| 2 | 1.603 | 0.3898 | 7.762 | 2.271 | 1.201 |
| 3 | 1.562 | 0.3719 | 7.698 | 2.23 | 1.17 |
| 4 | 1.539 | 0.3629 | 7.65 | 2.203 | 1.152 |
| 5 | 1.531 | 0.3598 | 7.623 | 2.191 | 1.145 |
| 6 | 1.538 | 0.3632 | 7.63 | 2.197 | 1.15 |
| 7 | 1.552 | 0.3695 | 7.653 | 2.212 | 1.161 |
| 8 | 1.568 | 0.3761 | 7.678 | 2.228 | 1.172 |
| 9 | 1.581 | 0.382 | 7.701 | 2.243 | 1.183 |
| 10 | 1.593 | 0.387 | 7.72 | 2.255 | 1.191 |

† Matsubara sum truncated at N = 4000 (needed more): a LOWER estimate of the float Tc.

### 4. Literature (LITERATURE values, quoted from memory; see confidence)

| material | Tc exp (K) | P (GPa) | λ (literature) | ω_log (literature) | confidence |
|---|---|---|---|---|---|
| Al | 1.175 (standard tabulated value (e.g. Kittel; CRC Handbook): 1.175 K) | 0.0 | 0.38-0.44 (DFPT, e.g. Savrasov & Savrasov, PRB 54, 16487 (1996): ~0.44); ~0.38-0.43 from specific heat / tunnelling (Grimvall, 'The electron-phonon interaction in metals', 1981) | ~270-320 K (uncertain; typical DFPT values) | Tc: high. lambda range: medium. omega_log: low/uncertain. |
| Pb | 7.19 (standard tabulated value: 7.19-7.20 K) | 0.0 | ~1.55 from tunnelling inversion (McMillan & Rowell, PRL 14, 108 (1965); Carbotte, RMP 62, 1027 (1990) table: 1.55); DFPT ~1.1-1.2 without spin-orbit coupling and ~1.5-1.6 with SOC (Heid, Bohnen, Sklyadneva, Chulkov, PRB 81, 174527 (2010)) | ~56 K (Carbotte 1990 table: 4.83 meV; medium confidence) | Tc: high. tunnelling lambda 1.55: high. DFPT no-SOC value: medium. omega_log: medium. |
| Nb | 9.25 (standard tabulated value: 9.25-9.26 K) | 0.0 | ~1.0-1.3 (DFPT; e.g. Savrasov & Savrasov 1996 report ~1.26; tunnelling/specific-heat estimates ~1.0) | ~150-200 K (uncertain) | Tc: high. lambda range: medium. omega_log: low/uncertain. |
| MgB2 | 39.0 (Nagamatsu et al., Nature 410, 63 (2001): 39 K) | 0.0 | isotropic (DOS-averaged) lambda ~0.7-0.9 (e.g. Kong, Dolgov, Jepsen, Andersen, PRB 64, 020501(R) (2001): ~0.87; Choi, Roundy, Sun, Cohen, Louie, Nature 418, 758 (2002): ~0.73 with anharmonic phonons) | ~60-75 meV, i.e. ~700-870 K (uncertain) | Tc: high. isotropic lambda range: medium. Size of the isotropic underestimate: medium (order of magnitude). omega_log: low. |
| H3S | 203.0 (Drozdov, Eremets, Troyan, Ksenofontov, Shylin, Nature 525, 73 (2015): ~203 K at ~155 GPa) | 155.0 | ~2.0-2.2 harmonic DFPT for Im-3m H3S near 150-200 GPa (Duan et al., Sci. Rep. 4, 6968 (2014)); anharmonicity lowers it (~1.8; Errea et al., PRL 114, 157004 (2015)) -- values uncertain | ~1300 K (harmonic; uncertain) | Tc/pressure: high. lambda, omega_log: medium/low. |

Note for Al: Weak coupling: Tc is exponentially sensitive to lambda and mu*; a 10% change in lambda moves Tc by a factor ~2. mu* ~0.10-0.13 commonly fitted.

### 5. Consistency checks

No serious finding (all rechecks OK, F7 bookkeeping consistent, every float Tc inside its proved bracket).

| check | result | detail |
|---|---|---|
| recheck[T300] | yes | stored REJECTED_PROVEN [None, 17.149707508706445] -> rechecked REJECTED_PROVEN [None, 17.149707508706445] |
| F7 combination[T300] | yes | verdict [None, 17.149707508706445] vs [min T_low, max T_up] = [None, 17.149707508706445] |
| recheck[Tlit_mu0.10] | yes | stored CERTIFIED_CANDIDATE [1.2968826847644568, 10.302066039359831] -> rechecked CERTIFIED_CANDIDATE [1.2968826847644568, 10.302066039359831] |
| F7 combination[Tlit_mu0.10] | yes | verdict [1.2968826847644568, 10.302066039359831] vs [min T_low, max T_up] = [1.2968826847644568, 10.302066039359831] |
| recheck[Tlit_mu0.15] | yes | stored UNDECIDED [0.5778994748817508, 10.302066039359831] -> rechecked UNDECIDED [0.5778994748817508, 10.302066039359831] |
| F7 combination[Tlit_mu0.15] | yes | verdict [0.5778994748817508, 10.302066039359831] vs [min T_low, max T_up] = [0.5778994748817508, 10.302066039359831] |
| T_exp <= combined T_up [Tlit_mu0.10] | yes | T_exp=1.175 K, T_up=10.3 K |
| T_exp <= combined T_up [Tlit_mu0.15] | yes | T_exp=1.175 K, T_up=10.3 K |
| T_low(mu*=0.15) <= T_exp | yes | T_low=0.5779 K, T_exp=1.175 K |
| float Tc(mu*=0.10) >= T_low [Tlit_mu0.10, a2F.dos2] | yes | float 2.270734 K vs proved T_low 1.354107 K |
| float Tc(mu*=0) <= T_up [Tlit_mu0.10, a2F.dos2] | yes | float 7.761999 K vs proved T_up 10.30207 K |
| float Tc(mu*=0.10) >= T_low [Tlit_mu0.10, a2F.dos3] | yes | float 2.22971 K vs proved T_low 1.323645 K |
| float Tc(mu*=0) <= T_up [Tlit_mu0.10, a2F.dos3] | yes | float 7.69849 K vs proved T_up 10.23299 K |
| float Tc(mu*=0.10) >= T_low [Tlit_mu0.10, a2F.dos4] | yes | float 2.202935 K vs proved T_low 1.304887 K |
| float Tc(mu*=0) <= T_up [Tlit_mu0.10, a2F.dos4] | yes | float 7.650298 K vs proved T_up 10.17674 K |
| float Tc(mu*=0.10) >= T_low [Tlit_mu0.10, a2F.dos5] | yes | float 2.190695 K vs proved T_low 1.296883 K |
| float Tc(mu*=0) <= T_up [Tlit_mu0.10, a2F.dos5] | yes | float 7.622733 K vs proved T_up 10.14259 K |
| float Tc(mu*=0.10) >= T_low [Tlit_mu0.10, a2F.dos6] | yes | float 2.197279 K vs proved T_low 1.301652 K |
| float Tc(mu*=0) <= T_up [Tlit_mu0.10, a2F.dos6] | yes | float 7.629537 K vs proved T_up 10.14898 K |
| float Tc(mu*=0.10) >= T_low [Tlit_mu0.10, a2F.dos7] | yes | float 2.211848 K vs proved T_low 1.313027 K |
| float Tc(mu*=0) <= T_up [Tlit_mu0.10, a2F.dos7] | yes | float 7.65275 K vs proved T_up 10.17472 K |
| float Tc(mu*=0.10) >= T_low [Tlit_mu0.10, a2F.dos8] | yes | float 2.228467 K vs proved T_low 1.324629 K |
| float Tc(mu*=0) <= T_up [Tlit_mu0.10, a2F.dos8] | yes | float 7.67836 K vs proved T_up 10.20413 K |
| float Tc(mu*=0.10) >= T_low [Tlit_mu0.10, a2F.dos9] | yes | float 2.24301 K vs proved T_low 1.334976 K |
| float Tc(mu*=0) <= T_up [Tlit_mu0.10, a2F.dos9] | yes | float 7.701098 K vs proved T_up 10.2289 K |
| float Tc(mu*=0.10) >= T_low [Tlit_mu0.10, a2F.dos10] | yes | float 2.2549 K vs proved T_low 1.343617 K |
| float Tc(mu*=0) <= T_up [Tlit_mu0.10, a2F.dos10] | yes | float 7.720006 K vs proved T_up 10.2504 K |
| float Tc(mu*=0.15) >= T_low [Tlit_mu0.15, a2F.dos2] | yes | float 1.200854 K vs proved T_low 0.6132242 K |
| float Tc(mu*=0) <= T_up [Tlit_mu0.15, a2F.dos2] | yes | float 7.761999 K vs proved T_up 10.30207 K |
| float Tc(mu*=0.15) >= T_low [Tlit_mu0.15, a2F.dos3] | yes | float 1.170288 K vs proved T_low 0.5935595 K |
| float Tc(mu*=0) <= T_up [Tlit_mu0.15, a2F.dos3] | yes | float 7.69849 K vs proved T_up 10.23299 K |
| float Tc(mu*=0.15) >= T_low [Tlit_mu0.15, a2F.dos4] | yes | float 1.151923 K vs proved T_low 0.5823845 K |
| float Tc(mu*=0) <= T_up [Tlit_mu0.15, a2F.dos4] | yes | float 7.650298 K vs proved T_up 10.17674 K |
| float Tc(mu*=0.15) >= T_low [Tlit_mu0.15, a2F.dos5] | yes | float 1.144534 K vs proved T_low 0.5778995 K |
| float Tc(mu*=0) <= T_up [Tlit_mu0.15, a2F.dos5] | yes | float 7.622733 K vs proved T_up 10.14259 K |
| float Tc(mu*=0.15) >= T_low [Tlit_mu0.15, a2F.dos6] | yes | float 1.149708 K vs proved T_low 0.5813594 K |
| float Tc(mu*=0) <= T_up [Tlit_mu0.15, a2F.dos6] | yes | float 7.629537 K vs proved T_up 10.14898 K |
| float Tc(mu*=0.15) >= T_low [Tlit_mu0.15, a2F.dos7] | yes | float 1.160831 K vs proved T_low 0.5883202 K |
| float Tc(mu*=0) <= T_up [Tlit_mu0.15, a2F.dos7] | yes | float 7.65275 K vs proved T_up 10.17472 K |
| float Tc(mu*=0.15) >= T_low [Tlit_mu0.15, a2F.dos8] | yes | float 1.172085 K vs proved T_low 0.5957939 K |
| float Tc(mu*=0) <= T_up [Tlit_mu0.15, a2F.dos8] | yes | float 7.67836 K vs proved T_up 10.20413 K |
| float Tc(mu*=0.15) >= T_low [Tlit_mu0.15, a2F.dos9] | yes | float 1.182786 K vs proved T_low 0.6024708 K |
| float Tc(mu*=0) <= T_up [Tlit_mu0.15, a2F.dos9] | yes | float 7.701098 K vs proved T_up 10.2289 K |
| float Tc(mu*=0.15) >= T_low [Tlit_mu0.15, a2F.dos10] | yes | float 1.191392 K vs proved T_low 0.6079147 K |
| float Tc(mu*=0) <= T_up [Tlit_mu0.15, a2F.dos10] | yes | float 7.720006 K vs proved T_up 10.2504 K |
| AD vs Eliashberg float [a2F.dos2, mu*=0.00] | yes | AD=7.56 K, Eliashberg=7.762 K, rel diff +2.7% |
| AD vs Eliashberg float [a2F.dos2, mu*=0.10] | info | AD=1.603 K, Eliashberg=2.271 K, rel diff +41.7% (info: different mu* conventions) |
| AD vs Eliashberg float [a2F.dos2, mu*=0.15] | info | AD=0.3898 K, Eliashberg=1.201 K, rel diff +208.1% (info: different mu* conventions) |
| AD vs Eliashberg float [a2F.dos3, mu*=0.00] | yes | AD=7.491 K, Eliashberg=7.698 K, rel diff +2.8% |
| AD vs Eliashberg float [a2F.dos3, mu*=0.10] | info | AD=1.562 K, Eliashberg=2.23 K, rel diff +42.8% (info: different mu* conventions) |
| AD vs Eliashberg float [a2F.dos3, mu*=0.15] | info | AD=0.3719 K, Eliashberg=1.17 K, rel diff +214.7% (info: different mu* conventions) |
| AD vs Eliashberg float [a2F.dos4, mu*=0.00] | yes | AD=7.443 K, Eliashberg=7.65 K, rel diff +2.8% |
| AD vs Eliashberg float [a2F.dos4, mu*=0.10] | info | AD=1.539 K, Eliashberg=2.203 K, rel diff +43.1% (info: different mu* conventions) |
| AD vs Eliashberg float [a2F.dos4, mu*=0.15] | info | AD=0.3629 K, Eliashberg=1.152 K, rel diff +217.5% (info: different mu* conventions) |
| AD vs Eliashberg float [a2F.dos5, mu*=0.00] | yes | AD=7.417 K, Eliashberg=7.623 K, rel diff +2.8% |
| AD vs Eliashberg float [a2F.dos5, mu*=0.10] | info | AD=1.531 K, Eliashberg=2.191 K, rel diff +43.1% (info: different mu* conventions) |
| AD vs Eliashberg float [a2F.dos5, mu*=0.15] | info | AD=0.3598 K, Eliashberg=1.145 K, rel diff +218.1% (info: different mu* conventions) |
| AD vs Eliashberg float [a2F.dos6, mu*=0.00] | yes | AD=7.425 K, Eliashberg=7.63 K, rel diff +2.8% |
| AD vs Eliashberg float [a2F.dos6, mu*=0.10] | info | AD=1.538 K, Eliashberg=2.197 K, rel diff +42.9% (info: different mu* conventions) |
| AD vs Eliashberg float [a2F.dos6, mu*=0.15] | info | AD=0.3632 K, Eliashberg=1.15 K, rel diff +216.5% (info: different mu* conventions) |
| AD vs Eliashberg float [a2F.dos7, mu*=0.00] | yes | AD=7.449 K, Eliashberg=7.653 K, rel diff +2.7% |
| AD vs Eliashberg float [a2F.dos7, mu*=0.10] | info | AD=1.552 K, Eliashberg=2.212 K, rel diff +42.5% (info: different mu* conventions) |
| AD vs Eliashberg float [a2F.dos7, mu*=0.15] | info | AD=0.3695 K, Eliashberg=1.161 K, rel diff +214.1% (info: different mu* conventions) |
| AD vs Eliashberg float [a2F.dos8, mu*=0.00] | yes | AD=7.475 K, Eliashberg=7.678 K, rel diff +2.7% |
| AD vs Eliashberg float [a2F.dos8, mu*=0.10] | info | AD=1.568 K, Eliashberg=2.228 K, rel diff +42.2% (info: different mu* conventions) |
| AD vs Eliashberg float [a2F.dos8, mu*=0.15] | info | AD=0.3761 K, Eliashberg=1.172 K, rel diff +211.6% (info: different mu* conventions) |
| AD vs Eliashberg float [a2F.dos9, mu*=0.00] | yes | AD=7.498 K, Eliashberg=7.701 K, rel diff +2.7% |
| AD vs Eliashberg float [a2F.dos9, mu*=0.10] | info | AD=1.581 K, Eliashberg=2.243 K, rel diff +41.9% (info: different mu* conventions) |
| AD vs Eliashberg float [a2F.dos9, mu*=0.15] | info | AD=0.382 K, Eliashberg=1.183 K, rel diff +209.6% (info: different mu* conventions) |
| AD vs Eliashberg float [a2F.dos10, mu*=0.00] | yes | AD=7.518 K, Eliashberg=7.72 K, rel diff +2.7% |
| AD vs Eliashberg float [a2F.dos10, mu*=0.10] | info | AD=1.593 K, Eliashberg=2.255 K, rel diff +41.6% (info: different mu* conventions) |
| AD vs Eliashberg float [a2F.dos10, mu*=0.15] | info | AD=0.387 K, Eliashberg=1.191 K, rel diff +207.9% (info: different mu* conventions) |
| independent float vs rtsc floatcore [a2F.dos2, mu*=0.00] | yes | 7.761999 K vs 7.761941 K |
| independent float vs rtsc floatcore [a2F.dos2, mu*=0.10] | yes | 2.270734 K vs 2.270734 K |
| independent float vs rtsc floatcore [a2F.dos2, mu*=0.15] | yes | 1.200854 K vs 1.200856 K |
| independent float vs rtsc floatcore [a2F.dos3, mu*=0.00] | yes | 7.69849 K vs 7.698435 K |
| independent float vs rtsc floatcore [a2F.dos3, mu*=0.10] | yes | 2.22971 K vs 2.229714 K |
| independent float vs rtsc floatcore [a2F.dos3, mu*=0.15] | yes | 1.170288 K vs 1.17029 K |
| independent float vs rtsc floatcore [a2F.dos4, mu*=0.00] | yes | 7.650298 K vs 7.650244 K |
| independent float vs rtsc floatcore [a2F.dos4, mu*=0.10] | yes | 2.202935 K vs 2.202939 K |
| independent float vs rtsc floatcore [a2F.dos4, mu*=0.15] | yes | 1.151923 K vs 1.151923 K |
| independent float vs rtsc floatcore [a2F.dos5, mu*=0.00] | yes | 7.622733 K vs 7.622677 K |
| independent float vs rtsc floatcore [a2F.dos5, mu*=0.10] | yes | 2.190695 K vs 2.190699 K |
| independent float vs rtsc floatcore [a2F.dos5, mu*=0.15] | yes | 1.144534 K vs 1.144536 K |
| independent float vs rtsc floatcore [a2F.dos6, mu*=0.00] | yes | 7.629537 K vs 7.629482 K |
| independent float vs rtsc floatcore [a2F.dos6, mu*=0.10] | yes | 2.197279 K vs 2.197284 K |
| independent float vs rtsc floatcore [a2F.dos6, mu*=0.15] | yes | 1.149708 K vs 1.14971 K |
| independent float vs rtsc floatcore [a2F.dos7, mu*=0.00] | yes | 7.65275 K vs 7.652697 K |
| independent float vs rtsc floatcore [a2F.dos7, mu*=0.10] | yes | 2.211848 K vs 2.211852 K |
| independent float vs rtsc floatcore [a2F.dos7, mu*=0.15] | yes | 1.160831 K vs 1.160831 K |
| independent float vs rtsc floatcore [a2F.dos8, mu*=0.00] | yes | 7.67836 K vs 7.678303 K |
| independent float vs rtsc floatcore [a2F.dos8, mu*=0.10] | yes | 2.228467 K vs 2.228471 K |
| independent float vs rtsc floatcore [a2F.dos8, mu*=0.15] | yes | 1.172085 K vs 1.172087 K |
| independent float vs rtsc floatcore [a2F.dos9, mu*=0.00] | yes | 7.701098 K vs 7.701043 K |
| independent float vs rtsc floatcore [a2F.dos9, mu*=0.10] | yes | 2.24301 K vs 2.243014 K |
| independent float vs rtsc floatcore [a2F.dos9, mu*=0.15] | yes | 1.182786 K vs 1.182786 K |
| independent float vs rtsc floatcore [a2F.dos10, mu*=0.00] | yes | 7.720006 K vs 7.719952 K |
| independent float vs rtsc floatcore [a2F.dos10, mu*=0.10] | yes | 2.2549 K vs 2.254904 K |
| independent float vs rtsc floatcore [a2F.dos10, mu*=0.15] | yes | 1.191392 K vs 1.191394 K |
| lambda vs QE lambda file [a2F.dos2] | yes | parsed 0.41519 vs file 0.4153 (-0.03%; histogram vs direct q-sum, rtsc tolerance 5%) |
| lambda vs QE lambda file [a2F.dos3] | yes | parsed 0.41274 vs file 0.4128 (-0.01%; histogram vs direct q-sum, rtsc tolerance 5%) |
| lambda vs QE lambda file [a2F.dos4] | yes | parsed 0.41159 vs file 0.4117 (-0.03%; histogram vs direct q-sum, rtsc tolerance 5%) |
| lambda vs QE lambda file [a2F.dos5] | yes | parsed 0.41128 vs file 0.4115 (-0.05%; histogram vs direct q-sum, rtsc tolerance 5%) |
| lambda vs QE lambda file [a2F.dos6] | yes | parsed 0.41182 vs file 0.412 (-0.04%; histogram vs direct q-sum, rtsc tolerance 5%) |
| lambda vs QE lambda file [a2F.dos7] | yes | parsed 0.4127 vs file 0.4129 (-0.05%; histogram vs direct q-sum, rtsc tolerance 5%) |
| lambda vs QE lambda file [a2F.dos8] | yes | parsed 0.41361 vs file 0.4139 (-0.07%; histogram vs direct q-sum, rtsc tolerance 5%) |
| lambda vs QE lambda file [a2F.dos9] | yes | parsed 0.41441 vs file 0.4147 (-0.07%; histogram vs direct q-sum, rtsc tolerance 5%) |
| lambda vs QE lambda file [a2F.dos10] | yes | parsed 0.41509 vs file 0.4154 (-0.08%; histogram vs direct q-sum, rtsc tolerance 5%) |
| smearing spread | info | lambda: 0.4113..0.4152 (0.9%); omega_log_K: 304.3..308 (1.2%); A_K: 67.07..67.53 (0.7%); omega_max_K: 471.5..471.5 (0.0%) |
| per-smearing brackets [Tlit_mu0.10] | info | [1.354, 10.3], [1.324, 10.23], [1.305, 10.18], [1.297, 10.14], [1.302, 10.15], [1.313, 10.17], [1.325, 10.2], [1.335, 10.23], [1.344, 10.25] |
| per-smearing brackets [Tlit_mu0.15] | info | [0.6132, 10.3], [0.5936, 10.23], [0.5824, 10.18], [0.5779, 10.14], [0.5814, 10.15], [0.5883, 10.17], [0.5958, 10.2], [0.6025, 10.23], [0.6079, 10.25] |
| M-ISO-SINGLEBAND present per verdict | yes | {"T300": true, "Tlit_mu0.10": true, "Tlit_mu0.15": true} |
| MULTIBAND_RISK on every REJECTED_PROVEN verdict | yes | statuses {'T300': 'REJECTED_PROVEN', 'Tlit_mu0.10': 'CERTIFIED_CANDIDATE', 'Tlit_mu0.15': 'UNDECIDED'}; MULTIBAND_RISK {'T300': True, 'Tlit_mu0.10': False, 'Tlit_mu0.15': False} (classify sets it only on single-band rejections) |
