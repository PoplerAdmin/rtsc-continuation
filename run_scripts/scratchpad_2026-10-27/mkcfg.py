import os
B=os.path.expanduser('~/rtsc_work/bench')
COMMON_TAIL='''
[certify]
T_target_K = 300.0
mu_hi = 0.15
omega_c_factor = 10.0
cutoff_rule = "block"
mode = "bracket"

[validity]
enabled = true
'''
LAUNCH='''[dfpt]
enabled = true
[dfpt.launcher]
kind = "mpirun"
nprocs = 8
npool = 8
mpi_args = ["--bind-to", "none", "-x", "OMP_NUM_THREADS"]
bin_dir = "/home/user/micromamba/envs/qe/bin"
mpirun = "/home/user/micromamba/envs/qe/bin/mpirun"
env = { OMP_NUM_THREADS = "1" }
'''
def cfg(name, header, hull_comment, magnetic, files, settings, extra=''):
    return f'''{header}
# Pseudopotential choice: QE 7.6 ph.x segfaults (adddvscf <- elphel) for ANY ultrasoft or PAW pseudo with
# electron_phonon set (tested 2026-10-27 on Al/Mg/Pb/Nb PAW and Mg/B/S GBRV-US), so only norm-conserving SG15 is used.
# Cutoffs from pw.x convergence tests (2026-10-27; CONVERGENCE TESTS, not physical results), see
# /mnt/g/ClaudeWorkspace/RTSC_Project/reports/step2_pseudopotentials_and_configs.md
[run]
name = "bench-{name}"
workdir = "run"
workers = 1
resume = true

[[sources]]
type = "local"
path = "structures"
recursive = false
backend = "auto"

[dedupe]
method = "auto"

[filters.hull]
enabled = false        # {hull_comment}
[filters.metallic]
enabled = true
max_band_gap_eV = 0.0
[filters.magnetic]
enabled = {magnetic}
[filters.dynamic]
enabled = true
imaginary_tol = 0.25
imaginary_tol_unit = "THz"

[budget]
max_materials = 1

[spectra]
use_attached = false
rel_envelope = 0.10

{LAUNCH}[dfpt.pseudos]
dir = "/home/user/pseudo/sg15"
files = {files}
[dfpt.settings]
{settings}{extra}{COMMON_TAIL}'''

configs = {
'Pb': cfg('Pb',
 '# BENCHMARK (step 2 of HANDOFF.md): fcc Pb (expt Tc 7.2 K), DFPT with QE 7.6, SG15 ONCV PBE 1.2 (norm-conserving, scalar-relativistic).\n'
 '# Scalar-relativistic DFPT underestimates lambda for Pb (SOC raises it from ~1.2 to ~1.5; Heid et al. PRB 81, 174527 (2010)).\n'
 '# SOC variant (NOT default; ~3-4x cost, untested in rtsc): files = { Pb = "Pb_ONCV_PBE_FR-1.0.upf" } and\n'
 '#   extra_pw_system = { noncolin = true, lspinorb = true }   (FR pseudo converged at the same 50 Ry in tests)',
 'single element: no hull to evaluate', 'true',
 '{ Pb = "Pb_ONCV_PBE-1.2.upf" }',
 '''do_relax = true
ecutwfc_Ry = 50.0      # SG15 Pb (14 e-, 5d semicore): 40 Ry already within 0.01 mRy / 0.05 kbar of 120 Ry
ecutrho_Ry = 200.0     # NC: 4*ecutwfc
smearing = "mv"
degauss_Ry = 0.02
k_scf = [16, 16, 16]
k_fit = [32, 32, 32]
q = [8, 8, 8]          # 29 irreducible q
k_matdyn = [32, 32, 32]
el_ph_sigma_Ry = 0.005
el_ph_nsigma = 10
tr2_ph = 1e-16
'''),
'Nb': cfg('Nb',
 '# BENCHMARK (step 2 of HANDOFF.md): bcc Nb (expt Tc 9.3 K), DFPT with QE 7.6, SG15 ONCV PBE 1.2 (norm-conserving, 4s4p semicore).',
 'single element: no hull to evaluate', 'true      # Nb is non-magnetic; the filter is a no-op safety check',
 '{ Nb = "Nb_ONCV_PBE-1.2.upf" }',
 '''do_relax = true
ecutwfc_Ry = 120.0     # SG15 Nb converges slowly: 120 Ry within 1.0 mRy of 200 Ry; pressure scatter +-2 kbar at k 8^3
ecutrho_Ry = 480.0     # NC: 4*ecutwfc
smearing = "mv"
degauss_Ry = 0.02
k_scf = [16, 16, 16]
k_fit = [32, 32, 32]
q = [8, 8, 8]          # 29 irreducible q
k_matdyn = [32, 32, 32]
el_ph_sigma_Ry = 0.005
el_ph_nsigma = 10
tr2_ph = 1e-16
'''),
'MgB2': cfg('MgB2',
 '# BENCHMARK (step 2 of HANDOFF.md): MgB2 P6/mmm (expt Tc 39 K; two-gap). Isotropic DFPT gives lambda ~0.7 and an\n'
 '# isotropic Tc well below 39 K: the expected single-band underestimate; the verdict must carry MULTIBAND_RISK.\n'
 '# SG15 ONCV PBE 1.2 for Mg (10 e-, 2s2p semicore) and B (3 e-).',
 'binary, but only one structure and no reference energies: hull not evaluated', 'true',
 '{ Mg = "Mg_ONCV_PBE-1.2.upf", B = "B_ONCV_PBE-1.2.upf" }',
 '''do_relax = true
ecutwfc_Ry = 90.0      # 90 Ry: 0.05 mRy/atom, 0.8 kbar vs 120 Ry (80 Ry: 0.15 mRy/atom, 1.7 kbar)
ecutrho_Ry = 360.0     # NC: 4*ecutwfc
smearing = "mv"
degauss_Ry = 0.02
k_scf = [18, 18, 12]   # ~isotropic k spacing (c/a = 1.14); commensurate with q
k_fit = [36, 36, 24]   # dense double-delta grid: the sigma-band FS sheets are small cylinders
q = [6, 6, 4]          # 21 irreducible q
k_matdyn = [36, 36, 24]
el_ph_sigma_Ry = 0.005
el_ph_nsigma = 10
tr2_ph = 1e-16
'''),
'H3S': cfg('H3S',
 '# BENCHMARK (step 2 of HANDOFF.md): Im-3m H3S at 200 GPa (expt Tc ~200 K at ~155 GPa). Primitive 4-atom bcc cell.\n'
 '# Harmonic PBE; anharmonicity (Errea et al. 2015/2016) lowers lambda substantially: validity flags expected.\n'
 '# SG15 ONCV PBE 1.2 for S (6 e-) and H. Core spheres overlap at this pressure (S-H 1.49 A = 2.82 bohr <\n'
 '# rc(S)+rc(H) = 3.27 bohr); cross-checked against GBRV/pslibrary US: pressure within 1.1%, see report.',
 'binary, but only one structure and no reference enthalpies: hull not evaluated', 'true',
 '{ S = "S_ONCV_PBE-1.2.upf", H = "H_ONCV_PBE-1.2.upf" }',
 '''do_relax = true
ecutwfc_Ry = 90.0      # 90 Ry: 0.03 mRy/atom, 0.3 kbar vs 140 Ry (at a = 3.09 A, P ~ 138 GPa)
ecutrho_Ry = 360.0     # NC: 4*ecutwfc
smearing = "mv"
degauss_Ry = 0.02
k_scf = [18, 18, 18]   # bcc primitive; commensurate with q
k_fit = [36, 36, 36]
q = [6, 6, 6]          # 16 irreducible q
k_matdyn = [36, 36, 36]
el_ph_sigma_Ry = 0.005
el_ph_nsigma = 10
tr2_ph = 1e-16
'''),
}
for n,t in configs.items():
    with open(f'{B}/{n}/config.toml','w') as f: f.write(t)
print('written')
