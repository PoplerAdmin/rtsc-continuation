set -e
lscpu | grep -E "^(Model name|Thread|Core|Socket|CPU\(s\))"
B=~/rtsc_work/bench/Al; mkdir -p $B/structures
. ~/rtsc_work/.venv/bin/activate
python - <<'PY'
from ase.build import bulk
from ase.io import write
import os
a = bulk("Al", "fcc", a=4.0495)   # experimental room-T lattice constant; vc-relaxed by the workflow
write(os.path.expanduser("~/rtsc_work/bench/Al/structures/Al_fcc.vasp"), a, format="vasp", direct=True)
PY
cat > $B/structures/manifest.toml <<'T'
# Benchmark input. Only facts about the input: target pressure.
["Al_fcc"]
material_id = "Al-fcc@0GPa"
pressure_GPa = 0.0
T
cat > $B/config.toml <<'T'
# BENCHMARK (step 2 of HANDOFF.md): fcc Al, real DFPT with QE 7.6 (conda-forge), SSSP 1.3.0 PBE precision.
[run]
name = "bench-Al"
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
enabled = false        # single element: no hull to evaluate
[filters.metallic]
enabled = true
max_band_gap_eV = 0.0
[filters.magnetic]
enabled = true
[filters.dynamic]
enabled = true
imaginary_tol = 0.25
imaginary_tol_unit = "THz"

[budget]
max_materials = 1

[spectra]
use_attached = false
rel_envelope = 0.10

[dfpt]
enabled = true
[dfpt.launcher]
kind = "mpirun"
nprocs = 16
npool = 16
mpi_args = ["--bind-to", "none", "-x", "OMP_NUM_THREADS"]
bin_dir = "/home/user/micromamba/envs/qe/bin"
env = { OMP_NUM_THREADS = "1" }
[dfpt.pseudos]
dir = "/home/user/pseudo/SSSP_1.3.0_PBE_precision"
sssp_json = "/home/user/pseudo/SSSP_1.3.0_PBE_precision.json"
[dfpt.settings]
do_relax = true
ecutwfc_Ry = 40.0      # SSSP precision recommends 30/240 for Al; 40/320 for margin
ecutrho_Ry = 320.0
smearing = "mv"
degauss_Ry = 0.02
k_scf = [16, 16, 16]
k_fit = [32, 32, 32]
q = [8, 8, 8]
k_matdyn = [32, 32, 32]
el_ph_sigma_Ry = 0.005
el_ph_nsigma = 10
tr2_ph = 1e-16

[certify]
T_target_K = 300.0
mu_hi = 0.15
omega_c_factor = 10.0
cutoff_rule = "block"
mode = "bracket"

[validity]
enabled = true
T
cd $B && rtsc run -c config.toml --plan-only 2>&1 | tail -15
