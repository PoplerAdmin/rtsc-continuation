export MAMBA_ROOT_PREFIX=~/micromamba; M=~/bin/micromamba
ls ~/micromamba/envs/*/lib/libelpa.so.* 2>/dev/null
for v in 2025.06.001 2026.02.002 2024.05.001 2024.03.001; do echo "== $v"; $M create -n dry2 --dry-run -y -c conda-forge --override-channels "elpa=$v=mpi_openmpi*" "openmpi=5.0.11" 2>&1 | grep -E 'elpa|error' | head -3; done
