export MAMBA_ROOT_PREFIX=~/micromamba; M=~/bin/micromamba
$M install -n qe76b0 -y -c conda-forge --override-channels "elpa=2026.02.002=mpi_openmpi*" 2>&1 | tail -4
ldd ~/micromamba/envs/qe76b0/bin/ph.x | grep -c 'not found'
