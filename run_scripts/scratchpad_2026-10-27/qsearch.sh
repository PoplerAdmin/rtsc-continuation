export MAMBA_ROOT_PREFIX=~/micromamba; M=~/bin/micromamba
$M search -c conda-forge qe 2>&1 | head -60
$M list -n qe 2>&1 | grep -E '^\s*(qe|openmpi|libgfortran|gfortran|fftw|openblas|libblas|liblapack|scalapack|elpa|hdf5|libxc|mpich|libopenblas)\b'
$M env list
which gdb; ls ~/micromamba/pkgs | head -50; du -sh ~/micromamba/pkgs
