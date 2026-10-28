set -e
mkdir -p ~/bin ~/rtsc_work/logs
cd ~
if [ ! -x ~/bin/micromamba ]; then
  curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest -o /tmp/mm.tar.bz2 && python3 -c "import tarfile; tarfile.open(\"/tmp/mm.tar.bz2\",\"r:bz2\").extract(\"bin/micromamba\", path=__import__(\"os\").path.expanduser(\"~\"))"
fi
~/bin/micromamba --version
export MAMBA_ROOT_PREFIX=~/micromamba
~/bin/micromamba create -y -n qe -c conda-forge "qe>=7.4=*mpi_openmpi*" openmpi > ~/rtsc_work/logs/qe_env.log 2>&1
~/bin/micromamba list -n qe | grep -E "^ *(qe|openmpi|openblas|fftw) "
ls ~/micromamba/envs/qe/bin | grep -E "^(pw|ph|q2r|matdyn|epw)\.x$"
