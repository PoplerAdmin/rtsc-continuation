export MAMBA_ROOT_PREFIX=~/micromamba; M=~/bin/micromamba
L=/mnt/g/ClaudeWorkspace/RTSC_Project/run_scripts/scratchpad_2026-10-27/mkenvs.log
echo "pkgs before: $(du -sm ~/micromamba/pkgs | cut -f1) MB" > $L
mk() { nice -n 10 $M create -y -n $1 -c conda-forge --override-channels "$2" >> $L.full 2>&1; echo "$1 $2 rc=$? pkgs=$(du -sm ~/micromamba/pkgs | cut -f1)MB" >> $L; }
mk qe75b "qe=7.5=h19104ac_2"
mk qe76mpich "qe=7.6=mpi_mpich_he3edec4_4"
mk qe74 "qe=7.4=hbc83ec2_0"
mk qe76b0 "qe=7.6=h2ec8f75_0"
mk qe75a "qe=7.5=hbc83ec2_0"
mk qe72 "qe=7.2=hac89879_4"
mk dbg "gdb"
echo ALLDONE >> $L
