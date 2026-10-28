export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
B=~/micromamba/envs/qe/bin
cd ~/rtsc_work/epw_Al/wvar_A/wannier
timeout 900 $B/mpirun -np 4 --bind-to none -x OMP_NUM_THREADS -x OPENBLAS_NUM_THREADS $B/epw.x -nk 4 -in epw.in > epw.rerun.out 2> epw.rerun.err
echo rc=$?
tail -n 5 epw.rerun.out; echo ---; tail -n 30 epw.rerun.err; echo ---; tail -n 6 pwscf.wout; ls ../out
