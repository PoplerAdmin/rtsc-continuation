export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
Q=$HOME/micromamba/envs/qe/bin
cd ~/rtsc_work/epw_Al
( time $Q/mpirun -np 4 --bind-to none -x OMP_NUM_THREADS -x OPENBLAS_NUM_THREADS $Q/pw.x -nk 4 -in nscf.in > nscf.out ) 2>&1 | tail -3
grep -E "Fermi energy|JOB DONE" nscf.out
