QE=~/micromamba/envs/qe/bin; export OMP_NUM_THREADS=1
nm $QE/ph.x 2>&1 | head -2; which gdb
for OCC in tetrahedra_lin tetrahedra_opt; do for NP in 1 4; do
D=~/rtsc_work/tetra_Al/linbug/${OCC#tetrahedra_}_np$NP; rm -rf $D; mkdir -p $D; cd $D
sed -e "s/tetrahedra_lin/$OCC/" -e "s/ecutwfc = 60.0/ecutwfc = 30.0/" -e "s/ecutrho = 240.0/ecutrho = 120.0/" -e "s/ 16 16 16 0 0 0/ 8 8 8 0 0 0/" ~/rtsc_work/tetra_Al/q6s_lin/scf.in > scf.in
sed -e "s/nq\([123]\) = 6/nq\1 = 2/" ~/rtsc_work/tetra_Al/q6s_lin/ph.in > ph.in
if [ $NP = 1 ]; then R="$QE/mpirun -np 1"; P=""; else R="$QE/mpirun -np 4 --bind-to none -x OMP_NUM_THREADS"; P="-nk 4"; fi
$R $QE/pw.x $P -in scf.in > scf.out 2>&1
timeout 900 $R $QE/ph.x $P -in ph.in > ph.out 2> ph.err; echo "$OCC np$NP ph rc=$? $(grep -c 'Calculation of q' ph.out) $(grep -c 'JOB DONE' ph.out)"
done; done
