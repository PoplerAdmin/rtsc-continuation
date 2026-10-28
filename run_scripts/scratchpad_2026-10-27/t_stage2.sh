# usage: NQ=6 NKS="16 24 32" bash t_stage2.sh
QE=~/micromamba/envs/qe/bin
export OMP_NUM_THREADS=1
NQ=${NQ:-6}
OCC=${OCC:-tetrahedra_lin}
W=~/rtsc_work/tetra_Al/q${NQ}s_${OCC#tetrahedra_}
cd $W
grep -q STAGE1_DONE timing.txt || { echo "stage1 not done"; exit 1; }
for NK in ${NKS:-16 24 32}; do
  D=$W/elph_k$NK
  if grep -q ELPH_DONE $D/timing.txt 2>/dev/null; then continue; fi
  rm -rf $D; mkdir -p $D; cp -r $W/tmp $D/tmp; cp $W/al.dyn* $D/
  cd $D
  cat > elph.in <<EOI
Al electron-phonon lambda_tetra
 &inputph
    prefix = 'al'
    outdir = './tmp/'
    tr2_ph = 1.0d-16
    fildyn = 'al.dyn'
    fildvscf = 'dv'
    ldisp = .true.
    lshift_q = .true.
    nq1 = $NQ
    nq2 = $NQ
    nq3 = $NQ
    amass(1) = 26.9815386
    electron_phonon = 'lambda_tetra'
    nk1 = $NK
    nk2 = $NK
    nk3 = $NK
 /
 &inputa2f
    nfreq = 500
 /
EOI
  t0=$(date +%s)
  $QE/mpirun -np 4 --bind-to none -x OMP_NUM_THREADS $QE/ph.x -nk 4 -in elph.in > elph.out 2> elph.err
  rc=$?
  t1=$(date +%s); echo "elph_wall_s $((t1-t0)) rc $rc" > timing.txt
  $QE/mpirun -np 1 $QE/alpha2f.x -in elph.in > a2F.out 2> a2F.err
  echo "alpha2f_rc $?" >> timing.txt
  rm -rf $D/tmp/_ph0/al.save $D/tmp/_ph0/*.wfc* 2>/dev/null
  echo ELPH_DONE >> timing.txt
  cd $W
done
echo STAGE2_DONE NKS="${NKS}" >> $W/timing.txt
