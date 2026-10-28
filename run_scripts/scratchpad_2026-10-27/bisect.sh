# Bisect the ph.x segfault on a cheap setup: q 2x2x2 only, same pseudos/flags.
Q=$HOME/micromamba/envs/qe/bin
export OMP_NUM_THREADS=1
run() { # name pseudo_dir pseudo ecut ecutrho kscf kfit np
  local n=$1 pd=$2 pp=$3 ec=$4 er=$5 k=$6 kf=$7 np=$8
  local d=~/rtsc_work/bisect/$n; rm -rf $d; mkdir -p $d; cd $d
  for kind in fit scf; do
    if [ $kind = fit ]; then kk=$kf; extra="la2F = .true."; else kk=$k; extra=""; fi
    cat > $kind.in <<P
&control
 calculation='scf', prefix='al', outdir='./tmp/', pseudo_dir='$pd'
/
&system
 ibrav=2, celldm(1)=7.65, nat=1, ntyp=1, ecutwfc=$ec, ecutrho=$er,
 occupations='smearing', smearing='mv', degauss=0.02 $( [ -n "$extra" ] && echo ", $extra")
/
&electrons
 conv_thr=1e-10
/
ATOMIC_SPECIES
 Al 26.98 $pp
ATOMIC_POSITIONS crystal
 Al 0 0 0
K_POINTS automatic
 $kk $kk $kk 0 0 0
P
  done
  cat > ph.in <<P
e-ph
&inputph
 prefix='al', outdir='./tmp/', tr2_ph=1e-16, fildyn='dyn', fildvscf='dvscf',
 electron_phonon='interpolated', el_ph_sigma=0.005, el_ph_nsigma=10,
 trans=.true., ldisp=.true., nq1=2, nq2=2, nq3=2, amass(1)=26.98
/
P
  if [ $np -gt 1 ]; then L="$Q/mpirun -np $np --bind-to none $Q/"; NK="-nk $np"; else L="$Q/"; NK=""; fi
  ${L}pw.x $NK -in fit.in > fit.out 2>&1 && ${L}pw.x $NK -in scf.in > scf.out 2>&1 && ${L}ph.x $NK -in ph.in > ph.out 2> ph.err
  echo "$n: ph rc=$? $(grep -c 'lambda(' ph.out) lambda-lines; $(grep -m1 -o 'Segmentation fault\|JOB DONE' ph.out ph.err | head -1)"
}
P=$HOME/pseudo/SSSP_1.3.0_PBE_precision; N=$HOME/pseudo/qe_std


run sg15_mpi8_k16 $HOME/pseudo/sg15 Al_ONCV_PBE-1.2.upf 40 160 8 16 8
