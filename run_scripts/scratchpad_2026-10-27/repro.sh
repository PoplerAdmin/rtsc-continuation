#!/bin/bash
# usage: repro.sh ENV CASE [VARIANT]   (serial unless VARIANT=nk2)
# CASE: al_paw_interp al_paw_simple al_paw_noeph al_sg15_interp mg_us_interp
# VARIANT: base nofildvscf nosym noinv nk2
ENV=$1; CASE=$2; VAR=${3:-base}
B=$HOME/micromamba/envs/$ENV/bin
export OMP_NUM_THREADS=1 GFORTRAN_ERROR_BACKTRACE=1
P=$HOME/pseudo/SSSP_1.3.0_PBE_precision; SG=$HOME/pseudo/sg15
D=$HOME/rtsc_work/qebug/$ENV/${CASE}_$VAR; rm -rf $D; mkdir -p $D; cd $D
EPH=interpolated
case $CASE in
 al_paw_*|al_sg15_*) SYS="ibrav=2, celldm(1)=7.65, nat=1, ntyp=1"; POS=$'ATOMIC_POSITIONS crystal\n Al 0.0 0.0 0.0'; K="8 8 8"; KF="16 16 16";;
 mg_us_*) SYS="ibrav=4, celldm(1)=6.05, celldm(3)=1.624, nat=2, ntyp=1"; POS=$'ATOMIC_POSITIONS crystal\n Mg 0.333333333 0.666666667 0.25\n Mg 0.666666667 0.333333333 0.75'; K="6 6 4"; KF="12 12 8";;
esac
case $CASE in
 al_paw_*) PD=$P; SPC=" Al 26.98 Al.pbe-n-kjpaw_psl.1.0.0.UPF"; EC="ecutwfc=30, ecutrho=240";;
 al_sg15_*) PD=$SG; SPC=" Al 26.98 Al_ONCV_PBE-1.2.upf"; EC="ecutwfc=30, ecutrho=120";;
 mg_us_*) PD=$P; SPC=" Mg 24.305 mg_pbe_v1.4.uspp.F.UPF"; EC="ecutwfc=30, ecutrho=240";;
esac
case $CASE in *_simple) EPH=simple;; *_noeph) EPH=none;; esac
SYMX=""; [ $VAR = nosym ] && SYMX=", nosym=.true."; [ $VAR = noinv ] && SYMX=", noinv=.true."
for kind in fit scf; do
  if [ $kind = fit ]; then kk=$KF; ex=", la2F=.true."; else kk=$K; ex=""; fi
  [ $EPH = none ] && [ $kind = fit ] && continue
  cat > $kind.in <<P
&control
 calculation='scf', prefix='t', outdir='./tmp/', pseudo_dir='$PD'
/
&system
 $SYS, $EC,
 occupations='smearing', smearing='mv', degauss=0.02 $ex $SYMX
/
&electrons
 conv_thr=1e-10
/
ATOMIC_SPECIES
$SPC
$POS
K_POINTS automatic
 $kk 0 0 0
P
done
if [ $EPH = none ]; then EP=""; else EP="electron_phonon='$EPH', el_ph_sigma=0.005, el_ph_nsigma=10,"; fi
FD="fildvscf='dvscf',"; [ $VAR = nofildvscf ] && FD=""
cat > ph.in <<P
ph
&inputph
 prefix='t', outdir='./tmp/', tr2_ph=1e-14, fildyn='dyn', $FD
 $EP
 trans=.true., ldisp=.true., nq1=1, nq2=1, nq3=1
/
P
if [ $VAR = nk2 ]; then L="nice -n 10 $B/mpirun -np 2 --bind-to none $B/"; [ -e $B/mpichversion ] && L="nice -n 10 $B/mpirun -np 2 $B/"; NK="-nk 2"; else L="nice -n 10 $B/"; NK=""; fi
ok=1
[ -f fit.in ] && { ${L}pw.x $NK -in fit.in > fit.out 2>&1 || ok=0; }
[ $ok = 1 ] && { ${L}pw.x $NK -in scf.in > scf.out 2>&1 || ok=0; }
rc=NA
[ $ok = 1 ] && { ${L}ph.x $NK -in ph.in > ph.out 2> ph.err; rc=$?; }
segv=$(cat ph.out ph.err 2>/dev/null | grep -c 'SIGSEGV\|Segmentation')
res=OTHER; [ "$rc" = 0 ] && grep -q 'JOB DONE' ph.out && res=OK; [ $segv -gt 0 ] && res=SEGV; [ $ok = 0 ] && res=PWFAIL
lam=$(grep -c 'lambda(' ph.out 2>/dev/null)
echo "$ENV $CASE $VAR -> $res (rc=$rc, lambda_lines=$lam) $(grep -m1 'Program PHONON' ph.out | xargs)"
[ $res = OTHER ] && grep -A3 'Error in routine' ph.out | head -5
true
