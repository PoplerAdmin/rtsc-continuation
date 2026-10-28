# ph.x e-ph TIMING tests (reduced k, q 2x2x2; 2 ranks; CONVERGENCE/TIMING TESTS ONLY)
Q=$HOME/micromamba/envs/qe/bin; export OMP_NUM_THREADS=1; S=$HOME/pseudo/sg15
R=~/rtsc_work/phtime; mkdir -p $R
L="nice -n 10 $Q/mpirun -np 2 --bind-to none -x OMP_NUM_THREADS"
t() { # name ec sys species pos k kf nq
  local n=$1 ec=$2 sys=$3 spc=$4 pos=$5 k=$6 kf=$7 nq=$8
  local d=$R/$n; rm -rf $d; mkdir -p $d; cd $d
  for kind in fit scf; do
    if [ $kind = fit ]; then kk="$kf"; ex=", la2F=.true."; else kk="$k"; ex=""; fi
    cat > $kind.in <<P
&control
 calculation='scf', prefix='t', outdir='./tmp/', pseudo_dir='$S'
/
&system
 $sys, ecutwfc=$ec, ecutrho=$((4*ec)), occupations='smearing', smearing='mv', degauss=0.02 $ex
/
&electrons
 conv_thr=1e-10
/
ATOMIC_SPECIES
$spc
$pos
K_POINTS automatic
 $kk 0 0 0
P
  done
  cat > ph.in <<P
ph
&inputph
 prefix='t', outdir='./tmp/', tr2_ph=1e-16, fildyn='dyn', fildvscf='dvscf',
 electron_phonon='interpolated', el_ph_sigma=0.005, el_ph_nsigma=10,
 trans=.true., ldisp=.true., $nq
/
P
  /usr/bin/time -f "%e s %M KB" $L $Q/pw.x -nk 2 -in fit.in > fit.out 2> fit.time
  /usr/bin/time -f "%e s %M KB" $L $Q/pw.x -nk 2 -in scf.in > scf.out 2> scf.time
  /usr/bin/time -f "%e s %M KB" $L $Q/ph.x -nk 2 -in ph.in > ph.out 2> ph.time
  echo "== $n fit: $(tail -1 fit.time) nk_fit=$(grep -m1 'number of k points' fit.out|awk '{print $5}') | scf: $(tail -1 scf.time) | ph: $(tail -1 ph.time) done=$(grep -c 'JOB DONE' ph.out) lambda=$(grep -c 'lambda(' ph.out)" | tee -a $R/results.txt
  grep -E 'Calculation of q|number of k points|PHONON.*WALL|Max static|dynamical RAM' ph.out | grep -v '^ *$' | tee -a $R/results.txt
  rm -rf tmp/_ph0/*.dwf* tmp/_ph0/*.bar* 2>/dev/null
}
FCC() { echo "ibrav=2, celldm(1)=$1, nat=1, ntyp=1"; }
BCC() { echo "ibrav=3, celldm(1)=$1, nat=$2, ntyp=$3"; }
case "$1" in
 al)   t al_cal 60 "$(FCC 7.65)" " Al 26.98 Al_ONCV_PBE-1.2.upf" $'ATOMIC_POSITIONS crystal\n Al 0 0 0' "8 8 8" "16 16 16" "nq1=2,nq2=2,nq3=2" ;;
 pb)   t pb 50 "$(FCC 9.354)" " Pb 207.2 Pb_ONCV_PBE-1.2.upf" $'ATOMIC_POSITIONS crystal\n Pb 0 0 0' "8 8 8" "16 16 16" "nq1=2,nq2=2,nq3=2" ;;
 nb)   t nb $2 "$(BCC 6.236 1 1)" " Nb 92.906 Nb_ONCV_PBE-1.2.upf" $'ATOMIC_POSITIONS crystal\n Nb 0 0 0' "8 8 8" "16 16 16" "nq1=2,nq2=2,nq3=2" ;;
 mgb2) t mgb2 $2 "ibrav=4, celldm(1)=5.832, celldm(3)=1.1419, nat=3, ntyp=2" $' Mg 24.305 Mg_ONCV_PBE-1.2.upf\n B 10.811 B_ONCV_PBE-1.2.upf' $'ATOMIC_POSITIONS crystal\n Mg 0 0 0\n B 0.333333333 0.666666667 0.5\n B 0.666666667 0.333333333 0.5' "6 6 4" "12 12 8" "nq1=2,nq2=2,nq3=2" ;;
 h3s)  t h3s $2 "$(BCC 5.631 4 2)" $' S 32.06 S_ONCV_PBE-1.2.upf\n H 1.008 H_ONCV_PBE-1.2.upf' $'ATOMIC_POSITIONS alat\n S 0 0 0\n H 0.5 0 0\n H 0 0.5 0\n H 0 0 0.5' "6 6 6" "12 12 12" "nq1=2,nq2=2,nq3=2" ;;
esac
