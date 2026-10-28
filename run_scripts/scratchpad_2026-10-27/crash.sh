# Crash characterisation: serial ph.x, Gamma-only (nq=1) unless NQ set. Usage: crash.sh lane
Q=$HOME/micromamba/envs/qe/bin; export OMP_NUM_THREADS=1
E=$HOME/pseudo/SSSP_1.3.0_PBE_efficiency; P=$HOME/pseudo/SSSP_1.3.0_PBE_precision; S=$HOME/pseudo/sg15
case_() { # name pd ec er k kf nq ephmode sysblock species positions
  local n=$1 pd=$2 ec=$3 er=$4 k=$5 kf=$6 nq=$7 eph=$8 sys=$9 spc=${10} pos=${11}
  local d=~/rtsc_work/crashchar/$n; rm -rf $d; mkdir -p $d; cd $d
  for kind in fit scf; do
    if [ $kind = fit ]; then kk="$kf"; extra=", la2F=.true."; else kk="$k"; extra=""; fi
    [ "$eph" = none ] && [ $kind = fit ] && continue
    cat > $kind.in <<P
&control
 calculation='scf', prefix='t', outdir='./tmp/', pseudo_dir='$pd'
/
&system
 $sys, ecutwfc=$ec, ecutrho=$er,
 occupations='smearing', smearing='mv', degauss=0.02 $extra
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
  if [ "$eph" = none ]; then EP=""; else EP="electron_phonon='$eph', el_ph_sigma=0.005, el_ph_nsigma=10,"; fi
  cat > ph.in <<P
ph
&inputph
 prefix='t', outdir='./tmp/', tr2_ph=1e-14, fildyn='dyn', fildvscf='dvscf',
 $EP
 trans=.true., ldisp=.true., nq1=$nq, nq2=$nq, nq3=$nq
/
P
  t0=$(date +%s)
  ok=1
  if [ -f fit.in ]; then nice -n 10 $Q/pw.x -in fit.in > fit.out 2>&1 || ok=0; fi
  [ $ok = 1 ] && { nice -n 10 $Q/pw.x -in scf.in > scf.out 2>&1 || ok=0; }
  [ $ok = 1 ] && { nice -n 10 $Q/ph.x -in ph.in > ph.out 2> ph.err; rc=$?; }
  t1=$(date +%s)
  echo "$n: pw_ok=$ok ph_rc=$rc lambda_lines=$(grep -c 'lambda(' ph.out) segv=$(grep -c SIGSEGV ph.err) done=$(grep -c 'JOB DONE' ph.out) t=$((t1-t0))s" | tee -a ~/rtsc_work/crashchar/summary.txt
}
FCC() { echo "ibrav=2, celldm(1)=$1, nat=1, ntyp=1"; }
ALPOS=$'ATOMIC_POSITIONS crystal\n X 0 0 0'
mk() { echo "${ALPOS/X/$1}"; }
MGSYS="ibrav=4, celldm(1)=6.05, celldm(3)=1.624, nat=2, ntyp=1"
MGPOS=$'ATOMIC_POSITIONS crystal\n Mg 0.333333333 0.666666667 0.25\n Mg 0.666666667 0.333333333 0.75'
MGB2SYS="ibrav=4, celldm(1)=5.833, celldm(3)=1.141, nat=3, ntyp=2"
MGB2POS=$'ATOMIC_POSITIONS crystal\n Mg 0 0 0\n B 0.333333333 0.666666667 0.5\n B 0.666666667 0.333333333 0.5'
mkdir -p ~/rtsc_work/crashchar
if [ "$1" = A ]; then
  case_ al_paw_noeph $P 30 240 "8 8 8" "16 16 16" 1 none "$(FCC 7.65)" " Al 26.98 Al.pbe-n-kjpaw_psl.1.0.0.UPF" "$(mk Al)"
  case_ al_paw_eph_simple $P 30 240 "8 8 8" "16 16 16" 1 simple "$(FCC 7.65)" " Al 26.98 Al.pbe-n-kjpaw_psl.1.0.0.UPF" "$(mk Al)"
  case_ al_paw_eph_interp $P 30 240 "8 8 8" "16 16 16" 1 interpolated "$(FCC 7.65)" " Al 26.98 Al.pbe-n-kjpaw_psl.1.0.0.UPF" "$(mk Al)"
  case_ mg_us_gbrv_interp $P 30 240 "6 6 4" "12 12 8" 1 interpolated "$MGSYS" " Mg 24.305 mg_pbe_v1.4.uspp.F.UPF" "$MGPOS"
  case_ mg_paw_interp $E 30 240 "6 6 4" "12 12 8" 1 interpolated "$MGSYS" " Mg 24.305 Mg.pbe-n-kjpaw_psl.0.3.0.UPF" "$MGPOS"
  case_ mgb2_us_gbrv_interp $P 35 280 "6 6 6" "12 12 12" 1 interpolated "$MGB2SYS" $' Mg 24.305 mg_pbe_v1.4.uspp.F.UPF\n B 10.811 B_pbe_v1.01.uspp.F.UPF' "$MGB2POS"
else
  case_ pb_paw_interp $P 30 240 "8 8 8" "16 16 16" 1 interpolated "$(FCC 9.35)" " Pb 207.2 Pb.pbe-dn-kjpaw_psl.0.2.2.UPF" "$(mk Pb)"
  case_ nb_paw_interp $P 35 280 "8 8 8" "16 16 16" 1 interpolated "ibrav=3, celldm(1)=6.24, nat=1, ntyp=1" " Nb 92.906 Nb.pbe-spn-kjpaw_psl.0.3.0.UPF" "$(mk Nb)"
  case_ s_us_gbrv_bcc_interp $P 30 240 "8 8 8" "16 16 16" 1 interpolated "ibrav=3, celldm(1)=4.5, nat=1, ntyp=1" " S 32.06 s_pbe_v1.4.uspp.F.UPF" "$(mk S)"
  case_ al_paw_eph_interp_q2 $P 30 240 "8 8 8" "16 16 16" 2 interpolated "$(FCC 7.65)" " Al 26.98 Al.pbe-n-kjpaw_psl.1.0.0.UPF" "$(mk Al)"
  case_ al_sg15_eph_interp_q2 $S 30 120 "8 8 8" "16 16 16" 2 interpolated "$(FCC 7.65)" " Al 26.98 Al_ONCV_PBE-1.2.upf" "$(mk Al)"
fi
