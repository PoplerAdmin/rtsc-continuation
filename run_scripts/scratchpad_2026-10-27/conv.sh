# pw.x ecutwfc convergence tests (CONVERGENCE TESTS ONLY, not physical results). 2 ranks, nice.
Q=$HOME/micromamba/envs/qe/bin; export OMP_NUM_THREADS=1
S=$HOME/pseudo/sg15; E=$HOME/pseudo/SSSP_1.3.0_PBE_efficiency
R=~/rtsc_work/conv; mkdir -p $R
pw() { # name ecut dual pd sys species positions kgrid extra_sys
  local n=$1 ec=$2 du=$3 pd=$4 sys=$5 spc=$6 pos=$7 k=$8 ex=$9
  local d=$R/$n; mkdir -p $d; cd $d
  cat > e$ec.in <<P
&control
 calculation='scf', prefix='c$ec', outdir='./tmp/', pseudo_dir='$pd', tstress=.true., tprnfor=.true.
/
&system
 $sys, ecutwfc=$ec, ecutrho=$((ec*du)), occupations='smearing', smearing='mv', degauss=0.02 $ex
/
&electrons
 conv_thr=1e-10, mixing_beta=0.5
/
ATOMIC_SPECIES
$spc
$pos
K_POINTS automatic
 $k 0 0 0
P
  local t0=$(date +%s.%N)
  nice -n 10 $Q/mpirun -np 2 --bind-to none $Q/pw.x -nk 2 -in e$ec.in > e$ec.out 2>&1
  local t1=$(date +%s.%N)
  local et=$(grep '^!' e$ec.out | awk '{print $5}'); local p=$(grep 'P=' e$ec.out | tail -1 | awk '{print $6}')
  local nat=$(grep 'number of atoms/cell' e$ec.out | awk '{print $5}')
  echo "$n ecut=$ec E=$et nat=$nat P=$p wall=$(echo "$t1-$t0"|bc) mem=$(grep -m1 'Estimated max dynamical RAM' e$ec.out | awk '{print $(NF-1),$NF}')" | tee -a $R/results.txt
  rm -rf tmp
}
PBSYS="ibrav=2, celldm(1)=9.354, nat=1, ntyp=1"; PBPOS=$'ATOMIC_POSITIONS crystal\n Pb 0 0 0'
NBSYS="ibrav=3, celldm(1)=6.236, nat=1, ntyp=1"; NBPOS=$'ATOMIC_POSITIONS crystal\n Nb 0 0 0'
MBSYS="ibrav=4, celldm(1)=5.832, celldm(3)=1.1419, nat=3, ntyp=2"
MBPOS=$'ATOMIC_POSITIONS crystal\n Mg 0 0 0\n B 0.333333333 0.666666667 0.5\n B 0.666666667 0.333333333 0.5'
HSYS="ibrav=3, celldm(1)=5.839, nat=4, ntyp=2"
HPOS=$'ATOMIC_POSITIONS alat\n S 0 0 0\n H 0.5 0 0\n H 0 0.5 0\n H 0 0 0.5'
case "$1" in
 pb)  for e in 30 40 50 60 70 80 100 120; do pw pb_sg15 $e 4 $S "$PBSYS" " Pb 207.2 Pb_ONCV_PBE-1.2.upf" "$PBPOS" "8 8 8"; done ;;
 pbfr) for e in 40 60 80 120; do pw pb_sg15fr_soc $e 4 $S "$PBSYS" " Pb 207.2 Pb_ONCV_PBE_FR-1.0.upf" "$PBPOS" "8 8 8" ", noncolin=.true., lspinorb=.true."; done ;;
 nb)  for e in ${NBLIST:-40 50 60 70 80 90 100 120}; do pw nb_sg15 $e 4 $S "$NBSYS" " Nb 92.906 Nb_ONCV_PBE-1.2.upf" "$NBPOS" "8 8 8"; done ;;
 mgb2) for e in 40 50 60 70 80 90 100 120; do pw mgb2_sg15 $e 4 $S "$MBSYS" $' Mg 24.305 Mg_ONCV_PBE-1.2.upf\n B 10.811 B_ONCV_PBE-1.2.upf' "$MBPOS" "8 8 6"; done ;;
 h3s) for e in 40 50 60 70 80 90 100 120 140; do pw h3s_sg15 $e 4 $S "$HSYS" $' S 32.06 S_ONCV_PBE-1.2.upf\n H 1.008 H_ONCV_PBE-1.2.upf' "$HPOS" "8 8 8"; done ;;
 h3str) # transferability: SG15 NC vs SSSP-eff US (GBRV S + pslibrary H) at two lattice constants
   for a in 5.839 5.631; do
     pw h3s_tr_sg15_a$a 100 4 $S "ibrav=3, celldm(1)=$a, nat=4, ntyp=2" $' S 32.06 S_ONCV_PBE-1.2.upf\n H 1.008 H_ONCV_PBE-1.2.upf' "$HPOS" "8 8 8"
     pw h3s_tr_us_a$a 80 10 $E "ibrav=3, celldm(1)=$a, nat=4, ntyp=2" $' S 32.06 s_pbe_v1.4.uspp.F.UPF\n H 1.008 H.pbe-rrkjus_psl.1.0.0.UPF' "$HPOS" "8 8 8"
   done ;;
esac
