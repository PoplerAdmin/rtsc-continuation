# Gamma-point phonons of H3S at a=5.631 bohr (~205 GPa): SG15 NC vs SSSP-eff US (phonon only, no e-ph). TRANSFERABILITY TEST ONLY.
Q=$HOME/micromamba/envs/qe/bin; export OMP_NUM_THREADS=1
L="nice -n 10 $Q/mpirun -np 2 --bind-to none -x OMP_NUM_THREADS"
for v in sg15 us; do
 d=~/rtsc_work/phtime/h3s_gamma_$v; rm -rf $d; mkdir -p $d; cd $d
 if [ $v = sg15 ]; then pd=$HOME/pseudo/sg15; ec=90; er=360; sp=$' S 32.06 S_ONCV_PBE-1.2.upf\n H 1.008 H_ONCV_PBE-1.2.upf'
 else pd=$HOME/pseudo/SSSP_1.3.0_PBE_efficiency; ec=80; er=800; sp=$' S 32.06 s_pbe_v1.4.uspp.F.UPF\n H 1.008 H.pbe-rrkjus_psl.1.0.0.UPF'; fi
 cat > scf.in <<P
&control
 calculation='scf', prefix='t', outdir='./tmp/', pseudo_dir='$pd'
/
&system
 ibrav=3, celldm(1)=5.631, nat=4, ntyp=2, ecutwfc=$ec, ecutrho=$er, occupations='smearing', smearing='mv', degauss=0.02
/
&electrons
 conv_thr=1e-11
/
ATOMIC_SPECIES
$sp
ATOMIC_POSITIONS alat
 S 0 0 0
 H 0.5 0 0
 H 0 0.5 0
 H 0 0 0.5
K_POINTS automatic
 12 12 12 0 0 0
P
 printf "ph\n&inputph\n prefix='t', outdir='./tmp/', tr2_ph=1e-16, fildyn='dyn'\n/\n0 0 0\n" > ph.in
 $L $Q/pw.x -nk 2 -in scf.in > scf.out 2>&1; $L $Q/ph.x -nk 2 -in ph.in > ph.out 2>&1
 echo "== $v"; grep 'freq (' ph.out | awk '{printf "%s ", $8}'; echo
done
