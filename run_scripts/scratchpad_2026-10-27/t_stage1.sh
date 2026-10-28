set -e
QE=~/micromamba/envs/qe/bin
export OMP_NUM_THREADS=1
NQ=${NQ:-6}
OCC=${OCC:-tetrahedra_lin}
W=~/rtsc_work/tetra_Al/q${NQ}s_${OCC#tetrahedra_}
mkdir -p $W && cd $W
cat > scf.in <<EOI
 &control
    calculation = 'scf'
    prefix = 'al'
    outdir = './tmp/'
    pseudo_dir = '/home/user/pseudo/sg15/'
 /
 &system
    ibrav = 0
    nat = 1
    ntyp = 1
    ecutwfc = 60.0
    ecutrho = 240.0
    occupations = '$OCC'
    nbnd = 8
 /
 &electrons
    conv_thr = 1.0d-12
    mixing_beta = 0.5
 /
ATOMIC_SPECIES
 Al  26.98153860 Al_ONCV_PBE-1.2.upf
CELL_PARAMETERS angstrom
     -0.00000000000000     2.01105392100000     2.01105392100000
      2.01105392100000    -0.00000000000000     2.01105392100000
      2.01105392100000     2.01105392100000    -0.00000000000000
ATOMIC_POSITIONS crystal
 Al      0.00000000000000     0.00000000000000     0.00000000000000
K_POINTS automatic
 16 16 16 0 0 0
EOI
cat > ph.in <<EOI
Al phonon tetra_opt shifted q
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
 /
EOI
t0=$(date +%s)
$QE/mpirun -np 4 --bind-to none -x OMP_NUM_THREADS $QE/pw.x -nk 4 -in scf.in > scf.out 2> scf.err
t1=$(date +%s); echo "scf_wall_s $((t1-t0))" > timing.txt
$QE/mpirun -np 4 --bind-to none -x OMP_NUM_THREADS $QE/ph.x -nk 4 -in ph.in > ph.out 2> ph.err
t2=$(date +%s); echo "ph_dfpt_wall_s $((t2-t1))" >> timing.txt
echo STAGE1_DONE >> timing.txt
