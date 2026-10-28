QE=~/micromamba/envs/qe/bin; export OMP_NUM_THREADS=1
W=~/rtsc_work/tetra_Al/q6s_opt; cd $W
[ -f scf.out ] && echo "scf_wall_s (see stage1_q6.log run) 112" > timing.txt
cp ph.out ph.out.part1
sed 's|    amass(1) = 26.9815386|    amass(1) = 26.9815386\n    recover = .true.|' ph.in > ph_recover.in
t1=$(date +%s)
$QE/mpirun -np 4 --bind-to none -x OMP_NUM_THREADS $QE/ph.x -nk 4 -in ph_recover.in > ph.out 2> ph.err || { echo "ph failed" >> timing.txt; exit 1; }
t2=$(date +%s); echo "ph_dfpt_wall_s_after_recover $((t2-t1)) (plus ~11.5 min before the restart)" >> timing.txt
echo STAGE1_DONE >> timing.txt
NQ=6 OCC=tetrahedra_opt NKS="16 24 32" bash /mnt/g/ClaudeWorkspace/RTSC_Project/run_scripts/scratchpad_2026-10-27/t_stage2.sh
