# Tetra stage 2 (electron_phonon='lambda_tetra'), run in place on q6s_opt/tmp, no copy.
# Derived from run_scripts/scratchpad_2026-10-27/t_stage2_inplace.sh with two changes:
#  - nothing is deleted here: prune_daemon.sh hash-then-deletes finished q points' wavefunctions,
#  - leftovers of an interrupted in-place run are moved aside, not removed.
# usage: NKS="24" setsid nohup bash t_stage2_inplace.sh &
QE=~/micromamba/envs/qe/bin; export OMP_NUM_THREADS=1
W=~/rtsc_work/tetra_Al/q6s_opt; cd $W
for NK in ${NKS:-24}; do
  D=$W/elph_k$NK; grep -q ELPH_DONE $D/timing.txt 2>/dev/null && continue
  mkdir -p $D
  sed "s/nk\([123]\) = 16/nk\1 = $NK/" $W/elph_k16/elph.in > $W/elph_k$NK.in
  cp $W/elph_k$NK.in $D/elph.in
  if ls $W/al.dyn*.elph.* $W/al.a2F.dat $W/al.McMillan.gp $W/elph_dir >/dev/null 2>&1; then
    I=$W/incomplete_$(date +%Y%m%dT%H%M%S); mkdir -p $I
    mv $W/al.dyn*.elph.* $W/al.a2F.dat $W/al.McMillan.gp $W/elph_dir $I/ 2>/dev/null
  fi
  t0=$(date +%s)
  $QE/mpirun -np 4 --bind-to none -x OMP_NUM_THREADS $QE/ph.x -nk 4 -in elph_k$NK.in > $D/elph.out 2> $D/elph.err
  rc=$?; t1=$(date +%s); echo "elph_wall_s $((t1-t0)) rc $rc" > $D/timing.txt
  $QE/mpirun -np 1 $QE/alpha2f.x -in elph_k$NK.in > $D/a2F.out 2> $D/a2F.err; echo "alpha2f_rc $?" >> $D/timing.txt
  mv $W/al.dyn*.elph.* $W/al.a2F.dat $W/al.McMillan.gp $D/ 2>/dev/null; mv $W/elph_dir $D/ 2>/dev/null
  echo ELPH_DONE >> $D/timing.txt
done
