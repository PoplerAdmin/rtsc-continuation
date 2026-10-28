# In-place e-ph runs (no tmp copy). usage: NKS="24 32" setsid nohup bash t_stage2_inplace.sh &
QE=~/micromamba/envs/qe/bin; export OMP_NUM_THREADS=1
W=~/rtsc_work/tetra_Al/q6s_opt; cd $W
rm -rf $W/elph_k24 $W/elph_k16/tmp
for NK in ${NKS:-24 32}; do
  D=$W/elph_k$NK; grep -q ELPH_DONE $D/timing.txt 2>/dev/null && continue
  mkdir -p $D
  sed "s/nk\([123]\) = 16/nk\1 = $NK/" $W/elph_k16/elph.in > $W/elph_k$NK.in
  cp $W/elph_k$NK.in $D/elph.in
  rm -f $W/al.dyn*.elph.* $W/al.a2F.dat $W/al.McMillan.gp; rm -rf $W/elph_dir
  # disk guard: drop nscf wavefunctions of q points already finished (dv files are kept)
  ( while true; do for f in $W/al.dyn*.elph.*; do [ -e "$f" ] || continue; n=${f##*.}; rm -f $W/tmp/_ph0/al.q_$n/al.wfc* $W/tmp/_ph0/al.q_$n/*.wfc* 2>/dev/null; done; sleep 60; done ) & G=$!
  t0=$(date +%s)
  $QE/mpirun -np 4 --bind-to none -x OMP_NUM_THREADS $QE/ph.x -nk 4 -in elph_k$NK.in > $D/elph.out 2> $D/elph.err
  rc=$?; kill $G; t1=$(date +%s); echo "elph_wall_s $((t1-t0)) rc $rc" > $D/timing.txt
  $QE/mpirun -np 1 $QE/alpha2f.x -in elph_k$NK.in > $D/a2F.out 2> $D/a2F.err; echo "alpha2f_rc $?" >> $D/timing.txt
  mv $W/al.dyn*.elph.* $W/al.a2F.dat $W/al.McMillan.gp $D/ 2>/dev/null; mv $W/elph_dir $D/ 2>/dev/null
  echo ELPH_DONE >> $D/timing.txt
done
