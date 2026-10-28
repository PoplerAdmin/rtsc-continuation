QE=~/micromamba/envs/qe/bin; export OMP_NUM_THREADS=1
cd ~/rtsc_work/tetra_Al/linbug/lin_np1
setarch -R $QE/ph.x -in ph.in > ph_bt.out 2> ph_bt.err; echo rc=$?
grep "^#" ph_bt.err
for a in $(grep -o "0x5555[0-9a-f]*" ph_bt.err); do off=$(printf "0x%x" $((a - 0x555555554000))); echo "$a $(addr2line -f -C -e $QE/ph.x $off | head -1)"; done
tail -5 ph_bt.out
