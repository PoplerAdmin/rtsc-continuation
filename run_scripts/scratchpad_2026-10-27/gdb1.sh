export OMP_NUM_THREADS=1
B=~/micromamba/envs/qe/bin; G=~/micromamba/envs/dbg/bin/gdb
D=~/rtsc_work/qebug/gdb_qe76; rm -rf $D; cp -r ~/rtsc_work/qebug/qe/al_paw_interp_base $D; cd $D; rm -rf tmp/_ph0 dyn* ph.out ph.err
nice -n 10 $G -q -batch -ex 'set pagination off' -ex 'handle SIGSEGV stop' -ex run -ex 'bt 20' -ex 'info registers rip rax rbx rcx rdx rsi rdi r8 r9 r10 r11 r12 r13 r14 r15' -ex 'x/12i $pc-30' -ex 'frame 1' -ex 'x/6i $pc-12' -ex 'info frame' --args $B/ph.x -in ph.in > gdb.out 2>&1
grep -v '^\s*$' gdb.out | grep -A80 'SIGSEGV' | head -90
