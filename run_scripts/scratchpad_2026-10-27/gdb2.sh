export OMP_NUM_THREADS=1
B=~/micromamba/envs/qe/bin; G=~/micromamba/envs/dbg/bin/gdb
cd ~/rtsc_work/qebug/gdb_qe76; rm -rf tmp/_ph0 dyn* ph.out
nm $B/ph.x | grep -i -E ' (adddvscf_|adddvscf__|__lrus_MOD_becpt|__lrus_MOD_becp1|elphel_)$'
nice -n 10 $G -q -batch -x /mnt/g/ClaudeWorkspace/RTSC_Project/run_scripts/scratchpad_2026-10-27/gdb.cmd --args $B/ph.x -in ph.in > gdb2.out 2>&1
grep -A30 '=== at adddvscf_' gdb2.out
