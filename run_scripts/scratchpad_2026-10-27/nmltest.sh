B=~/micromamba/envs/qe/bin
cd /mnt/g/ClaudeWorkspace/RTSC_Project/run_scripts/scratchpad_2026-10-27/wsl/run
for kv in "eps_acustic=2.0" "kmaps=.false." "eps_acoustic=2.0" "a2f_iso=.true." "dos_tetra=.true."; do
printf -- "--\n&inputepw\n prefix='x'\n $kv\n/\n" > t.in
echo "== $kv"; OMP_NUM_THREADS=1 timeout 60 $B/epw.x -i t.in 2>&1 | grep -E "Error|Bad line|xml data file|namelist" | head -3
done
