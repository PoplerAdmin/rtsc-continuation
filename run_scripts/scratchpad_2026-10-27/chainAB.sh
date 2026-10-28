export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
S=/mnt/g/ClaudeWorkspace/RTSC_Project/run_scripts/scratchpad_2026-10-27
while pgrep -x epw.x >/dev/null; do sleep 20; done
cd /mnt/g/ClaudeWorkspace/RTSC_Project/bench_tools
W=~/rtsc_work/epw_Al/wvar_A
python3 - <<'PY'
import json
m=json.load(open('/home/user/rtsc_work/epw_Al/material.json')); m['epw_dir']='/home/user/rtsc_work/epw_Al/wvar_A'
json.dump(m,open('/home/user/rtsc_work/epw_Al/wvar_A/material.json','w'))
PY
ls -la $W/out/pwscf.epmatwp
~/rtsc_work/.venv/bin/python -m ephproto.bandcheck run $W/material.json --wannier $W/wannier --epmatwp $W/out/pwscf.epmatwp --label A --np 4 --ef-grids 16 24 --kset random > ~/rtsc_work/epw_Al/g1_A.log 2>&1
echo "G1 A rc=$?"; tail -n 30 ~/rtsc_work/epw_Al/g1_A.log
echo "=== variant B (serial after A)"
cd ~/rtsc_work/epw_Al && rm -rf wvar_B/wannier
~/rtsc_work/.venv/bin/python $S/wvar.py B '{"scdm_mu": 3.5}' > wvar_B.log 2>&1
echo "B rc=$?"; tail -n 30 wvar_B.log
