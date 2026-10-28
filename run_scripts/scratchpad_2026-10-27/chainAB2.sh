export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
S=/mnt/g/ClaudeWorkspace/RTSC_Project/run_scripts/scratchpad_2026-10-27
cd ~/rtsc_work/epw_Al
for v in A B; do
  if [ $v = A ]; then o='{}'; else o='{"scdm_mu": 3.5}'; fi
  rm -rf wvar_$v/wannier wvar_$v/out/pwscf.epmatwp
  /usr/bin/time -f "wall %e s" ~/rtsc_work/.venv/bin/python $S/wvar.py $v "$o" > wvar_$v.log 2>&1
  echo "=== $v rc=$?"; tail -n 30 wvar_$v.log
done
