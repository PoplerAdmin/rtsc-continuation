cd ~/rtsc_work
export TORCH_DISABLE_NATIVE_JIT=1 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
for w in h3s mg; do
  echo "### START $w $(date -Is)"
  /usr/bin/time -v .venv/bin/python -I /mnt/g/ClaudeWorkspace/RTSC_Project/run_scripts/scratchpad_2026-10-27/ms_prescreen.py $w > /mnt/g/ClaudeWorkspace/RTSC_Project/run_scripts/scratchpad_2026-10-27/prescreen_$w.log 2> /mnt/g/ClaudeWorkspace/RTSC_Project/run_scripts/scratchpad_2026-10-27/prescreen_$w.err
  echo "### END $w rc=$? $(date -Is)"
done
