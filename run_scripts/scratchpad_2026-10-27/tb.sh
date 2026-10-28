cd ~/rtsc_work; ls rtsc_bundle/rtsc/rtsc/mlff; grep -rn "compile\|triton\|filter\|Frechet\|Unit" rtsc_bundle/rtsc/rtsc/mlff/*.py | grep -v re.compile | head -30
OMP_NUM_THREADS=2 .venv/bin/python -I /mnt/g/ClaudeWorkspace/RTSC_Project/run_scripts/scratchpad_2026-10-27/tb.py 2>&1 | tail -40
