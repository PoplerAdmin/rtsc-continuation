cd ~/rtsc_work
OMP_NUM_THREADS=2 .venv/bin/python -I /mnt/g/ClaudeWorkspace/RTSC_Project/run_scripts/scratchpad_2026-10-27/tb.py 2>&1 | head -30
T=.venv/lib/python3.12/site-packages/torch/_native
ls -R $T | head -50
grep -rn "environ\|getenv\|TORCH_" $T | head -30
