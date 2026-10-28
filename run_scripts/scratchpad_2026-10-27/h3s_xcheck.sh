cd ~/rtsc_work; export TORCH_DISABLE_NATIVE_JIT=1 OMP_NUM_THREADS=2
.venv/bin/python -I /mnt/g/ClaudeWorkspace/RTSC_Project/run_scripts/scratchpad_2026-10-27/h3s_xcheck.py 2>&1 | grep -v -i "warn\|INFO"
