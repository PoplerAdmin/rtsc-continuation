cd ~/rtsc_work
sed -n 1,40p .venv/lib/python3.12/site-packages/torch/_native/common_utils.py
echo "--- TORCHDYNAMO_DISABLE=1 TORCH_COMPILE_DISABLE=1"
TORCHDYNAMO_DISABLE=1 TORCH_COMPILE_DISABLE=1 OMP_NUM_THREADS=2 .venv/bin/python -I /mnt/g/ClaudeWorkspace/RTSC_Project/run_scripts/scratchpad_2026-10-27/tb.py 2>&1 | tail -2
echo "--- TORCH_DISABLE_NATIVE_JIT=1"
TORCH_DISABLE_NATIVE_JIT=1 OMP_NUM_THREADS=2 .venv/bin/python -I /mnt/g/ClaudeWorkspace/RTSC_Project/run_scripts/scratchpad_2026-10-27/tb.py 2>&1 | tail -3
