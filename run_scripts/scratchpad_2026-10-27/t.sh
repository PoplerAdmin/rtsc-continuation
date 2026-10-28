cd /mnt/g/ClaudeWorkspace/RTSC_Project/bench_tools && OMP_NUM_THREADS=1 nice -n 19 ~/rtsc_work/.venv/bin/python -m pytest -q tests 2>&1 | tail -3
