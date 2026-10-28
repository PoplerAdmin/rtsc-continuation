export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
cd /mnt/g/ClaudeWorkspace/RTSC_Project/bench_tools
if [ "$1" = "gen" ]; then taskset -c 0 ~/rtsc_work/.venv/bin/python tests/synthetic_epw.py; fi
taskset -c 0 ~/rtsc_work/.venv/bin/python -m pytest tests -q -p no:cacheprovider -x ${@:2} 2>&1 | tail -40
