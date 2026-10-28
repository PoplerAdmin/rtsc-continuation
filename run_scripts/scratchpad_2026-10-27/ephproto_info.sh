export OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
cd /mnt/g/ClaudeWorkspace/RTSC_Project/bench_tools
taskset -c 0 ~/rtsc_work/.venv/bin/python -m pytest tests -q -p no:cacheprovider -v 2>&1 | grep -c PASSED
for f in tests/test_qlist.py tests/test_phself.py tests/test_rtsc_tags.py; do echo -n "$f "; taskset -c 0 ~/rtsc_work/.venv/bin/python -m pytest $f -q -p no:cacheprovider 2>&1 | tail -1; done
taskset -c 0 ~/rtsc_work/.venv/bin/python - <<'PY'
import sys; sys.path[:0]=['.','tests']
import synthetic_epw as syn
from ephproto.phself import build_phself_spectrum
d=syn.fixture_dir('normal')
sp=build_phself_spectrum(d/syn.LAMBDA_FILE,d/syn.LINEWIDTH_FILE,d/syn.QFILE,d/'epw.in',d/'epw.out')
print(sp.summary()); print(sp.provenance['counts'])
PY
find tests ephproto -name __pycache__ -prune -exec rm -rf {} \; 2>/dev/null; true
