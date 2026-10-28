# Manual equivalent of verify_all.sh step 1: pip install --no-deps of a copy of rtsc/ into a
# throw-away venv (deps pre-installed there, since uv venvs can't be inherited), imports, entry point.
set -e
export PATH=$HOME/.local/bin:$PATH
W=~/rtsc_work/pkgcheck; rm -rf $W; mkdir -p $W/src
cp -a ~/rtsc_work/rtsc_bundle/rtsc/. $W/src/; rm -rf $W/src/*.egg-info
uv venv -q $W/venv -p 3.12 --seed
$W/venv/bin/python -m pip install -q numpy scipy python-flint==0.9.0
$W/venv/bin/python -m pip install --no-deps $W/src > $W/install.log 2>&1 && echo "pip install --no-deps: OK"
$W/venv/bin/python -c 'import numpy, scipy, flint, rtsc, rtsc.certify, rtsc.cli; print("imports OK; flint", flint.__version__)'
test -x $W/venv/bin/rtsc && echo "entry point OK: $W/venv/bin/rtsc"
