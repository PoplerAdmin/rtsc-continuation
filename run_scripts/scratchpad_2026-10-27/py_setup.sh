set -e
export PATH=$HOME/.local/bin:$PATH
cd ~/rtsc_work
uv python install 3.12 2>&1 | tail -2
uv venv .venv -p 3.12 --seed 2>&1 | tail -2
. .venv/bin/activate
python --version
uv pip install -e "rtsc_bundle/rtsc[test]" 2>&1 | tail -15
