cd ~/rtsc_work
export PYTHON=$HOME/rtsc_work/.venv/bin/python
export PATH=$HOME/rtsc_work/.venv/bin:$HOME/.elan/bin:$PATH
time ./rtsc_bundle/verify_all.sh --skip-python --skip-pytest --skip-demo --workdir ~/rtsc_work/verify_lean
