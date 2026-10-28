cd ~/rtsc_work
export PYTHON=$HOME/rtsc_work/.venv/bin/python
export PATH=$HOME/rtsc_work/.venv/bin:$PATH
time ./rtsc_bundle/verify_all.sh --cert --skip-python --skip-lean --workdir ~/rtsc_work/verify_py
