mkdir -p ~/pseudo/qe_std && cd ~/pseudo/qe_std
[ -s Al.pz-vbc.UPF ] || curl -sSL -o Al.pz-vbc.UPF https://pseudopotentials.quantum-espresso.org/upf_files/Al.pz-vbc.UPF
head -3 Al.pz-vbc.UPF; md5sum Al.pz-vbc.UPF
grep -rh "pz-vbc\|md5" ~/rtsc_work/rtsc_bundle/rtsc/tests/data_dfpt/qe67_al_smoke/*.in 2>/dev/null | head -3
cd ~/rtsc_work/rtsc_bundle/rtsc
. ~/rtsc_work/.venv/bin/activate
export RTSC_QE_BIN=$HOME/micromamba/envs/qe/bin PATH=$HOME/micromamba/envs/qe/bin:$PATH
export RTSC_RUN_QE=1 RTSC_QE_PSEUDO_DIR=$HOME/pseudo/qe_std RTSC_QE_NP=4 OMP_NUM_THREADS=1
time python -m pytest -q -p no:cacheprovider tests/test_dfpt_qe_real.py -rs 2>&1 | tail -20
