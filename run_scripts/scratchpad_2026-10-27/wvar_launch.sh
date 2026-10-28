cd ~/rtsc_work/epw_Al && . ~/rtsc_work/.venv/bin/activate
P=/mnt/g/ClaudeWorkspace/RTSC_Project/run_scripts/scratchpad_2026-10-27/wvar.py
nohup python -I $P A '{}' > wvar_A.log 2>&1 &
nohup python -I $P B '{"scdm_mu": 3.5}' > wvar_B.log 2>&1 &
sleep 20
pgrep -fa wvar.py; tail -n 3 wvar_A.log wvar_B.log
