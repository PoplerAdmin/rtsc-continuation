source ~/rtsc_work/.venv/bin/activate
python /mnt/g/ClaudeWorkspace/RTSC_Project/run_scripts/scratchpad_2026-10-27/mkcfg.py
for m in Pb Nb MgB2 H3S; do echo "===== $m"; cd ~/rtsc_work/bench/$m; rtsc run -c config.toml --plan-only 2>&1 | tail -25; echo "rc=$?"; done
