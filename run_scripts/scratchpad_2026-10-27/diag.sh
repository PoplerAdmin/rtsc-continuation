source ~/rtsc_work/.venv/bin/activate; python /mnt/g/ClaudeWorkspace/RTSC_Project/run_scripts/scratchpad_2026-10-27/diag.py
sed -n 1420,1440p ~/rtsc_work/epw_Al/wannier/pwscf.wout; grep -n -A3 "Energy windows\|Outer:\|Inner:" ~/rtsc_work/epw_Al/wannier/pwscf.wout | head -20
