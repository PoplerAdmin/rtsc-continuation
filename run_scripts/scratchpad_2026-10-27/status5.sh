cd ~/rtsc_work/epw_Al/wvar_A
ls -la . wannier; cat wannier/CRASH 2>/dev/null | head; head -40 ../wvar_A.log; dmesg 2>/dev/null | tail -5
sed -n 1,40p /mnt/g/ClaudeWorkspace/RTSC_Project/run_scripts/scratchpad_2026-10-27/wvar.py
