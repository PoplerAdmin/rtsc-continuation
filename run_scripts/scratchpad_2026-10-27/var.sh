R=/mnt/g/ClaudeWorkspace/RTSC_Project/run_scripts/scratchpad_2026-10-27/repro.sh
bash $R qe al_paw_simple nofildvscf 2>/dev/null
bash $R qe al_paw_interp nofildvscf 2>/dev/null
bash $R qe al_paw_interp nosym 2>/dev/null
bash $R qe al_paw_interp noinv 2>/dev/null
