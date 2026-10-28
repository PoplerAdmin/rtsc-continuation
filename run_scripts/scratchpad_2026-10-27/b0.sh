R=/mnt/g/ClaudeWorkspace/RTSC_Project/run_scripts/scratchpad_2026-10-27/repro.sh
for c in al_paw_interp al_paw_simple al_sg15_interp al_paw_noeph mg_us_interp; do bash $R qe76b0 $c 2>/dev/null; done
