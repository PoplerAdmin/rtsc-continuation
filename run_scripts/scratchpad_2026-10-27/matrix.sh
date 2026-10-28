R=/mnt/g/ClaudeWorkspace/RTSC_Project/run_scripts/scratchpad_2026-10-27/repro.sh
for e in qe76b0 qe76mpich qe75b qe75a qe74 qe72; do
 for c in al_paw_interp al_paw_simple al_sg15_interp al_paw_noeph mg_us_interp; do bash $R $e $c 2>/dev/null; done
done
