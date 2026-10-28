rm -f ~/rtsc_work/conv/results.txt
for m in pb nb mgb2 h3s h3str pbfr; do bash /mnt/g/ClaudeWorkspace/RTSC_Project/run_scripts/scratchpad_2026-10-27/conv.sh $m; done
echo ALLDONE >> ~/rtsc_work/conv/results.txt
