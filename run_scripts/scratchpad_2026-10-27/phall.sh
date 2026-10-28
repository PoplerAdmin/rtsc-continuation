P=/mnt/g/ClaudeWorkspace/RTSC_Project/run_scripts/scratchpad_2026-10-27
rm -f ~/rtsc_work/phtime/results.txt
bash $P/phtime.sh al; bash $P/phtime.sh pb; bash $P/phtime.sh nb 120; bash $P/phtime.sh mgb2 90; bash $P/phtime.sh h3s 90
bash $P/h3sgam.sh >> ~/rtsc_work/phtime/results.txt 2>&1
echo ALLDONE >> ~/rtsc_work/phtime/results.txt
