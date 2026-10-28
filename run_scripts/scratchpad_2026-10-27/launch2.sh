S=/mnt/g/ClaudeWorkspace/RTSC_Project/run_scripts/scratchpad_2026-10-27
setsid nohup bash $S/dfpt_queue2.sh > ~/rtsc_work/bench/queue2.out 2>&1 < /dev/null &
setsid nohup bash $S/chainAB2.sh > ~/rtsc_work/epw_Al/chainAB2.out 2>&1 < /dev/null &
sleep 60
grep -m1 -i recover ~/rtsc_work/bench/Nb/run/*/dfpt/ph.in; tail -n 3 ~/rtsc_work/bench/Nb/run_resume.log; pgrep -c ph.x; pgrep -c epw.x
