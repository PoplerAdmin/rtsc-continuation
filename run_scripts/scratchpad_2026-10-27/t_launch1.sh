cd ~/rtsc_work/tetra_Al && rm -rf q6s && NQ=6 setsid nohup bash /mnt/g/ClaudeWorkspace/RTSC_Project/run_scripts/scratchpad_2026-10-27/t_stage1.sh > stage1_q6.log 2>&1 < /dev/null & disown
sleep 5; ps aux | grep -c "[p]w.x -nk 4"
