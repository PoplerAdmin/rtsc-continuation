cd ~/rtsc_work/tetra_Al; ls q6s_opt; grep -c "Calculation of q" q6s_opt/ph.out; ls q6s_opt/al.dyn*
setsid nohup bash /mnt/g/ClaudeWorkspace/RTSC_Project/run_scripts/scratchpad_2026-10-27/t_resume_opt.sh > chain_opt.log 2>&1 < /dev/null & disown
sleep 30; tail -5 q6s_opt/ph.out; grep -n "recover\|Restart" q6s_opt/ph.out | head
