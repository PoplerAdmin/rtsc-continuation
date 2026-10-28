cd ~/rtsc_work/tetra_Al
NQ=6 OCC=tetrahedra_opt NKS="24 32" setsid nohup bash /mnt/g/ClaudeWorkspace/RTSC_Project/run_scripts/scratchpad_2026-10-27/t_stage2.sh > chain_opt_stage2b.log 2>&1 < /dev/null & disown
sleep 3; pgrep -af t_stage2
