kill 2379; kill 2500 2503 2504 2505 2506; sleep 3; pgrep -f "ph.x -nk 4" ; 
cd ~/rtsc_work/tetra_Al && mv q6s q6s_opt && sed -i "s|tetra_Al/q6s|tetra_Al/q6s_opt|" q6s_opt/ph.in 2>/dev/null; true
cat > chain_lin.sh <<'EOC'
S=/mnt/g/ClaudeWorkspace/RTSC_Project/run_scripts/scratchpad_2026-10-27
NQ=6 OCC=tetrahedra_lin bash $S/t_stage1.sh && NQ=6 OCC=tetrahedra_lin NKS="16 24 32" bash $S/t_stage2.sh
EOC
setsid nohup bash chain_lin.sh > chain_lin.log 2>&1 < /dev/null & disown
sleep 5; ls ~/rtsc_work/tetra_Al; pgrep -f "pw.x -nk 4" | wc -l
