# One-off cleanup after the C:-full crash (policy: reports/EVIDENCE_RETENTION.md). Hash-then-delete only.
P="python3 /mnt/g/ClaudeWorkspace/RTSC_Project/bench_tools/jobs/prune.py"
cd ~/rtsc_work/tetra_Al/q6s_opt
$P all elph_k16/tmp --reason "copy of q6s_opt/tmp used for the finished nk=16 lambda_tetra run; results kept in elph_k16/"
[ -d elph_k24 ] && mv elph_k24 elph_k24_crashed_20261027
$P all elph_k24_crashed_20261027/tmp --reason "scratch of nk=24 lambda_tetra run killed by the C:-full WSL crash; logs and partial outputs kept"
cd ~/rtsc_work/bench
A=$(ls -d Al/run/*/dfpt); $P wfc $A/tmp --reason "Al benchmark finished (verdict + certificate written); wavefunctions regenerable"
N=$(ls -d Nb/run/*/dfpt); $P qdone $N/tmp/_ph0 $N --reason "q point finished (dyn complete, next q started)"
$P inventory $A
df -h /mnt/g | tail -1; du -sh ~/rtsc_work/tetra_Al ~/rtsc_work/bench/Al ~/rtsc_work/bench/Nb
