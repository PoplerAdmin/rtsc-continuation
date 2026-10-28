P="python3 /mnt/g/ClaudeWorkspace/RTSC_Project/bench_tools/jobs/prune.py"
python3 --version
cd ~/rtsc_work
sum() { awk '$1=="DRY"{s+=$2;n++} END{printf "  %d files, %.2f GB\n", n, s/1e9}'; }
echo "== Al wfc (finished run)"; $P wfc bench/Al/run/*/dfpt/tmp --reason dry --dry-run | sum
echo "== Al wfc: anything matched that is not *.wfcN?"; $P wfc bench/Al/run/*/dfpt/tmp --reason dry --dry-run | awk '$1=="DRY"{print $3}' | grep -vE '\.wfc[0-9]+$' | head -3
N=$(ls -d bench/Nb/run/*/dfpt)
echo "== Nb qdone (interrupted at q18)"; $P qdone $N/tmp/_ph0 $N --reason dry --dry-run | awk '$1=="DRY"{print $3}' | grep -oE 'q_[0-9]+' | sort -t_ -k2 -n | uniq -c | tr '\n' ' '; echo
$P qdone $N/tmp/_ph0 $N --reason dry --dry-run | sum
echo "== Nb q dirs present / dyn files"; ls $N/tmp/_ph0 | grep -oE 'q_[0-9]+$' | sort -t_ -k2 -n | tail -3 | tr '\n' ' '; ls $N | grep -E '^dyn[0-9]+$' | sort -V | tail -3 | tr '\n' ' '; echo
echo "== tetra q6s_opt layout"; ls tetra_Al/q6s_opt/tmp/_ph0 | head -4; ls tetra_Al/q6s_opt/tmp/_ph0/*.q_2 2>/dev/null | head
echo "== tetra safe dirs, all-mode dry"; for d in tetra_Al/q6s_opt/elph_k16/tmp tetra_Al/q6s_opt/elph_k24 tetra_Al/q6s_lin/tmp; do echo " $d"; $P all $d --reason dry --dry-run | sum; done
echo "== manifest writes (dry must not create)"; ls bench/Al/run/*/dfpt/tmp/PRUNED_MANIFEST.tsv 2>&1 | tail -1
