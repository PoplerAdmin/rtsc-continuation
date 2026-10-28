cd ~/rtsc_work
pkill -f 'while pgrep -f wvar.py' ; echo killed-watcher
for v in A B; do echo "-- $v"; grep -B3 -A8 -i 'error' epw_Al/wvar_$v/wannier/epw.out | head -30; done
echo "== Pb"; ls bench/Pb/run/*/; grep -il -m1 error bench/Pb/run/*/*.json 2>/dev/null; find bench/Pb/run -name '*.log' | head; 
tail -n 30 $(ls -t bench/Pb/*.log bench/Pb/run/*/*.log 2>/dev/null | head -1)
echo "== Nb ph progress"; grep -c 'Calculation of q' bench/Nb/run/*/dfpt/ph.out; grep 'irreducible q' -A0 bench/Nb/run/*/dfpt/ph.out | head -2; grep -m1 'Dynamical matrices for' -A2 bench/Nb/run/*/dfpt/ph.out
echo "== tetra q6s_opt"; ls q6s_opt 2>/dev/null; ls tetra_Al/q6s_opt; tail -n 5 tetra_Al/q6s_opt/*/elph.out 2>/dev/null | tail -8
