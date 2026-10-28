journalctl --list-boots 2>/dev/null | tail -3
journalctl -b -1 --no-pager 2>/dev/null | grep -iE 'oom|out of memory|killed process|panic|hung' | tail -10
journalctl -b -1 --no-pager 2>/dev/null | tail -5
cd ~/rtsc_work
echo "== queue"; tail -n 4 bench/queue.log
echo "== Nb"; grep -c 'Calculation of q' bench/Nb/run/*/dfpt/ph.out; tail -n 3 bench/Nb/run/*/dfpt/ph.out
echo "== wvar"; for v in A B; do echo "-- $v"; tail -n 12 epw_Al/wvar_$v.log | cut -c1-200; done; ls -la epw_Al/wvar_A/out
echo "== tetra"; ls -la --time-style=+%H:%M tetra_Al/q6s_opt | tail -5; tail -n 3 tetra_Al/q6s_opt/ph.out
cat ~/.wslconfig 2>/dev/null; cat /mnt/c/Users/user/.wslconfig 2>/dev/null
