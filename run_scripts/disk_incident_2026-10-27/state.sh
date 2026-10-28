cd ~/rtsc_work
echo "== queue.log tail"; tail -n 8 bench/queue.log
echo "== Pb/Al run state"; for m in Al Pb; do ls bench/$m/run/*/ | tr '\n' ' '; echo; tail -n 3 bench/$m/run*.log 2>/dev/null; done
echo "== wvar logs"; for v in A B; do echo "-- $v"; tail -n 4 epw_Al/wvar_$v.log 2>/dev/null; ls epw_Al/wvar_$v 2>/dev/null | tr '\n' ' '; echo; done
echo "== elph_k24 contents (outside tmp)"; find tetra_Al/q6s_opt/elph_k24 -path '*/tmp' -prune -o -type f -printf '%s %p\n'
echo "== elph_k16 contents (outside tmp)"; find tetra_Al/q6s_opt/elph_k16 -path '*/tmp' -prune -o -type f -printf '%s %p\n' | head -20
echo "== stage-2 markers naming"; ls tetra_Al/q6s_opt/elph_k16 | grep -i elph | head -3
echo "== procs"; pgrep -fa 'ph.x|pw.x|epw.x|wannier|rtsc run' | head
