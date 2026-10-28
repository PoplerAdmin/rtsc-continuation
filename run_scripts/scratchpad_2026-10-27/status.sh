cd ~/rtsc_work
echo "== load"; uptime; nproc; free -g | head -2
echo "== top procs"; ps -eo pid,etime,pcpu,rss,args --sort=-pcpu | grep -v grep | egrep 'pw.x|ph.x|epw.x|wannier|python|mpirun' | head -25 | cut -c1-200
echo "== queue.log"; tail -n 15 bench/queue.log
echo "== wvar logs"; for v in A B; do echo "-- $v"; tail -n 8 epw_Al/wvar_$v.log; done
echo "== tetra"; ls -la tetra_Al 2>/dev/null | head -30
echo "== Pb ph progress"; ls bench_out 2>/dev/null; find bench -name '*.out' -newer bench/queue.log 2>/dev/null | head
