uptime; ps aux | grep -c "[e]pw.x"; ps aux | grep -c "[p]w.x "
cat ~/rtsc_work/epw_Al/wannier/DONE 2>/dev/null; tail -1 ~/rtsc_work/epw_Al/wannier/epw.out
for f in ~/rtsc_work/epw_Al/g1_*.log; do echo == $f; tail -12 $f; done 2>/dev/null
ls ~/rtsc_work/epw_Al/g1/*/ 2>/dev/null
