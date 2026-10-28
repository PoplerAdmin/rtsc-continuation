cd ~/rtsc_work/bench/Al/run/*/dfpt
ls -la; cat dfpt_state.json | head -40
grep -E 'Calculation of q|total cpu time|PHONON.*WALL|number of k points|Number of q|Dynamical matrices for' ph.out | awk '/Calculation of q/{print; next} /PHONON/{print; next} /Number|Dynamical|number/{print}' | head -60
grep -c 'Calculation of q' ph.out; tail -3 ph.out
grep -E 'PWSCF.*WALL|number of k points' scf.out ../dfpt/*.out 2>/dev/null | head
ls elph_dir | head; free -m; ps -eo pid,rss,etime,args | grep -E 'ph.x|pw.x' | grep -v grep | head
date
