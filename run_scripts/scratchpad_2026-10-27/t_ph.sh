R=~/rtsc_work/bench/Al/run/Al-fcc_0GPa-d2e010742b/dfpt
grep -n "MPI processes\|pools\|threads\|PHONON.*WALL\|number of k points\|Kohn-Sham states" $R/ph.out | head -30; grep "PHONON.*WALL" $R/ph.out | tail -2
grep -c "Calculation of q" $R/ph.out
cat $R/dfpt_state.json | head -60
ls ~/rtsc_work/.venv/bin | head; ps aux --sort=-%cpu | head -8
