cd ~/rtsc_work/epw_Al/wvar_A/wannier
grep -i -A5 'error\|exiting' pwscf.wout | head -30; tail -n 15 pwscf.wout; grep -iE 'dis_|num_bands|num_wann|froz|win_' pwscf.win
