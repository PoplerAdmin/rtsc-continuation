until grep -q "frozen window" ~/rtsc_work/epw_Al/wannier/DONE 2>/dev/null; do sleep 180; done; cat ~/rtsc_work/epw_Al/wannier/DONE; ls -la ~/rtsc_work/epw_Al/out/pwscf.epmatwp
