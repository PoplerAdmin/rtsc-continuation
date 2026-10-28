cd ~/rtsc_work/epw_Al
grep -m3 -i "number of electrons\|Fermi\|lattice parameter" nscf.out
cat wannier/epw.in; tail -5 wannier/epw.out; cat wannier/pwscf.win; cat discarded/wannier_nofrozen/epw.in discarded/wannier_nofrozen/DONE
grep -i "fermi" discarded/wannier_nofrozen/epw.out | head; ls discarded/v3_nofrozen/*; ls tmp/pwscf.save | head; nproc
