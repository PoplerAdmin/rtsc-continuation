cd ~/rtsc_work/tetra_Al/src
curl -sSfL -o tetra.f90 https://raw.githubusercontent.com/QEF/q-e/qe-7.6/PW/src/tetra.f90
grep -n "SUBROUTINE opt_tetra_init\|tetra_type\|wlsm(" tetra.f90 | head -40
curl -sSfL -o pw_restart_new.f90 https://raw.githubusercontent.com/QEF/q-e/qe-7.6/PW/src/pw_restart_new.f90
grep -n -i "tetra_type\|tetrahedra_lin\|tetrahedra_opt" pw_restart_new.f90 | head
curl -sSfL -o qexsd_copy.f90 https://raw.githubusercontent.com/QEF/q-e/qe-7.6/Modules/qexsd_copy.f90; grep -n -i "tetra" qexsd_copy.f90 | head
