cd ~/rtsc_work/tetra_Al/src
curl -sSfL -o dfpt_tetra_mod.f90 https://raw.githubusercontent.com/QEF/q-e/qe-7.6/LR_Modules/dfpt_tetra_mod.f90 || curl -sSfL -o dfpt_tetra_mod.f90 https://raw.githubusercontent.com/QEF/q-e/qe-7.6/PHonon/PH/dfpt_tetra_mod.f90
wc -l dfpt_tetra_mod.f90; grep -n "SUBROUTINE dfpt_tetra_main" -A80 dfpt_tetra_mod.f90 | head -120
grep -n "nntetra\|tetra_type" tetra.f90 | sed -n 1,40p
