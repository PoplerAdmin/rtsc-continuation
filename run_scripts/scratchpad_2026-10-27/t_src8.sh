cd ~/rtsc_work/tetra_Al/src
sed -n 60,135p elph_tetra_mod.f90; sed -n 155,265p elph_tetra_mod.f90; sed -n 292,316p elph_tetra_mod.f90
sed -n 405,440p phq_setup.f90
grep -n "SUBROUTINE" elph_tetra_mod.f90
