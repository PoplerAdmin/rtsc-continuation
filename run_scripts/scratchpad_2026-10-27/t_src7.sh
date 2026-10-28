cd ~/rtsc_work/tetra_Al/src
grep -n "SUBROUTINE\|errore\|nk1\|nbnd\|lshift\|wlsm\|elph_tetra ==\|ef\b\|WRITE(stdout" elph_tetra_mod.f90 | head -120
echo =====elphon around 360-420, 600-660
sed -n 395,420p elphon.f90; sed -n 610,640p elphon.f90
grep -n -i "tetra\|lshift\|nk1\|nscf" do_phonon.f90 phq_setup.f90 | head -40
