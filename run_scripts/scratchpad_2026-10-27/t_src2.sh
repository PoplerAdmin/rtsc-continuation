cd ~/rtsc_work/tetra_Al/src
sed -n 85,107p ~/rtsc_work/rtsc_bundle/HANDOFF.md
echo ======INPUTPH; grep -n -i -B2 -A40 "var electron_phonon" INPUT_PH.def | head -120
echo ======readin; grep -n -i "tetra" phq_readin.f90
