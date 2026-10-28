cd ~/rtsc_work/tetra_Al/src
sed -n 430,480p phq_readin.f90; sed -n 975,1000p phq_readin.f90
grep -n -i "lshift_q\|in_alpha2f" phq_readin.f90 INPUT_PH.def
grep -n -i -A25 "var lshift_q" INPUT_PH.def
echo ====alpha2f; cat alpha2f.f90
echo ====elphon tetra; grep -n -i "tetra" elphon.f90 | head -60
