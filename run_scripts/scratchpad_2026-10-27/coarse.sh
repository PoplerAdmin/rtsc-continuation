export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
source ~/rtsc_work/.venv/bin/activate
python /mnt/g/ClaudeWorkspace/RTSC_Project/run_scripts/scratchpad_2026-10-27/coarse.py frozen ~/rtsc_work/epw_Al/wannier ~/rtsc_work/epw_Al/out/pwscf.epmatwp
W=~/rtsc_work/epw_Al/wannier/pwscf.wout
grep -n "Final Spread\|<-- DIS\|dis_num_iter\|dis_conv_tol\|Omega I  \|WF centre and spread" $W | tail -8
grep -A6 "Final State" $W | head -8
grep -B2 -A3 "Final Omega_I\|Omega_I" $W | tail -6
