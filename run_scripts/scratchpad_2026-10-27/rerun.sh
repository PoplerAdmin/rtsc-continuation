cd ~/rtsc_work/mattergen_test
mv prescreen_MatterSim-v1.0.0-5M_0GPa.json prescreen_mgb_MatterSim-v1.0.0-5M_0GPa.json
mv prescreen_MatterSim-v1.0.0-5M_150GPa.json prescreen_mgb_MatterSim-v1.0.0-5M_150GPa.json
mkdir -p superseded && mv prescreen_MatterSim-v1.0.0-*GPa.json superseded/
cd ~/rtsc_work; export TORCH_DISABLE_NATIVE_JIT=1 OMP_NUM_THREADS=2
.venv/bin/python -I /mnt/g/ClaudeWorkspace/RTSC_Project/run_scripts/scratchpad_2026-10-27/ms_prescreen.py h3s 2>/dev/null | grep -v INFO > /mnt/g/ClaudeWorkspace/RTSC_Project/run_scripts/scratchpad_2026-10-27/prescreen_h3s_rerun.log; cat /mnt/g/ClaudeWorkspace/RTSC_Project/run_scripts/scratchpad_2026-10-27/prescreen_h3s_rerun.log
ls ~/rtsc_work/mattergen_test/
