cd ~/rtsc_work; export TORCH_DISABLE_NATIVE_JIT=1 OMP_NUM_THREADS=2
.venv/bin/python -I /mnt/g/ClaudeWorkspace/RTSC_Project/run_scripts/scratchpad_2026-10-27/alphaB.py 2>&1 | grep -v -i "warn\|INFO"
grep -o "supercell shortest length[^;]*" mattergen_test/prescreen_MatterSim-v1.0.0-5M_*GPa.json | sort | uniq -c | sort -rn | head -3
grep -c "supercell shortest" mattergen_test/prescreen_MatterSim-v1.0.0-5M_0GPa.json mattergen_test/prescreen_MatterSim-v1.0.0-5M_150GPa.json
