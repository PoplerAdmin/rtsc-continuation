cd ~/rtsc_work
TORCH_DISABLE_NATIVE_JIT=1 OMP_NUM_THREADS=2 .venv/bin/python -I /mnt/g/ClaudeWorkspace/RTSC_Project/run_scripts/scratchpad_2026-10-27/cmp_dev.py 2>&1 | grep -v INFO | tail -3
.venv/bin/python -I -c "import torch,triton,mattersim,ase,pymatgen,phonopy;print(torch.__version__,triton.__version__,ase.__version__,pymatgen.__version__ if hasattr(pymatgen,'__version__') else '',phonopy.__version__)" 2>&1|tail -1
.venv/bin/pip show mattersim 2>/dev/null | grep Version; nvidia-smi --query-gpu=name,driver_version,memory.used --format=csv
