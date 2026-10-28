# usage: g1_run.sh LABEL WANNIER_DIR EPMATWP
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
source ~/rtsc_work/.venv/bin/activate
cd /mnt/g/ClaudeWorkspace/RTSC_Project/bench_tools
python -m ephproto.bandcheck run ~/rtsc_work/epw_Al/material.json --wannier "$2" --epmatwp "$3" \
   --label "$1" --np 4 --ef-grids 16 24 > ~/rtsc_work/epw_Al/g1_$1.log 2>&1
echo "exit $?" >> ~/rtsc_work/epw_Al/g1_$1.log
