# usage: g1_both.sh LABEL WANNIER_DIR EPMATWP
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
source ~/rtsc_work/.venv/bin/activate
cd /mnt/g/ClaudeWorkspace/RTSC_Project/bench_tools
L=~/rtsc_work/epw_Al/g1_$1.log; : > $L
for ks in mp5 random; do
 echo "== kset $ks" >> $L
 python -m ephproto.bandcheck run ~/rtsc_work/epw_Al/material.json --wannier "$2" --epmatwp "$3" \
   --label "$1" --np 4 --ef-grids 16 24 --kset $ks --nrand 1000 >> $L 2>&1
 echo "exit $?" >> $L
done
