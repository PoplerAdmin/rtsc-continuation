for p in $(pgrep -u "$USER" -f "ph.x -nk 4|mpirun -np 4|t_stage1.sh"); do c=$(readlink /proc/$p/cwd); echo "$p $c $(tr '\0' ' ' < /proc/$p/cmdline | cut -c1-80)"; done
