B=~/micromamba/envs/qe/bin/epw.x
strings $B | grep -E "mp_mesh_q|full uniform (k|q)-mesh|irr wedge|Using uniform MP q-mesh" | sort -u
