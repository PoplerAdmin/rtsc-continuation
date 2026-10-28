~/rtsc_work/.venv/bin/python - <<'PY'
import spglib, numpy as np
a=2.011053921
lat=np.array([[0,a,a],[a,0,a],[a,a,0]])
cell=(lat,[[0,0,0]],[13])
for n in (4,6,8):
  for sh in ([0,0,0],[1,1,1]):
    m,g=spglib.get_ir_reciprocal_mesh([n]*3,cell,is_shift=sh)
    print(n,sh,len(np.unique(m)))
PY
nproc; cat /proc/cpuinfo | grep "model name" | head -1; free -g | head -2
