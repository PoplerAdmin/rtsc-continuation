source ~/rtsc_work/.venv/bin/activate; cd ~/rtsc_work/tetra_Al/src/ref
python - <<'PY'
import numpy as np
from rtsc.dfpt.units import to_kelvin_factor
t=np.loadtxt("aluminum.a2F.dat"); w,a=t[:,0],t[:,1]; dE=w[0]
print("n",len(w),"neg",(a<0).sum(),"min",a.min(), "lam_bins",2*np.sum(a*dE/w))
PY
rtsc certify aluminum.a2F.dat --units Ry --rule bins --target 300 --out /tmp/x.json 2>&1 | tail -4
