import sys, numpy as np
sys.path.insert(0,'/mnt/g/ClaudeWorkspace/RTSC_Project/bench_tools')
from ephproto.bandcheck import *
r='/home/user/rtsc_work/epw_Al/g1/'
kd,Ed=parse_pw_xml(r+'pw_bands/tmp/pwscf.save/data-file-schema.xml')
lab=sys.argv[1]
kw,Ew=parse_band_eig(r+lab+'/epw_bands/band.eig')
D=Ed[:,4:]; ef=11.109
print('DFT bands 1..8 at k1', Ed[0,:9].round(3)); print('Wan k1', Ew[0].round(3))
d=np.abs(D-ef).min(axis=1); print('min |E-EF| DFT per k: smallest 10 (meV)', np.sort(d)[:10]*1000)
d=np.abs(Ew-ef).min(axis=1); print('min |E-EF| Wan: smallest 10 (meV)', np.sort(d)[:10]*1000)
for h in (0.08,0.15,0.25,0.5):
  print(h, 'DFT states', int((np.abs(D-ef)<=h).sum()), 'Wan states', int((np.abs(Ew-ef)<=h).sum()))
# band-by-band diff below EF+2 for lowest bands
for j in range(4):
  m=D[:,j]<ef+2
  print('band',j+1,'n',m.sum(),'max|dE| meV',(np.abs(Ew[m,j]-D[m,j]).max()*1000 if m.any() else None))
