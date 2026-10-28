import sys, numpy as np
sys.path.insert(0,'/mnt/g/ClaudeWorkspace/RTSC_Project/bench_tools')
from ephproto.bandcheck import *
r='/home/user/rtsc_work/epw_Al/g1_random1000_s20261027/'
kd,Ed=parse_pw_xml(r+'pw_bands/tmp/pwscf.save/data-file-schema.xml')
for lab in ['frozen','discarded_nofrozen']:
  kw,Ew=parse_band_eig(r+lab+'/epw_bands/band.eig')
  D=np.sort(Ed[:,4:],1); W=np.sort(Ew,1); ef=11.109
  print('==',lab)
  for j in range(4):
    for lo,hi in [(-99,ef-5),(ef-5,ef-1),(ef-1,ef+1),(ef+1,ef+2),(ef+2,99)]:
      m=(D[:,j]>=lo)&(D[:,j]<hi)
      if m.any():
        d=(W[m,j]-D[m,j])*1000
        print(f' band{j+1} DFT in [{lo-ef:+.0f},{hi-ef:+.0f}] eV rel EF: n={m.sum():4d} mean {d.mean():8.1f} max|d| {np.abs(d).max():8.1f} meV')
  # number of DFT non-excluded bands below EF+2 per k
  nf=(D<ef+2).sum(1); print(' DFT states below EF+2 per k: max', nf.max(), 'hist', np.bincount(nf))
