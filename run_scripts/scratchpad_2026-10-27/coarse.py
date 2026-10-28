import sys, re, json, shutil, numpy as np
from pathlib import Path
sys.path.insert(0,'/mnt/g/ClaudeWorkspace/RTSC_Project/bench_tools')
from ephproto.bandcheck import _epw_dir, _sh, parse_band_eig, QE_BIN
from eph_converge import epw_input
m=json.load(open('/home/user/rtsc_work/epw_Al/material.json'))
txt=open('/home/user/rtsc_work/epw_Al/nscf.out').read()
# nscf k list (crystal) from nscf.in, eigenvalue blocks from nscf.out
kin=open('/home/user/rtsc_work/epw_Al/nscf.in').read().split('K_POINTS crystal')[1].split('\n')[2:514]
kc=np.array([[float(x) for x in l.split()[:3]] for l in kin])
blocks=re.findall(r'k =.*?bands \(ev\):\s*\n(.*?)\n\s*\n\s*occupation', txt, re.S)
E=np.array([[float(x) for x in re.findall(r'-?\d+\.\d+', b)] for b in blocks])
print('nscf blocks', E.shape)
sel=list(range(0,512,13))[:40]
lab=sys.argv[1]; wan=Path(sys.argv[2]); ep=Path(sys.argv[3])
d=Path('/home/user/rtsc_work/epw_Al/g1_coarsecheck')/lab
od=_epw_dir(m,d,wan,ep)
open(d/'kc.txt','w').write(f"{len(sel)} crystal\n"+''.join(f"{float(kc[i,0])!r} {float(kc[i,1])!r} {float(kc[i,2])!r} 1.0\n" for i in sel))
open(d/'gamma.txt','w').write("1 crystal\n0.0 0.0 0.0 1.0\n")
kw={"outdir":str(od)+"/","wannierize":False,"nbndsub":4,"epwwrite":False,"epwread":True,"elph":True,
    "band_plot":True,"filkf":"kc.txt","filqf":"gamma.txt","efermi_read":True,"fermi_energy":11.109,
    "bands_skipped":"exclude_bands = 1-4"}
open(d/'epw.in','w').write(epw_input(m,**kw))
_sh(f"{QE_BIN}/epw.x -nk 4 -in epw.in", d, d/'epw.out', 4)
k,W=parse_band_eig(d/'band.eig')
D=np.sort(E[sel][:,4:],1)[:,:4]
dd=W-D; mask=D<13.109
print(lab,'coarse-grid max|dE| inside frozen window (meV):', 1000*np.abs(dd[mask]).max(), ' n', mask.sum())
print(' first rows DFT', D[:3].round(4).tolist()); print(' W', W[:3].round(4).tolist())
