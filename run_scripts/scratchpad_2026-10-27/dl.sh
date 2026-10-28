python3 - <<'P'
import glob,json,os
E=os.path.expanduser("~/micromamba/envs/")
def recs(e): return {os.path.basename(f):json.load(open(f)).get("size",0) for f in glob.glob(E+e+"/conda-meta/*.json")}
old=recs("qe")
new={}
for e in "qe76mpich qe76b0 qe75b qe75a qe74 qe72 dbg".split(): new.update(recs(e))
d={k:v for k,v in new.items() if k not in old}
print(len(new),"distinct pkgs;",len(d),"not in qe env; size MB %.0f"%(sum(d.values())/1e6))
P
