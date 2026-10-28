export MAMBA_ROOT_PREFIX=~/micromamba; M=~/bin/micromamba
for e in qe qe76mpich qe76b0 qe75b qe75a qe74 qe72 dbg; do echo "## $e"; $M list -n $e 2>/dev/null | grep -E '^\s*(qe|openmpi|mpich|libgfortran5|elpa|hdf5|fftw|libopenblas|scalapack|gdb)\s'; done
md5sum ~/pseudo/SSSP_1.3.0_PBE_precision/Al.pbe-n-kjpaw_psl.1.0.0.UPF ~/pseudo/SSSP_1.3.0_PBE_precision/mg_pbe_v1.4.uspp.F.UPF ~/pseudo/sg15/Al_ONCV_PBE-1.2.upf
python3 - <<'P'
import os,json,glob,time
t0=time.mktime(time.strptime("2026-10-27 15:20","%Y-%m-%d %H:%M"))
tot=0;n=0
for d in glob.glob(os.path.expanduser("~/micromamba/pkgs/*/info/repodata_record.json")):
    pd=os.path.dirname(os.path.dirname(d))
    if os.path.getmtime(d)>=t0 or os.path.getctime(pd)>=t0:
        r=json.load(open(d)); tot+=r.get("size",0); n+=1
print("new packages:",n,"download MB: %.0f"%(tot/1e6))
P
du -sm ~/micromamba/envs/{qe76mpich,qe76b0,qe75b,qe75a,qe74,qe72,dbg} | tr '\n' ' '
