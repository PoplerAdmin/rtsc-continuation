export MAMBA_ROOT_PREFIX=~/micromamba; M=~/bin/micromamba
for spec in "qe=7.5=h19104ac_2" "qe=7.5=hbc83ec2_0" "qe=7.4=hbc83ec2_0" "qe=7.6=mpi_mpich_he3edec4_4" "qe=7.6=h2ec8f75_0" "qe=7.2=hac89879_4" "gdb"; do
  $M create -n drytest --dry-run -y -c conda-forge --override-channels "$spec" --json 2>/dev/null > /tmp/dry_$$.json
  python3 - "$spec" /tmp/dry_$$.json <<'P'
import json,sys,os
d=json.load(open(sys.argv[2]))
a=d.get("actions",{})
fetch=a.get("FETCH",[]); link=a.get("LINK",[])
cache=os.path.expanduser("~/micromamba/pkgs")
tot=0;new=[]
for p in link:
    fn=p.get("fn","")
    base=fn.replace(".conda","").replace(".tar.bz2","")
    if not os.path.isdir(os.path.join(cache,base)):
        tot+=p.get("size",0); new.append(base)
print(sys.argv[1], "link=",len(link),"to_download_MB=%.0f"%(tot/1e6), "n_new=",len(new))
print("   ", " ".join(n for n in new if any(k in n for k in ("qe-","openmpi","mpich","hdf5","elpa","gfortran","gdb","libgcc-"))))
P
done
