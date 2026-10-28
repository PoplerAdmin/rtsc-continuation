export MAMBA_ROOT_PREFIX=~/micromamba; M=~/bin/micromamba
$M repoquery search -c conda-forge "qe>=7.0" --json 2>/dev/null | python3 -c '
import json,sys
d=json.load(sys.stdin)
pk=d.get("result",{}).get("pkgs",[])
for p in sorted(pk,key=lambda p:(p["version"],p["build"])):
  print(p["version"],p["build"],p.get("size",0)//1024,"kB", p.get("timestamp",0), [x for x in p["depends"] if x.split()[0] in ("openmpi","mpich","libgfortran5","libgcc","libgcc-ng","elpa","hdf5","fftw")])
'
