set -x
B=~/micromamba/envs/qe/bin
ls $B | grep -iE "epw|wannier|alpha2f|lambda|pw.x|ph.x|q2r|matdyn"
W=/mnt/g/ClaudeWorkspace/RTSC_Project/run_scripts/scratchpad_2026-10-27/wsl
mkdir -p $W/run $W/src
cd $W/run
echo "" | timeout 60 $B/epw.x > epw_empty.out 2>&1; head -80 epw_empty.out
printf '&inputepw\n prefix="x"\n/\n' > e.in; timeout 60 $B/epw.x -i e.in > epw_min.out 2>&1; head -100 epw_min.out | tail -60
cp $B/epw_pp.py $W/src/epw_pp.py 2>/dev/null; wc -l $B/epw_pp.py
for f in ~/pseudo/sg15/*.upf; do echo "$f PSWFC=$(grep -c '<PP_PSWFC' $f) CHI=$(grep -c '<PP_CHI' $f) $(grep -o 'number_of_wfc="[^"]*"' $f)"; done
ls ~/pseudo
cd $W/src
curl -sL "https://api.github.com/repos/QEF/q-e/git/trees/qe-7.6?recursive=1" -o tree.json; ls -la tree.json
python3 -c "
import json;t=json.load(open('tree.json'));print(t.get('truncated'));
[print(x['path']) for x in t['tree'] if x['path'].startswith('EPW/') and (x['path'].endswith('.f90') or x['path'].endswith('.py') or 'bin' in x['path'])]" > epw_files.txt
wc -l epw_files.txt
