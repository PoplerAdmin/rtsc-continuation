cd ~/rtsc_work/tetra_Al/src
python3 -c "
import json;t=json.load(open('tree.json'))
for e in t['tree']:
  p=e['path']
  if 'PHonon/examples' in p and p.count('/')>=2 and p.count('/')<=3 and ('example1' in p or 'tetra' in p.lower()): print(p)
"
grep -rn "tetra" tree.json | head -0
