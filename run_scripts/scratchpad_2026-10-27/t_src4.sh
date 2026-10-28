cd ~/rtsc_work/tetra_Al/src
curl -sSfL -o elph_tetra_mod.f90 https://raw.githubusercontent.com/QEF/q-e/qe-7.6/PHonon/PH/elph_tetra_mod.f90 || echo FAIL
wc -l elph_tetra_mod.f90
curl -sSfL "https://api.github.com/repos/QEF/q-e/git/trees/qe-7.6?recursive=1" -o tree.json || echo treefail
python3 -c "
import json;t=json.load(open('tree.json'))
print(t.get('truncated'))
for e in t['tree']:
  p=e['path']
  if 'PHonon/examples' in p and p.count('/')<=3: print(p)
" | head -80
grep -rl "tetra" --include=*.in . 2>/dev/null | head
