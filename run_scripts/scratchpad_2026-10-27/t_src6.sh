cd ~/rtsc_work/tetra_Al/src
B=https://raw.githubusercontent.com/QEF/q-e/qe-7.6/PHonon/examples/tetra_example
curl -sSfL -o tetra_README $B/README; curl -sSfL -o tetra_run_example $B/run_example
cat tetra_README; sed -n '/cat > /,$p' tetra_run_example | grep -v "^\s*$" | head -250
python3 -c "
import json;t=json.load(open('tree.json'))
for e in t['tree']:
  p=e['path']
  if 'tetra_example/reference' in p: print(p)
"
