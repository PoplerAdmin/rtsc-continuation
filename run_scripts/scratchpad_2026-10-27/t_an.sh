source ~/rtsc_work/.venv/bin/activate; cd ~/rtsc_work/tetra_Al
grep -B2 -A12 "Calculation of alpha2F" q6s_opt/elph_k16/a2F.out | grep -v "^\s*$"
python /mnt/g/ClaudeWorkspace/RTSC_Project/run_scripts/scratchpad_2026-10-27/analyze.py q6s_opt/elph_k16 > analysis_k16.json; python -c "
import json;r=json.load(open('analysis_k16.json'))
for t in r['tetra']:
  print({k:v for k,v in t.items() if k not in ('nk_per_q',)})
for s in r['smear']: print(s)
print(r['shape_L1_vs_smear4'])
"
