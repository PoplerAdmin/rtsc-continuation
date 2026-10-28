cd ~/rtsc_work/bench/Pb
grep -nE 'ERROR|Error|error' run.log | head -15; echo ----; grep -E 'ERROR|REJECTED|UNDECIDED|CERTIFIED|EXCLUDED' run.log | head
python3 -c "import json;v=json.load(open([__import__('glob').glob('run/*/verdict.json')][0][0]));print({k:v[k] for k in list(v)[:8]})" 2>&1 | cut -c1-600
