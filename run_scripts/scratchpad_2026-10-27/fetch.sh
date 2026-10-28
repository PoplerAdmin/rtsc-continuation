W=/mnt/g/ClaudeWorkspace/RTSC_Project/run_scripts/scratchpad_2026-10-27/wsl/src
cd $W
python3 -c "
import json;t=json.load(open('tree.json'))
for x in t['tree']:
  p=x['path']
  if (p.startswith('EPW/examples') and (p.endswith('.in') )) or (p.startswith('PHonon/PH/') and ('alpha2f' in p or 'tetra' in p or 'elph' in p or 'phq_readin' in p)) or 'Doc/INPUT_EPW' in p or (p.startswith('EPW/') and 'INPUT' in p) or (p.startswith('PHonon/examples') and 'tetra' in p.lower()) or (p.startswith('test-suite/epw') and p.endswith('.in'))  or (p.startswith('PHonon/Doc') and 'INPUT_PH' in p): print(p)
" > list2.txt
wc -l list2.txt
mkdir -p raw
for p in EPW/src/readin.f90 EPW/src/input.f90 EPW/src/utilities/bcast_input.f90 EPW/src/supercond.f90 EPW/src/supercond_common.f90 EPW/src/supercond_driver.f90 EPW/src/supercond_iso.f90 EPW/src/io/io_supercond.f90 EPW/src/wannierization.f90 EPW/src/pw2wan.f90 EPW/src/wannier.f90 EPW/src/use_wannier.f90 EPW/src/printing.f90 EPW/src/bin/pp.py EPW/bin/pp.py EPW/src/wannier2bloch.f90 EPW/src/epw.f90 EPW/src/utilities/bzgrid.f90 EPW/src/init.f90 EPW/src/setups.f90 EPW/src/io/io.f90 EPW/src/close.f90 EPW/src/global_var.f90 EPW/src/utilities/kfold.f90 EPW/src/wann_common.f90; do
  f=raw/$(echo $p | tr / _); curl -sfL "https://raw.githubusercontent.com/QEF/q-e/qe-7.6/$p" -o $f || echo "FAIL $p"
done
ls -la raw | head -50
