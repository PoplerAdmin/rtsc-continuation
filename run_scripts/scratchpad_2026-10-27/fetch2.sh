W=/mnt/g/ClaudeWorkspace/RTSC_Project/run_scripts/scratchpad_2026-10-27/wsl/src
cd $W; mkdir -p raw2
for p in PHonon/PH/alpha2f.f90 PHonon/PH/elph_tetra_mod.f90 PHonon/PH/phq_readin.f90 PHonon/PH/elphon.f90 PHonon/examples/tetra_example/run_example PHonon/examples/tetra_example/README PHonon/examples/tetra_example/reference/al.elph.in PHonon/examples/tetra_example/reference/al.ph.in PHonon/examples/tetra_example/reference/al.scf.in PHonon/examples/tetra_example/reference/aluminum.a2F.dat PHonon/examples/tetra_example/reference/al.elph.out PHonon/examples/tetra_example/reference/al.dyn2.elph.2 PHonon/Doc/INPUT_PH.def EPW/examples/mgb2/epw/epw.in EPW/examples/mgb2/epw/nscf.in EPW/examples/pb/woSOC/epw/epw.in test-suite/epw_metal/epw1.in test-suite/epw_metal/pp.in test-suite/epw_super/pp.in test-suite/epw_super/epw1.in test-suite/epw_super/epw2.in test-suite/epw_super/epw3.in test-suite/epw_super/epw4.in test-suite/epw_super/epw5.in test-suite/epw_super/epw6.in test-suite/epw_super/epw7.in test-suite/epw_super/epw8.in test-suite/epw_super/nscf_epw.in test-suite/epw_base/epw1.in test-suite/epw_base/pp.in test-suite/epw_trev_uspp/epw1.in test-suite/epw_trev_paw/epw1.in EPW/src/wannier2bloch_opt.f90 EPW/src/supercond_aniso.f90 PW/src/input.f90; do
  f=raw2/$(echo $p | tr / _); curl -sfL "https://raw.githubusercontent.com/QEF/q-e/qe-7.6/$p" -o $f || echo "FAIL $p"
done
grep -il scdm raw2/* ../src/raw/* | head
python3 -c "
import json;t=json.load(open('tree.json'))
for x in t['tree']:
  p=x['path']
  if p.startswith('test-suite/benchmark') and 'epw_super' in p: print(p)
  if p.startswith('test-suite/epw_super') or p.startswith('test-suite/epw_metal'): print(p)
" | head -60
