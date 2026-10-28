S=/mnt/g/ClaudeWorkspace/RTSC_Project/run_scripts/scratchpad_2026-10-27
until grep -q 'qe72 mg_us_interp' $S/matrix.log; do sleep 5; done
bash $S/repro.sh qe al_paw_interp nk2 2>/dev/null
bash $S/repro.sh qe75b al_paw_interp nk2 2>/dev/null
# lambda comparison
for d in qe/al_sg15_interp_base qe75b/al_sg15_interp_base qe74/al_sg15_interp_base qe72/al_sg15_interp_base qe75b/al_paw_interp_base qe75a/al_paw_interp_base qe74/al_paw_interp_base qe75b/al_paw_interp_nk2 qe75b/mg_us_interp_base qe74/mg_us_interp_base; do
  echo "$d: $(grep -m1 -A0 'freq' ~/rtsc_work/qebug/$d/ph.out | head -1 | xargs) | $(grep 'lambda(' ~/rtsc_work/qebug/$d/ph.out | sed -n 5p | xargs)"; done
