R=~/rtsc_work/bench/Al/run/Al-fcc_0GPa-d2e010742b/dfpt
ls -la $R
for f in scf.in scf_fit.in ph.in matdyn.in; do echo "=== $f"; cat $R/$f; done
ls ~/rtsc_work/bench/Al/run/Al-fcc_0GPa-d2e010742b/
head -5 $R/a2F.dos1; cat $R/lambda | head -40
nproc; uptime
sed -n '/§3/,/§4/p' ~/rtsc_work/rtsc_bundle/HANDOFF.md | head -80
