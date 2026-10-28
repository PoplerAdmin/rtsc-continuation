source ~/rtsc_work/.venv/bin/activate
cd ~/rtsc_work/bench_out
mkdir -p demo_ana
( time nice -n 10 python /mnt/g/ClaudeWorkspace/RTSC_Project/bench_tools/analyze_benchmark.py demo_wd/SYNTHETIC-weak-einstein MgB2 --tc-lit 10 --out demo_ana --jobs 2 > demo_ana/run.log 2> demo_ana/run.err ) 2>> demo_ana/run.err
echo "demo exit $?"; grep -E "^> |\| T300|\| Tlit|run.*verdict itself|SERIOUS|^- \*\*" demo_ana/run.log | head -20; grep -iE "error|traceback" demo_ana/run.err | head; tail -3 demo_ana/run.err
( time nice -n 10 python /mnt/g/ClaudeWorkspace/RTSC_Project/bench_tools/analyze_benchmark.py /mnt/g/ClaudeWorkspace/RTSC_Project/rtsc_handoff_2026-10-27/rtsc_bundle/rtsc/tests/data_dfpt/qe67_al_smoke Al --out smoke_Al --jobs 2 --material-id Al-smoke > smoke_Al/run.log 2> smoke_Al/run.err ) 2>> smoke_Al/run.err
echo "smoke exit $?"; grep -E "NO\*\*|SERIOUS|†|no proved T_low" smoke_Al/run.log | head -8; tail -3 smoke_Al/run.err
