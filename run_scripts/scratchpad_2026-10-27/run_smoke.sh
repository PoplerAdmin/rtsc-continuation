source ~/rtsc_work/.venv/bin/activate
mkdir -p ~/rtsc_work/bench_out/smoke_Al
cd ~/rtsc_work/bench_out/smoke_Al
/usr/bin/time -v nice -n 10 python /mnt/g/ClaudeWorkspace/RTSC_Project/bench_tools/analyze_benchmark.py \
  /mnt/g/ClaudeWorkspace/RTSC_Project/rtsc_handoff_2026-10-27/rtsc_bundle/rtsc/tests/data_dfpt/qe67_al_smoke Al \
  --out ~/rtsc_work/bench_out/smoke_Al --jobs 2 --material-id Al-smoke > run.log 2> run.err
echo "exit $?" >> run.log
