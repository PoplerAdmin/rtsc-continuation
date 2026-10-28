for i in $(seq 1 110); do grep -q "^exit" ~/rtsc_work/bench_out/smoke_Al/run.log 2>/dev/null && break; sleep 5; done
cat ~/rtsc_work/bench_out/smoke_Al/run.log | tail -3; ls ~/rtsc_work/bench_out/smoke_Al
