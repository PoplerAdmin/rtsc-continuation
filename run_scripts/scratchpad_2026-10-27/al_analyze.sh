. ~/rtsc_work/.venv/bin/activate
D=$(ls -d ~/rtsc_work/bench/Al/run/Al-*)/dfpt
T=/mnt/g/ClaudeWorkspace/RTSC_Project/bench_tools/analyze_benchmark.py
# Variant A: smearings 2..10 only (sigma1 = 0.005 Ry excluded explicitly). Copy files, original untouched.
A=~/rtsc_work/bench_out/Al_inputs_s2-10; rm -rf $A; mkdir -p $A
for i in $(seq 2 10); do cp $D/a2F.dos$i $A/; done; cp $D/lambda $A/
( time nice python $T $A Al --out ~/rtsc_work/bench_out/Al_s2-10 --material-id "Al-fcc@0GPa[sigma2-10]" --jobs 4 ) > ~/rtsc_work/bench_out/Al_s2-10.log 2>&1 &
# Variant B: all 10 smearings, negative bins clipped (explicit opt-in)
( time nice python $T $D Al --out ~/rtsc_work/bench_out/Al_clip --material-id "Al-fcc@0GPa[clip]" --negative-policy clip --jobs 4 ) > ~/rtsc_work/bench_out/Al_clip.log 2>&1 &
wait
tail -4 ~/rtsc_work/bench_out/Al_s2-10.log ~/rtsc_work/bench_out/Al_clip.log
