source ~/rtsc_work/.venv/bin/activate
python /mnt/g/ClaudeWorkspace/RTSC_Project/run_scripts/scratchpad_2026-10-27/patch3.py && python -m py_compile /mnt/g/ClaudeWorkspace/RTSC_Project/bench_tools/analyze_benchmark.py && echo compiled
python - <<'P'
import json; d=json.load(open("/home/user/rtsc_work/bench_out/smoke_Al/Al_benchmark.json"))
for f in d["float_reference"]["per_smearing"][:3]:
    print(f["label"], {m:(round(v["Tc_K"],5), v.get("truncated"), v.get("N_matsubara"), v.get("truncation_shift_eta")) for m,v in f["eliashberg_float"].items()})
for t,r in d["proved"].items(): print(t, r["wall_certify_s"], [ (a["upper_search"] or {}).get("seconds") for a in r["per_alternative"]][:3], [(a["lower_search"] or {}).get("tc_float_lo") for a in r["per_alternative"]][:4])
P
ls ~/rtsc_work/bench/Al/run/Al-fcc_0GPa-*/; tail -3 ~/rtsc_work/bench/Al/run/Al-fcc_0GPa-*/dfpt/ph.out
