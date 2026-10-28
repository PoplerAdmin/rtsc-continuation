p="/mnt/g/ClaudeWorkspace/RTSC_Project/bench_tools/analyze_benchmark.py"; s=open(p,encoding="utf-8").read()
new=open("/mnt/g/ClaudeWorkspace/RTSC_Project/run_scripts/scratchpad_2026-10-27/newblock.py",encoding="utf-8").read()
a=s.index("def eliashberg_eta("); b=s.index("def rtsc_float_tc(")
s=s[:a]+new+s[b:]
s=s.replace("    os.nice(0)\n","")
open(p,"w",encoding="utf-8").write(s)
print("ok")
