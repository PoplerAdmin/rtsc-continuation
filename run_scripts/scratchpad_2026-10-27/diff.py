import re
R="/mnt/g/ClaudeWorkspace/RTSC_Project/run_scripts/scratchpad_2026-10-27/wsl/src/raw/"
src=open(R+"EPW_src_readin.f90").read()
m=re.search(r"NAMELIST / inputepw / &(.*?)\n\s*!\s*-{10}",src,re.S)
body=re.sub(r"!.*","",m.group(1))
keys=set(k.strip().lower() for k in body.replace("&"," ").replace("\n"," ").split(",") if k.strip())
txt=open("/mnt/g/ClaudeWorkspace/RTSC_Project/rtsc_handoff_2026-10-27/rtsc_bundle/rtsc/rtsc/dfpt/epw_inputs.py").read()
m2=re.search(r'EPW53_KEYS = frozenset\("""(.*?)"""',txt,re.S)
old=set(m2.group(1).split())
print(len(keys),len(old))
print("REMOVED:",sorted(old-keys))
print("NEW:",sorted(keys-old))
open(R+"../keys76.txt","w").write(" ".join(sorted(keys)))
