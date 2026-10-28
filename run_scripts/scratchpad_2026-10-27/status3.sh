cd ~/rtsc_work
tail -n 25 epw_Al/wvar_A/wannier/epw.out; ls epw_Al/wvar_A/wannier | head -30
echo "== Pb verdict"; python3 -c "
import json;d=json.load(open('bench/Pb/run/Pb-fcc_0GPa-fc9e9847a8/verdict.json'));print(json.dumps(d,indent=1)[:2500])"
