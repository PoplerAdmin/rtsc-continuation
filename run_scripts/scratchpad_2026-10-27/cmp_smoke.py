import glob, sys
from pathlib import Path
sys.path.insert(0, str(Path.home()/"rtsc_work/rtsc_bundle/rtsc"))
from rtsc.dfpt import collect_qe_eph, spectrum_lambda
ref = collect_qe_eph(Path.home()/"rtsc_work/rtsc_bundle/rtsc/tests/data_dfpt/qe67_al_smoke", "fixture")
cands = sorted(glob.glob("/tmp/pytest-of-*/pytest-current/*/work*"))
print("work dirs:", cands)
w = Path(cands[-1])
ws = [p for p in w.rglob("a2F.dos1")]
print("a2F found in", ws)
new = collect_qe_eph(ws[0].parent, "new")
print(" i   lam_QE7.6   lam_QE6.7   rel.diff   A_new(K)  A_ref(K)")
for i,(a,b) in enumerate(zip(new.spectra, ref.spectra),1):
    la, lb = spectrum_lambda(a), spectrum_lambda(b)
    print(f"{i:2d}  {la:9.4f}  {lb:9.4f}  {la/lb-1:+8.3%}  {sum(a.weight_K):8.2f}  {sum(b.weight_K):8.2f}")
