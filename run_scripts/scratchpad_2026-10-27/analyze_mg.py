import json, re, sys
from pathlib import Path
W = Path.home()/"rtsc_work/mattergen_test"; S = Path(sys.argv[1])
nat = {l.split("\t")[0]: (int(l.split("\t")[1]), l.split("\t")[2]) for l in (S/"inputs.tsv").read_text().splitlines()}
nat["MgB2-exp(KNOWN)"] = (3, "Mg1B2")
for P in ("0", "150"):
    d = json.loads((W/f"prescreen_MatterSim-v1.0.0-5M_{P}GPa.json").read_text())
    print(f"\n### 5M @ {P} GPa  wall {d['wall_s']:.0f} s")
    print("| label | N | formula | sg in | sg relaxed | E_hull ML (eV/at) | min freq (cm-1) | ML stable | lightHF | defer | relax steps/conv | H (eV/at) | V (A3/at) |")
    print("|"+"---|"*13)
    rows = []
    for r in d["rows"]:
        rel = next((n for n in r["notes"] if n.startswith("relax:")), "")
        g = lambda k: (re.search(k+r"=([-\d.e+]+|True|False)", rel) or [None, None])[1]
        e, v, steps, conv = g("energy_eV_per_atom"), g("volume_A3_per_atom"), g("steps"), g("converged")
        H = float(e) + float(P)*float(v)/160.21766208 if e else None
        dfr = [n.split(": ",1)[1].split(" (")[0].split(" eV")[0] for n in r["notes"] if n.startswith("ML_SUGGEST_DEFER")]
        other = [n for n in r["notes"] if "failed" in n or "NOT converged" in n]
        f = lambda x, fmt: (fmt % x) if x is not None else "-"
        n, form = nat[r["label"]]
        rows.append((r, H, float(v) if v else None, n, form))
        print(f"| {r['label']} | {n} | {form} | {r['sg_in']} | {r['sg_relaxed']} | {f(r['e_above_hull_ml_eV'],'%.3f')} | {f(r['min_freq_cm1'],'%.0f')} | {r['dynamically_stable_ml']} | {f(r['light_hf_fraction'],'%.2f')} | {'; '.join(dfr)} | {steps}/{conv} | {f(H,'%.4f')} | {f(float(v) if v else None,'%.3f')} |" + (" "+";".join(other) if other else ""))
        if r.get("a"): print(f"  relaxed a={r['a']:.4f} c={r['c']:.4f}")
    nd = sum(any(n.startswith("ML_SUGGEST_DEFER") for n in r["notes"]) for r,*_ in rows)
    nu = sum(r["dynamically_stable_ml"] is False for r,*_ in rows)
    nh = sum(any(n.startswith("ML_SUGGEST_DEFER: e_above") for n in r["notes"]) for r,*_ in rows)
    nf = sum(not r["ok"] for r,*_ in rows)
    print(f"\nTOTAL {len(rows)}; ML-unstable {nu}; DEFER any {nd}; DEFER by hull {nh}; failed {nf}")
    # duplicates
    groups = []
    for r, H, V, n, form in rows:
        if H is None: continue
        red = r["formula_in"]
        for gq in groups:
            r0, H0, V0 = gq[0]
            if r0["formula_in"] == red and r0["sg_relaxed"] == r["sg_relaxed"] and abs(H-H0) < 0.003 and abs(V-V0)/V0 < 0.01:
                gq.append((r, H, V)); break
        else: groups.append([(r, H, V)])
    for gq in groups:
        if len(gq) > 1: print("DUP:", gq[0][0]["formula_in"], gq[0][0]["sg_relaxed"], [x[0]["label"] for x in gq])
    print("REF notes sample:", [n for n in rows[1][0]["notes"] if n.startswith("ml_hull")])
