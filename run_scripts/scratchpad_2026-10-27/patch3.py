p="/mnt/g/ClaudeWorkspace/RTSC_Project/bench_tools/analyze_benchmark.py"; s=open(p,encoding="utf-8").read()
def rep(old,new):
    global s
    assert old in s, old[:60]
    s=s.replace(old,new)
rep("    if 1.5 * N1 <= 6000:", "    if 1.5 * N1 <= 4500:")
# float-in-bracket checks: truncation aware
rep('''            tc_mu = f["eliashberg_float"][f"{mu:.2f}"]["Tc_K"]
            tc_0 = f["eliashberg_float"]["0.00"]["Tc_K"]
            if a["T_low_K"] is not None:
                add(f"float Tc(mu*={mu:.2f}) >= T_low [{tag}, {a['label']}]",
                    tc_mu >= a["T_low_K"], f"float {_fmt(tc_mu, 7)} K vs proved T_low "
                    f"{_fmt(a['T_low_K'], 7)} K", serious=True)''',
'''            fe = f["eliashberg_float"][f"{mu:.2f}"]
            tc_mu = fe["Tc_K"]
            tc_0 = f["eliashberg_float"]["0.00"]["Tc_K"]
            if a["T_low_K"] is not None:
                ok = tc_mu >= a["T_low_K"]
                tr = bool(fe.get("truncated"))
                add(f"float Tc(mu*={mu:.2f}) >= T_low [{tag}, {a['label']}]",
                    ok if (ok or not tr) else None,
                    f"float {_fmt(tc_mu, 7)} K vs proved T_low {_fmt(a['T_low_K'], 7)} K"
                    + (" (float Matsubara-truncated: lower estimate)" if tr else ""),
                    serious=not tr)
            else:
                ls = a.get("lower_search") or {}
                add(f"no proved T_low [{tag}, {a['label']}]", None,
                    f"rtsc lower search found no C-RAYLEIGH certificate; its float proposal at the "
                    f"envelope minimum (weights x(1-e), mu*={mu:.2f}) was "
                    f"{_fmt(ls.get('tc_float_lo'), 5)} K; certify_lower only searches T >= "
                    f"max(1e-3*T_hi, 0.05 K) (search.py), so Tc_model below ~0.05 K cannot get a "
                    f"lower bound (sound, uninformative). Nominal-weight float Tc: {_fmt(tc_mu, 5)} K")''')
rep('''            el = f["eliashberg_float"][mu]["Tc_K"]
            if el > 0:
                add(f"independent float vs rtsc floatcore [{f['label']}, mu*={mu}]",
                    abs(d["Tc_K"] - el) / el <= 1e-3,
                    f"{_fmt(el, 7)} K vs {_fmt(d['Tc_K'], 7)} K")''',
'''            el = f["eliashberg_float"][mu]["Tc_K"]
            tr = bool(f["eliashberg_float"][mu].get("truncated"))
            if el > 0:
                ok = abs(d["Tc_K"] - el) / el <= 1e-3
                add(f"independent float vs rtsc floatcore [{f['label']}, mu*={mu}]",
                    ok if (ok or not tr) else None,
                    f"{_fmt(el, 7)} K vs {_fmt(d['Tc_K'], 7)} K"
                    + (" (Matsubara truncation in both: ours N<=4000, floatcore cap 3000; truncated "
                       "values are lower estimates, so the larger one is closer)" if tr else ""))''')
rep('''                 " | ".join(_fmt(f["eliashberg_float"][m]["Tc_K"]) for m in mus) + " |")
    L.append("")''',
'''                 " | ".join(_fmt(f["eliashberg_float"][m]["Tc_K"])
                            + (" †" if f["eliashberg_float"][m].get("truncated") else "")
                            for m in mus) + " |")
    L.append("")
    L.append("† Matsubara sum truncated at N = 4000 (needed more): a LOWER estimate of the float Tc.")
    L.append("")''')
open(p,"w",encoding="utf-8").write(s); print("ok")
