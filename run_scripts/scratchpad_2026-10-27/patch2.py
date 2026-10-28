p="/mnt/g/ClaudeWorkspace/RTSC_Project/bench_tools/analyze_benchmark.py"; s=open(p,encoding="utf-8").read()
old='''    for f in fl:
        for mu in ("0.10", "0.15"):
            ad = f["allen_dynes"].get(mu, {}).get("Tc_K")
            el = f["eliashberg_float"].get(mu, {}).get("Tc_K")
            if ad and el:
                rel = (el - ad) / ad
                add(f"AD vs Eliashberg float [{f['label']}, mu*={mu}]", abs(rel) <= 0.25,
                    f"AD={_fmt(ad)} K, Eliashberg={_fmt(el)} K, rel diff {rel:+.1%} "
                    "(different mu* conventions; ~10-25% expected)")'''
new='''    # mu* = 0 is the clean comparison (no cutoff-convention ambiguity): expect ~10-20 %.
    # For mu* > 0 the two use different mu* conventions (AD's fitted mu* vs our block cutoff at
    # omega_c = factor x max omega), and in weak coupling Tc is exponentially sensitive to that,
    # so the mu* > 0 differences are reported as information only.
    for f in fl:
        for mu in [m for m in f["allen_dynes"]]:
            ad = f["allen_dynes"].get(mu, {}).get("Tc_K")
            el = f["eliashberg_float"].get(mu, {}).get("Tc_K")
            if ad and el:
                rel = (el - ad) / ad
                clean = mu == "0.00"
                add(f"AD vs Eliashberg float [{f['label']}, mu*={mu}]",
                    (abs(rel) <= 0.20) if clean else None,
                    f"AD={_fmt(ad)} K, Eliashberg={_fmt(el)} K, rel diff {rel:+.1%}"
                    + ("" if clean else " (info: different mu* conventions)"))'''
assert old in s
s=s.replace(old,new)
open(p,"w",encoding="utf-8").write(s); print("ok")
