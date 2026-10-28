"""Float diagnostics (NOT proofs) for alpha2f.x tetra outputs vs QE smearing a2F.dos."""
import sys, re, json, math, glob
from pathlib import Path
import numpy as np
from rtsc.dfpt.units import to_kelvin_factor
from rtsc.dfpt.parsers import parse_qe_a2f_dos
RY_K = to_kelvin_factor("Ry")

def moments(w, a, dw):
    # bins rule: atom at w_i with mass a_i*dw_i  (w in Ry)
    m = a * dw
    lam = 2 * np.sum(m / w)
    wlog = math.exp(2 * np.sum(m / w * np.log(w)) / lam) if lam > 0 else float('nan')
    A = np.sum(m)
    return lam, wlog * RY_K, A * RY_K

def tetra(d):
    d = Path(d)
    t = np.loadtxt(d / "al.a2F.dat", comments="#")
    w, a = t[:, 0], t[:, 1]
    dE = w[0]
    assert np.allclose(np.diff(w), dE, rtol=1e-9, atol=1e-14)
    out = {"dir": str(d), "nfreq": len(w), "dE_Ry": dE, "wmax_K": w[-1] * RY_K}
    out["lam_bins"], out["wlog_bins_K"], out["A_bins_K"] = moments(w, a, np.full_like(w, dE))
    # trapezoid (with implicit a(0)=0 not included -- rtsc rule on the given points)
    h = np.diff(w); wt = np.zeros_like(w); wt[:-1] += h/2; wt[1:] += h/2
    out["lam_trap"], out["wlog_trap_K"], out["A_trap_K"] = moments(w, a, wt)
    neg = a < 0
    out["n_neg"] = int(neg.sum()); out["min_a2F"] = float(a.min())
    out["neg_mass_K"] = float(-np.sum(a[neg]) * dE * RY_K)
    out["neg_rel_A"] = out["neg_mass_K"] / out["A_bins_K"]
    out["neg_freqs_K"] = [float(x) for x in (w[neg] * RY_K)[:20]]
    s = (d / "a2F.out").read_text()
    f = lambda pat: [float(x) for x in re.findall(pat, s)]
    lams = f(r"lambda :\s+(\S+)")
    out["qe_lambda_from_a2F"], out["qe_lambda_direct"] = lams[0], lams[1]
    out["qe_wln_from_a2F_K"] = f(r"omega_ln \[Ry\] :\s+(\S+)")[0] * RY_K
    out["qe_wln_direct_K"] = f(r"omega_ln \[K\] :\s+(\S+)")[0]
    out["qe_lambda_from_a2F_str"] = re.findall(r"lambda :\s+(\S+)", s)[0]
    e = (d / "elph.out").read_text()
    out["dosEf_states_spin_Ry"] = sorted(set(float(x) for x in re.findall(r"DOS =\s*(\S+) states/spin", e)))
    out["nq_done"] = len(re.findall(r"Tetrahedron method", e))
    out["nk_per_q"] = [int(x) for x in re.findall(r"number of k points=\s*(\d+) \(tetra", e)]
    return out, (w, a)

def smear(path):
    tab = parse_qe_a2f_dos(path)
    w, a = np.asarray(tab.E_Ry), np.asarray(tab.a2f_total)
    lam, wl, A = moments(np.where(w > 0, w, np.nan)[w > 0], a[w > 0], np.full((w > 0).sum(), tab.delta_Ry))
    return {"file": path, "lam": lam, "wlog_K": wl, "A_K": A, "qe_lam": tab.lambda_reported,
            "n_neg": int((a[w > 0] < 0).sum()), "min": float(a.min())}, (w, a, tab.delta_Ry)

def rebin(w, a, dw, edges):
    # histogram of the bins-rule mass onto common edges (K)
    m, _ = np.histogram(w * RY_K, bins=edges, weights=a * dw * RY_K)
    return m

if __name__ == "__main__":
    res = {"tetra": [], "smear": []}
    edges = np.arange(0, 520, 10.0)
    hist = {}
    for d in sys.argv[1:]:
        try:
            o, (w, a) = tetra(d)
        except Exception as ex:
            print("skip", d, ex); continue
        res["tetra"].append(o)
        hist[d] = rebin(w, a, np.full_like(w, w[0]), edges)
    S = "/home/user/rtsc_work/bench/Al/run/Al-fcc_0GPa-d2e010742b/dfpt/a2F.dos"
    for i in range(1, 11):
        o, (w, a, dl) = smear(S + str(i)); res["smear"].append(o)
        hist["smear%d" % i] = rebin(w[w > 0], a[w > 0], dl, edges)
    # shape: L1 distance of normalised 10 K histograms vs smearing sigma=0.02 (index 4)
    ref = hist["smear4"] / hist["smear4"].sum()
    res["shape_L1_vs_smear4"] = {k: float(np.abs(v / v.sum() - ref).sum()) for k, v in hist.items()}
    res["hist_edges_K"] = edges.tolist()
    res["hist"] = {k: v.tolist() for k, v in hist.items()}
    print(json.dumps(res, indent=1))
