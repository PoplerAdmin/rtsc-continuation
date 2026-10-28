"""Build the SYNTHETIC oracle test suite and compute oracle Tc values.

Every spectrum here is SYNTHETIC (invented for testing); none is data for a real compound.

    python make_suite.py --list
    python make_suite.py --case NAME        # computes and writes cases/NAME.json
    python make_suite.py --index            # collects cases/*.json into suite_index.json

Per case: oracle Tc at mu* in {0, 0.1, 0.15}, omega_c = 10 * max omega, block cutoff
N_c(T) = #{n >= 0 : (2n+1) pi T <= omega_c}.  Envelope cases (rel_envelope > 0) also get Tc at
the all-low corner, all-high corner and 4 random corners of the mass box.
"""
from __future__ import annotations

import argparse
import glob
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..")))

from oracle import A0, Model, allen_dynes, find_tc  # noqa: E402
from rtsc.spectrum import Spectrum  # noqa: E402  (shared data type, ARCH section 2)

MUS = (0.0, 0.1, 0.15)
CASEDIR = os.path.join(HERE, "cases")
LAM_OPT = 1.1436            # near-extremal Einstein lambda (PROOF.md, float exploration)


# ------------------------------------------------------------------------------------------
# spectrum builders (all return (omega list, weight list), kelvin)
# ------------------------------------------------------------------------------------------
def einstein(lam, om):
    return [om], [lam * om / 2.0]


def from_density(edges_or_grid, dens, lam_target):
    """Midpoint-rule atoms from a density on a grid, normalised to lambda = lam_target."""
    x = np.asarray(edges_or_grid, float)
    d = np.asarray(dens, float)
    dx = np.gradient(x)
    w = d * dx
    lam = 2 * np.sum(w / x)
    w = w * lam_target / lam
    keep = w > 0
    return list(x[keep]), list(w[keep])


def debye(lam, wD, n=200):
    x = (np.arange(n) + 0.5) * wD / n
    return from_density(x, x ** 2, lam)


def join(*parts):
    o, w = [], []
    for a, b in parts:
        o += list(a)
        w += list(b)
    return o, w


def gauss_band(lam, c, s, lo, hi, n):
    x = np.linspace(lo, hi, n)
    return from_density(x, np.exp(-0.5 * ((x - c) / s) ** 2), lam)


def hydride(lam_ac, lam_h, n_ac=120, n_h=180, h_lo=1500.0, h_hi=2500.0, h_c=2000.0, h_s=250.0):
    x = (np.arange(n_ac) + 0.5) * 600.0 / n_ac
    ac = from_density(x, x ** 2 * np.exp(-((x / 450.0) ** 4)), lam_ac)
    return join(ac, gauss_band(lam_h, h_c, h_s, h_lo, h_hi, n_h))


def scale(spec, s):
    o, w = spec
    return [s * a for a in o], [s * b for b in w]


def dfpt_like(n, wmax, lam, seed, peaks):
    rng = np.random.default_rng(seed)
    x = np.linspace(wmax / n, wmax, n)
    d = np.zeros(n)
    for (c, s, h) in peaks:
        d += h * np.exp(-0.5 * ((x - c) / s) ** 2)
    d += 0.02 * (x / wmax) ** 2                                     # Debye-like floor
    d *= np.clip(1 + 0.25 * rng.standard_normal(n), 0.0, None)      # smearing-like noise
    d[x > 0.98 * wmax] = 0.0
    return from_density(x, d, lam)


# ------------------------------------------------------------------------------------------
def build_cases():
    C = {}

    def add(name, spec, desc, e=0.0, tags=()):
        C[name] = dict(spec=spec, desc=desc, e=e, tags=list(tags))

    for lam in (0.2, 0.3, 0.5, 0.75, 1.0, LAM_OPT, 1.5, 2.0, 3.0, 4.0, 5.0):
        add(f"einstein_w1000_lam{lam:g}", einstein(lam, 1000.0),
            f"single Einstein mode at 1000 K, lambda={lam:g}", tags=["einstein"])
    add("twopeak_a", join(einstein(0.5, 200.0), einstein(0.7, 1200.0)),
        "two Einstein peaks 200 K (lam .5) + 1200 K (lam .7)", tags=["twopeak"])
    add("twopeak_b", join(einstein(1.0, 80.0), einstein(0.5, 600.0)),
        "two Einstein peaks 80 K (lam 1) + 600 K (lam .5)", tags=["twopeak"])
    add("twopeak_c", join(einstein(0.4, 400.0), einstein(1.2, 1600.0)),
        "two Einstein peaks 400 K (lam .4) + 1600 K (lam 1.2)", tags=["twopeak"])
    add("debye_wD400_lam0.6", debye(0.6, 400.0), "Debye omega^2 up to 400 K, lam .6, 200 atoms",
        tags=["debye"])
    add("debye_wD1000_lam1.2", debye(1.2, 1000.0), "Debye omega^2 up to 1000 K, lam 1.2",
        tags=["debye"])
    add("debye_wD2000_lam2", debye(2.0, 2000.0), "Debye omega^2 up to 2000 K, lam 2",
        tags=["debye"])
    add("hydride_lam1.6", hydride(0.6, 1.0), "acoustic band <600 K (lam .6) + H band 1500-2500 K "
        "(lam 1.0)", tags=["hydride"])
    add("hydride_lam2.2", hydride(0.6, 1.6), "acoustic (lam .6) + H band 1500-2500 K (lam 1.6)",
        tags=["hydride"])
    add("hydride_lam2.8", hydride(0.6, 2.2), "acoustic (lam .6) + H band 1500-2500 K (lam 2.2)",
        tags=["hydride"])
    add("soft_atom2K", join(einstein(1.0, 800.0), einstein(0.5, 2.0)),
        "800 K mode (lam 1) + soft atom at 2 K (lam .5)", tags=["soft"])
    add("soft_atom0.5K", join(einstein(0.8, 500.0), einstein(2.0, 0.5)),
        "500 K mode (lam .8) + soft atom at 0.5 K (lam 2)", tags=["soft"])
    add("soft_band1-20K", join(gauss_band(1.5, 8.0, 5.0, 1.0, 20.0, 40), einstein(0.5, 400.0)),
        "soft band 1-20 K (lam 1.5) + 400 K mode (lam .5)", tags=["soft"])
    xw = np.geomspace(10.0, 3000.0, 300)
    add("wide_flat_10-3000", from_density(xw, np.ones_like(xw), 1.5),
        "flat a2F on 10..3000 K (log grid, 300 atoms), lam 1.5", tags=["wide"])
    add("wide_twohump_10-3000", from_density(xw, np.exp(-0.5 * ((xw - 60) / 30) ** 2)
                                             + 0.5 * np.exp(-0.5 * ((xw - 2200) / 300) ** 2), 3.0),
        "two humps (60 K, 2200 K) on 10..3000 K log grid, lam 3", tags=["wide"])
    add("dfpt_3000atoms", dfpt_like(3000, 2500.0, 1.8, 11, [(250, 80, 1.0), (700, 120, .6),
                                                           (1900, 250, .5)]),
        "DFPT-like noisy grid, 3000 atoms to 2500 K, lam 1.8", tags=["dfpt", "many_atoms"])
    add("dfpt_6000atoms", dfpt_like(6000, 1200.0, 0.9, 12, [(150, 50, 1.0), (450, 80, .7),
                                                           (900, 100, .3)]),
        "DFPT-like noisy grid, 6000 atoms to 1200 K, lam 0.9", tags=["dfpt", "many_atoms"])
    # adversarial near the universal bound (Tc/A ~ 0.2308679 at mu* = 0)
    for A in (1299.3, 1299.6, 1302.0):
        add(f"adv_opt_einstein_A{A:g}", einstein(LAM_OPT, 2 * A / LAM_OPT),
            f"near-extremal Einstein (lam {LAM_OPT}) with A = {A} K: Tc(mu*=0) ~ A/4.3314815",
            tags=["adversarial", "near_bound"])
    A = 1300.0
    Om = 2 * A / LAM_OPT
    add("adv_opt_twopeak_A1300", ([Om * 0.98, Om * 1.02], [A / 2, A / 2]),
        "near-extremal split pair (Omega*(1 -/+ 0.02)), A = 1300 K", tags=["adversarial",
                                                                        "near_bound"])
    add("adv_opt_threepeak_A1300", ([Om * 0.9, Om, Om * 1.1], [A * .2, A * .6, A * .2]),
        "near-extremal triple, A = 1300 K", tags=["adversarial", "near_bound"])
    # T*-straddling cases (scaled by exact similarity so that Tc is a chosen distance from 300)
    C["adv_hydride_mu0.15_300p3"] = dict(spec=("SCALE_TO", "hydride_lam2.2", 0.15, 300.3),
                                         desc="hydride-like scaled: Tc(mu*=.15) = 300.3 K",
                                         e=0.0, tags=["adversarial", "T*_straddle"])
    C["adv_hydride_mu0.15_299p7"] = dict(spec=("SCALE_TO", "hydride_lam2.2", 0.15, 299.7),
                                         desc="hydride-like scaled: Tc(mu*=.15) = 299.7 K",
                                         e=0.0, tags=["adversarial", "T*_straddle"])
    C["adv_einstein2_mu0_300p5"] = dict(spec=("SCALE_TO", "einstein_w1000_lam2", 0.0, 300.5),
                                        desc="Einstein lam 2 scaled: Tc(mu*=0) = 300.5 K",
                                        e=0.0, tags=["adversarial", "T*_straddle"])
    C["adv_einstein2_mu0_299p5"] = dict(spec=("SCALE_TO", "einstein_w1000_lam2", 0.0, 299.5),
                                        desc="Einstein lam 2 scaled: Tc(mu*=0) = 299.5 K",
                                        e=0.0, tags=["adversarial", "T*_straddle"])
    # weak coupling (Tc ~ 1-5 K)
    add("weak_einstein_w100_lam0.45", einstein(0.45, 100.0), "Einstein 100 K, lam .45",
        tags=["weak"])
    add("weak_debye_wD300_lam0.45", debye(0.45, 300.0), "Debye 300 K, lam .45", tags=["weak"])
    add("weak_debye_wD150_lam0.6", debye(0.6, 150.0), "Debye 150 K, lam .6", tags=["weak"])
    # strong coupling lambda ~ 4
    add("strong_einstein_w500_lam4", einstein(4.0, 500.0), "Einstein 500 K, lam 4",
        tags=["strong"])
    add("strong_debye_wD1500_lam4", debye(4.0, 1500.0), "Debye 1500 K, lam 4", tags=["strong"])
    add("strong_twopeak_lam4", join(einstein(3.0, 150.0), einstein(1.0, 900.0)),
        "150 K (lam 3) + 900 K (lam 1)", tags=["strong"])
    # envelope cases
    add("env_hydride_lam2.2_e0.05", hydride(0.6, 1.6), "hydride-like lam 2.2, rel_envelope .05",
        e=0.05, tags=["envelope"])
    add("env_einstein_lam1_e0.1", einstein(1.0, 1000.0), "Einstein 1000 K lam 1, rel_envelope .1",
        e=0.1, tags=["envelope"])
    return C


def to_spectrum(name, spec, desc, e):
    o, w = spec
    return Spectrum(tuple(float(a) for a in o), tuple(float(b) for b in w), e,
                    label=f"SYNTHETIC: {name} -- {desc}",
                    provenance={"synthetic": True, "note": "SYNTHETIC test spectrum, not data "
                                "for any real compound", "generator": "tests/oracle/make_suite.py"})


def run_tc(o, w, mu, omega_c, A_ceiling=None):
    m = Model.from_spectrum(o, w, mu, omega_c)
    t = time.time()
    r = find_tc(m, A_ceiling=A_ceiling, max_jumps=int(os.environ.get("ORACLE_MAX_JUMPS", "700")))
    d = r.to_dict()
    d["seconds"] = round(time.time() - t, 2)
    return d


def resolve(C, name):
    spec = C[name]["spec"]
    if isinstance(spec, tuple) and spec and spec[0] == "SCALE_TO":
        _, base, mu, target = spec
        o, w = C[base]["spec"]
        wc = 10 * max(o)
        r = run_tc(o, w, mu, wc)
        s = target / r["Tc"]          # exact similarity: (o, w, omega_c) -> s*(...), Tc -> s*Tc
        return scale((o, w), s), {"scaled_from": base, "scale": s, "target_mu": mu,
                                  "target_Tc": target}
    return spec, {}


def compute_case(name):
    C = build_cases()
    c = C[name]
    spec, extra = resolve(C, name)
    S = to_spectrum(name, spec, c["desc"], c["e"])
    S.provenance.update(extra)
    o, w = list(S.omega_K), list(S.weight_K)
    omega_c = 10.0 * max(o)
    out = {"name": name, "synthetic": True, "description": c["desc"], "tags": c["tags"],
           "spectrum": S.to_dict(), "omega_c_K": omega_c, "omega_c_rule": "10*max(omega_K)",
           "N_c_rule": "#{n>=0 : (2n+1)*pi*T <= omega_c}",
           "A_K": S.A(), "lambda": S.lambda_(), "omega_log_K": S.omega_log(),
           "ceiling_K": S.A() / A0, "n_atoms": len(o), "oracle": {}, "allen_dynes": {}}
    for mu in MUS:
        out["oracle"][str(mu)] = run_tc(o, w, mu, omega_c)
        out["allen_dynes"][str(mu)] = allen_dynes(S.lambda_(), S.omega_log(), S.omega2(), mu)
    if c["e"] > 0:
        e = c["e"]
        rng = np.random.default_rng(1234)
        corners = {"all_low": [1 - e] * len(w), "all_high": [1 + e] * len(w)}
        for i in range(4):
            corners[f"random{i}"] = list(1 + e * rng.choice([-1.0, 1.0], size=len(w)))
        env = {}
        for mu in MUS:
            res = {}
            for cn, f in corners.items():
                ww = [a * b for a, b in zip(w, f)]
                res[cn] = run_tc(o, ww, mu, omega_c, A_ceiling=sum(ww))
            Tcs = [r["Tc"] for r in res.values() if r["Tc"] is not None]
            env[str(mu)] = {"corners": res,
                            "Tc_min": min(Tcs) if Tcs else None,
                            "Tc_max": max(Tcs) if Tcs else None,
                            "err_max": max(r["err"] for r in res.values() if r["Tc"])}
        out["envelope"] = env
    os.makedirs(CASEDIR, exist_ok=True)
    with open(os.path.join(CASEDIR, name + ".json"), "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    return out


def index():
    rows = []
    for f in sorted(glob.glob(os.path.join(CASEDIR, "*.json"))):
        d = json.load(open(f))
        row = {"name": d["name"], "A_K": d["A_K"], "lambda": d["lambda"],
               "n_atoms": d["n_atoms"], "ceiling_K": d["ceiling_K"]}
        for mu, r in d["oracle"].items():
            row[f"Tc_mu{mu}"] = r["Tc"]
            row[f"err_mu{mu}"] = r["err"]
            row[f"status_mu{mu}"] = r["status"]
            row[f"floor_mu{mu}"] = r["T_floor"]
            row[f"reentrant_mu{mu}"] = r["reentrant"]
            row[f"AD_mu{mu}"] = d["allen_dynes"][mu]
        rows.append(row)
    json.dump(rows, open(os.path.join(HERE, "suite_index.json"), "w"), indent=1)
    return rows


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--case")
    ap.add_argument("--index", action="store_true")
    a = ap.parse_args()
    if a.list:
        print("\n".join(build_cases()))
    if a.case:
        t = time.time()
        d = compute_case(a.case)
        s = " ".join(f"mu{mu}:Tc={r['Tc'] if r['Tc'] is None else round(r['Tc'], 6)}"
                     f"({r['status']},{r['seconds']}s)" for mu, r in d["oracle"].items())
        print(f"{a.case}: A={d['A_K']:.4g} lam={d['lambda']:.4g} {s} total {time.time()-t:.1f}s",
              flush=True)
    if a.index:
        for r in index():
            print(r["name"], {k: r[k] for k in r if k.startswith("Tc_")})
