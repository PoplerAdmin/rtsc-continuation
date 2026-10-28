import re
p='/mnt/g/ClaudeWorkspace/RTSC_Project/bench_tools/ephproto/bandcheck.py'
s=open(p).read()
s=s.replace('''def write_kpoint_files(outdir: str | Path, n=(5, 5, 5), coarse=(8, 8, 8)) -> dict:
    d = Path(outdir)
    d.mkdir(parents=True, exist_ok=True)
    keep, drop = validation_kpoints(n, coarse)''','''def random_kpoints(nrand: int, seed: int, coarse: Sequence[int] = (8, 8, 8)):
    """SUPPLEMENTARY (post hoc, not the pre-registered rule): nrand uniform crystal points in
    [0,1)^3 from numpy default_rng(seed), rounded to 1e-6 (filkf/band.eig precision).  Points within
    1e-4 of a coarse-grid point (in any component combination) are dropped."""
    rng = np.random.default_rng(seed)
    pts = np.round(rng.random((nrand, 3)), 6)
    c = np.asarray(coarse, float)
    near = np.all(np.abs(pts * c - np.round(pts * c)) < 1e-4 * c, axis=1)
    keep = [tuple(Fraction(repr(float(x))) for x in p) for p in pts[~near]]
    drop = [tuple(Fraction(repr(float(x))) for x in p) for p in pts[near]]
    return keep, drop


def write_kpoint_files(outdir: str | Path, n=(5, 5, 5), coarse=(8, 8, 8),
                       kset: str = "mp5", nrand: int = 1000, seed: int = 20261027) -> dict:
    d = Path(outdir)
    d.mkdir(parents=True, exist_ok=True)
    if kset == "mp5":
        keep, drop = validation_kpoints(n, coarse)
    elif kset == "random":
        keep, drop = random_kpoints(nrand, seed, coarse)
    else:
        raise BandCheckError(f"unknown kset {kset!r}")''')
s=s.replace('''    return {"n": list(n), "coarse": list(coarse), "N_k": len(keep),''','''    return {"kset": kset, "n": list(n) if kset == "mp5" else None,
            "nrand": nrand if kset == "random" else None, "seed": seed if kset == "random" else None,
            "coarse": list(coarse), "N_k": len(keep),''')
s=s.replace('''              ef_grids: Sequence[int] = (16,)) -> dict:
    root = Path(m["epw_dir"]) / "g1"
    kmeta = write_kpoint_files(root, coarse=m["nk"])''','''              ef_grids: Sequence[int] = (16,), kset: str = "mp5", nrand: int = 1000,
              seed: int = 20261027) -> dict:
    root = Path(m["epw_dir"]) / ("g1" if kset == "mp5" else f"g1_{kset}{nrand}_s{seed}")
    kmeta = write_kpoint_files(root, coarse=m["nk"], kset=kset, nrand=nrand, seed=seed)''')
# E_F runs do not depend on k set: reuse the mp5 tree's runs if present
s=s.replace('''    d = root / label / f"epw_ef_K{K[0]}"''','''    d = Path(m["epw_dir"]) / "g1" / label / f"epw_ef_K{K[0]}"      # independent of the k set''')
s=s.replace('''    r.add_argument("--ef-grids", type=int, nargs="+", default=[16])''','''    r.add_argument("--ef-grids", type=int, nargs="+", default=[16])
    r.add_argument("--kset", choices=["mp5", "random"], default="mp5")
    r.add_argument("--nrand", type=int, default=1000)
    r.add_argument("--seed", type=int, default=20261027)''')
s=s.replace('''        res = stage_run(m, Path(args.wannier), Path(args.epmatwp), args.label, args.np,
                        tuple(args.ef_grids))''','''        res = stage_run(m, Path(args.wannier), Path(args.epmatwp), args.label, args.np,
                        tuple(args.ef_grids), args.kset, args.nrand, args.seed)''')
# record distinct energies in each window (degeneracy / symmetry diagnostic) and inconclusive flag
s=s.replace('''               "pass_dE": bool(a.size > 0 and mx <= w["tol_eV"]),''','''               "n_k_with_states": int(sum(1 for ik in range(len(D)) if
                                      np.any((D[ik] >= lo) & (D[ik] <= hi)) or
                                      np.any((Wn[ik] >= lo) & (Wn[ik] <= hi)))),
               "inconclusive_no_states": bool(a.size == 0),
               "pass_dE": bool(a.size > 0 and mx <= w["tol_eV"]),''')
open(p,'w').write(s)
print(s.count('kset'))
