"""Hash-then-delete pruning of regenerable DFT intermediates (policy: reports/EVIDENCE_RETENTION.md).

Every file is SHA-256 hashed and recorded BEFORE it is deleted. The record is appended to
<root>/PRUNED_MANIFEST.tsv and to a central copy on the Windows side, outside the WSL disk image.
Nothing is deleted unless its manifest line has been written and flushed to both places.

    python prune.py wfc  ROOT --reason TEXT [--dry-run]           # only QE wavefunction files
    python prune.py all  ROOT --reason TEXT [--dry-run]           # every file under ROOT (whole dir)
    python prune.py qdone PH0DIR DYNDIR --reason TEXT [--dry-run] # wfc of finished q dirs of a live ph.x
    python prune.py qmark PH0DIR MARKDIR --mark "al.dyn{n}.elph.{n}" --reason TEXT  # same, progress from markers
    python prune.py inventory ROOT                                # hash every KEPT file; deletes nothing
    python prune.py verify MANIFEST [BASE]                        # re-hash listed files ("~" = BASE); exit 1 on any mismatch

Wavefunction files: QE per-process '<prefix>.wfc<N>' and 'wfc<N>.dat'/'wfcup|dw<N>.dat'.
Never touched: dvscf, dyn, elph*, phsave, *.xml, charge density, anything inside <prefix>.save/, anything else.
"""
from __future__ import annotations

import argparse
import hashlib
import os
import re
import sys
import time
from pathlib import Path

CENTRAL = Path("/mnt/g/ClaudeWorkspace/RTSC_Project/reports/data/pruned_manifests")
WFC_RE = re.compile(r"(^.+\.wfc\d+$)|(^wfc(up|dw)?\d+\.(dat|hdf5)$)")
HEADER = "utc_time\tmachine\tpath\tbytes\tsha256\tmtime_utc\treason\n"
# Manifests are published: paths are recorded relative to the home directory and the machine by a
# neutral label (RTSC_MACHINE_LABEL), never by hostname or username.
MACHINE = os.environ.get("RTSC_MACHINE_LABEL", "run-machine-1")


def disp(p: Path) -> str:
    s, home = str(p.resolve()), str(Path.home())
    return "~" + s[len(home):] if s == home or s.startswith(home + "/") else s


def tag_of(root: Path) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", disp(root).lstrip("~/"))


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 24), b""):
            h.update(chunk)
    return h.hexdigest()


def manifests(root: Path) -> list[Path]:
    CENTRAL.mkdir(parents=True, exist_ok=True)
    tag = tag_of(root)
    return [root / "PRUNED_MANIFEST.tsv", CENTRAL / f"{tag}.tsv"]


def record(paths: list[Path], line: str) -> None:
    for m in paths:
        new = not m.exists()
        with open(m, "a") as f:
            if new:
                f.write(HEADER)
            f.write(line)
            f.flush()
            os.fsync(f.fileno())


def prune(files: list[Path], root: Path, reason: str, dry: bool) -> int:
    mans = manifests(root)
    total = 0
    for p in sorted(files):
        if p.is_symlink() or not p.is_file() or p.name == "PRUNED_MANIFEST.tsv":
            continue
        st = p.stat()
        digest = sha256(p)
        line = "\t".join([time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), MACHINE,
                          disp(p), str(st.st_size), digest,
                          time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(st.st_mtime)), reason]) + "\n"
        if dry:
            print("DRY", st.st_size, p)
            continue
        record(mans, line)
        p.unlink()
        total += st.st_size
    print(f"pruned {total/1e9:.2f} GB under {disp(root)}; manifests: {', '.join(map(disp, mans))}")
    return total


def _dyn_complete(dyn: Path) -> bool:
    return dyn.is_file() and dyn.stat().st_size > 0 and bool(dyn.read_text(errors="replace").strip())


def inventory(root: Path) -> None:
    """Hash every file kept under root into <root>/RUN_MANIFEST.tsv (+ central copy); deletes nothing."""
    tag = tag_of(root)
    CENTRAL.mkdir(parents=True, exist_ok=True)
    rows = [HEADER.replace("\treason", "\tstatus")]
    stamp = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    for p in sorted(root.rglob("*")):
        if p.is_file() and not p.is_symlink() and p.name not in ("RUN_MANIFEST.tsv", "PRUNED_MANIFEST.tsv"):
            st = p.stat()
            rows.append("\t".join([stamp, MACHINE, disp(p), str(st.st_size), sha256(p),
                                   time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(st.st_mtime)), "kept"]) + "\n")
    for m in (root / "RUN_MANIFEST.tsv", CENTRAL / f"{tag}.RUN.tsv"):
        m.write_text("".join(rows))
    print(f"inventoried {len(rows) - 1} files under {disp(root)}")


def verify(manifest: Path, base: Path) -> int:
    """Re-hash the files a manifest lists ('~' resolved to base) and compare. Returns #non-matching."""
    counts = {"MATCH": 0, "MISMATCH": 0, "MISSING": 0}
    for row in manifest.read_text().splitlines()[1:]:
        f = row.split("\t")
        p = base / f[2][2:] if f[2].startswith("~/") else Path(f[2])
        st = "MISSING" if not p.is_file() else ("MATCH" if sha256(p) == f[4] else "MISMATCH")
        counts[st] += 1
        if st != "MATCH":
            print(st, f[2])
    print(" ".join(f"{k}={v}" for k, v in counts.items()))
    return counts["MISMATCH"] + counts["MISSING"]


def finished_q_dirs(ph0: Path, dyndir: Path) -> list[Path]:
    """q dirs whose dyn file is complete AND whose successor's dyn file exists (ph.x has moved on).

    ph.x creates every _ph0/<prefix>.q_N directory up front, so directory presence says nothing
    about progress; the dyn files do. The last q is never pruned here (use 'wfc' once the run ends).
    """
    out = []
    for d in ph0.iterdir():
        m = re.search(r"\.q_(\d+)$", d.name)
        if d.is_dir() and m:
            n = int(m.group(1))
            if _dyn_complete(dyndir / f"dyn{n}") and (dyndir / f"dyn{n + 1}").exists():
                out.append(d)
    return out


def marked_q_dirs(ph0: Path, markdir: Path, mark: str) -> list[Path]:
    """q dirs N for which markdir/mark(N) AND markdir/mark(N+1) exist (q N done, ph.x moved on).

    For ph.x runs whose dyn files pre-exist (e.g. electron_phonon='lambda_tetra' re-reading a
    finished DFPT), progress is shown by per-q output such as 'al.dyn{n}.elph.{n}'.
    """
    out = []
    for d in ph0.iterdir():
        m = re.search(r"\.q_(\d+)$", d.name)
        if d.is_dir() and m:
            n = int(m.group(1))
            if (markdir / mark.format(n=n)).exists() and (markdir / mark.format(n=n + 1)).exists():
                out.append(d)
    return out


def wfc_files(root: Path) -> list[Path]:
    """Wavefunction files under root, excluding <prefix>.save/ (scf restart data is kept)."""
    return [p for p in root.rglob("*") if WFC_RE.match(p.name)
            and not any(part.endswith(".save") for part in p.relative_to(root).parts[:-1])]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["wfc", "all", "qdone", "qmark", "inventory", "verify"])
    ap.add_argument("root")
    ap.add_argument("dyndir", nargs="?", help="qdone: dir with dynN files; qmark: dir with marker files")
    ap.add_argument("--mark", help="qmark: marker filename with {n}, e.g. al.dyn{n}.elph.{n}")
    ap.add_argument("--reason", default="")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    if a.mode == "verify":  # root = manifest file, dyndir = directory that "~" stands for
        sys.exit(1 if verify(Path(a.root), Path(a.dyndir or Path.home())) else 0)
    root = Path(a.root)
    if not root.is_dir():
        sys.exit(f"not a directory: {root}")
    if a.mode == "inventory":
        return inventory(root)
    if not a.reason:
        sys.exit("--reason is required for pruning modes")
    if a.mode == "all":
        files = [p for p in root.rglob("*")]
    elif a.mode == "wfc":
        files = wfc_files(root)
    elif a.mode == "qdone":
        if not a.dyndir:
            sys.exit("qdone needs DYNDIR")
        files = [p for d in finished_q_dirs(root, Path(a.dyndir)) for p in wfc_files(d)]
    else:
        if not (a.dyndir and a.mark):
            sys.exit("qmark needs MARKDIR and --mark")
        files = [p for d in marked_q_dirs(root, Path(a.dyndir), a.mark) for p in wfc_files(d)]
    prune(files, root, a.reason, a.dry_run)
    if a.mode == "all" and not a.dry_run:
        for d in sorted((p for p in root.rglob("*") if p.is_dir()), key=lambda p: -len(p.parts)):
            try:
                d.rmdir()
            except OSError:
                pass


if __name__ == "__main__":
    main()
