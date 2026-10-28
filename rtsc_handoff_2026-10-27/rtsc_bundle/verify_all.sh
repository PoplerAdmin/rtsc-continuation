#!/usr/bin/env bash
# verify_all.sh: check the rtsc bundle end to end.
#
#   1  python  : package installs (pip, no network, into a throw-away venv) and imports
#   2  pytest  : the rtsc test suite (real-QE test runs only with RTSC_RUN_QE=1)
#   3  demo    : `rtsc demo` on SYNTHETIC spectra; every input has exactly one status (S4)
#   4  recheck : `rtsc recheck` re-verifies every classified demo verdict (REJECTED / CERTIFIED /
#                UNDECIDED), plus `rtsc certify` on a bundled QE a2F test file
#   5  lean    : `lake build`; no sorry/admit/axiom/native_decide in the sources; #print axioms
#                for the Check*.lean files and for EVERY declaration (lean_checks/AllAxioms.lean)
#                shows only propext, Classical.choice, Quot.sound
#   6  cert    : (optional, minutes) the constant certificate certificates/wE1/run_all.sh
#
# Usage: ./verify_all.sh [options]
#   --cert            also run the wE1 certificate for a0 = 4.33148145 (about 6 min/core, longer
#                     on a loaded machine)
#   --cert-all        --cert plus the independent route-B certificate (wB: upper bound at
#                     a = 4.3314814, sharpness at a = 4.33148149 as used in PROOF.md Prop. 7)
#                     and the multiband Proposition N check (wMB/cert_neg.py)
#   --workdir DIR     scratch/output directory (default: ./rtsc_verify); the report is
#                     DIR/verify_report.txt
#   --lean-dir DIR    build/check the Lean project in DIR instead of ./lean. DIR's sources must
#                     be identical to ./lean (checked with diff). Use it to reuse an existing
#                     built project.
#   --lean-seed DIR   if ./lean/.lake does not exist, copy DIR/.lake into ./lean first (reuse a
#                     Mathlib build for exactly the same lake-manifest.json; saves ~1 h)
#   --no-cache-get    never try `lake exe cache get` (build Mathlib from source if needed)
#   --skip-python | --skip-pytest | --skip-demo | --skip-lean
#   --pytest-args "ARGS"   extra arguments for pytest (e.g. "-x -k certify")
#   --nice N          niceness for heavy steps (default 10)
#   -h | --help
set -uo pipefail

BUNDLE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK="$PWD/rtsc_verify"
DO_CERT=0; DO_CERT_ALL=0; SKIP_PY=0; SKIP_PYTEST=0; SKIP_DEMO=0; SKIP_LEAN=0
LEAN_DIR=""; LEAN_SEED=""; CACHE_GET=1; PYTEST_ARGS=""; NICE=10

while [ $# -gt 0 ]; do
  case "$1" in
    --cert) DO_CERT=1;;
    --cert-all) DO_CERT=1; DO_CERT_ALL=1;;
    --workdir) WORK="$2"; shift;;
    --lean-dir) LEAN_DIR="$2"; shift;;
    --lean-seed) LEAN_SEED="$2"; shift;;
    --no-cache-get) CACHE_GET=0;;
    --skip-python) SKIP_PY=1;;
    --skip-pytest) SKIP_PYTEST=1;;
    --skip-demo) SKIP_DEMO=1;;
    --skip-lean) SKIP_LEAN=1;;
    --pytest-args) PYTEST_ARGS="$2"; shift;;
    --nice) NICE="$2"; shift;;
    -h|--help) sed -n '2,32p' "${BASH_SOURCE[0]}"; exit 0;;
    *) echo "unknown option: $1 (see --help)" >&2; exit 2;;
  esac
  shift
done

mkdir -p "$WORK"
WORK="$(cd "$WORK" && pwd)"
REPORT="$WORK/verify_report.txt"
LOGDIR="$WORK/logs"; mkdir -p "$LOGDIR"
: > "$REPORT"
export PYTHONDONTWRITEBYTECODE=1          # keep the bundle free of __pycache__
PY="${PYTHON:-python3}"
NICECMD=(nice -n "$NICE")
declare -a RESULTS=()
FAILED=0

say()    { echo "[verify] $*" | tee -a "$REPORT"; }
result() { # result STEP STATUS DETAIL
  RESULTS+=("$(printf '%-9s %-5s %s' "$1" "$2" "$3")")
  say "== $1: $2  $3"
  [ "$2" = "FAIL" ] && FAILED=1
  return 0
}

say "bundle  : $BUNDLE"
say "workdir : $WORK"
say "date    : $(date -u +%Y-%m-%dT%H:%M:%SZ)   host: $(uname -srm)   cores: $(nproc 2>/dev/null || echo ?)"
say "python  : $($PY --version 2>&1)"
[ -f "$BUNDLE/VERSIONS.txt" ] && say "versions: $(head -1 "$BUNDLE/VERSIONS.txt")"

# ---------------------------------------------------------------------------------------------
# 1. Python package: install into a venv that sees the system site-packages. Tries an offline
#    build first, then an isolated build (fetches setuptools from the index). Dependencies are
#    never installed here (--no-deps). Falls back to PYTHONPATH if both fail.
# ---------------------------------------------------------------------------------------------
RTSC=()
if [ "$SKIP_PY" -eq 0 ]; then
  VENV="$WORK/venv"; PKGSRC="$WORK/pkgsrc"
  rm -rf "$VENV" "$PKGSRC"
  # build from a copy: setuptools writes build/ and *.egg-info into the source tree
  how=""
  if mkdir -p "$PKGSRC" && cp -a "$BUNDLE/rtsc/." "$PKGSRC/" && \
     "$PY" -m venv --system-site-packages "$VENV" >"$LOGDIR/venv.log" 2>&1; then
    # (a) offline: uses the setuptools/wheel already installed
    if "$VENV/bin/python" -m pip install --no-deps --no-build-isolation --no-index \
         "$PKGSRC" >>"$LOGDIR/venv.log" 2>&1; then
      how="offline build (--no-build-isolation --no-index)"
    else
      # (b) isolated build: pip fetches setuptools>=68 + wheel from the package index
      echo "--- offline build failed; retrying with an isolated build environment" >>"$LOGDIR/venv.log"
      rm -rf "$PKGSRC/build" "$PKGSRC"/*.egg-info
      if "$VENV/bin/python" -m pip install --no-deps "$PKGSRC" >>"$LOGDIR/venv.log" 2>&1; then
        how="isolated build (offline build failed with the system setuptools, see venv.log)"
      fi
    fi
  fi
  if [ -n "$how" ]; then
    if "$VENV/bin/python" -c 'import numpy, scipy, flint, rtsc, rtsc.certify, rtsc.cli' \
         >>"$LOGDIR/venv.log" 2>&1 && [ -x "$VENV/bin/rtsc" ]; then
      RTSC=("$VENV/bin/rtsc"); PY_RUN="$VENV/bin/python"
      result python PASS "pip install --no-deps OK via $how; core imports OK; entry point $VENV/bin/rtsc"
    else
      result python FAIL "installed, but core imports failed (numpy/scipy/python-flint missing?) - see $LOGDIR/venv.log"
    fi
  else
    say "venv/pip install failed (see $LOGDIR/venv.log); falling back to PYTHONPATH=$BUNDLE/rtsc"
    if PYTHONPATH="$BUNDLE/rtsc" "$PY" -c 'import numpy, scipy, flint, rtsc.cli' >>"$LOGDIR/venv.log" 2>&1; then
      result python WARN "pip install failed (packaging problem, or no setuptools and no index access); PYTHONPATH import works"
    else
      result python FAIL "rtsc not importable - see $LOGDIR/venv.log"
    fi
  fi
fi
if [ ${#RTSC[@]} -eq 0 ]; then
  export PYTHONPATH="$BUNDLE/rtsc${PYTHONPATH:+:$PYTHONPATH}"
  RTSC=("$PY" -m rtsc.cli); PY_RUN="$PY"
fi

# ---------------------------------------------------------------------------------------------
# 2. pytest
# ---------------------------------------------------------------------------------------------
if [ "$SKIP_PYTEST" -eq 0 ]; then
  if "$PY_RUN" -c 'import pytest' 2>/dev/null; then
    say "running pytest (log: $LOGDIR/pytest.log) ..."
    # shellcheck disable=SC2086
    ( cd "$BUNDLE/rtsc" && "${NICECMD[@]}" "$PY_RUN" -m pytest -q -p no:cacheprovider -rs $PYTEST_ARGS ) \
      > "$LOGDIR/pytest.log" 2>&1
    rc=$?
    summary=$(grep -E '^(=+ )?[0-9]+ (passed|failed|error)|passed|failed' "$LOGDIR/pytest.log" | tail -1)
    if [ $rc -eq 0 ]; then result pytest PASS "$summary"
    else result pytest FAIL "exit $rc: $summary (see $LOGDIR/pytest.log)"; fi
  else
    result pytest FAIL "pytest not installed (pip install 'rtsc[test]')"
  fi
else
  result pytest SKIP "--skip-pytest"
fi

# ---------------------------------------------------------------------------------------------
# 3. demo  +  4. recheck
# ---------------------------------------------------------------------------------------------
if [ "$SKIP_DEMO" -eq 0 ]; then
  DEMO="$WORK/demo"
  rm -rf "$DEMO"
  "${NICECMD[@]}" "${RTSC[@]}" demo --workdir "$DEMO" --no-resume > "$LOGDIR/demo.log" 2>&1
  rc=$?
  if [ $rc -ne 0 ] || [ ! -f "$DEMO/run_summary.json" ]; then
    result demo FAIL "rtsc demo exit $rc (see $LOGDIR/demo.log)"
  else
    acct=$("$PY_RUN" - "$DEMO" <<'EOF'
import json, sys, pathlib
d = pathlib.Path(sys.argv[1])
s = json.loads((d / "run_summary.json").read_text())
n, counts = s["n_inputs"], s["counts"]
dirs = [p for p in d.iterdir() if (p / "verdict.json").is_file()]
st = [json.loads((p / "verdict.json").read_text())["status"] for p in dirs]
ok = (sum(counts.values()) == n == len(dirs)
      and all(st.count(k) == v for k, v in counts.items())
      and "SYNTHETIC" in (d / "summary.md").read_text())
print(("OK" if ok else "BAD"), f"inputs={n}", " ".join(f"{k}={v}" for k, v in counts.items()))
EOF
)
    case "$acct" in
      OK*) result demo PASS "${acct#OK } (counts add up, one verdict per input, SYNTHETIC banner present)";;
      *)   result demo FAIL "accounting check failed: $acct";;
    esac
  fi

  # recheck every classified verdict (status that rtsc.certify produced from spectra)
  n_ok=0; n_bad=0; seen_rej=0; seen_cert=0; bad_list=""
  : > "$LOGDIR/recheck.log"
  for v in "$DEMO"/*/verdict.json; do
    [ -f "$v" ] || continue
    st=$("$PY_RUN" -c 'import json,sys; print(json.load(open(sys.argv[1]))["status"])' "$v")
    case "$st" in REJECTED_PROVEN|CERTIFIED_CANDIDATE|UNDECIDED) ;; *) continue;; esac
    out=$("${NICECMD[@]}" "${RTSC[@]}" recheck "$v" 2>&1); rc=$?
    echo "$out" >> "$LOGDIR/recheck.log"
    if [ $rc -eq 0 ] && grep -q "RECHECK OK" <<<"$out"; then
      n_ok=$((n_ok+1))
      [ "$st" = REJECTED_PROVEN ] && seen_rej=1
      [ "$st" = CERTIFIED_CANDIDATE ] && seen_cert=1
    else
      n_bad=$((n_bad+1)); bad_list="$bad_list $(basename "$(dirname "$v")")"
    fi
  done
  # rtsc certify on real (deliberately unconverged, smoke-test) QE output shipped with the tests
  A2F="$BUNDLE/rtsc/tests/data_dfpt/qe67_al_smoke/a2F.dos5"
  cert_note="certify: test file missing"
  if [ -f "$A2F" ]; then
    out=$("${NICECMD[@]}" "${RTSC[@]}" certify "$A2F" --format qe_a2f_dos --mode decide \
            --out "$WORK/certify_al_smoke.json" 2>&1); rc=$?
    echo "$out" >> "$LOGDIR/recheck.log"
    if [ $rc -eq 0 ] && grep -q "REJECTED_PROVEN" <<<"$out"; then
      cert_note="rtsc certify on QE a2F.dos5 (Al smoke test) -> REJECTED_PROVEN as expected"
    else
      cert_note="rtsc certify on QE a2F.dos5 unexpected (exit $rc)"; n_bad=$((n_bad+1))
    fi
  fi
  if [ $n_bad -eq 0 ] && [ $seen_rej -eq 1 ] && [ $seen_cert -eq 1 ]; then
    result recheck PASS "$n_ok demo verdicts re-verified (incl. REJECTED_PROVEN and CERTIFIED_CANDIDATE); $cert_note"
  else
    result recheck FAIL "ok=$n_ok bad=$n_bad [$bad_list] rejected_seen=$seen_rej certified_seen=$seen_cert; $cert_note (see $LOGDIR/recheck.log)"
  fi
else
  result demo SKIP "--skip-demo"; result recheck SKIP "--skip-demo"
fi

# ---------------------------------------------------------------------------------------------
# 5. Lean
# ---------------------------------------------------------------------------------------------
lean_sources() { # list of files that define the Lean project (relative paths)
  ( cd "$1" && find . -path ./.lake -prune -o -type f \
      \( -name '*.lean' -o -name lakefile.toml -o -name lakefile.lean -o -name lean-toolchain \
         -o -name lake-manifest.json \) -print | sort )
}
check_axiom_output() { # FILE -> prints "OK n" or "BAD reason"
  "$PY_RUN" - "$1" "${2:-0}" <<'EOF'
import re, sys
allowed = {"propext", "Classical.choice", "Quot.sound"}
text = open(sys.argv[1]).read()
expected = int(sys.argv[2])
lines = [l for l in text.splitlines() if l.strip()]
n, bad = 0, []
for l in lines:
    m = re.match(r"^'(.+)' depends on axioms: \[(.*)\]\s*$", l)
    if m:
        axs = {a.strip() for a in m.group(2).split(",") if a.strip()}
        n += 1
        if not axs <= allowed: bad.append(l)
    elif re.match(r"^'(.+)' does not depend on any axioms\s*$", l):
        n += 1
    else:
        bad.append("unexpected output: " + l)
if expected and n != expected:
    bad.append(f"expected {expected} #print axioms results, got {n}")
print("OK %d" % n if not bad else "BAD " + " | ".join(bad[:5]))
EOF
}

if [ "$SKIP_LEAN" -eq 0 ]; then
  export PATH="$HOME/.elan/bin:$PATH"
  LDIR="${LEAN_DIR:-$BUNDLE/lean}"
  LDIR="$(cd "$LDIR" 2>/dev/null && pwd)" || LDIR=""
  lean_ok=1; lean_detail=""
  if [ -z "$LDIR" ]; then
    lean_ok=0; lean_detail="Lean directory not found"
  elif ! command -v lake >/dev/null 2>&1; then
    lean_ok=0; lean_detail="lake not found: install elan (README section 5)"
  fi
  if [ $lean_ok -eq 1 ] && [ "$LDIR" != "$BUNDLE/lean" ]; then
    # the external project must contain exactly the bundled sources
    if diff <(lean_sources "$BUNDLE/lean") <(lean_sources "$LDIR") >"$LOGDIR/lean_diff.log" 2>&1 && \
       ( cd "$BUNDLE/lean" && lean_sources . | while read -r f; do
           cmp -s "$f" "$LDIR/$f" || { echo "differs: $f"; exit 1; }; done ) >>"$LOGDIR/lean_diff.log" 2>&1; then
      say "Lean: using $LDIR (sources identical to the bundle's lean/)"
    else
      lean_ok=0; lean_detail="--lean-dir $LDIR has different sources than the bundle (see $LOGDIR/lean_diff.log)"
    fi
  fi
  if [ $lean_ok -eq 1 ] && [ -n "$LEAN_SEED" ] && [ ! -d "$LDIR/.lake" ]; then
    if cmp -s "$LEAN_SEED/lake-manifest.json" "$LDIR/lake-manifest.json"; then
      say "Lean: seeding $LDIR/.lake from $LEAN_SEED/.lake (cp -a) ..."
      cp -a "$LEAN_SEED/.lake" "$LDIR/.lake" || { lean_ok=0; lean_detail="seed copy failed"; }
    else
      say "Lean: --lean-seed ignored, lake-manifest.json differs"
    fi
  fi
  if [ $lean_ok -eq 1 ]; then
    cd "$LDIR" || exit 1
    say "Lean: $(lake --version 2>&1 | head -1); toolchain file: $(cat lean-toolchain)"
    if [ ! -d .lake/packages/mathlib/.lake/build/lib ] && [ $CACHE_GET -eq 1 ]; then
      say "Lean: no Mathlib build found; trying 'lake exe cache get' (needs network access to the Mathlib cache) ..."
      if "${NICECMD[@]}" lake exe cache get > "$LOGDIR/lake_cache_get.log" 2>&1; then
        say "Lean: Mathlib cache downloaded"
      else
        say "Lean: 'lake exe cache get' failed (see $LOGDIR/lake_cache_get.log); lake build will compile the needed Mathlib modules from source (~1 h on 2 cores)"
      fi
    fi
    t0=$(date +%s)
    "${NICECMD[@]}" lake build > "$LOGDIR/lake_build.log" 2>&1; rc=$?
    t1=$(date +%s)
    if [ $rc -ne 0 ]; then
      lean_ok=0; lean_detail="lake build failed (exit $rc, see $LOGDIR/lake_build.log)"
    else
      lean_detail="lake build OK ($((t1-t0)) s)"
      # (a) source scan: no sorry/admit, no new axioms, no unchecked evaluation shortcuts
      lean_sources . | grep '\.lean$' | grep -v '^\./lean_checks/' > "$LOGDIR/lean_files.txt"
      : > "$LOGDIR/lean_source_scan.log"
      while read -r f; do
        grep -nHwE 'sorry|admit' "$f"
        grep -nHE '^[[:space:]]*(private[[:space:]]+|protected[[:space:]]+)?axiom[[:space:]]' "$f"
        grep -nHE 'native_decide|ofReduceBool|implemented_by|@\[extern' "$f"
      done < "$LOGDIR/lean_files.txt" >> "$LOGDIR/lean_source_scan.log"
      nfiles=$(wc -l < "$LOGDIR/lean_files.txt")
      if [ -s "$LOGDIR/lean_source_scan.log" ]; then
        lean_ok=0; lean_detail="$lean_detail; source scan found: $(head -3 "$LOGDIR/lean_source_scan.log" | tr '\n' ' ')"
      else
        lean_detail="$lean_detail; $nfiles .lean files contain no sorry/admit/axiom/native_decide/implemented_by/extern"
      fi
      # (b) the hand-written #print axioms files
      for chk in Check.lean CheckPipeline.lean CheckMultiBand.lean CheckDiscovery.lean; do
        [ -f "$chk" ] || { lean_ok=0; lean_detail="$lean_detail; $chk missing"; continue; }
        expect=$(grep -c '^#print axioms' "$chk")
        "${NICECMD[@]}" lake env lean "$chk" > "$LOGDIR/${chk%.lean}_axioms.txt" 2>&1; rc=$?
        res=$(check_axiom_output "$LOGDIR/${chk%.lean}_axioms.txt" "$expect")
        if [ $rc -eq 0 ] && [[ "$res" == OK* ]]; then
          lean_detail="$lean_detail; $chk: ${res#OK } theorems OK"
        else
          lean_ok=0; lean_detail="$lean_detail; $chk: exit $rc $res"
        fi
      done
      # (c) every declaration of TcR1, PipelineLogic, MultiBand (not a hand-picked list)
      "${NICECMD[@]}" lake env lean "$BUNDLE/lean_checks/AllAxioms.lean" > "$LOGDIR/all_axioms.txt" 2>&1; rc=$?
      summ=$(grep '^AXIOM-AUDIT' "$LOGDIR/all_axioms.txt" | tail -1)
      used=$(grep '^AXIOMS ' "$LOGDIR/all_axioms.txt" | sed 's/.*: \[\(.*\)\].*/\1/' | tr ',' '\n' \
             | sed 's/^ *//;s/ *$//' | grep -v '^$' | sort -u | tr '\n' ' ')
      extra=$(echo "$used" | tr ' ' '\n' | grep -vxE 'propext|Classical\.choice|Quot\.sound|' || true)
      if [ $rc -eq 0 ] && [[ "$summ" =~ decls=([0-9]+)\ bad=0$ ]] && [ "${BASH_REMATCH[1]}" -gt 0 ] && [ -z "$extra" ]; then
        lean_detail="$lean_detail; AllAxioms: ${BASH_REMATCH[1]} declarations, axioms used = { $used}"
      else
        lean_ok=0; lean_detail="$lean_detail; AllAxioms: exit $rc, '$summ', axioms used = { $used} (see $LOGDIR/all_axioms.txt)"
      fi
    fi
    cd - >/dev/null || true
  fi
  if [ $lean_ok -eq 1 ]; then result lean PASS "$lean_detail"; else result lean FAIL "$lean_detail"; fi
else
  result lean SKIP "--skip-lean"
fi

# ---------------------------------------------------------------------------------------------
# 6. certificates (optional). Run on a copy so the shipped logs stay untouched.
# ---------------------------------------------------------------------------------------------
if [ "$DO_CERT" -eq 1 ]; then
  CW="$WORK/cert_wE1"; rm -rf "$CW"; mkdir -p "$CW"
  cp "$BUNDLE"/certificates/wE1/*.py "$BUNDLE"/certificates/wE1/run_all.sh "$CW"/
  say "certificate wE1: running run_all.sh (a0=4.33148145, eta=1e-9, sharpness a1=4.3314815) in $CW ..."
  t0=$(date +%s)
  ( cd "$CW" && PATH="$(dirname "$(command -v "$PY_RUN")"):$PATH" bash ./run_all.sh ) > "$LOGDIR/cert_wE1.log" 2>&1
  t1=$(date +%s)
  if grep -q "^FINAL VERDICT: CERTIFIED" "$LOGDIR/cert_wE1.log"; then
    result cert PASS "wE1: $(grep '^VERDICT cells' "$LOGDIR/cert_wE1.log" | tail -1 | sed 's/VERDICT cells: //'); FINAL VERDICT CERTIFIED ($((t1-t0)) s)"
  else
    result cert FAIL "wE1: no 'FINAL VERDICT: CERTIFIED' (see $LOGDIR/cert_wE1.log)"
  fi
  if [ "$DO_CERT_ALL" -eq 1 ]; then
    BW="$WORK/cert_wB"; rm -rf "$BW"; mkdir -p "$BW"; cp "$BUNDLE"/certificates/wB/*.py "$BW"/
    MW="$WORK/cert_wMB"; rm -rf "$MW"; mkdir -p "$MW"; cp "$BUNDLE"/certificates/wMB/*.py "$MW"/
    ok=1; det=""
    ( cd "$BW" && "${NICECMD[@]}" "$PY_RUN" cert_ends.py 43314814 10000000 13/10 40 ) > "$LOGDIR/cert_wB_ends.log" 2>&1 || ok=0
    if grep -q "False" "$LOGDIR/cert_wB_ends.log" || ! grep -q "Lemma H" "$LOGDIR/cert_wB_ends.log"; then ok=0; det="$det ends:FAIL"; else det="$det ends:OK"; fi
    # sharpness at a = 4.33148149 (PROOF.md Prop. 7: c* >= 1/4.33148149)
    ( cd "$BW" && "${NICECMD[@]}" "$PY_RUN" cert_lower.py 433148149 100000000 75748661 10000000 300 ) > "$LOGDIR/cert_wB_lower.log" 2>&1 || ok=0
    if grep -q "CERTIFIED Q(x) > 0: True" "$LOGDIR/cert_wB_lower.log"; then det="$det lower:OK"; else ok=0; det="$det lower:FAIL"; fi
    ( cd "$BW" && "${NICECMD[@]}" "$PY_RUN" cert_upper.py 43314814 10000000 13/10 40 80 1000 ) > "$LOGDIR/cert_wB_upper.log" 2>&1 || ok=0
    if grep -q "^CERTIFIED:" "$LOGDIR/cert_wB_upper.log" && ! grep -q "FAIL" "$LOGDIR/cert_wB_upper.log"; then det="$det upper:OK"; else ok=0; det="$det upper:FAIL"; fi
    ( cd "$MW" && "${NICECMD[@]}" "$PY_RUN" cert_neg.py ) > "$LOGDIR/cert_wMB_neg.log" 2>&1 || ok=0
    if grep -q "VERDICT Proposition N scalar inequalities: PASS" "$LOGDIR/cert_wMB_neg.log"; then det="$det wMB-PropN:OK"; else ok=0; det="$det wMB-PropN:FAIL"; fi
    if [ $ok -eq 1 ]; then result cert-B PASS "wB (upper a=4.3314814 on [1.3,40] + ends; sharpness a=4.33148149) and wMB Prop. N:$det"; else result cert-B FAIL "$det (logs in $LOGDIR)"; fi
  fi
else
  result cert SKIP "not requested (use --cert; about 6 core-minutes)"
fi

# ---------------------------------------------------------------------------------------------
{
  echo
  echo "================ verify_all.sh summary ================"
  for r in "${RESULTS[@]}"; do echo "$r"; done
  if [ $FAILED -eq 0 ]; then echo "OVERALL: PASS"; else echo "OVERALL: FAIL"; fi
  echo "report: $REPORT   logs: $LOGDIR"
} | tee -a "$REPORT"
exit $FAILED
