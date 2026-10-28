#!/bin/bash
# Re-certify a finished benchmark material with the branch rtsc code, reusing its QE outputs.
# usage: validate_material.sh <Material> <T_target_K> <mu_hi> [tag]
#   Creates ~/rtsc_work/validate/<Material>/<tag>/ with:
#     config.toml  = bench config + smearing_admission = "mode_positive" + the given [certify] T*, mu_hi
#     structures/  = copy of the bench structures
#     run/<id>/dfpt = copy of the finished DFPT directory WITHOUT scratch (only tmp/<prefix>.a2Fsave,
#                     which the scf_fit completion check needs), so every QE step is adopted, not rerun.
#   Then runs `rtsc run` with PYTHONPATH pointing at the repository's rtsc (this branch).
set -e
M=$1; T=$2; MU=$3; TAG=${4:-T${2}_mu${3}}
B=~/rtsc_work/bench/$M
R=/mnt/g/ClaudeWorkspace/RTSC_Project/rtsc_handoff_2026-10-27/rtsc_bundle/rtsc
V=~/rtsc_work/validate/$M/$TAG
SRC=$(ls -d $B/run/*/dfpt); ID=$(basename "$(dirname "$SRC")")
rm -rf "$V"; mkdir -p "$V/run/$ID/dfpt/tmp"
cp -r $B/structures "$V/"
rsync -a --exclude tmp "$SRC/" "$V/run/$ID/dfpt/"
cp "$SRC"/tmp/*.a2Fsave "$V/run/$ID/dfpt/tmp/"
sed -e 's/^\[dfpt.settings\]$/[dfpt.settings]\nsmearing_admission = "mode_positive"/' \
    -e "s/^T_target_K = .*/T_target_K = $T/" -e "s/^mu_hi = .*/mu_hi = $MU/" \
    -e 's/^name = "\(.*\)"/name = "\1-validate"/' $B/config.toml > "$V/config.toml"
grep -q '^smearing_admission = "mode_positive"' "$V/config.toml"
grep -q "^T_target_K = $T\$" "$V/config.toml"
cd "$V"
t0=$(date +%s.%N)
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONPATH=$R nice -n 10 ~/rtsc_work/.venv/bin/python -c "import sys, rtsc; assert rtsc.__file__.startswith('$R'), rtsc.__file__; from rtsc.cli import main; sys.argv=['rtsc','run','-c','config.toml']; main()" > run.log 2>&1 || true
t1=$(date +%s.%N)
echo "wall_s $(echo "$t1 - $t0" | bc)" >> run.log
python3 - "$V/run/$ID" <<'EOF'
import json, sys, pathlib
d = pathlib.Path(sys.argv[1])
v = json.loads((d / "verdict.json").read_text())
st = json.loads((d / "dfpt" / "dfpt_state.json").read_text())["steps"]
print(f"{v['material_id']}: {v['status']}  T_low={v.get('T_low_K')}  T_up={v.get('T_up_K')}")
print("  rules:", v.get("rule_ids"), " flags:", v.get("flags"))
print("  steps:", {k: ("adopted" if s.get("adopted") else s.get("status")) for k, s in st.items()})
for r in v.get("reasons", [])[:3]:
    print("  reason:", str(r)[:300])
EOF
