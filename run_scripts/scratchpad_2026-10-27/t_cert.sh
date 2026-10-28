# usage: t_cert.sh <elph dir>
source ~/rtsc_work/.venv/bin/activate
S=/mnt/g/ClaudeWorkspace/RTSC_Project/run_scripts/scratchpad_2026-10-27
D=$1; cd $D; mkdir -p cert
L=$(grep -m1 "lambda :" a2F.out | awk '{print $3}')
echo "expected_lambda (alpha2f.x stdout, from a2F) = $L"
F=al.a2F.dat; ID=Al-fcc-tetraopt-$(basename $(dirname $D))-$(basename $D)
python $S/dump_spectra.py $F Ry bins error 0.10 cert/spectra.json $L
run() { tag=$1; shift
  echo "=== certify $tag: $*"
  rtsc certify $F --format two_column --units Ry --rule bins --expected-lambda $L --material-id $ID "$@" --out cert/verdict_$tag.json 2> cert/certify_$tag.err | tail -n +2 | cut -c1-600
  cat cert/certify_$tag.err | head -5
  echo "--- recheck $tag"; rtsc recheck cert/verdict_$tag.json --spectra cert/spectra.json 2>&1 | grep -v "^warning: no --config"
}
run T300_mu015 --target 300 --mu-hi 0.15 --mode bracket
run T1.2_mu010 --target 1.2 --mu-hi 0.10 --mode bracket
run T1.2_mu015 --target 1.2 --mu-hi 0.15 --mode bracket
