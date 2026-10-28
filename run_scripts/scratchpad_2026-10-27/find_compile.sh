cd ~/rtsc_work
SP=$(.venv/bin/python -I -c 'import mattersim,os;print(os.path.dirname(mattersim.__file__))')
echo $SP
grep -rn "compile\|triton\|jit" $SP --include=*.py | grep -v "re.compile" | head -40
echo ---- rtsc
grep -rn "compile\|triton\|TORCH" ~/rtsc_work/rtsc_bundle/rtsc/src/rtsc/mlff 2>/dev/null | grep -v "re.compile" | head; ls ~/rtsc_work/rtsc_bundle/rtsc
echo ---- env
env | grep -i "torch\|triton\|CC" ; which gcc cc clang
ls ~/rtsc_work/mattergen_test
