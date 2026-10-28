B=~/micromamba/envs/qe/bin
for e in Al Pb Nb Mg B H S; do
 f=$(ls ~/pseudo/sg15/${e}_ONCV_PBE-1.2.upf 2>/dev/null || ls ~/pseudo/sg15/${e}_ONCV_PBE-1.*.upf | tail -1)
 echo "== $f"; grep -n -A2 '<PP_PSWFC' $f | head -4; grep -o 'number_of_wfc="[^"]*"' $f; grep -c "PP_CHI" $f
done
echo "== SSSP efficiency"
for e in Al Pb Nb Mg B H S; do f=$(ls ~/pseudo/SSSP_1.3.0_PBE_efficiency/ | grep -i "^${e}[._-]" | head -1); echo "$e $f CHI=$(grep -c '<PP_CHI' ~/pseudo/SSSP_1.3.0_PBE_efficiency/$f) $(grep -o 'number_of_wfc="[^"]*"' ~/pseudo/SSSP_1.3.0_PBE_efficiency/$f | head -1)"; done
echo "== epw_pp.py"; cat $B/epw_pp.py
echo "== w90"; $B/wannier90.x -v 2>&1 | head -5; strings $B/epw.x | grep -i "wannier90 v\|Version 3\|w90" | head -5; ldd $B/epw.x | grep -i wann
ls ~/micromamba/envs/qe/lib | grep -i wann
ls $B | grep -i -E "\.py$|kmesh|pp"
micromamba list -n qe 2>/dev/null | grep -i -E "qe|quantum|wannier" ; ~/micromamba/bin/micromamba list -n qe 2>/dev/null | grep -i -E "qe|quantum|wannier"
