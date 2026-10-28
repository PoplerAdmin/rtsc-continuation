cd ~/pseudo/sg15
for f in Pb_ONCV_PBE-1.2.upf Pb_ONCV_PBE_FR-1.0.upf Nb_ONCV_PBE-1.2.upf Mg_ONCV_PBE-1.2.upf B_ONCV_PBE-1.2.upf H_ONCV_PBE-1.2.upf S_ONCV_PBE-1.2.upf; do
 echo "=== $f"; grep -m1 -o 'z_valence="[^"]*"' $f; grep -m1 -o 'relativistic="[^"]*"' $f
 awk '/# lmax/,/# Model core/' $f | head -12; grep -A2 -m1 'icmod' $f
done
cat ~/rtsc_work/crashchar/summary.txt
