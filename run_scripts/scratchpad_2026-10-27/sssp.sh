set -e
cd ~/pseudo && rm -f sssp_eff.json sssp_eff.tar.gz
B="https://archive.materialscloud.org/records/rcyfm-68h65/files"
for f in SSSP_1.3.0_PBE_efficiency SSSP_1.3.0_PBE_precision; do
  curl -sSL -o $f.tar.gz "$B/$f.tar.gz?download=1"
  curl -sSL -o $f.json "$B/$f.json?download=1"
done
cat > md5.txt <<'M'
a58f1b3373f330179fd0832c48bb9a52  SSSP_1.3.0_PBE_efficiency.tar.gz
3153c4b20fc90a44fba0236627525644  SSSP_1.3.0_PBE_efficiency.json
fde94756886f32ada7bf597547557eb5  SSSP_1.3.0_PBE_precision.tar.gz
1692c5c9ce89e1c7c783f8f0eee0cbfa  SSSP_1.3.0_PBE_precision.json
M
md5sum -c md5.txt
for f in SSSP_1.3.0_PBE_efficiency SSSP_1.3.0_PBE_precision; do mkdir -p $f && tar -xzf $f.tar.gz -C $f; done
ls SSSP_1.3.0_PBE_precision | head; ls SSSP_1.3.0_PBE_precision | wc -l
python3 -c "
import json
for v in ['precision','efficiency']:
  d=json.load(open('SSSP_1.3.0_PBE_%s.json'%v))
  print(v, {e:(d[e]['filename'],d[e]['cutoff_wfc'],d[e]['cutoff_rho']) for e in ['Al','Pb','Nb','Mg','B','H','S']})"
