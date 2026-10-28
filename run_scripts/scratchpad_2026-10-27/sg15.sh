set -e
mkdir -p ~/pseudo/sg15 && cd ~/pseudo/sg15
[ -s sg15.tar.gz ] || curl -sSL -o sg15.tar.gz http://www.quantum-simulation.org/potentials/sg15_oncv/sg15_oncv_upf_2020-02-06.tar.gz
sha256sum sg15.tar.gz; tar -xzf sg15.tar.gz
ls | grep -E "^(Al|Pb|Nb|Mg|B|H|S)_ONCV_PBE" 
md5sum Al_ONCV_PBE-1.2.upf Pb_ONCV_PBE-1.2.upf Nb_ONCV_PBE-1.2.upf 2>/dev/null
grep -m3 -iE "z_valence|core_correction|pseudo_type" Al_ONCV_PBE-1.2.upf
