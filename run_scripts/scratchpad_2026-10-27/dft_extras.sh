set -e
export PATH=$HOME/.local/bin:$PATH
cd ~/rtsc_work && . .venv/bin/activate
uv pip install -e "rtsc_bundle/rtsc[dft,mp,test]" 2>&1 | tail -4
python -c "import pymatgen.core, ase, phonopy, spglib, yaml; print('dft extras OK', phonopy.__version__)"
# SSSP 1.3 efficiency (PBE) from Materials Cloud
mkdir -p ~/pseudo && cd ~/pseudo
curl -sSL -o sssp_eff.tar.gz "https://archive.materialscloud.org/record/file?record_id=1732&filename=SSSP_1.3.0_PBE_efficiency.tar.gz" || true
curl -sSL -o sssp_eff.json "https://archive.materialscloud.org/record/file?record_id=1732&filename=SSSP_1.3.0_PBE_efficiency.json" || true
ls -la ~/pseudo; file ~/pseudo/*
