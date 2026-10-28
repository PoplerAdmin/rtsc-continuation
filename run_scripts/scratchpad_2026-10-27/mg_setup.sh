set -e
export PATH=$HOME/.local/bin:$PATH
cd ~
[ -d mattergen ] || git clone -q https://github.com/microsoft/mattergen
cd mattergen && git log -1 --format='%H %cd' && git lfs version 2>/dev/null || echo "no git-lfs (checkpoints via Hugging Face)"
uv venv -q .venv --python 3.10
. .venv/bin/activate
( time uv pip install -e . ) > ~/rtsc_work/logs/mattergen_install.log 2>&1 || { tail -30 ~/rtsc_work/logs/mattergen_install.log; exit 1; }
tail -3 ~/rtsc_work/logs/mattergen_install.log
python -c "import torch; print('torch', torch.__version__, 'cuda', torch.version.cuda, 'avail', torch.cuda.is_available(), torch.cuda.get_device_name(0) if torch.cuda.is_available() else '-', torch.cuda.get_arch_list() if torch.cuda.is_available() else '')"
python -c "import torch; x=torch.randn(2048,2048,device='cuda'); print('matmul ok', float((x@x).sum())!=0)"
which mattergen-generate
