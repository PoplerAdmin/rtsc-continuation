set -e
curl -sSf https://raw.githubusercontent.com/leanprover/elan/master/elan-init.sh -o /tmp/elan-init.sh
sh /tmp/elan-init.sh -y --default-toolchain none --no-modify-path > ~/rtsc_work/elan_install.log 2>&1
export PATH=$HOME/.elan/bin:$PATH
cd ~/rtsc_work/rtsc_bundle/lean
cat lean-toolchain
lake --version
( time lake exe cache get ) > ~/rtsc_work/lake_cache_get.log 2>&1
tail -5 ~/rtsc_work/lake_cache_get.log
