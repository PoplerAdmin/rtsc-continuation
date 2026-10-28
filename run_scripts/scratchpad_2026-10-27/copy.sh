D=/mnt/g/ClaudeWorkspace/RTSC_Project/reports/data; mkdir -p $D
cd ~/rtsc_work/mattergen_test && cp prescreen_h3s_*.json prescreen_mgb_*.json h3s_supercell_xcheck.json $D/ && cp -r hull_cache $D/ && ls -R $D
