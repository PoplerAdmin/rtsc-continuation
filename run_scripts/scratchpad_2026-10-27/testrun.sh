export QE_BIN=~/micromamba/envs/qe/bin SSSP_DIR=~/pseudo/SSSP_1.3.0_PBE_precision SG15_DIR=~/pseudo/sg15 WORK=~/rtsc_work/qebug/runsh
nice -n 10 bash /mnt/g/ClaudeWorkspace/RTSC_Project/reports/qe_bug_repro/run.sh all 2>/dev/null
QE_BIN=~/micromamba/envs/qe75b/bin nice -n 10 bash /mnt/g/ClaudeWorkspace/RTSC_Project/reports/qe_bug_repro/run.sh al_paw 2>/dev/null
