export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd ~/rtsc_work/epw_Al/wvar_A/wannier
mkdir -p dbg && cp pwscf.win pwscf.mmn pwscf.nnkp dbg/ 2>/dev/null
ls pwscf.amn pwscf.eig 2>&1
grep -iE 'scdm|auto_proj|wdata|nbndsub|bands_skipped' epw.in
