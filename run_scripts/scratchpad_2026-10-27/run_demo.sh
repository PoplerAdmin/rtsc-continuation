source ~/rtsc_work/.venv/bin/activate
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd ~/rtsc_work/bench_out
( time nice -n 10 rtsc recheck smoke_Al/Al_verdict_Tlit_mu0.10.json --spectra smoke_Al/Al_spectra.json ) 2>&1 | tail -6
rm -rf demo_wd; nice -n 10 rtsc demo --workdir ~/rtsc_work/bench_out/demo_wd > demo.log 2>&1; tail -8 demo.log
ls demo_wd
