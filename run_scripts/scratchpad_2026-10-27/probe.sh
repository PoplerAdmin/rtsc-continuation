source ~/rtsc_work/.venv/bin/activate
python -c "import rtsc,flint,numpy,scipy;print(rtsc.__file__,flint.__version__,numpy.__version__,scipy.__version__)"
nproc; uptime
ls ~/rtsc_work/bench/ 2>&1; ls -la ~/rtsc_work/bench/Al/run/ 2>&1
for d in ~/rtsc_work/bench/Al/run/Al-fcc_0GPa-*/; do echo "== $d"; ls -la $d; ls $d/dfpt 2>/dev/null | head -40; done
find ~/rtsc_work/bench -maxdepth 3 -name "*.toml" | head
