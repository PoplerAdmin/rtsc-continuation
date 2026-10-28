cd /mnt/g/ClaudeWorkspace/RTSC_Project/run_scripts/scratchpad_2026-10-27/wsl/src/raw
for p in EPW/src/spectral.f90 EPW/src/selfen.f90 EPW/src/io/io_selfen.f90 EPW/src/utilities/utilities.f90; do
  f=$(echo $p | tr / _); curl -sfL "https://raw.githubusercontent.com/QEF/q-e/qe-7.6/$p" -o $f || echo "FAIL $p"
done
ls -la EPW_src_spectral.f90 EPW_src_selfen.f90 EPW_src_io_io_selfen.f90
