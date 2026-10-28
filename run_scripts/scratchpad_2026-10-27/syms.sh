B=~/micromamba/envs/qe/bin
file -L $B/ph.x; nm $B/ph.x 2>&1 | head -2; nm -D $B/ph.x | wc -l; nm -D --defined-only $B/ph.x | grep -i -E 'adddvscf|elphel|elphon' | head
readelf -S $B/ph.x | grep -E 'debug|symtab'
ls ~/micromamba/envs/qe/lib | grep -i -E '\.debug|qe' | head
~/micromamba/envs/dbg/bin/gdb --version | head -1
ls ~/micromamba/pkgs/*.conda 2>/dev/null | wc -l; du -cm ~/micromamba/pkgs/*.conda ~/micromamba/pkgs/*.tar.bz2 2>/dev/null | tail -1
