#!/bin/bash
# After moving the distro: filesystem health, sizes of the scratch dirs to clean, host free space.
echo "== host G: free"; df -h /mnt/g | tail -1
echo "== ext4 errors this boot"; dmesg 2>/dev/null | grep -iE 'ext4|I/O error|EXT4-fs' | tail -10
echo "== previous boot ext4 errors"; journalctl -b -1 -k --no-pager 2>/dev/null | grep -iE 'EXT4-fs (error|warning)|I/O error' | tail -10
cd ~/rtsc_work
echo "== top-level usage"; du -sh * 2>/dev/null | sort -h | tail -12
echo "== cleanup candidates"
for d in tetra_Al/q6s_opt/elph_k16/tmp tetra_Al/q6s_opt/elph_k24 tetra_Al/q6s_lin/tmp; do du -sh "$d" 2>/dev/null || echo "missing $d"; done
echo "== Nb dyn file tails"; for f in bench/Nb/run/*/dfpt/dyn*; do tail -c 200 "$f" | grep -q . || echo "EMPTY TAIL $f"; done; ls bench/Nb/run/*/dfpt/ | grep -c '^dyn'
