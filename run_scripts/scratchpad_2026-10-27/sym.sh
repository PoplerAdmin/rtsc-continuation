Q=$HOME/micromamba/envs/qe/bin
nm $Q/ph.x | grep -E ' (_start|main|MAIN__)$'
f=~/rtsc_work/crashchar/al_paw_eph_simple/ph.err
start=$(nm $Q/ph.x | awk '$3=="_start"{print $1}')
a10=$(grep '^#10' $f | awk '{print $2}')
python3 - "$start" "$a10" $(grep -E '^#[1-7] ' $f | awk '{print $2}') <<'PY' > /tmp/offs
import sys
st=int(sys.argv[1],16); a10=int(sys.argv[2],16)
# _start return address is _start+0x29 typically; find base as page-aligned
base=(a10-st)&~0xfff
for a in sys.argv[3:]:
    print(hex(int(a,16)-base-1))
PY
cat /tmp/offs
which addr2line || ls $Q/*addr2line* 2>/dev/null
# fallback: nearest preceding symbol via nm
nm -n $Q/ph.x | grep -iE ' [tT] ' > /tmp/syms
python3 - /tmp/syms $(cat /tmp/offs) <<'PY'
import sys,bisect
syms=[(int(l.split()[0],16),l.split()[2]) for l in open(sys.argv[1])]
addrs=[s[0] for s in syms]
for o in sys.argv[2:]:
    o=int(o,16); i=bisect.bisect_right(addrs,o)-1; print(hex(o), syms[i][1], hex(o-syms[i][0]))
PY
