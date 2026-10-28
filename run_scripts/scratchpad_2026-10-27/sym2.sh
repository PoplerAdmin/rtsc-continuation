Q=$HOME/micromamba/envs/qe/bin
nm -n $Q/ph.x | grep -iE ' [tT] ' > /tmp/syms
for f in ~/rtsc_work/crashchar/{v_twostep,mg_us_gbrv_interp,nb_paw_interp}/ph.err; do echo "== $f"
python3 - /tmp/syms $f <<'PY'
import sys,bisect,re
syms=[(int(l.split()[0],16),l.split()[2]) for l in open(sys.argv[1])]
addrs=[s[0] for s in syms]
fr=[int(m.group(1),16) for m in (re.match(r'#\d+\s+0x([0-9a-f]+) in',l) for l in open(sys.argv[2])) if m]
st=[s for s in syms if s[1]=='_start'][0][0]
# last frame before 0xffff.. is _start
a=[x for x in fr if x<0x7f0000000000 and x>0x500000000000]
base=(a[-1]-st)&~0xfff
for x in a:
    o=x-base-1; i=bisect.bisect_right(addrs,o)-1; print(hex(o), syms[i][1], hex(o-syms[i][0]))
PY
done
addr2line -f -e $Q/ph.x 0x287160 0x160929
