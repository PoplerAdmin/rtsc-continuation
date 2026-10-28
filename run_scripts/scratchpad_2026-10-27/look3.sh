cd ~/rtsc_work/bisect; ls
for d in */; do echo "== $d"; tail -5 $d/ph.out; tail -15 $d/ph.err; done 2>/dev/null | head -150
