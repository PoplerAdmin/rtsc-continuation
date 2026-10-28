cd /mnt/g/ClaudeWorkspace/RTSC_Project/run_scripts/scratchpad_2026-10-27/wsl/src
mkdir -p ts
grep -E "^(test-suite/epw|EPW/examples)" list2.txt | grep -E "epw[0-9a-z_]*\.in$" > tslist.txt
while read p; do f=ts/$(echo $p | tr / _); [ -s $f ] || curl -sL "https://raw.githubusercontent.com/QEF/q-e/qe-7.6/$p" -o $f; done < tslist.txt
grep -il "scdm" ts/* 
echo ---
grep -l "a2f *= *.true." ts/* 
echo ---
grep -il "delta_approx" ts/*
