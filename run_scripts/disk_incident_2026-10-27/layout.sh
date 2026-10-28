cd ~/rtsc_work/bench/Al/run/*/dfpt/tmp
ls | head; echo; ls _ph0 | head; echo; d=$(ls -d _ph0/*.q_2 | head -1); echo "$d"; ls -la $d | head -20; ls -la $d/*.save 2>/dev/null | head -10
echo "== size by type (Al tmp)"; find . -type f -printf '%s %f\n' | awk '{n=$2; sub(/[0-9]+(\.dat)?$/,"",n); s[n]+=$1} END{for(k in s) printf "%8.2f GB %s\n", s[k]/1e9, k}' | sort -rn | head -15
