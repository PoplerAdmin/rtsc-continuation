# Real (non-dry) test on a synthetic tree in /tmp; central manifests redirected to /tmp too.
set -e; set +e
J=/mnt/g/ClaudeWorkspace/RTSC_Project/bench_tools/jobs
T=$(mktemp -d -p "$HOME"); C=$T/central
sed "s#^CENTRAL = .*#CENTRAL = Path('$C')#" $J/prune.py > $T/prune.py
P="python3 $T/prune.py"
D=$T/dfpt; mkdir -p $D/tmp/_ph0/x.save
for n in 1 2 3 4; do q=$D/tmp/_ph0/x.q_$n; mkdir -p $q/x.save; head -c 100000 /dev/urandom > $q/x.wfc1; head -c 5000 /dev/urandom > $q/x.dvscf1; echo w > $q/x.save/wfc1.hdf5; done
echo "freq" > $D/dyn1; echo "freq" > $D/dyn2; : > $D/dyn3     # q1,q2 done; q3 live (empty dyn)
sha_q1=$(sha256sum $D/tmp/_ph0/x.q_1/x.wfc1 | cut -d' ' -f1)
$P qdone $D/tmp/_ph0 $D --reason test
echo "-- remaining wfc:"; find $D -name 'x.wfc1' | sort
echo "-- kept dvscf/save count: $(find $D -name 'x.dvscf1' | wc -l) $(find $D -name 'wfc1.hdf5' | wc -l)"
echo "-- manifest:"; cut -f3-5,7 $D/tmp/_ph0/PRUNED_MANIFEST.tsv
grep -q "$sha_q1" $D/tmp/_ph0/PRUNED_MANIFEST.tsv && echo "q1 hash matches pre-delete sha256sum: OK"
cmp $D/tmp/_ph0/PRUNED_MANIFEST.tsv $C/*.tsv && echo "central copy identical: OK"
echo "-- qmark:"; M=$T/marks; mkdir $M; touch $M/x.dyn3.elph.3 $M/x.dyn4.elph.4; $P qmark $D/tmp/_ph0 $M --mark 'x.dyn{n}.elph.{n}' --reason test2; find $D -name 'x.wfc1' | sort
echo "-- inventory:"; $P inventory $D; wc -l < $D/RUN_MANIFEST.tsv
echo "-- missing reason:"; $P wfc $D || true
echo "-- verify (intact):"; $P verify $D/RUN_MANIFEST.tsv; echo "exit $?"
echo "x" >> $D/dyn1; echo "-- verify (after tampering dyn1):"; $P verify $D/RUN_MANIFEST.tsv; echo "exit $?"
rm -rf $T
