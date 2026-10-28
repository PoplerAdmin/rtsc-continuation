# Same header rename on the per-run copies inside WSL (manifests written before the rename).
for m in $(find ~/rtsc_work -name 'PRUNED_MANIFEST.tsv' -o -name 'RUN_MANIFEST.tsv'); do sed -i '1s/^utc_time\thost\t/utc_time\tmachine\t/' "$m"; head -1 "$m" | cut -f2; done | sort | uniq -c
python3 /mnt/g/ClaudeWorkspace/RTSC_Project/bench_tools/jobs/prune.py verify ~/rtsc_work/bench/Al/run/*/dfpt/RUN_MANIFEST.tsv
