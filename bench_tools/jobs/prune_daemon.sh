#!/bin/bash
# Every 120 s: hash-then-delete the wavefunctions of finished q points in live ph.x runs
# (policy: reports/EVIDENCE_RETENTION.md). Only wavefunction files are touched; see prune.py.
P="python3 /mnt/g/ClaudeWorkspace/RTSC_Project/bench_tools/jobs/prune.py"
LOG=${LOG:-$HOME/rtsc_work/prune_daemon.log}
# Runs whose scratch is held until a decision is made about them (Pb: verdict ERROR, unresolved).
HOLD=${HOLD:-/Pb/}
echo "$(date) start (hold: $HOLD)" >> "$LOG"
while true; do
  for d in "$HOME"/rtsc_work/bench/*/run/*/dfpt; do
    [ -d "$d/tmp/_ph0" ] || continue
    case "$d" in *$HOLD*) continue ;; esac
    $P qdone "$d/tmp/_ph0" "$d" --reason "q point finished (dyn complete, next q started)" 2>&1 \
      | grep -v '^pruned 0.00 GB' | sed "s/^/$(date +%F\ %T) /" >> "$LOG"
  done
  W=$HOME/rtsc_work/tetra_Al/q6s_opt
  [ -d "$W/tmp/_ph0" ] && $P qmark "$W/tmp/_ph0" "$W" --mark 'al.dyn{n}.elph.{n}' \
      --reason "lambda_tetra q point finished (next q written)" 2>&1 \
      | grep -v '^pruned 0.00 GB' | sed "s/^/$(date +%F\ %T) /" >> "$LOG"
  sleep 120
done
