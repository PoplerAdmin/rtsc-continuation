#!/bin/bash
# Host-disk guard for WSL jobs. The WSL disk image lives on G:, so the free space inside WSL
# (df /) is NOT the limit; the host drive's free space is (df /mnt/g).
# Below LOW_GB: SIGSTOP all QE/EPW/Wannier processes (paused, not killed; work is kept).
# Above HIGH_GB: SIGCONT them. Every transition is logged.
LOW_GB=${LOW_GB:-50}
HIGH_GB=${HIGH_GB:-60}
LOG=${LOG:-$HOME/rtsc_work/disk_guard.log}
PROCS='pw.x|ph.x|epw.x|wannier90.x|q2r.x|matdyn.x|lambda.x'
paused=0
echo "$(date) start: low=${LOW_GB}G high=${HIGH_GB}G" >> "$LOG"
while true; do
  free_gb=$(df -BG --output=avail /mnt/g | tail -1 | tr -dc '0-9')
  if [ "$paused" = 0 ] && [ "$free_gb" -lt "$LOW_GB" ]; then
    pgrep -f -- "$PROCS" | xargs -r kill -STOP
    paused=1; echo "$(date) PAUSED jobs: G: free ${free_gb}G < ${LOW_GB}G" >> "$LOG"
  elif [ "$paused" = 1 ] && [ "$free_gb" -ge "$HIGH_GB" ]; then
    pgrep -f -- "$PROCS" | xargs -r kill -CONT
    paused=0; echo "$(date) RESUMED jobs: G: free ${free_gb}G" >> "$LOG"
  fi
  sleep 30
done
