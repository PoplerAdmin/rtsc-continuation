#!/bin/bash
# End-to-end certificate.  usage: ./run_all.sh [A0 ETA A1 N K]
# defaults = the tight certificate documented in CERT.md
set -u
cd "$(dirname "$0")"
A0=${1:-4.33148145}; ETA=${2:-1e-9}; A1=${3:-4.3314815}; N=${4:-150}; K=${5:-1500}
LOG=log_${A0}.txt
{
echo "=== certificate run $(date -u)  a0=$A0 eta=$ETA a1=$A1 N=$N K=$K"
nice python3 ends.py "$A0" "$ETA" 1 400;            r1=$?
nice python3 run_cells.py "$A0" "$ETA" "$N" "$K" 1 400 "cells_${A0}.log"; r2=$?
nice python3 sharp.py "$A1" 300;                    r3=$?
if [ $r1 -eq 0 ] && [ $r2 -eq 0 ] && [ $r3 -eq 0 ]; then
  echo "FINAL VERDICT: CERTIFIED  R(a0=$A0, eta=$ETA) for all w>0 (both sectors);  rho+(a1=$A1, w) > 0 for some w"
else
  echo "FINAL VERDICT: NOT CERTIFIED (ends=$r1 cells=$r2 sharp=$r3)"
fi
} 2>&1 | tee "$LOG"
