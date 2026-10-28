#!/bin/bash
# Sequential validation runs: <Material> <T*> <mu_hi> triples as arguments.
# Experimental Tc used as T*: Al 1.175 K, Pb 7.19 K, Nb 9.25 K, MgB2 39 K, H3S 203 K (~155 GPa).
J=/mnt/g/ClaudeWorkspace/RTSC_Project/bench_tools/jobs
LOG=~/rtsc_work/validate/chain.log; mkdir -p ~/rtsc_work/validate
while [ $# -ge 3 ]; do
  echo "$(date '+%F %T') START $1 T*=$2 mu_hi=$3" >> $LOG
  bash $J/validate_material.sh "$1" "$2" "$3" >> $LOG 2>&1
  echo "$(date '+%F %T') END $1 T*=$2 mu_hi=$3 ($(grep wall_s ~/rtsc_work/validate/$1/T$2_mu$3/run.log))" >> $LOG
  shift 3
done
