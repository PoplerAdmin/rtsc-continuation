#!/bin/bash
# Submit a sharded rtsc run plus the collect job that runs when all shards have ended.
#   env/slurm/submit.sh config.toml 32 [max_concurrent]
set -euo pipefail
CONFIG="$(realpath "${1:?usage: submit.sh config.toml NSHARDS [max_concurrent]}")"
N="${2:?NSHARDS}"
MAXC="${3:-8}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p rtsc_logs
jid=$(sbatch --parsable --array="0-$((N-1))%${MAXC}" \
        --export=ALL,RTSC_CONFIG="$CONFIG",RTSC_NSHARDS="$N" "$HERE/rtsc_array.sbatch")
# afterany: collect also runs when some shards failed; their materials are reported as ERROR
cid=$(sbatch --parsable --dependency="afterany:${jid}" \
        --export=ALL,RTSC_CONFIG="$CONFIG" "$HERE/rtsc_collect.sbatch")
echo "array job $jid ($N shards, <= $MAXC at once); collect job $cid"
