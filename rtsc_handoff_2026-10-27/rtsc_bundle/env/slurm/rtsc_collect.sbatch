#!/bin/bash
# Merge the per-material records of a sharded run and write the reports (candidates.csv,
# certificates/, summary.md). Cheap: one core, minutes. Materials without a finished record
# are reported as ERROR ("no finished record"), never dropped (S4).
#SBATCH --job-name=rtsc-collect
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=rtsc_logs/%x_%j.out
set -euo pipefail
CONFIG="${RTSC_CONFIG:?set RTSC_CONFIG to the absolute path of config.toml}"
# (activate the same environment as in rtsc_array.sbatch)
rtsc collect -c "$CONFIG"
