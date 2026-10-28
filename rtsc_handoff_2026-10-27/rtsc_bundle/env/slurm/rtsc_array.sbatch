#!/bin/bash
# Slurm job-array template for `rtsc run --shard K/N`.  Submit with env/slurm/submit.sh, which
# also queues the `rtsc collect` job; or by hand:
#   sbatch --array=0-31%8 --export=ALL,RTSC_CONFIG=$PWD/config.toml,RTSC_NSHARDS=32 rtsc_array.sbatch
#
# Shard K computes positions K, K+N, K+2N, ... of the global priority order. Stages 1-5
# (sources, dedupe, filters, ordering, budget) are cheap and deterministic, so every task repeats
# them. With [ml] enabled, run `rtsc run -c config.toml --plan-only` once before the array.
# Re-running a shard is safe: finished materials are skipped and QE steps restart from their
# own checkpoints.
#
# Resources: rtsc starts one `srun -n <nprocs>` job step per QE call ([dfpt.launcher]
# kind = "srun"). With [run] workers = W, W materials run at the same time, so request
# ntasks = W * nprocs. Add "--exact" to mpi_args so that concurrent steps do not each try to
# take the whole allocation (Slurm >= 21.08; older versions use "--exclusive").
#
#SBATCH --job-name=rtsc-shard
#SBATCH --array=0-31%8               # N = 32 shards, at most 8 running at once (edit)
#SBATCH --nodes=1
#SBATCH --ntasks=64                  # = workers x [dfpt.launcher] nprocs   (edit)
#SBATCH --cpus-per-task=1            # pure MPI; OMP_NUM_THREADS=1 (rtsc sets omp_threads = 1)
#SBATCH --mem-per-cpu=3G             # hydrides with dense k grids: 2-4 GB per rank
#SBATCH --time=48:00:00              # one shard = several materials; see README "sizing"
#SBATCH --output=rtsc_logs/%x_%A_%a.out
# #SBATCH --partition=...  --account=...  --qos=...   (site specific)

set -euo pipefail

CONFIG="${RTSC_CONFIG:?set RTSC_CONFIG to the absolute path of config.toml}"
NSHARDS="${RTSC_NSHARDS:-${SLURM_ARRAY_TASK_COUNT:-32}}"
K="${SLURM_ARRAY_TASK_ID:?not running inside a job array}"
if [ -n "${SLURM_ARRAY_TASK_COUNT:-}" ] && [ "$SLURM_ARRAY_TASK_COUNT" != "$NSHARDS" ]; then
  echo "RTSC_NSHARDS=$NSHARDS but the array has $SLURM_ARRAY_TASK_COUNT tasks" >&2; exit 2
fi

# --- software environment (edit for your site) -------------------------------------------------
# Option A: the conda environment from env/environment.yml (QE 7.x + EPW from conda-forge)
#   source "$HOME/miniforge3/etc/profile.d/conda.sh"; conda activate rtsc
# Option B: site modules for QE, plus a Python venv with `pip install './rtsc[dft,mp]'`
#   module load quantum-espresso/7.4 ; source "$HOME/venvs/rtsc/bin/activate"
export OMP_NUM_THREADS=1
# export MP_API_KEY=...       # only if the config has a materials_project source; better: keep
#                             # it in ~/.config and `source` it here, never in the config file
# Scratch: QE writes large *.wfc/_ph0 files under <workdir>/<id>/dfpt. Put [run] workdir on a
# parallel file system with enough space (10-100 GB per material in flight).

echo "rtsc shard $K/$NSHARDS on $(hostname) at $(date -u)  config=$CONFIG  ntasks=${SLURM_NTASKS:-?}"
rtsc run -c "$CONFIG" --shard "$K/$NSHARDS"
