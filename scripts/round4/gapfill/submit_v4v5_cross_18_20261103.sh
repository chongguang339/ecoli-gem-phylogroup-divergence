#!/bin/bash
#SBATCH --job-name=v4v5_cross18
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=52
#SBATCH --mem=96G
#SBATCH --time=18:00:00
#SBATCH --output=/home/u2024110181382/Project/04.pipeline/revision5_v4v5_cross_18_20261103/logs/cross18_%j.out
set -euo pipefail
ROOT=/home/u2024110181382/Project/04.pipeline/revision5_v4v5_cross_18_20261103
export PATH=/home/u2024110181382/miniconda3/envs/gapseq/bin:$PATH
export OPENBLAS_NUM_THREADS=4
sha256sum "$ROOT/scripts/run_v4v5_cross_18_20261103.py" "$ROOT/scripts/submit_v4v5_cross_18_20261103.sh" > "$ROOT/logs/script_sha256.txt"
python3 "$ROOT/scripts/run_v4v5_cross_18_20261103.py"
