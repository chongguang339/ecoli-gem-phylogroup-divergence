#!/bin/bash
#SBATCH --job-name=v1v2_matched_20
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --output=/home/u2024110181382/Project/04.pipeline/revision5_v1v2_matched_20_20261102/logs/matched20_%j.out
set -euo pipefail
ROOT=/home/u2024110181382/Project/04.pipeline/revision5_v1v2_matched_20_20261102
mkdir -p "$ROOT/logs" "$ROOT/scripts"
export PATH=/home/u2024110181382/miniconda3/envs/gapseq/bin:$PATH
export OPENBLAS_NUM_THREADS=4
sha256sum "$ROOT/scripts/run_v1v2_matched_20_20261102.py" "$ROOT/scripts/submit_v1v2_matched_20_20261102.sh" > "$ROOT/logs/script_sha256.txt"
python3 "$ROOT/scripts/run_v1v2_matched_20_20261102.py"
