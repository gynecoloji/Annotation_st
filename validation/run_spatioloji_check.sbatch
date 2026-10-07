#!/usr/bin/env bash
#SBATCH --job-name=annotst_sjcheck
#SBATCH --account=2016502
#SBATCH --partition=scavenger
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=02:00:00
#SBATCH --output=validation/spatioloji_check/slurm_%j.out
# Real-data check of the spatioloji_s bridge. Submit from the package root:
#   sbatch validation/run_spatioloji_check.sbatch
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
export MPLBACKEND=Agg TF_CPP_MIN_LOG_LEVEL=3
/users/jiwang1/.conda/envs/spatioloji/bin/python validation/spatioloji_check.py
