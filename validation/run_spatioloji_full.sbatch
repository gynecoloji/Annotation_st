#!/usr/bin/env bash
#SBATCH --job-name=annotst_sjfull
#SBATCH --account=2016502
#SBATCH --partition=scavenger
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=220G
#SBATCH --time=06:00:00
#SBATCH --output=validation/spatioloji_full/slurm_%j.out
# Full in-memory annotation of the OV section through the spatioloji_s bridge, compared with the
# validated run. Submit from the package root:  sbatch validation/run_spatioloji_full.sbatch
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
export MPLBACKEND=Agg TF_CPP_MIN_LOG_LEVEL=3
export OMP_NUM_THREADS=32 MKL_NUM_THREADS=32 OPENBLAS_NUM_THREADS=32 NUMBA_NUM_THREADS=32
/users/jiwang1/.conda/envs/spatioloji/bin/python validation/spatioloji_full_run.py
