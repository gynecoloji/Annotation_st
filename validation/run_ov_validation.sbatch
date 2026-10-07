#!/usr/bin/env bash
#SBATCH --job-name=annotst_validate
#SBATCH --account=2016502
#SBATCH --partition=scavenger
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --time=08:00:00
#SBATCH --requeue
#SBATCH --output=validation/ov_full/slurm_%j.out
# ----------------------------------------------------------------------------
# Real-data validation of annotation_st: run the whole pipeline on the Xenium 5K OV
# section into validation/ov_full/ and compare the per-cell labels with the accepted
# run (Grant/st_datasets/04-analysis/Analysis_01/01_annot_Xenium_OV/2026-09-08_run).
# Submit from the package root:   sbatch validation/run_ov_validation.sbatch
# Safe to requeue: finished stages are skipped (stage_state.json).
#
# Thread counts mirror the reference jobs because Leiden partitions can shift with the
# number of threads: QC + lineage came from a 32-CPU job (1178155), every later stage
# from 16-CPU jobs.
# ----------------------------------------------------------------------------
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
BIN=/users/jiwang1/.conda/envs/spatioloji/bin
CFG=configs/xenium_ov_5k.yaml
OUT=$PWD/validation/ov_full
REF=/easley/scratch/projects/amitra/amitra2016502/Grant/st_datasets/04-analysis/Analysis_01/01_annot_Xenium_OV/2026-09-08_run/results/final_annotation_by_cell.tsv.gz
export MPLBACKEND=Agg TF_CPP_MIN_LOG_LEVEL=3
threads() { export OMP_NUM_THREADS=$1 MKL_NUM_THREADS=$1 OPENBLAS_NUM_THREADS=$1 NUMBA_NUM_THREADS=$1; }
COMMON=(-c "$CFG" --outdir "$OUT" --checkpoint-dir "$OUT/checkpoints")

echo "[$(date)] job $SLURM_JOB_ID on $(hostname)"
$BIN/annot-st validate "${COMMON[@]}"
threads 32
$BIN/annot-st run "${COMMON[@]}" --n-jobs 32 --stages qc lineage
threads 16
$BIN/annot-st run "${COMMON[@]}" --n-jobs 16 --stages caf csc immune merge spatial
$BIN/annot-st status "${COMMON[@]}"
$BIN/annot-st compare --new "$OUT/results/final_annotation_by_cell.tsv.gz" --ref "$REF" -o "$OUT/results/compare"
echo "[$(date)] done"
