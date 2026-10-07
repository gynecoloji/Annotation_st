#!/bin/bash
#SBATCH --job-name=annot_ov
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G
#SBATCH --time=08:00:00
#SBATCH --output=/easley/scratch/projects/amitra/amitra2016502/annotation_st/_logs/annot_ov.%j.out
#SBATCH --error=/easley/scratch/projects/amitra/amitra2016502/annotation_st/_logs/annot_ov.%j.err
set -euo pipefail
A=/easley/scratch/projects/amitra/amitra2016502/annotation_st
PY=/users/jiwang1/.conda/envs/spatioloji/bin/python
cd "$A"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export NUMBA_NUM_THREADS=${SLURM_CPUS_PER_TASK}
echo "start : $(date -Is)  on $(hostname)  cpus=${SLURM_CPUS_PER_TASK}"
"$PY" -m annotation_st.cli run -c configs/xenium_ov_5k_rerun.yaml
echo "pipeline done: $(date -Is)"
"$PY" - <<'PY'
import pandas as pd
f = "/easley/scratch/projects/amitra/amitra2016502/annotation_st_runs/ov/results/final_annotation_by_cell.tsv.gz"
d = pd.read_csv(f, sep="\t", index_col=0)
print(f"\n{len(d):,} cells annotated; columns: {list(d.columns)}")
for c in ("lineage", "cell_type", "cell_type_fine"):
    if c in d:
        v = d[c].value_counts()
        print(f"\n== {c} ({len(v)} labels)")
        print(pd.DataFrame({"n": v, "pct": (100 * v / len(d)).round(2)}).to_string())
PY
echo "=== regression vs the accepted 2026-09-08 labels ==="
"$PY" - <<'PY'
from annotation_st.compare import compare_labels, format_summary, load_labels
new = load_labels("/easley/scratch/projects/amitra/amitra2016502/annotation_st_runs/ov/results/final_annotation_by_cell.tsv.gz")
ref = load_labels("/easley/scratch/projects/amitra/amitra2016502/annotation_st/validation/ov_full/results/final_annotation_by_cell.tsv.gz")
print(format_summary(compare_labels(new, ref)))
PY
echo "finish: $(date -Is)"
