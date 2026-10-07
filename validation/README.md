# Real-data validation — Xenium 5K human ovarian adenocarcinoma

`run_ov_validation.sbatch` runs the whole package on the OV section into `ov_full/` and
compares the per-cell labels with the accepted run
(`Grant/st_datasets/04-analysis/Analysis_01/01_annot_Xenium_OV/2026-09-08_run/results/final_annotation_by_cell.tsv.gz`).

## Run of 2026-10-03 (job 1256130)

Scavenger partition, 32 CPUs / 120 GB requested, 1 h 34 min, peak memory 79 GB. QC and
lineage ran with 32 threads, every later stage with 16, mirroring the reference jobs.
Config: `configs/xenium_ov_5k.yaml` with the output folders redirected here; the resolved
copy is `ov_full/config_used.yaml`. Figures were on (74 files).

### Comparison with the accepted labels (1,083,808 cells, all shared)

| Column | Agreement | ARI | Cells that differ |
|---|---|---|---|
| `leiden` | 1.0000 | 1.0000 | 0 |
| `lineage` | 1.0000 | 1.0000 | 0 |
| `caf_subtype_cell` | 1.0000 | 1.0000 | 0 |
| `csc_status` | 1.0000 | 1.0000 | 0 |
| `csc_OCSC_z` (numeric) | r = 1.0 | — | max abs diff 0.0 |
| `cell_type` | 0.9963 | 0.9982 | 3,971 |
| `immune_compartment` | 0.9890 | 0.9913 | 11,944 |
| `immune_subtype` | 0.9189 | 0.9940 | 87,931 |
| `cell_type_fine` | 0.9189 | 0.9915 | 87,931 |

Full tables: `ov_full/results/compare/`.

### Reading the result

- **QC, lineage, CAF and CSC are reproduced exactly.** Filter counts (5,879 / 895 / 69 /
  23,132 / 46,628; 1,083,808 kept), 21 Leiden clusters, every lineage and cell-type count,
  the CAF per-cell counts (myCAF 42,068; apCAF 37,595; iCAF 32,159; FB_unassigned 37,160),
  30,061 CSC_like cells and the OCSC z-score are identical. The spatial graph matches too
  (edge cut 38.5 µm, 3.2 M edges; myCAF neighbour log2 ratio 1.86).
- **All 835,106 cells outside the immune lineages carry identical fine labels.** Every one
  of the 87,931 differing cells is in T_NK, Plasma or Macrophage_Mono.
- **The immune differences are a scheme difference, not a porting error.** The accepted
  labels were regenerated on 2026-09-15 with a revised immune hand-off (converged
  multi-round passes, `Immune_mixed`, larger coarse maps; job 1203856, 49 fine labels). The
  package implements the scheme documented on 2026-09-13 (single hop, 43 fine labels), and
  it reproduces that scheme exactly: all 34 immune subtype counts of the earlier final run
  (job 1180895) match to the cell, including the hand-off of 1,001 B and 889 plasma cells.
  Labels present only in the reference: `Immune_mixed`, `cDC2`, `cDC_unspecified`,
  `CD4_naive_memory`, `Monocyte_nonclassical`, `Bcell_spillover_Tumor`.

### Open item

Port the revised hand-off (`01_xenium_immune_subcluster.py --max-rounds 4`, hand-off
entries of the form `label: (target compartment, identity or none)`, the `coarse` maps of
all three compartments) into `stages/immune.py`, `markers.py` and the config validator,
then rerun this job; the immune columns should then agree as the others do.

### Housekeeping

`ov_full/checkpoints/` holds 20 GB of h5ad files. They can be deleted once the result is
no longer needed; `annot-st status` will then report the stages as not run.

## spatioloji_s compatibility check (2026-10-04, job 1258411)

`run_spatioloji_check.sbatch` loads the same OV section with spatioloji_s 0.5.1's own Xenium
loader and checks the bridge (`annotation_st.spatioloji_bridge`) against it. Three minutes,
peak memory 47 GB. Result in `spatioloji_check/report.json`:

| Check | Result |
|---|---|
| Cells from `sj.spatioloji.from_xenium` + `bridge.to_anndata` vs annotation_st's reader | same 1,157,659 cells |
| Genes after separating the 5,028 control features | same 5,001 genes |
| Count matrix | identical |
| Per-cell control sums (negative probes, negative codewords, genomic controls, unassigned and deprecated codewords) | all identical |
| Coordinates | identical |
| Validated labels written into `sp.cell_meta` with `bridge.add_annotations` | identical for all 1,083,808 annotated cells; the 73,851 QC-failed cells are `NA` |
| `bridge.palette` vs the run's `cell_type_colors.csv` | identical |
| `sj.visualization.xenium_plot_spatial(sp, "cell_type_fine", color_dict=palette)` | drawn: `spatioloji_check/spatioloji_xenium_plot_spatial_cell_type_fine.png` |

Because the bridge's input is identical to the reader's, a pipeline run that starts from a
spatioloji object gives the labels validated above.

## Full in-memory run through the spatioloji bridge (2026-10-04, job 1258421)

`run_spatioloji_full.sbatch` annotates the whole OV section in memory,
`sj.spatioloji.from_xenium(...)` then `bridge.annotate(sp, cfg, qc=True)`, and compares the labels
written to `sp.cell_meta` with the validated command-line run (`ov_full/`). Thread counts mirror
that run (32 for QC and lineage, 16 afterwards). 39 minutes, peak memory 59 GB; the command-line
run took 94 minutes and 79 GB with figures and checkpoints.

| Column | Agreement | Cells that differ |
|---|---|---|
| `leiden`, `lineage`, `cell_type`, `cell_type_fine` | 1.0000 | 0 |
| `caf_subtype_cell`, `csc_status` | 1.0000 | 0 |
| `immune_compartment`, `immune_subtype` | 1.0000 | 0 |
| `csc_OCSC_z` (numeric) | r = 1.0 | max abs diff 0.0 |

All 1,083,808 QC-passed cells are shared and identically labelled; the other 73,851 cells of
the object are `NA`. Tables in `spatioloji_full/compare/`, per-cell labels in
`spatioloji_full/annotation_by_cell.tsv.gz`, the spatioloji plot in
`spatioloji_full/spatioloji_cell_type_fine.png`.
