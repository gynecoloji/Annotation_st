# annotation_st

[![CI](https://github.com/gynecoloji/Annotation_st/actions/workflows/ci.yml/badge.svg)](https://github.com/gynecoloji/Annotation_st/actions/workflows/ci.yml)
[![Release](https://img.shields.io/github/v/release/gynecoloji/annotation_st)](https://github.com/gynecoloji/Annotation_st/releases/latest)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue)](pyproject.toml)

Config-driven cell annotation for Xenium spatial transcriptomics. The package turns the
step-wise annotation documented in
`Grant/st_datasets/01-documentation/annotation_methods_01.md` (Xenium 5K human ovarian
adenocarcinoma, analysis `01_annot_Xenium_OV`) into one installable tool: every threshold,
marker set and compartment definition lives in a YAML file, each step is a library function
on an `AnnData`, and the `annot-st` command chains the steps with checkpoints.

```
read Xenium ─► qc ─► lineage ─► caf ┐
                              ├─► csc ├─► merge ─► spatial
                              └─► immune ┘
```

| Stage | What it does | Decision rule (defaults = documented run) |
|---|---|---|
| `qc` | load the input and compute per-cell QC metrics | **no cell or gene is removed**; `qc_filter_summary.tsv` reports what the thresholds (transcripts ≥ 30, genes ≥ 10, control fraction ≤ 0.02, area within 1–99 pct, nucleus present; genes in ≥ 10 cells) would remove. `qc.apply_filters: true` opts in to filtering |
| `lineage` | normalise, PCA(50) → kNN(15) → Leiden 1.0, score 16 lineage sets | cluster mean z > 0.5 else `Unassigned`; proliferation z > 1.0 → `Epithelial_tumor_cycling` |
| `caf` | fibroblasts re-embedded, iCAF / myCAF / apCAF scored within fibroblasts | per cell: best z > 0 and margin ≥ 0.25 else `FB_unassigned`; sub-cluster label also reported |
| `csc` | OCSC score within tumour cells | z > 1.5 and ≥ 2 core markers detected → `CSC_like` |
| `immune` | T/NK → B → myeloid; identity per sub-cluster, then state per parent | core identities at mean z ≥ −0.5, minority / cross at > 0.5, spillover fall-through, `unresolved`; states at > 0.25; B/plasma handed from T/NK to the B compartment |
| `merge` | `cell_type_fine` (finest label per cell) and updated coarse `cell_type` | CAF subtype / CSC status / immune subtype / coarse label |
| `spatial` | Delaunay graph (99th-pct edge cut), neighbourhood enrichment, CSC neighbourhood composition | 100 label permutations, `n_jobs=1` |

Ligand–receptor communication (analysis 02) and the pptx report are not part of this
package, and neither is QC filtering: filter cells and genes upstream (for example with
spatioloji_s), then annotate. The `qc` stage only measures and reports. The one exception is
the shipped `configs/xenium_ov_5k.yaml`, which sets `qc.apply_filters: true` because the run it
reproduces did filter.

## Install

```bash
# on easley the pipeline environment already has every dependency
/users/jiwang1/.conda/envs/spatioloji/bin/python -m pip install -e . --no-deps --no-build-isolation
# elsewhere
pip install -e .            # pulls scanpy, squidpy, anndata, pyyaml, ...
```

## Quick start

```bash
# 1. write a config; its markers section opts in to the default marker library, one line per section
annot-st init-config -o my_run.yaml --xenium-dir /path/to/xenium_outs
#    or start from configs/xenium_ov_5k.yaml, which reproduces the 2026-09-08 OV run exactly
#    (its gene sets are written out; it points at the OV data and the project folders)

# 2. edit paths / parameters, then run everything
annot-st run -c my_run.yaml

# 3. or run stage by stage, re-run one stage, skip figures, smoke-test on a subset
annot-st qc      -c my_run.yaml --subsample 50000
annot-st lineage -c my_run.yaml
annot-st run     -c my_run.yaml --stages immune merge spatial --force
annot-st run     -c my_run.yaml --no-figures
```

An up-to-date stage is skipped, so a job can be requeued safely; see "Trusting a run
folder" below for what "up to date" means. `--n-jobs` defaults to `SLURM_CPUS_PER_TASK` when
set.

### Output layout

| Config key | Default | Content |
|---|---|---|
| `project.outdir` | `annot_run` | `results/` (tables), `figures/` (png + pdf), `logs/`, `qc/`, `params.json`, `analysis.log`, `config_used.yaml`, `stage_state.json` |
| `project.checkpoint_dir` | `<outdir>/checkpoints` | `<prefix>_01_qc.h5ad`, `_02_lineage`, `_03_caf`, `_04_csc`, `_05_immune_<compartment>`, `_06_final` |

`params.json` records every parameter and the gene sets actually used (after panel
intersection) under one key per stage (`qc_filter`, `cluster_annotate`, `fibroblast_caf`,
`cancer_csc`, `immune_subcluster`, `merge`, `spatial_squidpy`). Result-table names match the
original scripts (`leiden_lineage_scores.tsv`, `caf_subtype_by_cell.tsv.gz`,
`immune_<comp>_subcluster_scores.tsv`, `final_annotation_by_cell.tsv.gz`,
`nhood_enrichment_zscore.tsv`, ...).

## Using it with spatioloji_s

`annotation_st` works directly on the `spatioloji` object of the
[spatioloji_s](https://github.com/gynecoloji/spatioloji_s) package, and puts its results where
spatioloji's spatial, communication and plotting functions look for them.

**Annotate an object in memory.** No files, no checkpoints; the labels are the ones the
command-line pipeline gives for the same data and config.

```python
import spatioloji_s as sj
import annotation_st as ast
from annotation_st import spatioloji_bridge as bridge

sp = sj.spatioloji.from_xenium("/path/to/xenium_outs")
cfg = ast.load_config("my_run.yaml")              # or ast.preset_config() for the default library

out = bridge.annotate(sp, cfg)                    # filter with sj.xenium_qc before, if you want filtering
sp.cell_meta[["cell_type", "cell_type_fine", "csc_status"]].head()

# spatioloji analyses use the new columns like any other
sj.visualization.plot_umap(sp, color_by="cell_type_fine", colors=out.palettes["cell_type_fine"])
sj.visualization.xenium_plot_spatial(sp, "cell_type_fine", color_dict=out.palettes["cell_type_fine"])
config = sj.ccc.CCCConfig(group_col="cell_type_fine", layer="log_normalized")
```

What `annotate` does:

- **Columns written to `sp.cell_meta`:** `lineage`, `cell_type`, `cell_type_fine`, `caf_subtype_cell`,
  `caf_subtype_cluster`, `csc_status`, `csc_OCSC_z`, `csc_n_core_detected`, `immune_compartment`,
  `immune_subtype`, `proliferating`, and `annot_leiden` (the lineage clusters; a `leiden` column of
  your own is never touched). `prefix="ast_"` namespaces them; `overwrite=False` refuses to replace.
- **No cell is removed.** `qc=True` only computes the QC metrics; with `cfg.qc.apply_filters`
  also set, cells that fail get the label `NA`.
- **Label columns are categorical**, with categories in legend-hierarchy order, so spatioloji
  legends read tumour, CAF, stroma, T, NK, B, myeloid, QC.
- **Colours:** `out.palettes[column]`, or `bridge.palette(sp, column, cfg)` at any time, is an
  ordered `{label: colour}` mapping with exactly the categories of that column. It is correct
  for every way spatioloji plots take colours: by label (`colors=`, `color_map=`, `color_dict=`)
  and by position (`palette=` in the embedding plots). Colours are the same as in annotation_st's own figures.
- **Input details handled for you:** spatioloji's Xenium loader keeps control probes and
  codewords in the expression matrix. `bridge.to_anndata(sp)` separates them, sums them into the
  QC columns, and stores coordinates as (x, y). `sp.to_anndata()` stores four coordinate
  columns, which would break a spatial neighbour graph.
- `log_layer="log_normalized"` reuses a layer you normalised with spatioloji instead of
  normalising the raw counts again.

**Use a saved object as pipeline input.** Give `project.spatioloji` a pickle (`sp.to_pickle`) or a
components folder (`sp.save_components`) instead of `project.xenium_dir`:

```yaml
project: {spatioloji: sample.pkl, outdir: annot_run}
```

```bash
annot-st run -c my_run.yaml
annot-st to-spatioloji -c my_run.yaml -o sample_annotated.pkl    # labels back into the object
```

`bridge.add_annotations(sp, "annot_run/results/final_annotation_by_cell.tsv.gz")` does the last
step from Python. The same in-memory entry point exists for a plain AnnData:
`annotation_st.api.annotate_adata(adata, cfg, qc=True, spatial=True)`.

## Trusting a run folder

**Check the config before spending compute.** `annot-st validate -c my_run.yaml` (also run
automatically at the start of `annot-st run`) reports every problem at once: wrong value
types and ranges, immune-compartment cross-references (a core identity or state parent
without a marker set, a hand-off to a compartment that does not exist), lineage names that
no lineage set defines (they would silently select zero cells), unreadable marker files and
a missing Xenium folder. `--stages` restricts the check to what those stages need.

**Stale results are detected, not reused.** `stage_state.json` records for every stage the
inputs that determine its result (parameters and the resolved gene lists, so editing a
marker file counts), a run id, and the run ids of the upstream stages it consumed. A stage
is up to date only if its outputs exist, its inputs are unchanged and no upstream stage was
re-run since. Otherwise `annot-st run` re-runs it and everything downstream:

```bash
annot-st status -c my_run.yaml          # up to date / stale (with the reason) / not run
# stage    status      reason
# csc      stale       inputs changed: csc.z_thresh
# merge    stale       upstream stage 'csc' is out of date (inputs changed: csc.z_thresh)
annot-st run -c my_run.yaml             # re-runs csc, merge, spatial only
```

Running a single stage on top of an out-of-date upstream stage is refused
(`--allow-stale` overrides). A stage that skips itself (too few cells) is recorded as
such, and `merge` only uses refinements that the state file lists as current, so a
checkpoint left over from an earlier configuration cannot leak into the final labels.
Settings that only affect figures or speed (`project.figures`, `project.n_jobs`,
`*.umap_cells`, `lineage.skip_umap`) do not invalidate anything. Leiden partitions can
shift with the thread count; `n_jobs` is recorded per stage for that reason.

**The run folder keeps its own configuration.** `config_used.yaml` is the resolved config
of the latest invocation: absolute paths, CLI overrides applied, marker files inlined. It
can be passed back to `annot-st run -c` to reproduce the run without the original files.

**Compare against an accepted annotation.**

```bash
annot-st compare --new run/results/final_annotation_by_cell.tsv.gz --ref reference.tsv.gz -o cmp/
annot-st compare -c my_run.yaml --ref reference.tsv.gz --min-agreement 0.99   # exit 1 below the bar
```

Cells are matched by id. For every label column both tables share the report gives the
agreement, the adjusted Rand index (robust to renumbered clusters), per-label recall /
precision and the confusion table; numeric columns get a correlation and the largest
difference. Tables can be `.tsv[.gz]`, `.csv`, `.parquet` or `.h5ad` (`obs`).
`validation/run_ov_validation.sbatch` runs the whole pipeline on the OV section and
compares it with the accepted run; see `validation/README.md` for the result.

## Colours: one scheme for every figure and every dataset

Every UMAP and tissue map that shows cell-type labels uses one label -> colour map: the
rough-type maps (`umap_cell_type`, `spatial_cell_type`), the CAF, CSC and immune maps
(`umap_fb_caf_subtype`, `spatial_fb_caf_subtype`, `umap_tumor_csc`, `spatial_tumor_csc`,
`umap_<compartment>_subtype`, `spatial_<compartment>_subtype`), the CSC neighbour bar
charts and the final map. A label therefore has the same colour wherever it appears. The
colours are also stored in each checkpoint (`uns['<column>_colors']`). Panels that show
cluster numbers rather than cell types keep scanpy's default colours.

## Final annotation map: legend hierarchy and colours

`figures/spatial_cell_type_fine` shows every cell in the colour of its fine type. Its legend
is hierarchical: bold group headers (tumour epithelium, fibroblast / CAF, other stroma, T
cell, NK cell, B lineage, three myeloid groups, unresolved / QC), each followed by its
labels with cell counts.

Colours are computed, not hand-picked, in three tiers:

| Tier | Labels | Colours |
|---|---|---|
| grey | unassigned, spillover and unresolved classes (the `Unresolved / QC` group) | greys with next to no hue; the classes seen most often get the best-separated greys |
| bright | cell types with a specific fine type (`Treg`, `myCAF`, `TREM2_TAM`, `Plasma`, ...) | light, saturated colours |
| deep | lineage or identity labels without a specific fine type (`CD8_T_unspecified`, `Macrophage_unpolarized`, `FB_unassigned`, `T_NK`, ...; list in `plot.generic_labels`) | dark, still clearly coloured tones |

Inside a legend group the labels are placed by farthest-point selection in OKLab colour
space, so subtypes that sit together in the tissue are clearly different (at least 15
OKLab units x 100 apart in every cell-type group of the default hierarchy; the tests
enforce this). Rough types (the lineage names) are also kept at least that far from each
other, and no colour is used twice. A label's colour depends only on the config, never on
which labels happen to be present. The deep tier exists because bright colours alone
cannot keep the 13 to 14 labels of the largest groups apart; setting
`plot.generic_labels: []` makes every cell type bright at the cost of that separation.

**Every colour is recorded.** `results/cell_type_colors.csv` has one row per label:

| Column | Meaning |
|---|---|
| `label`, `hex` | cell-type label and its colour |
| `level` | `fine`, `coarse` (rough type) or `fine+coarse` |
| `legend_group` | the big group in the legend hierarchy |
| `coarse_type` | for a fine type, the rough type most of its cells carry |
| `n_cells_fine`, `n_cells_coarse` | cells with that label at each level |

The same colours are stored in the checkpoints (`uns['cell_type_colors']`,
`uns['cell_type_fine_colors']`), so scanpy and other viewers pick them up.

**Reproduce or change colours** through the `plot` section; nothing has to be recomputed:

```yaml
plot:
  palette_file: cell_type_colors.csv     # a saved table (label, hex): those labels keep their colours
  colors: {Treg: "#d62728"}              # single overrides; other labels are placed around them
  generic_labels: [T_NK, CD8_T_unspecified, ...]   # labels drawn in deep tones instead of bright ones
  legend_groups:                         # group title -> labels, both in display order
    Tumour epithelium: [Tumor_nonCSC, CSC_like]
    ...
```

```bash
annot-st figures -c my_run.yaml          # redraw all label maps from the checkpoints (minutes, no recomputation)
annot-st figures -c my_run.yaml --which cell_type_fine immune
```

### Using the scheme on other datasets

The scheme is a rule set, not a fixed list, so it carries over:

- **Another section run through this package** gets the same colours automatically: a
  label's colour depends only on the config. To freeze the colours for a whole project
  regardless of later config edits, export them once and point every dataset at the file:

  ```bash
  annot-st palette -c my_run.yaml -o project_palette.csv      # label, hex, tier, legend_group
  ```
  ```yaml
  plot: {palette_file: project_palette.csv}
  ```
  `configs/cell_type_palette.csv` is this export for the default scheme.
- **New labels** (another tissue, another panel) need no setup. Labels the hierarchy does
  not list go under "Other" and are coloured by the same rules: bright and mutually
  distinct, deep tones if the name ends in `_unspecified` / `_unpolarized` / `_unassigned`,
  grey if it is an unassigned / spillover / unresolved class. Adding them to
  `plot.legend_groups` gives them a proper legend group.
- **A dataset annotated elsewhere** (any h5ad with label columns):

  ```bash
  annot-st palette --palette-file project_palette.csv --h5ad liver.h5ad --keys cell_type subtype \
      --write-h5ad -o liver_colors.csv --map figures/liver_map --map-key subtype
  ```
  Known labels keep their colour, new ones are added to the scheme; `--write-h5ad` stores
  `uns['<key>_colors']` in the file (scanpy, squidpy and cellxgene then use them), `--map`
  draws the hierarchical-legend tissue map. From Python:
  `annotation_st.pipeline.label_palette(cfg, labels)`,
  `annotation_st.plotting.apply_palette(adata, key, palette)` and
  `annotation_st.plotting.build_fine_map(xy, labels, palette)`.

Pointing `palette_file` at a copy of a run's `cell_type_colors.csv` gives exactly that run's
colours in any other run; the project's older `fine_label_colors.tsv` can be used the same
way. Labels that no group lists appear under "Other" with their own distinct colours.

## Configuration

One YAML file with sections `project`, `qc`, `lineage`, `caf`, `csc`, `immune`,
`spatial`, `plot`, `markers`. A file may be partial: missing parameters take their defaults,
unknown keys raise an error.

### Marker sets are opt-in

Nothing is scored unless the `markers` section lists it. A config without that section has
no lineage sets, no CSC set, no immune compartments. Each set is written in one of three ways:

| You write | Meaning |
|---|---|
| `default` | the set of that name from the package's default library (every gene referenced, see below) |
| `[GENE1, GENE2]` | your own gene list |
| `path/to/file.txt` | your own list, one gene per line (relative to the config file) |

```yaml
markers:
  lineage:                              # rough types: only these four are scored
    Epithelial_tumor: default
    Fibroblast: default
    T_NK: default
    Hepatocyte: [ALB, APOA1, TTR]
  proliferation: default
  csc: default                          # leave the line out and the CSC stage is skipped
  csc_core: default
  compartments:                         # fine types
    TNK: default                        # a whole compartment from the library
    Myeloid:                            # or set by set, state by state
      lineages: [Macrophage_Mono]
      core: [Macrophage]
      unresolved: Myeloid_unresolved
      unspecified: {Macrophage: Macrophage_unpolarized}
      level1: {Macrophage: default, Monocyte: default, MyCell: markers/mycell.txt}
      level2: {Macrophage: {TREM2_TAM: default, M1_macrophage: default}}
```

Shortcuts: a section can be `default` as a whole (`lineage: default`, `compartments: default`,
`caf: default`, `spillover: default`), and `markers: default` or `markers: {preset: default, csc: [...]}`
opts in to the whole library, with your own keys replacing the corresponding part.
`annot-st init-config` writes the one-line-per-section form; `--markers explicit` writes the
gene lists out and `--markers none` leaves everything to you.

What happens when something is not opted in:

- no lineage sets: the lineage stage cannot run, and `annot-st run` / `validate` say so before any compute;
- no CAF sets, no CSC set, or no immune compartments: that stage is skipped and recorded as skipped,
  and the merged labels simply keep the rough type for those cells.

Where every set came from (default library, inline, or which file) is written to
`results/marker_sources.json`; the gene lists actually used, after panel intersection, are in
`results/*_marker_sets_used.json` and `config_used.yaml`.

### The default library and its references

```bash
annot-st markers                              # sets, gene counts, reference counts
annot-st markers --set TNK.level2.CD4_T/Treg  # each gene with its references
annot-st markers --references refs.tsv --dropped dropped.tsv
```

The library holds 69 sets and 550 set-gene pairs: 16 lineages, proliferation, the
ovarian cancer stem-cell set, the three spillover sets, identity and state sets of the T/NK, B and
myeloid compartments, and iCAF / myCAF / apCAF. **Every gene of every set has at least one
reference** (1371 references, 414 distinct PMIDs; a test enforces it). A reference is one of:

| Evidence | What it is |
|---|---|
| `curated` | a CellMarker 2.0 record (Hu et al., *Nucleic Acids Res* 2023, PMID 36300619) for that cell type and gene; the record's own PMID is cited |
| `primary_text` | the gene named with the cell type in a sentence of one of the set's primary references (full text from PubMed Central); the sentence is stored |
| `literature_text` | the same kind of sentence from another open-access paper, read and accepted by hand (`tools/reviewed_literature_evidence.tsv`); the sentence is stored |

Every PMID was resolved in PubMed when the table was built. Candidate genes for which none of the
three was found are **not** in the library; `annot-st markers --dropped` lists all 91 of them
with the reason. Examples: `WDFY4` is not in the mature-DC set because the literature names it for
cDC1, and `IDO1` is not in the M1 set because the sentences found point the other way. The table
is rebuilt with `tools/build_marker_references.py` (needs network access).

Limits to keep in mind: a `curated` reference means a curated database links that paper to the
gene for that cell type; the papers behind those records were not re-read here. `literature_text`
sentences include review articles. A reference shows that a gene is used as a marker for the type,
not that it is specific on your panel or in your tissue.

The gene sets of the validated OV run are kept unchanged in `configs/xenium_ov_5k.yaml` (written
out, so it does not depend on the library); 14 of its set-gene pairs outside the CAF lists are not in
the default library.

### Using your own cancer stem cell markers

The CSC stage accepts a custom list in either of two ways (a file wins over an inline list,
which wins over the built-in OCSC set):

```yaml
csc:
  marker_file: my_stemness_genes.txt        # one gene per line, '#' comments; relative to the config file
  core_marker_file: my_core_genes.txt       # optional: genes for the ">= min_core detected" rule
  z_thresh: 1.5
  min_core: 2
# or inline
markers:
  csc: [CD44, PROM1, ALDH1A1, CD24, LGR5]
  csc_core: [CD44, PROM1]
```

If no core list is given, the built-in core markers that occur in your set are used; if none
of them does, every gene of your set counts as core. The set actually used, its source and the
genes dropped because they are not on the panel are recorded in
`results/csc_marker_sets_used.json` and in `params.json["cancer_csc"]["marker_source"]`. From
Python: `call_csc(tum, cfg.csc, csc_genes=[...], core_genes=[...])`.

## Library use

```python
import scanpy as sc
from annotation_st import load_config
from annotation_st.stages.lineage import annotate_lineage
from annotation_st.stages.caf import annotate_caf
from annotation_st.stages.immune import run_immune

cfg = load_config("my_run.yaml")
adata = sc.read_h5ad("xenium_ov_01_qc.h5ad")               # raw counts
res = annotate_lineage(adata, cfg.lineage, cfg.markers)    # adds obs['leiden','lineage','cell_type']
res.means                                                   # cluster × lineage mean-z table

fb = adata[adata.obs["lineage"] == "Fibroblast"].copy()
caf = annotate_caf(fb, cfg.caf, {"iCAF": [...], "myCAF": [...], "apCAF": [...]})

immune = run_immune(adata, cfg.immune, cfg.markers)        # {'TNK': ImmuneResult, 'Bcell': ..., 'Myeloid': ...}
```

Stage functions never write files or draw figures; `annotation_st.pipeline` adds the
checkpoints, tables, figures and `params.json` around them.

## Layout

```
src/annotation_st/
  config.py     dataclasses + YAML loading (defaults = documented OV run)
  markers.py    gene sets of the OV run (library candidates), marker-file readers, panel intersection
  library.py    the default marker library and its reference table
  data/marker_library/   library sets, references, dropped candidates, build metadata
  io.py         Xenium reader, run-folder layout, params.json / logs
  scoring.py    score_sets, cluster_assign, per_cell_assign, embed
  state.py      per-stage run records and staleness checks (stage_state.json)
  compare.py    label comparison against a reference table
  plotting.py   png+pdf figure helpers, computed colour palette, hierarchical legend
  stages/       qc, lineage, caf, csc, immune, merge, spatial (pure functions)
  pipeline.py   stage runners with checkpoints
  api.py        annotate_adata: the whole annotation on an AnnData in memory
  spatioloji_bridge.py   spatioloji_s objects in, labels and colours back
  cli.py        annot-st
configs/xenium_ov_5k.yaml   full config for the OV section
validation/                 real-data validation job and its report
tools/                      build_marker_references.py, find_literature_evidence.py, reviewed sentences
tests/                      synthetic Xenium-like data with planted populations
docs/superpowers/           design spec and implementation plan
```

## Tests

```bash
MPLBACKEND=Agg /users/jiwang1/.conda/envs/spatioloji/bin/python -m pytest tests -q            # ~3 min
MPLBACKEND=Agg /users/jiwang1/.conda/envs/spatioloji/bin/python -m pytest tests -q -m slow    # + end-to-end run with figures
```

## Mapping to the original scripts

| Original (`02-scripts/Script_01/`) | Package |
|---|---|
| `01_xenium_common.py` | `markers.py`, `config.py`, `io.py`, `scoring.py`, `plotting.py` |
| `01_xenium_qc_filter.py` | `stages/qc.py` + `pipeline.run_qc` |
| `01_xenium_cluster_annotate.py` | `stages/lineage.py` + `pipeline.run_lineage` |
| `01_xenium_fibroblast_caf.py` | `stages/caf.py` + `pipeline.run_caf` |
| `01_xenium_cancer_csc.py` | `stages/csc.py` + `pipeline.run_csc` |
| `01_xenium_immune_subcluster.py` | `stages/immune.py` + `pipeline.run_immune_stage` |
| `01_xenium_spatial_squidpy.py` | `stages/merge.py`, `stages/spatial.py` + `pipeline.run_merge`, `run_spatial` |
| `01_run_pipeline.sbatch` | `annot-st run -c config.yaml` |
