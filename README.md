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
| `lineage` | normalise, PCA(50) → kNN(15) → Leiden 1.0, score the lineage sets | cluster mean z > 0.5 else `Unassigned`; proliferation z > 1.0 → `Epithelial_tumor_cycling` |
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

## Contents

- [Install](#install)
- [Quick start](#quick-start)
- [Step by step](#step-by-step)
  1. [Input data](#1-input-data)
  2. [Write the config](#2-write-the-config)
  3. [Choose marker sets](#3-choose-marker-sets)
  4. [Validate](#4-validate)
  5. [Run](#5-run)
  6. [Read the outputs](#6-read-the-outputs)
  7. [Figures and colours](#7-figures-and-colours)
  8. [Compare with a reference annotation](#8-compare-with-a-reference-annotation)
  9. [Use it from Python](#9-use-it-from-python)
- [Command reference](#command-reference)
- [Troubleshooting](#troubleshooting)
- [How a run folder stays trustworthy](#how-a-run-folder-stays-trustworthy)
- [The default marker library and its references](#the-default-marker-library-and-its-references)
- [Package layout](#package-layout), [Tests](#tests), [Mapping to the original scripts](#mapping-to-the-original-scripts)

## Install

Python 3.10 or newer. The dependencies (scanpy, squidpy, anndata, pyyaml, matplotlib,
pyarrow, h5py, igraph, leidenalg) are declared in `pyproject.toml`.

```bash
git clone https://github.com/gynecoloji/Annotation_st.git
cd Annotation_st
pip install -e .                      # pulls every dependency
pip install -e .[spatioloji]          # also spatioloji_s, for the bridge
pip install -e .[test]                # also pytest

# on easley the pipeline environment already has every dependency
/users/jiwang1/.conda/envs/spatioloji/bin/python -m pip install -e . --no-deps --no-build-isolation

annot-st --version                    # check the install
```

## Quick start

```bash
# 1. write a config that points at your Xenium folder; its markers section opts in to the
#    default marker library (one editable line per section)
annot-st init-config -o my_run.yaml --xenium-dir /path/to/xenium_outs

# 2. edit paths and parameters (see "Write the config"), then check it
annot-st validate -c my_run.yaml

# 3. smoke-test on 50 000 random cells into a scratch folder (minutes)
annot-st run -c my_run.yaml --subsample 50000 --outdir smoke --no-figures

# 4. the real run (about 1.5 h for a full 5K section; see "Run" for a SLURM script)
annot-st run -c my_run.yaml

# 5. look at the results
annot-st status -c my_run.yaml
column -t annot_run/results/final_cell_type_counts.tsv | head
open annot_run/figures/spatial_cell_type_fine.png
```

To reproduce the documented ovarian cancer run instead, start from
`configs/xenium_ov_5k.yaml` (its gene sets are written out, it points at the OV data and it
opts in to QC filtering as that run did).

## Step by step

### 1. Input data

The pipeline accepts three kinds of input; all of them end up as an `AnnData` with raw counts.

| Input | How to point at it | Requirements |
|---|---|---|
| 10x Xenium output folder | `project.xenium_dir` (or `--xenium-dir`) | `cell_feature_matrix.h5` and `cells.parquet`. Every column of `cells.parquet` is kept in `obs`; centroids (µm) go to `obsm['spatial']`; control probes and codewords are summed into per-cell QC columns and removed from the matrix |
| saved spatioloji_s object | `project.spatioloji` | a pickle (`sp.to_pickle`) or a components folder (`sp.save_components`). Only one of `xenium_dir` / `spatioloji` may be set |
| object in memory | Python API, see [Use it from Python](#9-use-it-from-python) | `AnnData` with raw counts in `X` and, for the spatial stage, `obsm['spatial']`; or a `spatioloji` object |

The package does not filter. If cells or genes should be removed, do that before (for example
`sj.xenium_qc` in spatioloji_s) and hand the filtered object over. `qc.subsample` /
`--subsample N` keeps N random cells and exists for smoke tests only.

### 2. Write the config

```bash
annot-st init-config -o my_run.yaml --xenium-dir /path/to/xenium_outs    # library sets by name
annot-st init-config -o my_run.yaml --markers explicit                   # gene lists written out
annot-st init-config -o my_run.yaml --markers none                       # no marker sets at all
```

One YAML file with sections `project`, `qc`, `lineage`, `caf`, `csc`, `immune`, `spatial`,
`plot`, `markers`. A file may be partial: missing parameters take the defaults listed below,
unknown keys raise a `ConfigError`. Relative paths (marker files, palette file) are resolved
against the folder of the config file. Every `project` path, `project.seed` and
`qc.subsample` can be overridden on the command line (`--outdir`, `--xenium-dir`,
`--prefix`, `--seed`, `--subsample`, ...), and `--n-jobs` defaults to `SLURM_CPUS_PER_TASK`
when that variable is set.

**project**

| Key | Default | Meaning |
|---|---|---|
| `name` | `sample` | free text, written to `params.json` |
| `xenium_dir` | `null` | Xenium output folder (input) |
| `spatioloji` | `null` | saved spatioloji_s object (alternative input) |
| `outdir` | `annot_run` | run folder: results, figures, logs, state |
| `checkpoint_dir` | `null` → `<outdir>/checkpoints` | where the per-stage h5ad files go (can be a scratch disk) |
| `prefix` | `sample` | prefix of the checkpoint files |
| `seed` | `0` | PCA / Leiden / UMAP / permutation seed |
| `n_jobs` | `8` | threads for scanpy; recorded per stage because Leiden partitions can shift with it |
| `figures` | `true` | draw figures (`--no-figures` switches them off for one invocation) |

**qc**

| Key | Default | Meaning |
|---|---|---|
| `subsample` | `0` | keep N random cells (0 = all); smoke tests |
| `min_transcripts`, `min_genes` | `30`, `10` | per-cell thresholds |
| `max_control_frac` | `0.02` | control-probe fraction of total counts |
| `area_pct` | `[1, 99]` | cell-area percentile window |
| `require_nucleus` | `true` | flag cells without a nucleus |
| `min_cells_per_gene` | `10` | gene threshold |
| `apply_filters` | `false` | **only report** what the thresholds would remove (default), or remove it |

**lineage**

| Key | Default | Meaning |
|---|---|---|
| `n_pcs`, `n_neighbors`, `resolution` | `50`, `15`, `1.0` | embedding and Leiden on all cells |
| `min_z` | `0.5` | cluster mean z a lineage set needs to label the cluster |
| `cycling_z`, `cycling_lineage`, `cycling_label` | `1.0`, `Epithelial_tumor`, `Epithelial_tumor_cycling` | proliferation call inside one lineage |
| `unassigned` | `Unassigned` | label of clusters below `min_z` |
| `umap_cells`, `skip_umap` | `300000`, `false` | UMAP on a random subset (0 = all), or none |

**caf**

| Key | Default | Meaning |
|---|---|---|
| `include_lineages` | `[Fibroblast]` | rough types re-embedded and scored |
| `n_pcs`, `n_neighbors`, `resolution` | `30`, `15`, `0.5` | sub-clustering |
| `min_z`, `margin` | `0.0`, `0.25` | per-cell call: best z above `min_z` and ahead of the runner-up by `margin` |
| `min_cells` | `50` | fewer cells → the stage is skipped and recorded as skipped |
| `unassigned` | `FB_unassigned` | label when no set wins |
| `marker_files` | `{}` | `{iCAF: file.txt, ...}`; an alternative to `markers.caf` |

**csc**

| Key | Default | Meaning |
|---|---|---|
| `tumor_lineages` | `[Epithelial_tumor]` | cells scored |
| `z_thresh`, `min_core` | `1.5`, `2` | CSC call: score z above threshold and at least this many core markers detected |
| `n_pcs`, `n_neighbors`, `resolution`, `umap_cells`, `min_cells` | `30`, `15`, `0.4`, `250000`, `50` | tumour sub-embedding (for figures and tables) |
| `positive_label`, `negative_label` | `CSC_like`, `Tumor_nonCSC` | labels |
| `marker_file`, `core_marker_file` | `null` | your own gene lists; see [custom CSC markers](#your-own-cancer-stem-cell-markers) |

**immune**

| Key | Default | Meaning |
|---|---|---|
| `compartments` | `null` | which compartments to run (`null` = every compartment in `markers.compartments`, in the order they hand cells on: T/NK → B → myeloid) |
| `n_pcs`, `n_neighbors`, `resolution` | `30`, `15`, `1.5` | level-1 sub-clustering per compartment |
| `min_z`, `low_identity` | `0.5`, `-0.5` | identity thresholds: core identities qualify at mean z ≥ `low_identity`, minority / cross identities at > `min_z` |
| `resolution2`, `min_z2` | `1.0`, `0.25` | level-2 state sub-clustering inside each parent identity |
| `spill_flag` | `1.0` | per-cell `immune_spillover` flag threshold |
| `umap_cells`, `min_cells`, `min_cells_level2` | `200000`, `100`, `50` | subset for UMAP; minimum cells for a compartment / a level-2 parent |

**spatial**

| Key | Default | Meaning |
|---|---|---|
| `n_perms` | `100` | label permutations for neighbourhood enrichment |
| `radius`, `edge_cut_pct` | `0.0`, `99.0` | Delaunay graph with the longest 1 % of edges removed, or a fixed-radius graph when `radius > 0` |
| `fine_key` | `cell_type_fine` | label column analysed |
| `composition_groups` | `[CSC_like, Tumor_nonCSC]` | groups whose neighbourhood composition is tabulated |

**plot** is described under [Figures and colours](#7-figures-and-colours); **markers** next.

### 3. Choose marker sets

Nothing is scored unless the `markers` section lists it. A config without that section has
no lineage sets, no CSC set, no immune compartments. Each set is written in one of three ways:

| You write | Meaning |
|---|---|
| `default` | the set of that name from the package's default library (every gene referenced, see [below](#the-default-marker-library-and-its-references)) |
| `[GENE1, GENE2]` | your own gene list |
| `path/to/file.txt` | your own list, one gene per line, `#` comments allowed (relative to the config file) |

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
  spillover: default                    # tumour / fibroblast / endothelial contamination sets for the immune stage
  caf: default                          # iCAF, myCAF, apCAF
  compartments:                         # fine types
    TNK: default                        # a whole compartment from the library
    Myeloid:                            # or set by set, state by state
      lineages: [Macrophage_Mono]       # rough types that feed this compartment
      core: [Macrophage]                # identities accepted at the low threshold
      unresolved: Myeloid_unresolved    # label when nothing qualifies
      unspecified: {Macrophage: Macrophage_unpolarized}   # label when no level-2 state wins
      level1: {Macrophage: default, Monocyte: default, MyCell: markers/mycell.txt}
      level2: {Macrophage: {TREM2_TAM: default, M1_macrophage: default}}
```

Shortcuts: a section can be `default` as a whole (`lineage: default`, `compartments: default`,
`caf: default`, `spillover: default`), and `markers: default` or `markers: {preset: default, csc: [...]}`
opts in to the whole library, with your own keys replacing the corresponding part.

The default library offers 16 lineages (`Epithelial_tumor`, `Fibroblast`, `Pericyte_SMC`,
`Endothelial`, `Lymphatic_EC`, `T_NK`, `B_cell`, `Plasma`, `Macrophage_Mono`, `Dendritic`,
`Neutrophil`, `Mast`, `Mesothelial`, `Ovarian_stroma`, `Adipocyte`, `Schwann`), a
proliferation set, the ovarian CSC set with its core, three spillover sets, three CAF sets and
three immune compartments (`TNK`, `Bcell`, `Myeloid`) with their level-1 identities and
level-2 states. Lineage names are shared vocabulary: `caf.include_lineages`,
`csc.tumor_lineages` and each compartment's `lineages` must name lineages that
`markers.lineage` defines, and `annot-st validate` checks that.

```bash
annot-st markers                                 # every set with gene and reference counts
annot-st markers --set lineage/Fibroblast        # the genes of one set with their references
annot-st markers --set TNK.level2.CD4_T/Treg
```

What happens when something is not opted in:

- no lineage sets: the lineage stage cannot run, and `annot-st run` / `validate` say so before any compute;
- no CAF sets, no CSC set, or no immune compartments: that stage is skipped and recorded as skipped,
  and the merged labels simply keep the rough type for those cells.

Where every set came from (default library, inline, or which file) is written to
`results/marker_sources.json`; the gene lists actually used, after intersection with the panel,
are in `results/*_marker_sets_used.json` and `config_used.yaml`.

#### Your own cancer stem cell markers

The CSC stage accepts a custom list in either of two ways (a file wins over an inline list):

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

If no core list is given, the library's core markers that occur in your set are used; if none
of them does, every gene of your set counts as core. The set actually used, its source and the
genes dropped because they are not on the panel are recorded in
`results/csc_marker_sets_used.json` and in `params.json["cancer_csc"]["marker_source"]`.

### 4. Validate

```bash
annot-st validate -c my_run.yaml                       # the whole pipeline
annot-st validate -c my_run.yaml --stages qc lineage   # only what these stages need
```

The check also runs at the start of every `annot-st run`. It reports every problem at once:
wrong value types and ranges, immune-compartment cross-references (a core identity or state
parent without a marker set, a hand-off to a compartment that does not exist), lineage names
that no lineage set defines (they would silently select zero cells), unreadable marker files
and a missing input folder. Exit status 1 with a list of problems means nothing was run.

### 5. Run

```bash
annot-st run -c my_run.yaml                                  # everything, in order
annot-st run -c my_run.yaml --stages immune merge spatial    # a subset, in order
annot-st run -c my_run.yaml --no-figures                     # compute only; draw later with `annot-st figures`
annot-st csc -c my_run.yaml --force                          # one stage, even if it is up to date
annot-st status -c my_run.yaml                               # what is up to date, stale or not run
```

A stage that is up to date is skipped, so a job can be requeued safely and a changed
parameter re-runs only the stages it affects (see [trust](#how-a-run-folder-stays-trustworthy)).
Running a single stage on top of an out-of-date upstream stage is refused;
`--allow-stale` overrides that for a quick look.

**Resources.** The full OV section (1.16 M cells × 5 001 genes, QC filtering on) ran on one
node with 32 CPUs in 1 h 34 min and peaked at 79 GB of memory:

| Stage | Wall time |
|---|---|
| qc | 6 min |
| lineage | 37 min |
| caf | 5 min |
| csc | 17 min |
| immune | 8 min |
| merge | 9 min |
| spatial | 11 min |
| `annot-st figures` afterwards (all label maps, 1 M points each) | 47 min |

A SLURM script for a run of that size:

```bash
#!/bin/bash
#SBATCH --job-name=annot_st
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --time=06:00:00
#SBATCH --output=logs/annot_st_%j.out

source activate spatioloji                    # an environment with annotation_st installed
annot-st validate -c my_run.yaml || exit 1
annot-st run -c my_run.yaml                   # --n-jobs defaults to SLURM_CPUS_PER_TASK
```

Smaller sections need proportionally less. For a first look use `--subsample 50000
--outdir smoke --no-figures`, which finishes in minutes.

### 6. Read the outputs

```
<outdir>/
  results/            tables (below)
  figures/            every figure as png and pdf
  qc/                 qc_filter_summary.tsv, qc_metrics_summary_{pre,post}filter.tsv, QC histograms and map
  logs/<stage>.log    one log per stage
  analysis.log        all stages, appended
  params.json         every parameter and the gene sets used, one key per stage
  config_used.yaml    the resolved config of the latest invocation (absolute paths, overrides applied, marker files inlined)
  stage_state.json    per-stage run records (inputs fingerprint, run id, upstream run ids, outputs)
  checkpoints/        <prefix>_01_qc.h5ad, _02_lineage, _03_caf, _04_csc, _05_immune_<compartment>, _06_final
```

**The main result** is `results/final_annotation_by_cell.tsv.gz`, one row per cell:

| Column | Content |
|---|---|
| `leiden` | lineage-stage cluster |
| `lineage` | rough type from the cluster's lineage scores |
| `cell_type` | rough type after merge (cycling tumour, immune hand-offs and spillover classes applied) |
| `cell_type_fine` | the finest label: CAF subtype for fibroblasts, CSC status for tumour cells, immune subtype for immune cells, otherwise `cell_type` |
| `caf_subtype_cell` | iCAF / myCAF / apCAF / `FB_unassigned` (fibroblasts only) |
| `csc_status`, `csc_OCSC_z` | `CSC_like` / `Tumor_nonCSC` and the stemness score (tumour cells only) |
| `immune_compartment`, `immune_subtype` | compartment that labelled the cell and its fine immune label |

`checkpoints/<prefix>_06_final.h5ad` holds the same cells with every column, the scores
(`score_<set>_z`, `caf_<set>_z`, `csc_OCSC_z`), the per-cell `proliferating` and
`immune_spillover` flags, the embeddings and the colours (`uns['<column>_colors']`).

**Tables per stage** (`results/`, names as in the original scripts):

| Stage | Tables |
|---|---|
| lineage | `leiden_lineage_scores.tsv` (cluster × set mean z), `lineage_assignment.tsv`, `leiden_top_markers.tsv`, `cell_type_counts.tsv`, `marker_sets_used.json`, `marker_sources.json` |
| caf | `fb_subcluster_caf_scores.tsv`, `caf_subtype_counts.tsv`, `caf_subtype_by_cell.tsv.gz`, `caf_marker_sets_used.json` |
| csc | `csc_counts.tsv`, `csc_marker_positivity.tsv`, `csc_vs_nonCSC_markers.tsv`, `csc_by_cell.tsv.gz`, `csc_marker_sets_used.json` |
| immune | per compartment `immune_<comp>_subcluster_scores.tsv`, `immune_<comp>_level2_scores.tsv`, `immune_<comp>_subtype_counts.tsv`, `immune_<comp>_top_markers.tsv`, `immune_handoff_<from>_to_<to>.tsv`; `immune_subtype_by_cell.tsv.gz`, `immune_marker_sets_used.json` |
| merge | `final_cell_type_counts.tsv`, `final_cell_type_coarse_counts.tsv`, `final_annotation_by_cell.tsv.gz`, `cell_type_colors.csv` |
| spatial | `nhood_enrichment_zscore.tsv`, `nhood_enrichment_count.tsv`, `csc_neighbourhood_composition.tsv` |

**Figures** (`figures/`, png + pdf):

| Stage | Figures |
|---|---|
| qc | `qc/qc_hist_prefilter`, `qc_hist_postfilter`, `qc_spatial_prefilter` |
| lineage | `umap_leiden`, `umap_cell_type`, `spatial_cell_type`, `dotplot_lineage_markers`, `heatmap_cluster_lineage_scores` |
| caf | `umap_fb_caf_subtype`, `spatial_fb_caf_subtype`, `umap_fb_caf_scores`, `dotplot_fb_caf_markers`, `violin_fb_caf_scores`, `heatmap_fb_subcluster_caf_scores` |
| csc | `umap_tumor_csc`, `spatial_tumor_csc`, `dotplot_tumor_csc_markers`, `tumor_csc_score_hist` |
| immune | per compartment `umap_<comp>_subtype`, `spatial_<comp>_subtype`, `dotplot_<comp>_markers`, `heatmap_<comp>_subcluster_scores`, `heatmap_<comp>_level2_<parent>` |
| merge | `spatial_cell_type_fine` (every cell, hierarchical legend with counts) |
| spatial | `nhood_enrichment_fine`, `csc_neighbour_caf_fraction`, `csc_neighbour_other_fraction` |

**What the labels mean.** Labels that are not a named cell type carry information too:

| Label | Meaning |
|---|---|
| `Unassigned` | lineage cluster in which no lineage set reached `lineage.min_z` |
| `Epithelial_tumor_cycling` | tumour cell whose proliferation score exceeds `lineage.cycling_z` |
| `FB_unassigned` | fibroblast without a clear CAF subtype (no set above `caf.min_z`, or runner-up within `caf.margin`) |
| `CSC_like` / `Tumor_nonCSC` | tumour cell above / below the stemness rule |
| `<Compartment>_spillover_<Lineage>` (e.g. `TNK_spillover_Tumor`) | immune sub-cluster whose contamination set scored higher than every identity: most likely segmentation spillover from the neighbouring lineage |
| `<Identity>_unspecified`, `Macrophage_unpolarized` | identity assigned, no level-2 state reached `immune.min_z2` |
| `T_NK_unresolved`, `B_lineage_unresolved`, `Myeloid_unresolved` | sub-cluster in which no identity qualified |

### 7. Figures and colours

Every UMAP and tissue map that shows cell-type labels uses one label → colour map: the
rough-type maps, the CAF, CSC and immune maps, the CSC neighbour bar charts and the final
map. A label therefore has the same colour wherever it appears, and the colours are also
stored in each checkpoint (`uns['<column>_colors']`), so scanpy, squidpy and cellxgene pick
them up. Panels that show cluster numbers keep scanpy's default colours.

`figures/spatial_cell_type_fine` shows every cell in the colour of its fine type. Its legend
is hierarchical: bold group headers (tumour epithelium, fibroblast / CAF, other stroma, T
cell, NK cell, B lineage, three myeloid groups, unresolved / QC), each followed by its labels
with cell counts. Colours are computed, not hand-picked, in three tiers:

| Tier | Labels | Colours |
|---|---|---|
| grey | unassigned, spillover and unresolved classes (the `Unresolved / QC` group) | greys with next to no hue; the classes seen most often get the best-separated greys |
| bright | cell types with a specific fine type (`Treg`, `myCAF`, `TREM2_TAM`, `Plasma`, ...) | light, saturated colours |
| deep | lineage or identity labels without a specific fine type (`CD8_T_unspecified`, `Macrophage_unpolarized`, `FB_unassigned`, `T_NK`, ...; list in `plot.generic_labels`) | dark, still clearly coloured tones |

Inside a legend group the labels are placed by farthest-point selection in OKLab colour
space, so subtypes that sit together in the tissue are clearly different (at least 15
OKLab units × 100 apart in every cell-type group of the default hierarchy; the tests
enforce this). Rough types are also kept at least that far from each other, and no colour is
used twice. A label's colour depends only on the config, never on which labels happen to be
present. Setting `plot.generic_labels: []` makes every cell type bright at the cost of that
separation in the largest groups.

**Every colour is recorded** in `results/cell_type_colors.csv`, one row per label:
`label`, `hex`, `level` (`fine`, `coarse` or `fine+coarse`), `legend_group`, `coarse_type`
(for a fine type, the rough type most of its cells carry), `n_cells_fine`, `n_cells_coarse`.

**Reproduce or change colours** through the `plot` section; nothing has to be recomputed:

```yaml
plot:
  palette_file: cell_type_colors.csv     # a saved table (label, hex): those labels keep their colours
  colors: {Treg: "#d62728"}              # single overrides; other labels are placed around them
  generic_labels: [T_NK, CD8_T_unspecified, ...]   # labels drawn in deep tones instead of bright ones
  legend_groups:                         # group title -> labels, both in display order
    Tumour epithelium: [Tumor_nonCSC, CSC_like]
    ...
  qc_group: Unresolved / QC              # the group drawn in greys
```

```bash
annot-st figures -c my_run.yaml                                   # redraw all label maps from the checkpoints
annot-st figures -c my_run.yaml --which cell_type_fine immune     # some of them
```

`figures` redraws `umap_cell_type`, `spatial_cell_type`, the CAF, CSC and immune subtype
maps and `spatial_cell_type_fine` from the checkpoints without recomputing anything; dot
plots and heatmaps need the expression data and are drawn by the stage itself.

**Using the scheme on other datasets.** The scheme is a rule set, not a fixed list:

- **Another section run through this package** gets the same colours automatically. To
  freeze the colours for a whole project regardless of later config edits, export them once
  and point every dataset at the file:

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
  `uns['<key>_colors']` in the file, `--map` draws the hierarchical-legend tissue map.

### 8. Compare with a reference annotation

```bash
annot-st compare -c my_run.yaml --ref accepted_labels.tsv.gz                 # report in <outdir>/results/compare/
annot-st compare --new run/results/final_annotation_by_cell.tsv.gz --ref reference.h5ad -o cmp/
annot-st compare -c my_run.yaml --ref reference.tsv.gz --min-agreement 0.99  # exit 1 below the bar (CI gate)
```

Cells are matched by id. For every label column both tables share the report gives the
agreement, the adjusted Rand index (robust to renumbered clusters), per-label recall /
precision and the confusion table; numeric columns get a correlation and the largest
difference. `--columns` restricts the comparison. Tables can be `.tsv[.gz]`, `.csv`,
`.parquet` or `.h5ad` (`obs`).

### 9. Use it from Python

**Whole annotation on an AnnData in memory** (no files, no checkpoints; the labels are the
ones the command line gives for the same data and config):

```python
import annotation_st as ast
from annotation_st.api import annotate_adata
from annotation_st.io import read_xenium

cfg = ast.load_config("my_run.yaml")         # or ast.preset_config() for the whole default library
adata = read_xenium("/path/to/xenium_outs")  # raw counts in X, centroids in obsm['spatial']

res = annotate_adata(adata, cfg, qc=True, spatial=True)
res.adata.obs[["cell_type", "cell_type_fine", "csc_status"]]
res.palette                                  # label -> hex, for every rough and fine type present
res.stages                                   # {'lineage': 'completed', 'caf': 'skipped: no CAF marker sets opted in', ...}
res.tables["final_cell_type_counts"]         # the result tables as DataFrames
```

Arguments: `qc=True` computes the QC metrics first (nothing is removed unless
`cfg.qc.apply_filters`), `spatial=True` adds the neighbourhood statistics, `normalize=False`
says `X` is already log-normalised (raw counts then expected in `layers['counts']`),
`umap=True` computes the lineage UMAP, `copy=False` works in place.

**With spatioloji_s.** `annotation_st` works directly on the `spatioloji` object of the
[spatioloji_s](https://github.com/gynecoloji/spatioloji_s) package and puts its results where
spatioloji's spatial, communication and plotting functions look for them.

```python
import spatioloji_s as sj
import annotation_st as ast
from annotation_st import spatioloji_bridge as bridge

sp = sj.spatioloji.from_xenium("/path/to/xenium_outs")
cfg = ast.load_config("my_run.yaml")

out = bridge.annotate(sp, cfg)                    # filter with sj.xenium_qc before, if you want filtering
sp.cell_meta[["cell_type", "cell_type_fine", "csc_status"]].head()

# spatioloji analyses use the new columns like any other
sj.visualization.plot_umap(sp, color_by="cell_type_fine", colors=out.palettes["cell_type_fine"])
sj.visualization.xenium_plot_spatial(sp, "cell_type_fine", color_dict=out.palettes["cell_type_fine"])
config = sj.ccc.CCCConfig(group_col="cell_type_fine", layer="log_normalized")
```

What `bridge.annotate` does:

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
  and by position (`palette=` in the embedding plots).
- **Input details handled for you:** spatioloji's Xenium loader keeps control probes and
  codewords in the expression matrix. `bridge.to_anndata(sp)` separates them, sums them into the
  QC columns, and stores coordinates as (x, y).
- `log_layer="log_normalized"` reuses a layer you normalised with spatioloji instead of
  normalising the raw counts again.

A saved object can also be the input of the command-line pipeline, and a finished run can be
written back into one:

```yaml
project: {spatioloji: sample.pkl, outdir: annot_run}
```

```bash
annot-st run -c my_run.yaml
annot-st to-spatioloji -c my_run.yaml -o sample_annotated.pkl    # labels back into the object
```

`bridge.add_annotations(sp, "annot_run/results/final_annotation_by_cell.tsv.gz")` does the
last step from Python.

**Single stages.** The stage functions are pure: they take an `AnnData` and a config section,
add columns, and return a result object with the tables. They never write files or draw figures.

```python
import scanpy as sc
from annotation_st import load_config
from annotation_st.stages.lineage import annotate_lineage
from annotation_st.stages.caf import annotate_caf
from annotation_st.stages.csc import call_csc
from annotation_st.stages.immune import run_immune

cfg = load_config("my_run.yaml")
adata = sc.read_h5ad("sample_01_qc.h5ad")                   # raw counts
res = annotate_lineage(adata, cfg.lineage, cfg.markers)     # adds obs['leiden','lineage','cell_type']
res.means                                                   # cluster × lineage mean-z table

fb = adata[adata.obs["lineage"] == "Fibroblast"].copy()
caf = annotate_caf(fb, cfg.caf, {"iCAF": [...], "myCAF": [...], "apCAF": [...]})

tum = adata[adata.obs["lineage"] == "Epithelial_tumor"].copy()
csc = call_csc(tum, cfg.csc, csc_genes=[...], core_genes=[...])

immune = run_immune(adata, cfg.immune, cfg.markers)         # {'TNK': ImmuneResult, 'Bcell': ..., 'Myeloid': ...}
```

Colour helpers: `annotation_st.pipeline.label_palette(cfg, labels)` gives the scheme's colours
for any labels, `annotation_st.plotting.apply_palette(adata, key, palette)` stores them in
`uns`, and `annotation_st.plotting.build_fine_map(xy, labels, palette)` draws the
hierarchical-legend tissue map.

## Command reference

| Command | Purpose | Main options |
|---|---|---|
| `annot-st init-config` | write a config for editing | `-o FILE`, `--xenium-dir DIR`, `--markers default\|explicit\|none`, `--force` |
| `annot-st markers` | show the default marker library | `--set NAME`, `--references refs.tsv`, `--dropped dropped.tsv` |
| `annot-st validate -c CFG` | check the config without running | `--stages ...` |
| `annot-st run -c CFG` | run all or selected stages in order | `--stages ...`, `--force`, `--allow-stale`, `--no-figures`, `--subsample N`, `--n-jobs N` |
| `annot-st qc\|lineage\|caf\|csc\|immune\|merge\|spatial -c CFG` | run one stage | same as `run` |
| `annot-st status -c CFG` | up to date / stale (with reason) / not run, per stage | |
| `annot-st figures -c CFG` | redraw the label maps from the checkpoints | `--which cell_type caf csc immune cell_type_fine` |
| `annot-st palette` | export the colour scheme or apply it to another h5ad | `-c CFG`, `--palette-file`, `-o`, `--h5ad`, `--keys`, `--write-h5ad`, `--map`, `--map-key` |
| `annot-st to-spatioloji -c CFG -o OUT.pkl` | write a run's labels into a spatioloji_s object | `--input OBJ`, `--label-prefix` |
| `annot-st compare --ref REF` | compare per-cell labels with a reference | `-c CFG` or `--new TABLE`, `--columns`, `-o DIR`, `--min-agreement` |

Every command that takes `-c` also accepts `--outdir`, `--checkpoint-dir`, `--xenium-dir`,
`--prefix`, `--seed` and `--subsample` as overrides of the config. Exit status is 0 on success
and 1 on a configuration or input problem.

## Troubleshooting

| Symptom | Cause and fix |
|---|---|
| `validate` says the lineage stage has no marker sets | the `markers` section is missing or has no `lineage` entry; add `lineage: default` or your own sets |
| `ConfigError: unknown key ...` | a typo or a key from another section; keys are checked strictly |
| a lineage name is "not defined by any lineage set" | `caf.include_lineages`, `csc.tumor_lineages` or a compartment's `lineages` names a type that `markers.lineage` does not score; it would select zero cells |
| stage recorded as `skipped` | too few cells for it (`min_cells`) or no marker sets opted in; the merged labels keep the rough type for those cells. Check `analysis.log` |
| a single-stage command refuses to run on a stale upstream | an upstream stage changed since; run `annot-st run` to refresh in order, or pass `--allow-stale` for a quick look |
| labels differ slightly between two machines | Leiden partitions can shift with the thread count; `n_jobs` is recorded per stage in `params.json`. Use the same `--n-jobs` to reproduce |
| genes of a marker set are missing from the panel | expected; the set is intersected with the panel and the dropped genes are listed in `results/*_marker_sets_used.json` |
| the full run runs out of memory | a 1 M-cell section peaked at 79 GB; request more memory, or put `checkpoint_dir` on a fast disk and run stage by stage |
| figures take very long | the full-section maps draw a million points each; run with `--no-figures` and draw later with `annot-st figures`, or lower `umap_cells` |
| `annot-st: command not found` | the install went into another environment; run `python -m annotation_st.cli ...` from the right interpreter |

## How a run folder stays trustworthy

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

A stage that skips itself (too few cells) is recorded as such, and `merge` only uses
refinements that the state file lists as current, so a checkpoint left over from an earlier
configuration cannot leak into the final labels. Settings that only affect figures or speed
(`project.figures`, `project.n_jobs`, `*.umap_cells`, `lineage.skip_umap`) do not invalidate
anything.

**The run folder keeps its own configuration.** `config_used.yaml` is the resolved config
of the latest invocation: absolute paths, CLI overrides applied, marker files inlined. It
can be passed back to `annot-st run -c` to reproduce the run without the original files.

**Validated on real data.** The package was run on the full Xenium 5K OV section against the
accepted annotation: QC, lineage, CAF, CSC and OCSC scores identical on all 1,083,808 cells,
immune labels differing only where the reference scripts moved to a revised hand-off scheme;
the in-memory and spatioloji_s routes reproduce the command-line labels exactly. The
validation jobs, reports and figures are kept outside the repository.

## The default marker library and its references

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

## Package layout

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
configs/xenium_ov_5k.yaml              full config for the OV section
configs/default_library_example.yaml   a config that opts in to the library
configs/cell_type_palette.csv          the default colour scheme, exported
tools/                                 build_marker_references.py, find_literature_evidence.py, reviewed sentences
tests/                                 synthetic Xenium-like data with planted populations
```

## Tests

```bash
MPLBACKEND=Agg python -m pytest tests -q            # ~3 min
MPLBACKEND=Agg python -m pytest tests -q -m slow    # + end-to-end run with figures
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

## Citing

See `CITATION.cff`. Contributions follow `CONTRIBUTING.md` (Conventional Commits; releases
are cut by release-please).
