# annotation_st — design (2026-09-13)

Python package that packages the cell-annotation workflow documented in
`Grant/st_datasets/01-documentation/annotation_methods_01.md` (Xenium 5K human ovarian
adenocarcinoma, analysis `01_annot_Xenium_OV`) so it can be rerun, re-parameterised and
applied to other Xenium sections without editing scripts.

## Goals

1. Reproduce the documented decision path (Steps 0–6) from a single config file whose
   defaults equal the parameters of the final run (job 1180895).
2. Make every threshold, marker set and compartment definition data, not code
   (YAML), so another panel/tissue only needs a new config.
3. Expose each step as a library function on an `AnnData` (no file I/O), plus a CLI
   that chains them with checkpoints.
4. Tests on small synthetic data for every decision rule.

## Non-goals (v1)

- Analysis 02 (ligand–receptor communication) and the pptx report. These are
  downstream of annotation and are left for a later `ccc` stage.
- Any change to the biology of the rules. Where the source scripts and the methods
  document differ, the scripts win (they produced the run).

## Layout

```
annotation_st/                      (this folder)
  pyproject.toml                    package metadata; console script `annot-st`
  README.md
  configs/xenium_ov_5k.yaml         full config reproducing the 2026-09-08 run
  src/annotation_st/
    __init__.py
    config.py      dataclasses (QCConfig, LineageConfig, CafConfig, CscConfig,
                   ImmuneConfig, SpatialConfig, ProjectConfig, MarkerConfig,
                   PipelineConfig) + load_config(yaml) / default_config()
    markers.py     default marker sets (lineage, proliferation, CSC, spillover,
                   immune compartments), read_marker_file, filter_to_panel
    io.py          read_xenium(dir, subsample), RunPaths, setup_logging,
                   update_params, append_log
    scoring.py     score_sets, cluster_assign, per_cell_assign, embed
    plotting.py    save_fig, spatial_plot, umap_plot, dotplot, score_heatmap, grouped_bar
    stages/
      qc.py        compute_qc_metrics, qc_filter
      lineage.py   annotate_lineage
      caf.py       annotate_caf
      csc.py       call_csc
      immune.py    annotate_compartment, run_immune (ordering + hand-off)
      merge.py     merge_labels
      spatial.py   build_spatial_graph, nhood_enrichment, neighbourhood_composition
    pipeline.py    stage registry; run_stage(name, cfg, paths, force) reads the
                   checkpoint, calls the stage, writes h5ad/tables/figures/params
    cli.py         annot-st init-config | run | qc | lineage | caf | csc | immune |
                   merge | spatial
  tests/           synthetic AnnData fixtures + one test module per source module
```

## Data flow and checkpoints

```
read_xenium ─► qc ─► lineage ─► {caf, csc, immune} ─► merge ─► spatial
              01      02          03   04   05_*        06        (in place on 06)
```

Checkpoints (`<checkpoint_dir>/<prefix>_NN_<stage>.h5ad`, gzip) mirror the current
`xenium_ov_NN_*.h5ad` files. `lineage` writes the full log-normalised object; `caf`,
`csc` and `immune` write subsets; `merge` writes the final full object and `spatial`
adds `obsp`/`uns` to it. `pipeline.run` skips a stage whose checkpoint exists unless
`--force`.

Output directories are configurable independently so the package fits the project
convention (h5ad → `03-data/.../Processed`, everything else → a run folder):

| Key | Default | Content |
|---|---|---|
| `project.outdir` | `./annot_run` | `results/`, `figures/`, `logs/`, `qc/`, `params.json`, `analysis.log` |
| `project.checkpoint_dir` | `<outdir>/checkpoints` | h5ad checkpoints |
| `project.prefix` | `sample` | checkpoint file prefix |

## Stage contracts

Each core function takes an `AnnData` (already subset where relevant) and a config
dataclass and returns a small result object holding the modified `AnnData` and the
tables the run folder needs. No core function writes files or draws figures.

| Stage | Core function | Adds to obs | Tables returned |
|---|---|---|---|
| qc | `compute_qc_metrics(adata)`, `qc_filter(adata, cfg) -> (adata, flag_summary)` | gene_counts, control_frac, nucleus_ratio, counts_per_area | qc_filter_summary, pre/post metric summaries |
| lineage | `annotate_lineage(adata, cfg, markers) -> LineageResult` | leiden, score_*, lineage, proliferating, cell_type | cluster×set mean z with assigned/n_cells, cell_type_counts, lineage_assignment, top markers |
| caf | `annotate_caf(fb, cfg, caf_sets) -> CafResult` | fb_leiden, caf_*, caf_subtype_cell, caf_subtype_cluster, caf_best_z, caf_margin | subcluster scores, subtype counts |
| csc | `call_csc(tum, cfg, markers) -> CscResult` | csc_OCSC, csc_OCSC_z, csc_n_core_detected, csc_status, tumor_leiden | counts, marker positivity, DE table |
| immune | `annotate_compartment(sub, name, comp_cfg, cfg, prior) -> ImmuneResult`; `run_immune(full, cfg) -> dict[name, ImmuneResult]` | `<name>_leiden`, imm1_*, imm2_*, immune_level1, immune_spillover, immune_subtype, immune_state_cell | level-1 and level-2 score tables, counts, hand-off table, top markers |
| merge | `merge_labels(full, fb_obs, tum_obs, immune) -> full` | cell_type_fine, updated cell_type, immune_subtype, immune_compartment, caf_*, csc_* | fine and coarse counts |
| spatial | `build_spatial_graph(adata, cfg)`, `nhood_enrichment(adata, key, cfg)`, `neighbourhood_composition(adata, key, groups)` | obsp spatial_*, uns nhood | z-score and count matrices, CSC composition |

Decision rules are copied verbatim from the scripts:

- lineage: cluster mean z of best set > `min_z` (0.5) else `Unassigned`; cycling =
  proliferation z > `cycling_z` (1.0) within `Epithelial_tumor`.
- caf: per cell best z > `min_z` (0) and best − second ≥ `margin` (0.25) else
  `FB_unassigned`; sub-cluster label by mean z.
- csc: z > `z_thresh` (1.5) and ≥ `min_core` (2) core markers with raw count > 0.
- immune level 1: core identities qualify at mean z ≥ `low_identity` (−0.5), minority /
  cross at > `min_z` (0.5); otherwise spillover set > `min_z` → `<comp>_spillover_<X>`;
  otherwise `<comp>_unresolved`; `immune_spillover` flag at > `spill_flag` (1.0);
  hand-off labels move to the receiving compartment and keep their identity there.
- immune level 2: re-cluster within parent at `resolution2` (1.0); state mean z >
  `min_z2` (0.25) else the parent's default label; cross labels map to the other
  compartment's coarse label.
- merge: fine label = CAF subtype / CSC status / immune subtype / coarse label;
  coarse label updated only through each compartment's `coarse` map.
- spatial: Delaunay graph, edges > 99th percentile removed, permutation enrichment
  with `n_jobs=1`, non-finite z set to 0, CSC vs non-CSC neighbour fractions with
  log2 ratio (pseudo-count 1e-4).

## Configuration

One YAML file, sections `project`, `qc`, `lineage`, `caf`, `csc`, `immune`,
`spatial`, `markers`. `markers` holds every gene list (lineage sets, proliferation,
csc + csc_core, spillover sets, and per compartment: lineages, core, level1, level2,
unspecified, cross, handoff, coarse). CAF sets can be inline lists or
`caf.marker_files: {iCAF: path.txt, ...}`. `annot-st init-config` writes the shipped
OV config for editing. Unknown keys raise; missing keys take defaults so a config can
be partial.

## CLI

```
annot-st init-config [-o config.yaml]
annot-st run   -c config.yaml [--stages qc lineage ...] [--force] [--no-figures]
annot-st qc|lineage|caf|csc|immune|merge|spatial -c config.yaml [--force] [--no-figures]
```

Common overrides: `--outdir`, `--checkpoint-dir`, `--xenium-dir`, `--n-jobs`,
`--seed`, `--subsample N` (qc only, for smoke tests).

## Error handling

- Missing input checkpoint → clear error naming the stage to run first.
- Compartment/lineage with fewer cells than `min_cells` (50 for caf/csc, 100 for
  immune) → stage logs a warning and returns `None`; merge tolerates missing pieces.
- Marker sets with zero panel genes are skipped and recorded under
  `dropped_not_on_panel` in `results/*_marker_sets_used.json`.
- Every stage records its parameters and versions in `params.json[<stage>]`.

## Testing

`tests/conftest.py` builds a synthetic `AnnData` (~3,000 cells, ~250 genes) with
planted populations: tumour (with a CSC-like subset), fibroblasts split iCAF/myCAF,
CD8 T, CD4 T (incl. Treg), NK, plasma, B, macrophages, endothelium, spatial
coordinates clustered by population. Tests assert each rule on this data (majority
recovery, threshold edges, hand-off, unresolved/spillover fall-through, merge
precedence, graph edge cut, composition log2 ratio) and a CLI smoke run of all stages
with `--no-figures`. Run with the `spatioloji` env
(`/users/jiwang1/.conda/envs/spatioloji/bin/python -m pytest`).
