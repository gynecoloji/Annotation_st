# annotation_st Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A pip-installable package `annotation_st` (CLI `annot-st`) that reruns the Xenium OV annotation pipeline (Steps 0–6) from one YAML config, with every stage usable as a library function on an AnnData.

**Architecture:** Pure stage functions in `annotation_st.stages.*` take an AnnData + config dataclass and return a result object (AnnData + tables). `pipeline.py` wraps them with checkpoint I/O, figures and `params.json`. `config.py` turns YAML into dataclasses whose defaults equal the documented OV run; `markers.py` holds the default gene sets. `cli.py` is a thin argparse layer.

**Tech Stack:** Python 3.12, scanpy 1.11.5, squidpy 1.8.1, anndata 0.12.7, PyYAML, pytest. Development interpreter: `/users/jiwang1/.conda/envs/spatioloji/bin/python` (alias `$PY` below).

**Spec:** `docs/superpowers/specs/2026-09-13-annotation-st-design.md`

## Global Constraints

- All package files live under `/easley/scratch/projects/amitra/amitra2016502/annotation_st` (src layout `src/annotation_st`).
- Decision rules and default values are copied from `Grant/st_datasets/02-scripts/Script_01/01_xenium_*.py`; when the methods doc and the scripts disagree, the scripts win.
- Core stage functions do no file I/O and draw no figures.
- Not a git repository; no commit steps. Each task ends with `$PY -m pytest tests -q` green.
- Install once in editable mode: `$PY -m pip install -e . --no-deps --no-build-isolation`.
- Set `MPLBACKEND=Agg` and `NUMBA_NUM_THREADS=4` when running tests.

---

### Task 1: Package skeleton, config dataclasses, YAML loading

**Files:**
- Create: `pyproject.toml`, `README.md`, `src/annotation_st/__init__.py`, `src/annotation_st/config.py`
- Test: `tests/test_config.py`

**Interfaces:**
- Produces:
  - `@dataclass ProjectConfig(name="sample", xenium_dir=None, outdir="annot_run", checkpoint_dir=None, prefix="sample", seed=0, n_jobs=8, figures=True)`
  - `@dataclass QCConfig(subsample=0, min_transcripts=30, min_genes=10, max_control_frac=0.02, area_pct=(1.0, 99.0), require_nucleus=True, min_cells_per_gene=10)`
  - `@dataclass LineageConfig(n_pcs=50, n_neighbors=15, resolution=1.0, min_z=0.5, cycling_z=1.0, umap_cells=300000, skip_umap=False, cycling_lineage="Epithelial_tumor", cycling_label="Epithelial_tumor_cycling")`
  - `@dataclass CafConfig(include_lineages=("Fibroblast",), n_pcs=30, n_neighbors=15, resolution=0.5, min_z=0.0, margin=0.25, min_cells=50, unassigned="FB_unassigned", marker_files={})`
  - `@dataclass CscConfig(tumor_lineages=("Epithelial_tumor",), z_thresh=1.5, min_core=2, n_pcs=30, n_neighbors=15, resolution=0.4, umap_cells=250000, min_cells=50, positive_label="CSC_like", negative_label="Tumor_nonCSC")`
  - `@dataclass ImmuneConfig(compartments=None (=all), n_pcs=30, n_neighbors=15, resolution=1.5, min_z=0.5, low_identity=-0.5, resolution2=1.0, min_z2=0.25, spill_flag=1.0, umap_cells=200000, min_cells=100, min_cells_level2=50)`
  - `@dataclass SpatialConfig(n_perms=100, radius=0.0, edge_cut_pct=99.0, fine_key="cell_type_fine", composition_groups=("CSC_like","Tumor_nonCSC"))`
  - `@dataclass CompartmentConfig(lineages, level1: dict[str,list[str]], level2: dict[str,dict[str,list[str]]], core: list[str], unresolved: str, cross: dict[str,str], unspecified: dict[str,str], handoff: dict[str,str]={}, coarse: dict[str,str]={})`
  - `@dataclass MarkerConfig(lineage: dict, proliferation: list, csc: list, csc_core: list, spillover: dict, compartments: dict[str, CompartmentConfig], caf: dict[str,list[str]]={})`
  - `@dataclass PipelineConfig(project, qc, lineage, caf, csc, immune, spatial, markers)`
  - `default_config() -> PipelineConfig` (uses `markers.DEFAULTS`), `load_config(path) -> PipelineConfig`, `config_to_dict(cfg) -> dict`, `dump_config(cfg, path)`.
  - Unknown keys in any section raise `ValueError("unknown key 'x' in section 'qc'")`.

- [ ] **Step 1: Write failing tests**

```python
# tests/test_config.py
import textwrap, pytest
from annotation_st.config import default_config, load_config, dump_config, PipelineConfig

def test_defaults_match_documented_run():
    cfg = default_config()
    assert cfg.qc.min_transcripts == 30 and cfg.qc.max_control_frac == 0.02
    assert cfg.lineage.resolution == 1.0 and cfg.lineage.min_z == 0.5
    assert cfg.caf.margin == 0.25 and cfg.csc.z_thresh == 1.5 and cfg.csc.min_core == 2
    assert cfg.immune.low_identity == -0.5 and cfg.immune.min_z2 == 0.25
    assert cfg.spatial.n_perms == 100
    assert set(cfg.markers.compartments) == {"TNK", "Bcell", "Myeloid"}
    assert cfg.markers.compartments["TNK"].handoff == {"B_cell": "Bcell", "Plasma": "Bcell"}

def test_partial_yaml_overrides_only_given_keys(tmp_path):
    p = tmp_path / "c.yaml"
    p.write_text(textwrap.dedent("""
        project: {name: demo, outdir: out}
        qc: {min_transcripts: 300}
        markers:
          lineage: {Foo: [A, B], Bar: [C]}
    """))
    cfg = load_config(p)
    assert cfg.qc.min_transcripts == 300 and cfg.qc.min_genes == 10
    assert cfg.markers.lineage == {"Foo": ["A", "B"], "Bar": ["C"]}
    assert cfg.markers.proliferation  # untouched default

def test_unknown_key_raises(tmp_path):
    p = tmp_path / "c.yaml"; p.write_text("qc: {min_trnscripts: 3}\n")
    with pytest.raises(ValueError, match="unknown key"):
        load_config(p)

def test_roundtrip(tmp_path):
    cfg = default_config(); p = tmp_path / "c.yaml"
    dump_config(cfg, p)
    assert load_config(p) == cfg
```

- [ ] **Step 2: Run** `$PY -m pytest tests/test_config.py -q` → fails (ImportError).
- [ ] **Step 3: Implement** `pyproject.toml` (setuptools, `[project.scripts] annot-st = "annotation_st.cli:main"`, deps scanpy, squidpy, anndata, pyyaml, pandas, numpy, scipy, matplotlib, h5py, pyarrow), `config.py` with the dataclasses above, a generic `_from_dict(cls, d, section)` that checks unknown keys, converts lists→tuples for tuple fields, and builds `CompartmentConfig` objects. `markers.py` must exist for `default_config()` — create it in this task with `DEFAULTS` (see Task 2 for the full content; here copy the dicts from `01_xenium_common.py`).
- [ ] **Step 4: Install editable, run tests** → PASS.

### Task 2: Marker utilities

**Files:**
- Modify: `src/annotation_st/markers.py`
- Test: `tests/test_markers.py`

**Interfaces:**
- Produces: `LINEAGE_MARKERS, PROLIF_MARKERS, CSC_MARKERS, CSC_CORE, SPILLOVER_SETS, TNK_LEVEL1, TNK_LEVEL2, MYELOID_LEVEL1, MYELOID_LEVEL2, B_LEVEL1, B_LEVEL2, IMMUNE_COMPARTMENTS` (verbatim from `01_xenium_common.py`), `DEFAULTS: dict` (the `markers` YAML section), `read_marker_file(path) -> list[str]`, `load_marker_files(mapping: dict[str, path]) -> dict[str, list[str]]`, `filter_to_panel(sets, var_names) -> (kept, dropped)`.

- [ ] **Step 1: Tests**

```python
def test_read_marker_file_dedups_and_skips_comments(tmp_path):
    f = tmp_path / "x.txt"; f.write_text("A\n# c\nB\n\nA\n")
    assert read_marker_file(f) == ["A", "B"]

def test_filter_to_panel():
    kept, dropped = filter_to_panel({"S": ["A", "Z"]}, ["A", "B"])
    assert kept == {"S": ["A"]} and dropped == {"S": ["Z"]}

def test_defaults_have_16_lineages_and_prolif():
    assert len(LINEAGE_MARKERS) == 16 and PROLIF_MARKERS[0] == "MKI67"
```

- [ ] **Step 2–4:** run → fail; implement; run → pass.

### Task 3: Scoring primitives

**Files:**
- Create: `src/annotation_st/scoring.py`, `tests/conftest.py`
- Test: `tests/test_scoring.py`

**Interfaces:**
- Produces:
  - `score_sets(adata, sets, prefix, seed=0, ctrl_size=50) -> list[str]` (obs `prefix+name`, `prefix+name+"_z"`; z across `adata`).
  - `cluster_assign(adata, cluster_key, sets, prefix, min_z=0.0, unassigned="Unassigned") -> (means: DataFrame, mapping: dict)`.
  - `per_cell_assign(Z: ndarray, names, min_z, margin, unassigned) -> (call: ndarray[str], best, margin_arr)`.
  - `embed(adata, n_pcs, n_neighbors, resolution, seed, key, umap=True, umap_cells=0) -> None` (PCA→kNN→Leiden igraph 2 iterations; UMAP on a random subset with NaN elsewhere when `umap_cells` < n_obs; stores `adata.uns[key+"_umap_subset_idx"]`).
  - `conftest.synthetic_adata(seed=0) -> AnnData` fixture `synth`: ~3,000 cells, ~250 genes = union of all default marker genes that matter + filler; populations with counts drawn Poisson(λ=0.3 background, λ=4 for that population's markers); `obs["truth"]`; `layers["counts"]`, X log-normalised; `obsm["spatial"]` with each population in its own blob; `obs` QC columns (`cell_area`, `nucleus_area`, `nucleus_count`, ctrl counts) for qc tests. Populations: Epithelial_tumor (600, of which 100 `CSC` with CD44/PROM1/POU5F1/CXCR4 λ=6), Fibroblast (500: 250 iCAF-high, 250 myCAF-high), CD8_T (250: 120 with cytotoxic genes), CD4_T (250: 100 Treg), NK (150), Plasma (200), B_cell (120), Macrophage (400: 150 TREM2), Endothelial (200), Ovarian_stroma (150).

- [ ] **Step 1: Tests**

```python
def test_per_cell_assign_margin_rule():
    Z = np.array([[1.0, 0.5], [1.0, 0.9], [-0.1, -0.5]])
    call, best, marg = per_cell_assign(Z, ["a", "b"], min_z=0.0, margin=0.25, unassigned="U")
    assert list(call) == ["a", "U", "U"]

def test_cluster_assign_threshold(synth):
    sets, _ = filter_to_panel(LINEAGE_MARKERS, synth.var_names)
    scored = score_sets(synth, sets, "score_")
    means, mapping = cluster_assign(synth, "truth_coarse", scored, "score_", min_z=0.5)
    assert mapping["Fibroblast"] == "Fibroblast" and mapping["Epithelial_tumor"] == "Epithelial_tumor"

def test_embed_umap_subset(synth):
    embed(synth, 20, 15, 0.5, 0, "leiden", umap_cells=500)
    assert np.isnan(synth.obsm["X_umap"][:, 0]).sum() == synth.n_obs - 500
```

- [ ] **Step 2–4:** fail → implement → pass.

### Task 4: QC stage

**Files:**
- Create: `src/annotation_st/io.py` (`read_xenium(xenium_dir, subsample=0, seed=0) -> AnnData` copied from `read_xenium_h5` + parquet join; `RunPaths(cfg) ` with `.outdir .results .figures .logs .qc .checkpoints .params .analysis_log`, `.checkpoint(stage) -> Path` using `{prefix}_{NN}_{stage}.h5ad` (qc=01, lineage=02, caf=03, csc=04, immune_<name>=05, final=06); `setup_logging(paths, stage)`, `update_params(paths, stage, dict)`, `append_log(paths, msg)`, `now()`), `src/annotation_st/stages/__init__.py`, `src/annotation_st/stages/qc.py`
- Test: `tests/test_qc.py`, `tests/test_io.py`

**Interfaces:**
- Produces: `compute_qc_metrics(adata) -> None` (adds gene_counts, control_counts, control_frac, nucleus_ratio, counts_per_area, n_genes_by_counts); `qc_filter(adata, cfg: QCConfig) -> QCResult(adata, flag_summary: DataFrame, area_bounds: tuple, summary_pre: DataFrame, summary_post: DataFrame)`; copies `X` to `layers["counts"]`.

- [ ] **Step 1: Tests**

```python
def test_qc_filter_removes_planted_failures(synth_raw):
    # synth_raw: counts in X, 20 cells with gene_counts<30, 10 with nucleus_count 0
    compute_qc_metrics(synth_raw)
    res = qc_filter(synth_raw, QCConfig(area_pct=(0, 100)))
    assert res.flag_summary.loc["low_transcripts", "n_cells_flagged"] >= 20
    assert res.flag_summary.loc["no_nucleus", "n_cells_flagged"] == 10
    assert "counts" in res.adata.layers and res.adata.n_obs == res.flag_summary.loc["KEPT", "n_cells_flagged"]

def test_runpaths_checkpoint_names(tmp_path):
    cfg = default_config(); cfg.project.outdir = str(tmp_path); cfg.project.prefix = "ov"
    p = RunPaths(cfg)
    assert p.checkpoint("lineage").name == "ov_02_lineage.h5ad"
    assert p.checkpoint("immune_TNK").name == "ov_05_immune_TNK.h5ad"
```

- [ ] **Step 2–4:** fail → implement → pass.

### Task 5: Lineage stage

**Files:**
- Create: `src/annotation_st/stages/lineage.py`
- Test: `tests/test_lineage.py`

**Interfaces:**
- Produces: `annotate_lineage(adata, cfg: LineageConfig, markers: MarkerConfig, seed=0, normalize=True) -> LineageResult(adata, means, mapping, counts, sets_used, sets_dropped, prolif_used, top_markers: DataFrame)`. Normalises (`normalize_total` + `log1p`, sets `raw`) when `normalize`, calls `embed(..., key="leiden")`, scores, assigns, cycling rule, `rank_genes_groups` top-15 per cluster.

- [ ] **Step 1: Tests**

```python
def test_lineage_recovers_planted_populations(synth_raw):
    res = annotate_lineage(synth_raw, LineageConfig(n_pcs=20, umap_cells=0, skip_umap=True), default_config().markers)
    ct = pd.crosstab(res.adata.obs["truth_coarse"], res.adata.obs["lineage"])
    for lin in ["Epithelial_tumor", "Fibroblast", "T_NK", "Macrophage_Mono", "Plasma", "Endothelial"]:
        assert ct.loc[lin, lin] / ct.loc[lin].sum() > 0.8
    assert (res.adata.obs["cell_type"] == "Epithelial_tumor_cycling").sum() > 0
    assert set(res.means.columns) >= {"assigned", "n_cells", "mean_prolif_z"}
```

- [ ] **Step 2–4:** fail → implement → pass.

### Task 6: CAF and CSC stages

**Files:**
- Create: `src/annotation_st/stages/caf.py`, `src/annotation_st/stages/csc.py`
- Test: `tests/test_caf.py`, `tests/test_csc.py`

**Interfaces:**
- `annotate_caf(fb, cfg: CafConfig, caf_sets: dict, seed=0) -> CafResult(adata, means, counts, sets_used, sets_dropped) | None` (None if `fb.n_obs < cfg.min_cells`). Adds `fb_leiden`, `caf_<s>`, `caf_<s>_z`, `caf_subtype_cell`, `caf_best_z`, `caf_margin`, `caf_subtype_cluster`.
- `call_csc(tum, cfg: CscConfig, markers: MarkerConfig, seed=0) -> CscResult(adata, counts, positivity, de: DataFrame|None, sets_used, sets_dropped) | None`. Adds `csc_OCSC`, `csc_OCSC_z`, `csc_n_core_detected`, `csc_status`, `tumor_leiden`, `X_umap`.

- [ ] **Step 1: Tests**

```python
def test_caf_per_cell_and_cluster_labels(synth_fb):
    res = annotate_caf(synth_fb, CafConfig(n_pcs=15), {"iCAF": [...], "myCAF": [...], "apCAF": [...]})
    ct = pd.crosstab(res.adata.obs["truth"], res.adata.obs["caf_subtype_cell"])
    assert ct.loc["iCAF", "iCAF"] > ct.loc["iCAF"].drop("iCAF").max()
    assert ct.loc["myCAF", "myCAF"] > ct.loc["myCAF"].drop("myCAF").max()
    assert set(res.counts.columns) == {"per_cell", "per_subcluster", "frac_per_cell"}

def test_caf_too_few_cells_returns_none(synth_fb):
    assert annotate_caf(synth_fb[:10].copy(), CafConfig(), {...}) is None

def test_csc_rule(synth_tum):
    res = call_csc(synth_tum, CscConfig(n_pcs=15, umap_cells=0), default_config().markers)
    obs = res.adata.obs
    assert ((obs["csc_status"] == "CSC_like") == ((obs["csc_OCSC_z"] > 1.5) & (obs["csc_n_core_detected"] >= 2))).all()
    ct = pd.crosstab(obs["truth"], obs["csc_status"])
    assert ct.loc["CSC", "CSC_like"] / ct.loc["CSC"].sum() > 0.5
    assert ct.loc["Epithelial_tumor", "CSC_like"] / ct.loc["Epithelial_tumor"].sum() < 0.1
```

- [ ] **Step 2–4:** fail → implement → pass.

### Task 7: Immune stage (level 1, hand-off, level 2)

**Files:**
- Create: `src/annotation_st/stages/immune.py`
- Test: `tests/test_immune.py`

**Interfaces:**
- `assign_level1(means1: DataFrame, identities, spills, comp: CompartmentConfig, cfg: ImmuneConfig, name) -> (level1: dict, flag: dict)` — pure function on the cluster×set mean-z table (unit-testable without scanpy).
- `annotate_compartment(sub, name, comp: CompartmentConfig, cfg: ImmuneConfig, spillover_sets, seed=0, prior: dict|None=None) -> ImmuneResult(adata, means1, level2_tables: dict, counts, handoff: dict[target, dict[cell, label]], sets_used, top_markers) | None`.
- `run_immune(full, cfg: ImmuneConfig, markers: MarkerConfig, seed=0, compartments=None) -> dict[str, ImmuneResult]` — order: compartments with `handoff` first; passes `extra_cells`/`prior`.

- [ ] **Step 1: Tests**

```python
def test_assign_level1_rules():
    comp = CompartmentConfig(lineages=["T_NK"], level1={}, level2={}, core=["CD8_T"], unresolved="TNK_unresolved",
                             cross={"cross_Myeloid": "Macrophage_Mono"}, unspecified={})
    cfg = ImmuneConfig()
    means = pd.DataFrame({"CD8_T": [-0.2, -0.9, -0.9, -0.9], "B_cell": [0.1, 0.8, 0.1, 0.1],
                          "cross_Myeloid": [0, 0, 0, 0], "spillover_Tumor": [0.2, 1.2, 0.7, 0.1]},
                         index=["c0", "c1", "c2", "c3"])
    lv, flag = assign_level1(means, ["CD8_T", "B_cell", "cross_Myeloid"], ["spillover_Tumor"], comp, cfg, "TNK")
    assert lv == {"c0": "CD8_T", "c1": "B_cell", "c2": "TNK_spillover_Tumor", "c3": "TNK_unresolved"}
    assert flag == {"c0": False, "c1": True, "c2": False, "c3": False}

def test_run_immune_handoff_and_states(synth_lineage):
    res = run_immune(synth_lineage, ImmuneConfig(n_pcs=15, umap_cells=0, resolution=1.0), default_config().markers)
    tnk, b = res["TNK"], res["Bcell"]
    assert "Bcell" in tnk.handoff and len(tnk.handoff["Bcell"]) > 50      # planted B cells sat in T_NK lineage
    moved = list(tnk.handoff["Bcell"])
    assert set(b.adata.obs.loc[moved, "immune_level1"]) <= {"B_cell", "Plasma"}
    assert not tnk.adata.obs_names.isin(moved).any()
    sub = tnk.adata.obs
    assert (sub.loc[sub["truth"] == "Treg", "immune_subtype"] == "Treg").mean() > 0.5
    assert (sub.loc[sub["truth"] == "CD8_cyto", "immune_subtype"] == "CD8_cytotoxic").mean() > 0.5
```

(`synth_lineage` fixture = `synth` with `obs["lineage"]` set from truth, B cells deliberately placed in `T_NK`.)

- [ ] **Step 2–4:** fail → implement → pass.

### Task 8: Merge and spatial stages

**Files:**
- Create: `src/annotation_st/stages/merge.py`, `src/annotation_st/stages/spatial.py`
- Test: `tests/test_merge.py`, `tests/test_spatial.py`

**Interfaces:**
- `merge_labels(full, fb_obs: DataFrame|None, tum_obs: DataFrame|None, immune: dict[str, DataFrame], compartments: dict[str, CompartmentConfig]) -> MergeResult(adata, fine_counts, coarse_counts)`.
- `build_spatial_graph(adata, cfg: SpatialConfig) -> float` (returns edge cut in µm; radius graph when `cfg.radius > 0`).
- `nhood_enrichment(adata, key, cfg, seed) -> (z: DataFrame, count: DataFrame)` (n_jobs=1, NaN→0).
- `neighbourhood_composition(adata, key, groups) -> DataFrame` with `log2_ratio_<g0>_vs_<g1>` when exactly two groups.

- [ ] **Step 1: Tests**

```python
def test_merge_precedence_and_coarse_update():
    full = tiny_adata(cell_type=["Fibroblast","Epithelial_tumor_cycling","T_NK","T_NK","Plasma"])
    fb = pd.DataFrame({"caf_subtype_cell": ["myCAF"], "caf_subtype_cluster": ["myCAF"]}, index=["c0"])
    tum = pd.DataFrame({"csc_status": ["CSC_like"], "csc_OCSC_z": [2.0], "csc_n_core_detected": [3]}, index=["c1"])
    imm = {"TNK": pd.DataFrame({"immune_subtype": ["Treg"], "immune_level1": ["CD4_T"]}, index=["c2"]),
           "Bcell": pd.DataFrame({"immune_subtype": ["Memory_B", "Plasma"], "immune_level1": ["B_cell", "Plasma"]}, index=["c3", "c4"])}
    res = merge_labels(full, fb, tum, imm, default_config().markers.compartments)
    assert list(res.adata.obs["cell_type_fine"]) == ["myCAF", "CSC_like", "Treg", "Memory_B", "Plasma"]
    assert list(res.adata.obs["cell_type"]) == ["Fibroblast", "Epithelial_tumor_cycling", "T_NK", "B_cell", "Plasma"]
    assert list(res.adata.obs["immune_compartment"]) == ["NA", "NA", "TNK", "Bcell", "Bcell"]

def test_graph_edge_cut_and_composition(synth):
    cut = build_spatial_graph(synth, SpatialConfig())
    D = synth.obsp["spatial_distances"]
    assert D.data.max() <= cut
    synth.obs["lab"] = synth.obs["truth_coarse"]
    comp = neighbourhood_composition(synth, "lab", ("Epithelial_tumor", "Fibroblast"))
    assert "log2_ratio_Epithelial_tumor_vs_Fibroblast" in comp.columns
    assert comp.loc["Epithelial_tumor", "Epithelial_tumor"] > comp.loc["Epithelial_tumor", "Fibroblast"]
    z, c = nhood_enrichment(synth, "lab", SpatialConfig(n_perms=20), seed=0)
    assert np.isfinite(z.values).all() and z.loc["Epithelial_tumor", "Epithelial_tumor"] > 0
```

- [ ] **Step 2–4:** fail → implement → pass.

### Task 9: Plotting helpers, pipeline runner, CLI, shipped config, README

**Files:**
- Create: `src/annotation_st/plotting.py` (port of the figure helpers), `src/annotation_st/pipeline.py`, `src/annotation_st/cli.py`, `configs/xenium_ov_5k.yaml`
- Modify: `README.md`
- Test: `tests/test_cli.py`

**Interfaces:**
- `pipeline.STAGES = ["qc", "lineage", "caf", "csc", "immune", "merge", "spatial"]`; `run_stage(name, cfg, paths, force=False, log=None) -> None`; `run(cfg, stages=None, force=False) -> None`; `write_synthetic_xenium(dir, adata)` lives in `tests/conftest.py` (writes `cell_feature_matrix.h5` + `cells.parquet` in 10x layout) so the CLI can be smoke-tested end to end.
- `cli.main(argv=None) -> int`.

- [ ] **Step 1: Tests**

```python
def test_cli_full_run_on_synthetic_xenium(tmp_path, synth_raw):
    xdir = tmp_path / "xenium"; write_synthetic_xenium(xdir, synth_raw)
    cfgp = tmp_path / "cfg.yaml"
    assert main(["init-config", "-o", str(cfgp)]) == 0
    cfg = load_config(cfgp); cfg.project.xenium_dir = str(xdir); cfg.project.outdir = str(tmp_path / "run")
    cfg.project.figures = False; cfg.qc.area_pct = (0, 100)
    cfg.lineage.n_pcs = 20; cfg.lineage.skip_umap = True; cfg.caf.n_pcs = 15; cfg.csc.n_pcs = 15; cfg.csc.umap_cells = 0
    cfg.immune.n_pcs = 15; cfg.immune.umap_cells = 0; cfg.spatial.n_perms = 10
    cfg.markers.caf = {"iCAF": [...], "myCAF": [...], "apCAF": [...]}
    dump_config(cfg, cfgp)
    assert main(["run", "-c", str(cfgp)]) == 0
    run = tmp_path / "run"
    assert (run / "checkpoints" / "sample_06_final.h5ad").exists()
    assert (run / "results" / "final_cell_type_counts.tsv").exists()
    params = json.loads((run / "params.json").read_text())
    assert set(params) >= {"qc_filter", "cluster_annotate", "fibroblast_caf", "cancer_csc", "immune_subcluster", "merge", "spatial_squidpy"}
    # second run skips everything
    assert main(["run", "-c", str(cfgp)]) == 0
```

- [ ] **Step 2–4:** fail → implement → pass. Result-table file names match the original scripts (see spec table and script docstrings) so downstream docs remain valid.
- [ ] **Step 5:** `configs/xenium_ov_5k.yaml` = `dump_config(default_config())` with `project.xenium_dir` set to `/easley/scratch/projects/amitra/amitra2016502/Xenium_5k_Human_OV_FFPE_outs` and `caf.marker_files` pointing at the three Reference_Data txt files. README documents install, config, CLI, library use, output layout, and the mapping to the original scripts.

### Task 10: Full test run and verification

- [ ] `MPLBACKEND=Agg $PY -m pytest tests -q` all green.
- [ ] `annot-st --help` and `annot-st init-config -o /tmp/x.yaml` work from the spatioloji env.
- [ ] Figures path exercised once: rerun the CLI smoke test with `figures: true` on the synthetic data (manual, not in the suite) and confirm png+pdf files appear.
