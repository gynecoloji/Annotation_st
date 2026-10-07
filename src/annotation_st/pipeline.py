"""Stage runners: checkpoint I/O, result tables, figures and ``params.json`` around the
pure stage functions in :mod:`annotation_st.stages`.

Result-table file names match the original ``01_xenium_*.py`` scripts so existing
documentation stays valid.
"""
from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

from . import plotting as P
from .config import ConfigError, PipelineConfig, config_to_dict, validate_config
from .io import RunPaths, append_log, now, read_xenium, setup_logging, update_params, versions, write_json
from .stages.caf import annotate_caf, resolve_caf_sets
from .stages.csc import call_csc, resolve_csc_markers
from .stages.immune import compartment_order, run_immune
from .stages.lineage import annotate_lineage
from .stages.merge import merge_labels
from .stages.qc import QC_METRICS, compute_qc_metrics, qc_filter
from .stages.spatial import build_spatial_graph, n_edges, neighbourhood_composition, nhood_enrichment
from .state import DEPENDS, StageState, selected_compartments, stage_status

STAGES = ["qc", "lineage", "caf", "csc", "immune", "merge", "spatial"]


class StageInputMissing(FileNotFoundError):
    pass


def _need(path: Path, produced_by: str) -> Path:
    if not Path(path).exists():
        raise StageInputMissing(f"{path} not found — run stage '{produced_by}' first")
    return Path(path)


def _read(path: Path):
    import scanpy as sc
    return sc.read_h5ad(path)


def _read_obs(path: Path) -> pd.DataFrame:
    import anndata as ad
    return ad.read_h5ad(path, backed="r").obs.copy()


def input_path(cfg: PipelineConfig, config_dir=None) -> Path | None:
    """The input of the qc stage: the Xenium folder, or the saved spatioloji object (relative to ``config_dir``)."""
    if cfg.project.spatioloji:
        f = Path(cfg.project.spatioloji)
        return f if f.is_absolute() or config_dir is None else Path(config_dir) / f
    return Path(cfg.project.xenium_dir) if cfg.project.xenium_dir else None


def caf_sets_for(cfg: PipelineConfig, base: str | Path | None = None) -> dict[str, list[str]]:
    """Inline ``markers.caf`` if given, else the files in ``caf.marker_files`` (relative to ``base``)."""
    return resolve_caf_sets(cfg.caf, cfg.markers.caf, base=base)


def _selected_compartments(cfg: PipelineConfig) -> list[str]:
    return selected_compartments(cfg)


def _done(*outputs) -> dict:
    return {"status": "completed", "outputs": [str(o) for o in outputs], "note": ""}


def _skipped(note: str) -> dict:
    return {"status": "skipped", "outputs": [], "note": note}


# ----------------------------------------------------------------------------
# colours and the maps that show them (shared by the stage runners and `annot-st figures`)
# ----------------------------------------------------------------------------
def colour_columns(cfg: PipelineConfig, adata, keys, config_dir=None, aliases: dict[str, str] | None = None) -> dict:
    """Give the label columns ``keys`` of ``adata`` the run's colours (``uns[key + '_colors']``),
    so every UMAP / tissue map of every stage shows a label in the colour it has in the final
    map. ``aliases`` maps helper labels to the label whose colour they share
    (``cross_Myeloid`` -> ``Macrophage_Mono``). Returns the palette."""
    keys = [k for k in keys if k in adata.obs]
    for k in keys:
        if not isinstance(adata.obs[k].dtype, pd.CategoricalDtype):
            adata.obs[k] = adata.obs[k].astype(str).astype("category")
    aliases = aliases or {}
    labels = [str(c) for k in keys for c in adata.obs[k].cat.categories if str(c) not in aliases]
    pal = label_palette(cfg, labels + list(aliases.values()), config_dir)
    for alias, target in aliases.items():
        pal[alias] = pal[target]
    for k in keys:
        P.apply_palette(adata, k, pal)
    return pal


def draw_lineage_maps(adata, paths: RunPaths) -> None:
    P.umap_plot(adata, "cell_type", paths.figures / "umap_cell_type", title="Lineage annotation")
    P.spatial_plot(adata, "cell_type", paths.figures / "spatial_cell_type", size=1, title="Lineage annotation (tissue)")


def draw_caf_maps(fb, paths: RunPaths) -> None:
    cols = [c for c in ("fb_leiden", "caf_subtype_cell", "caf_subtype_cluster") if c in fb.obs]
    P.umap_plot(fb, cols, paths.figures / "umap_fb_caf_subtype", ncols=3, wspace=0.4)
    P.spatial_plot(fb, "caf_subtype_cell", paths.figures / "spatial_fb_caf_subtype", size=3,
                   title="CAF subtype (fibroblasts only)")


def draw_csc_maps(tum, paths: RunPaths) -> None:
    cols = [c for c in ("tumor_leiden", "csc_status", "csc_OCSC_z") if c in tum.obs]
    P.umap_plot(tum, cols, paths.figures / "umap_tumor_csc", ncols=3, cmap="magma", vmax=3, wspace=0.4)
    P.spatial_plot(tum, "csc_status", paths.figures / "spatial_tumor_csc", size=2, title="CSC-like tumour cells (tissue)")


def draw_immune_maps(sub, name: str, paths: RunPaths) -> None:
    cols = [c for c in (f"{name}_leiden", "immune_level1", "immune_subtype") if c in sub.obs]
    P.umap_plot(sub, cols, paths.figures / f"umap_{name}_subtype", ncols=3, wspace=0.45)
    P.spatial_plot(sub, "immune_subtype", paths.figures / f"spatial_{name}_subtype", size=3,
                   title=f"{name} functional subtypes (tissue)")


IMMUNE_KEYS = ("immune_subtype", "immune_level1")


# ----------------------------------------------------------------------------
# stage runners
# ----------------------------------------------------------------------------
def run_qc(cfg: PipelineConfig, paths: RunPaths, log: logging.Logger, config_dir=None, state=None) -> dict:
    seed = cfg.project.seed
    source = input_path(cfg, config_dir)
    if source is None:
        raise ValueError("no input: set project.xenium_dir or project.spatioloji")
    log.info("reading %s", source)
    if cfg.project.spatioloji:
        from .spatioloji_bridge import read_spatioloji
        adata = read_spatioloji(source, subsample=cfg.qc.subsample, seed=seed)
    else:
        adata = read_xenium(source, subsample=cfg.qc.subsample, seed=seed)
    n_in = adata.n_obs
    log.info("matrix: %d cells x %d genes", adata.n_obs, adata.n_vars)
    compute_qc_metrics(adata)
    if cfg.project.figures:
        P.qc_histograms(adata, paths.qc / "qc_hist_prefilter", f"per-cell QC before filtering (n={adata.n_obs:,})",
                        [("gene_counts", True), ("n_genes_by_counts", True), ("cell_area", True),
                         ("nucleus_area", True), ("control_frac", False), ("counts_per_area", True)])
        P.qc_spatial_maps(adata, paths.qc / "qc_spatial_prefilter", seed=seed)
    res = qc_filter(adata, cfg.qc)
    if not cfg.qc.apply_filters:
        log.info("QC filtering is off (qc.apply_filters: false): every cell and gene is kept; "
                 "qc_filter_summary.tsv reports what the thresholds would remove")
    res.summary_pre.to_csv(paths.qc / "qc_metrics_summary_prefilter.tsv", sep="\t")
    res.summary_post.to_csv(paths.qc / "qc_metrics_summary_postfilter.tsv", sep="\t")
    res.flag_summary.to_csv(paths.qc / "qc_filter_summary.tsv", sep="\t")
    log.info("filter summary:\n%s", res.flag_summary.to_string())
    out = res.adata
    if cfg.project.figures:
        P.qc_histograms(out, paths.qc / "qc_hist_postfilter", f"after filtering (n={out.n_obs:,})",
                        [("gene_counts", True), ("n_genes_by_counts", True), ("cell_area", True)], colour="#55A868")
    out.write_h5ad(paths.checkpoint("qc"), compression="gzip")
    log.info("kept %d cells, %d genes -> %s", out.n_obs, out.n_vars, paths.checkpoint("qc"))
    update_params(paths, "qc_filter", {
        "timestamp": now(), "input": str(source), "input_type": "spatioloji" if cfg.project.spatioloji else "xenium",
        "xenium_dir": str(cfg.project.xenium_dir), "subsample": cfg.qc.subsample, "seed": seed,
        "apply_filters": cfg.qc.apply_filters,
        "n_cells_input": int(n_in), "n_cells_kept": int(out.n_obs), "n_genes_kept": int(out.n_vars),
        "min_transcripts": cfg.qc.min_transcripts, "min_genes": cfg.qc.min_genes,
        "max_control_frac": cfg.qc.max_control_frac, "area_pct": list(cfg.qc.area_pct),
        "area_bounds_um2": list(res.area_bounds), "require_nucleus": cfg.qc.require_nucleus,
        "min_cells_per_gene": cfg.qc.min_cells_per_gene, "versions": versions()})
    append_log(paths, f"qc: {out.n_obs} cells kept of {n_in} -> {paths.checkpoint('qc').name}")
    return _done(paths.checkpoint("qc"))


def run_lineage(cfg: PipelineConfig, paths: RunPaths, log: logging.Logger, config_dir=None, state=None) -> dict:
    from dataclasses import replace

    adata = _read(_need(paths.checkpoint("qc"), "qc"))
    log.info("loaded %d cells x %d genes", adata.n_obs, adata.n_vars)
    lcfg = replace(cfg.lineage, skip_umap=cfg.lineage.skip_umap or not cfg.project.figures)  # UMAP is for figures only
    res = annotate_lineage(adata, lcfg, cfg.markers, seed=cfg.project.seed)
    log.info("leiden: %d clusters\n%s", len(res.means), res.counts.to_string())
    write_json(paths.results / "marker_sets_used.json",
               {"lineage": res.sets_used, "lineage_dropped_not_on_panel": res.sets_dropped, "proliferation": res.prolif_used})
    res.means.to_csv(paths.results / "leiden_lineage_scores.tsv", sep="\t")
    res.counts.to_csv(paths.results / "cell_type_counts.tsv", sep="\t")
    pd.DataFrame({"cluster": list(res.mapping), "lineage": list(res.mapping.values())}).to_csv(
        paths.results / "lineage_assignment.tsv", sep="\t", index=False)
    if res.top_markers is not None:
        res.top_markers.to_csv(paths.results / "leiden_top_markers.tsv", sep="\t", index=False)
    pal = label_palette(cfg, adata.obs["cell_type"].cat.categories, config_dir)
    P.apply_palette(adata, "cell_type", pal)             # used by the figures below and kept in the h5ad
    if not (paths.results / COLOR_TABLE).exists() or state is None or not state.completed("merge"):
        write_color_table(cfg, paths, pal, coarse=adata.obs["cell_type"])
    if cfg.project.figures:
        P.umap_plot(adata, "leiden", paths.figures / "umap_leiden", title="Leiden clusters")
        draw_lineage_maps(adata, paths)
        P.dotplot_save(adata, {k: v[:6] for k, v in res.sets_used.items()}, "cell_type",
                       paths.figures / "dotplot_lineage_markers", use_raw=True, standard_scale="var")
        m = res.means.copy(); m["assigned"] = m["assigned"].astype(str)
        P.score_heatmap(m, res.scored, paths.figures / "heatmap_cluster_lineage_scores",
                        title="Leiden cluster × lineage marker score")
    adata.write_h5ad(paths.checkpoint("lineage"), compression="gzip")
    log.info("wrote %s", paths.checkpoint("lineage"))
    update_params(paths, "cluster_annotate", {
        "timestamp": now(), "normalization": "normalize_total(median)+log1p", "n_pcs": cfg.lineage.n_pcs,
        "n_neighbors": cfg.lineage.n_neighbors, "leiden_resolution": cfg.lineage.resolution, "leiden_flavor": "igraph",
        "n_clusters": int(len(res.means)), "umap": not cfg.lineage.skip_umap, "umap_cells": cfg.lineage.umap_cells,
        "min_z": cfg.lineage.min_z, "cycling_z": cfg.lineage.cycling_z, "lineage_marker_sets": res.sets_used,
        "proliferation_markers": res.prolif_used, "cell_type_counts": res.counts["n_cells"].to_dict(),
        "seed": cfg.project.seed})
    append_log(paths, f"lineage: {len(res.means)} Leiden clusters, {len(res.counts)} cell types -> {paths.checkpoint('lineage').name}")
    return _done(paths.checkpoint("lineage"))


def run_caf(cfg: PipelineConfig, paths: RunPaths, log: logging.Logger, config_dir=None, state=None) -> dict:
    full = _read(_need(paths.checkpoint("lineage"), "lineage"))
    fb = full[full.obs["lineage"].isin(list(cfg.caf.include_lineages)).values].copy()
    del full
    log.info("fibroblast compartment: %d cells", fb.n_obs)
    raw_sets = caf_sets_for(cfg, base=config_dir)
    if not raw_sets:
        log.info("no CAF marker sets opted in (markers.caf or caf.marker_files); CAF stage skipped")
        return _skipped("no CAF marker sets opted in")
    res = annotate_caf(fb, cfg.caf, raw_sets, seed=cfg.project.seed, umap=cfg.project.figures)
    if res is None:
        log.warning("fewer than %d fibroblasts; CAF stage skipped", cfg.caf.min_cells)
        return _skipped(f"fewer than {cfg.caf.min_cells} cells in lineages {list(cfg.caf.include_lineages)}")
    write_json(paths.results / "caf_marker_sets_used.json",
               {"n_in_file": {k: len(v) for k, v in raw_sets.items()}, "used_on_panel": res.sets_used,
                "dropped_not_on_panel": res.sets_dropped})
    for k in res.sets_used:
        log.info("%s: %d/%d markers on panel", k, len(res.sets_used[k]), len(raw_sets[k]))
    log.info("CAF subtype counts:\n%s", res.counts.to_string())
    res.means.to_csv(paths.results / "fb_subcluster_caf_scores.tsv", sep="\t")
    res.counts.to_csv(paths.results / "caf_subtype_counts.tsv", sep="\t")
    fb.obs[["fb_leiden", "caf_subtype_cell", "caf_subtype_cluster", "caf_best_z", "caf_margin"]
           + [f"caf_{s}" for s in res.scored]].to_csv(paths.results / "caf_subtype_by_cell.tsv.gz", sep="\t")
    colour_columns(cfg, fb, ["caf_subtype_cell", "caf_subtype_cluster"], config_dir)
    if cfg.project.figures:
        draw_caf_maps(fb, paths)
        P.umap_plot(fb, [f"caf_{s}_z" for s in res.scored], paths.figures / "umap_fb_caf_scores", ncols=3, cmap="magma", vmax=3)
        P.dotplot_save(fb, {k: v[:12] for k, v in res.sets_used.items()}, "caf_subtype_cell",
                       paths.figures / "dotplot_fb_caf_markers", use_raw=False, standard_scale="var")
        P.violin_save(fb, [f"caf_{s}_z" for s in res.scored], "caf_subtype_cell", paths.figures / "violin_fb_caf_scores")
        P.score_heatmap(res.means, res.scored, paths.figures / "heatmap_fb_subcluster_caf_scores",
                        title="Fibroblast sub-cluster × CAF signature", vmin=-1.5, vmax=1.5)
    fb.write_h5ad(paths.checkpoint("caf"), compression="gzip")
    log.info("wrote %s", paths.checkpoint("caf"))
    update_params(paths, "fibroblast_caf", {
        "timestamp": now(), "fb_lineages": list(cfg.caf.include_lineages), "n_fibroblasts": int(fb.n_obs),
        "n_pcs": cfg.caf.n_pcs, "n_neighbors": cfg.caf.n_neighbors, "leiden_resolution": cfg.caf.resolution,
        "n_subclusters": int(fb.obs["fb_leiden"].nunique()), "min_z": cfg.caf.min_z, "margin": cfg.caf.margin,
        "caf_markers_on_panel": res.sets_used, "caf_markers_dropped": res.sets_dropped,
        "caf_subtype_counts_per_cell": res.counts["per_cell"].to_dict(), "seed": cfg.project.seed})
    append_log(paths, f"caf: {fb.n_obs} fibroblasts, {res.counts['per_cell'].to_dict()} -> {paths.checkpoint('caf').name}")
    return _done(paths.checkpoint("caf"))


def run_csc(cfg: PipelineConfig, paths: RunPaths, log: logging.Logger, config_dir=None, state=None) -> dict:
    full = _read(_need(paths.checkpoint("lineage"), "lineage"))
    tum = full[full.obs["lineage"].isin(list(cfg.csc.tumor_lineages)).values].copy()
    del full
    log.info("tumour epithelial compartment: %d cells", tum.n_obs)
    csc_genes, core_genes, source = resolve_csc_markers(cfg.csc, cfg.markers, base=config_dir)
    if not csc_genes:
        log.info("no CSC marker set opted in (markers.csc or csc.marker_file); CSC stage skipped")
        return _skipped("no CSC marker set opted in")
    log.info("CSC markers (%d, from %s); core (%d, from %s)", len(csc_genes), source["csc"], len(core_genes), source["core"])
    if source.get("core_not_in_set"):
        log.warning("core markers not in the CSC set (they still count for the detection rule): %s", source["core_not_in_set"])
    res = call_csc(tum, cfg.csc, seed=cfg.project.seed, umap=cfg.project.figures, csc_genes=csc_genes, core_genes=core_genes)
    if res is None:
        log.warning("fewer than %d tumour cells; CSC stage skipped", cfg.csc.min_cells)
        return _skipped(f"fewer than {cfg.csc.min_cells} cells in lineages {list(cfg.csc.tumor_lineages)}")
    write_json(paths.results / "csc_marker_sets_used.json",
               {"source": source, "requested": {"OCSC": csc_genes, "OCSC_core": core_genes},
                "used_on_panel": res.sets_used, "dropped_not_on_panel": res.sets_dropped})
    log.info("CSC status:\n%s", res.counts.to_string())
    res.counts.to_csv(paths.results / "csc_counts.tsv", sep="\t")
    res.positivity.to_csv(paths.results / "csc_marker_positivity.tsv", sep="\t")
    if res.de is not None:
        res.de.to_csv(paths.results / "csc_vs_nonCSC_markers.tsv", sep="\t", index=False)
    tum.obs[["csc_OCSC", "csc_OCSC_z", "csc_n_core_detected", "csc_status"]].to_csv(paths.results / "csc_by_cell.tsv.gz", sep="\t")
    colour_columns(cfg, tum, ["csc_status"], config_dir)
    if cfg.project.figures:
        draw_csc_maps(tum, paths)
        P.dotplot_save(tum, res.sets_used["OCSC"], "csc_status", paths.figures / "dotplot_tumor_csc_markers", use_raw=False, standard_scale="var")
        P.histogram_with_threshold(tum.obs["csc_OCSC_z"], cfg.csc.z_thresh, paths.figures / "tumor_csc_score_hist",
                                   "OCSC score (z)", "tumour cells", f"z > {cfg.csc.z_thresh}")
    tum.write_h5ad(paths.checkpoint("csc"), compression="gzip")
    log.info("wrote %s", paths.checkpoint("csc"))
    update_params(paths, "cancer_csc", {
        "timestamp": now(), "tumor_lineages": list(cfg.csc.tumor_lineages), "n_tumor_cells": int(tum.n_obs),
        "ocsc_markers_on_panel": res.sets_used["OCSC"], "ocsc_markers_dropped": res.sets_dropped["OCSC"],
        "core_markers": res.sets_used["OCSC_core"], "marker_source": source,
        "z_thresh": cfg.csc.z_thresh, "min_core": cfg.csc.min_core,
        "n_csc_like": res.n_csc, "frac_csc_like": res.n_csc / tum.n_obs, "umap_cells": cfg.csc.umap_cells,
        "seed": cfg.project.seed})
    append_log(paths, f"csc: {res.n_csc}/{tum.n_obs} tumour cells CSC-like -> {paths.checkpoint('csc').name}")
    return _done(paths.checkpoint("csc"))


def run_immune_stage(cfg: PipelineConfig, paths: RunPaths, log: logging.Logger, config_dir=None, state=None) -> dict:
    if not _selected_compartments(cfg):
        log.info("no immune compartment opted in (markers.compartments); immune stage skipped")
        return _skipped("no immune compartment opted in")
    full = _read(_need(paths.checkpoint("lineage"), "lineage"))
    results = run_immune(full, cfg.immune, cfg.markers, seed=cfg.project.seed, umap=cfg.project.figures, log=log)
    del full
    used, per_cell, summary, written = {}, [], {}, []
    for name, res in results.items():
        sub = res.adata
        res.means1.to_csv(paths.results / f"immune_{name}_subcluster_scores.tsv", sep="\t")
        if res.level2_tables:
            pd.concat(res.level2_tables.values()).to_csv(paths.results / f"immune_{name}_level2_scores.tsv", sep="\t")
        res.counts.to_csv(paths.results / f"immune_{name}_subtype_counts.tsv", sep="\t")
        if res.top_markers is not None:
            res.top_markers.to_csv(paths.results / f"immune_{name}_top_markers.tsv", sep="\t", index=False)
        for target, ids in res.handoff.items():
            pd.Series(ids, name="prior").rename_axis("cell_id").to_csv(
                paths.results / f"immune_handoff_{name}_to_{target}.tsv", sep="\t")
        log.info("[%s] final subtypes:\n%s", name, res.counts.to_string())
        colour_columns(cfg, sub, IMMUNE_KEYS, config_dir, aliases=cfg.markers.compartments[name].cross)
        if cfg.project.figures:
            draw_immune_maps(sub, name, paths)
            l1 = res.sets_used["level1"]["used_on_panel"]
            dot_sets = {k: v[:6] for k, v in l1.items() if not k.startswith("spillover_")}
            for parent, info in res.sets_used["level2"].items():
                for st, genes in info["used_on_panel"].items():
                    if genes:
                        dot_sets[st] = genes[:6]
            P.dotplot_save(sub, dot_sets, "immune_subtype", paths.figures / f"dotplot_{name}_markers", use_raw=False, standard_scale="var")
            m1 = res.means1.copy(); m1["assigned"] = m1["assigned"].astype(str)
            P.score_heatmap(m1, res.scored1, paths.figures / f"heatmap_{name}_subcluster_scores",
                            title=f"{name} sub-cluster × identity / spillover signature", vmin=-1.5, vmax=1.5)
            for parent, m2 in res.level2_tables.items():
                cols = [c for c in m2.columns if c not in ("parent", "assigned", "n_cells")]
                mm = m2.copy(); mm["assigned"] = mm["assigned"].astype(str)
                P.score_heatmap(mm, cols, paths.figures / f"heatmap_{name}_level2_{parent}",
                                title=f"{name}: {parent} sub-cluster × state", vmin=-1.5, vmax=1.5)
        h5 = paths.checkpoint(f"immune_{name}")
        sub.write_h5ad(h5, compression="gzip")
        written.append(h5)
        log.info("[%s] wrote %s", name, h5)
        used[name] = res.sets_used
        summary[name] = res.counts["n_cells"].to_dict()
        pc = sub.obs[[f"{name}_leiden", "immune_level1", "immune_subtype", "immune_state_cell", "immune_spillover"]].copy()
        pc.columns = ["subcluster", "immune_level1", "immune_subtype", "immune_state_cell", "immune_spillover"]
        pc.insert(0, "compartment", name)
        per_cell.append(pc)
    write_json(paths.results / "immune_marker_sets_used.json", used)
    if per_cell:
        pd.concat(per_cell).to_csv(paths.results / "immune_subtype_by_cell.tsv.gz", sep="\t")
    ic = cfg.immune
    update_params(paths, "immune_subcluster", {
        "timestamp": now(), "compartments": {k: cfg.markers.compartments[k].lineages for k in results},
        "n_pcs": ic.n_pcs, "n_neighbors": ic.n_neighbors, "leiden_resolution": ic.resolution, "min_z_spillover": ic.min_z,
        "low_identity": ic.low_identity, "resolution_level2": ic.resolution2, "min_z_level2": ic.min_z2,
        "spill_flag": ic.spill_flag, "umap_cells": ic.umap_cells, "markers_on_panel": used, "subtype_counts": summary,
        "seed": cfg.project.seed})
    append_log(paths, "immune: " + "; ".join(f"{k}: {len(v)} subtypes" for k, v in summary.items()))
    if not written:
        return _skipped(f"no compartment had at least {cfg.immune.min_cells} cells")
    return _done(*written)


def run_merge(cfg: PipelineConfig, paths: RunPaths, log: logging.Logger, config_dir=None, state=None) -> dict:
    adata = _read(_need(paths.checkpoint("lineage"), "lineage"))
    def usable(stage: str, h5: Path) -> bool:
        # A checkpoint left behind by an earlier configuration (the stage has since skipped
        # itself, or never ran under this state file) must not leak into the merged labels.
        if not h5.exists():
            return False
        return True if state is None else (state.completed(stage) and str(h5) in state.outputs(stage))

    fb = _read_obs(paths.checkpoint("caf")) if usable("caf", paths.checkpoint("caf")) else None
    tum = _read_obs(paths.checkpoint("csc")) if usable("csc", paths.checkpoint("csc")) else None
    immune = {}
    for name, h5 in paths.immune_checkpoints(_selected_compartments(cfg)).items():
        if usable("immune", h5):
            immune[name] = _read_obs(h5)
        else:
            log.warning("no current immune result for %s; keeping coarse labels", name)
    for what, tab in (("caf", fb), ("csc", tum)):
        if tab is None:
            log.warning("no current %s result; keeping coarse labels", what)
    res = merge_labels(adata, fb, tum, immune, cfg.markers.compartments)
    res.fine_counts.to_csv(paths.results / "final_cell_type_counts.tsv", sep="\t")
    res.coarse_counts.to_csv(paths.results / "final_cell_type_coarse_counts.tsv", sep="\t")
    log.info("final cell types:\n%s", res.fine_counts.to_string())
    cols = [c for c in ["leiden", "lineage", "cell_type", "cell_type_fine", "caf_subtype_cell", "immune_compartment",
                        "immune_subtype", "csc_status", "csc_OCSC_z"] if c in adata.obs]
    adata.obs[cols].to_csv(paths.results / "final_annotation_by_cell.tsv.gz", sep="\t")
    pal = label_palette(cfg, list(adata.obs["cell_type_fine"].cat.categories) + list(adata.obs["cell_type"].cat.categories),
                        config_dir)
    P.apply_palette(adata, "cell_type", pal)
    P.apply_palette(adata, "cell_type_fine", pal)
    write_color_table(cfg, paths, pal, coarse=adata.obs["cell_type"], fine=adata.obs["cell_type_fine"])
    if cfg.project.figures:
        draw_fine_map(cfg, paths, adata.obsm["spatial"], adata.obs["cell_type_fine"].to_numpy(), pal, log)
    adata.write_h5ad(paths.checkpoint("final"), compression="gzip")
    log.info("wrote %s", paths.checkpoint("final"))
    update_params(paths, "merge", {"timestamp": now(), "n_fine_labels": int(len(res.fine_counts)),
                                   "n_coarse_labels": int(len(res.coarse_counts)),
                                   "final_cell_type_counts": res.fine_counts["n_cells"].to_dict(),
                                   "inputs": {"caf": fb is not None, "csc": tum is not None, "immune": list(immune)}})
    append_log(paths, f"merge: {len(res.fine_counts)} fine / {len(res.coarse_counts)} coarse labels -> {paths.checkpoint('final').name}")
    return _done(paths.checkpoint("final"))


def run_spatial(cfg: PipelineConfig, paths: RunPaths, log: logging.Logger, config_dir=None, state=None) -> dict:
    adata = _read(_need(paths.checkpoint("final"), "merge"))
    sp = cfg.spatial
    key = sp.fine_key
    cut = build_spatial_graph(adata, sp)
    log.info("spatial graph: %d edges (edge cut at %.1f um)", n_edges(adata), cut)
    z, count = nhood_enrichment(adata, key, sp, seed=cfg.project.seed)
    z.to_csv(paths.results / "nhood_enrichment_zscore.tsv", sep="\t")
    count.to_csv(paths.results / "nhood_enrichment_count.tsv", sep="\t")
    groups = list(sp.composition_groups)
    ncomp = neighbourhood_composition(adata, key, groups)
    ncomp.to_csv(paths.results / "csc_neighbourhood_composition.tsv", sep="\t")
    log.info("neighbourhood composition:\n%s", ncomp.round(4).to_string())
    if cfg.project.figures:
        P.nhood_heatmap(adata, key, paths.figures / "nhood_enrichment_fine")
        if ncomp.shape[1] >= 2 and len(groups) == 2:
            caf_labels = [c for c in ["iCAF", "myCAF", "apCAF", cfg.caf.unassigned] if c in ncomp.index]
            lab_col = dict(zip(map(str, adata.obs[key].cat.categories), adata.uns.get(f"{key}_colors", [])))
            series = [lab_col[g] for g in groups] if all(g in lab_col for g in groups) else None
            if caf_labels:
                P.grouped_bar(ncomp.loc[caf_labels, groups], paths.figures / "csc_neighbour_caf_fraction",
                              ylabel="mean fraction of neighbours", title="CAF subtypes around tumour cells", colors=series)
            others = [c for c in ncomp.index if c not in caf_labels and c not in groups and ncomp.loc[c, groups].max() > 0.002]
            if others:
                P.grouped_bar(ncomp.loc[others, groups], paths.figures / "csc_neighbour_other_fraction",
                              ylabel="mean fraction of neighbours", title="Other neighbours of tumour cells", colors=series)
    adata.write_h5ad(paths.checkpoint("final"), compression="gzip")
    log.info("wrote %s", paths.checkpoint("final"))
    update_params(paths, "spatial_squidpy", {
        "timestamp": now(), "graph": "delaunay" if sp.radius <= 0 else f"radius_{sp.radius}um",
        "edge_cut_um": cut, "edge_cut_pct": sp.edge_cut_pct, "n_edges": n_edges(adata), "n_perms": sp.n_perms,
        "label_key": key, "n_labels": int(len(z)), "composition_groups": groups, "seed": cfg.project.seed,
        "versions": versions()})
    append_log(paths, f"spatial: nhood enrichment on {len(z)} '{key}' labels")
    return _done(paths.checkpoint("final"), paths.results / "nhood_enrichment_zscore.tsv")


def pinned_colors(cfg: PipelineConfig, config_dir=None) -> dict[str, str]:
    """Fixed label colours: ``plot.palette_file`` (relative to ``config_dir``) overlaid by ``plot.colors``."""
    pinned: dict[str, str] = {}
    if cfg.plot.palette_file:
        f = Path(cfg.plot.palette_file)
        if not f.is_absolute() and config_dir is not None:
            f = Path(config_dir) / f
        pinned.update(P.read_palette_file(f))
    pinned.update(cfg.plot.colors)
    return pinned


COLOR_TABLE = "cell_type_colors.csv"


def label_palette(cfg: PipelineConfig, labels=(), config_dir=None) -> dict[str, str]:
    """The run's label -> colour map for rough and fine cell types alike. Fine types are
    mutually distinct inside their legend group, rough types (the lineage names) among
    themselves. Depends only on the config (legend groups, lineage names, pinned colours),
    so every stage and every later run with the same config gets the same colours."""
    coarse = list(cfg.markers.lineage) + [cfg.lineage.cycling_label, cfg.lineage.unassigned]
    return P.fine_palette(labels, cfg.plot.legend_groups, pinned_colors(cfg, config_dir), cfg.plot.qc_group,
                          also_distinct=[coarse], generic=cfg.plot.generic_labels)


def write_color_table(cfg: PipelineConfig, paths: RunPaths, pal: dict[str, str], coarse=None, fine=None) -> Path:
    """``results/cell_type_colors.csv``: the colour of every rough / fine cell type of this run
    (see :func:`annotation_st.plotting.color_table`). Point ``plot.palette_file`` at a copy of
    it to reproduce exactly these colours in another run."""
    out = paths.results / COLOR_TABLE
    P.color_table(pal, coarse=coarse, fine=fine, groups=cfg.plot.legend_groups, qc_group=cfg.plot.qc_group,
                  generic=cfg.plot.generic_labels).to_csv(out, index=False)
    return out


def draw_fine_map(cfg: PipelineConfig, paths: RunPaths, xy, labels, pal: dict[str, str],
                  log: logging.Logger | None = None) -> Path:
    """``figures/spatial_cell_type_fine`` (png + pdf): tissue map with the hierarchical legend
    of ``cfg.plot`` in the colours of ``pal``."""
    labels = np.asarray(labels).astype(str)
    groups, qc = cfg.plot.legend_groups, cfg.plot.qc_group
    fig = P.build_fine_map(xy, labels, pal, groups, qc_group=qc)
    stem = paths.figures / "spatial_cell_type_fine"
    P.save_fig(fig, stem)
    if log:
        present = np.unique(labels)
        log.info("drew %s (.png, .pdf): %d labels in %d legend groups", stem, len(present),
                 len(P.legend_layout(present, groups, qc)))
    return stem


FIGURES = ["cell_type", "caf", "csc", "immune", "cell_type_fine"]


def _light(h5: Path):
    """Labels, colours and coordinates of a checkpoint without its expression matrix."""
    import anndata as ad

    a = ad.read_h5ad(h5, backed="r")
    try:
        light = ad.AnnData(obs=a.obs.copy())
        for k in ("spatial", "X_umap"):
            if k in a.obsm:
                light.obsm[k] = np.asarray(a.obsm[k])
    finally:
        a.file.close()
    return light


def redraw_figures(cfg: PipelineConfig, which=None, config_dir=None) -> list[Path]:
    """Redraw the label maps (UMAPs and tissue maps of rough types, CAF subtypes, CSC status,
    immune subtypes and the final annotation) from existing checkpoints, in the colours of the
    current config, without re-running any stage. Figures that need expression data (dot
    plots, score heat maps) are not touched. With ``which=None`` everything available is
    drawn; an explicitly requested figure whose checkpoint is missing is an error."""
    explicit = which is not None
    which = list(which) if explicit else list(FIGURES)
    unknown = set(which) - set(FIGURES)
    if unknown:
        raise KeyError(f"unknown figure(s) {sorted(unknown)} (choose from {FIGURES})")
    paths = RunPaths(cfg)
    log = setup_logging(paths, "figures")
    state = StageState(paths)
    drawn: list[Path] = []

    def source(stage: str, h5: Path, what: str):
        """The light checkpoint to draw ``what`` from, or None when it is not there."""
        if not h5.exists():
            if explicit:
                raise StageInputMissing(f"{h5} not found — run stage '{stage}' first")
            return None
        status, reason = stage_status(stage, cfg, paths, state, config_dir)
        if status != "up_to_date":
            log.warning("%s: stage '%s' is not up to date (%s); drawing its labels as they are", what, stage, reason or status)
        return _light(h5)

    if "cell_type" in which and (a := source("lineage", paths.checkpoint("lineage"), "cell_type")) is not None:
        colour_columns(cfg, a, ["cell_type"], config_dir)
        draw_lineage_maps(a, paths)
        drawn += [paths.figures / "umap_cell_type", paths.figures / "spatial_cell_type"]
    if "caf" in which and (a := source("caf", paths.checkpoint("caf"), "caf")) is not None:
        colour_columns(cfg, a, ["caf_subtype_cell", "caf_subtype_cluster"], config_dir)
        draw_caf_maps(a, paths)
        drawn += [paths.figures / "umap_fb_caf_subtype", paths.figures / "spatial_fb_caf_subtype"]
    if "csc" in which and (a := source("csc", paths.checkpoint("csc"), "csc")) is not None:
        colour_columns(cfg, a, ["csc_status"], config_dir)
        draw_csc_maps(a, paths)
        drawn += [paths.figures / "umap_tumor_csc", paths.figures / "spatial_tumor_csc"]
    if "immune" in which:
        for name, h5 in paths.immune_checkpoints(_selected_compartments(cfg)).items():
            if not h5.exists():
                continue
            a = source("immune", h5, f"immune {name}")
            colour_columns(cfg, a, IMMUNE_KEYS, config_dir, aliases=cfg.markers.compartments[name].cross)
            draw_immune_maps(a, name, paths)
            drawn += [paths.figures / f"umap_{name}_subtype", paths.figures / f"spatial_{name}_subtype"]
    if "cell_type_fine" in which and (a := source("merge", paths.checkpoint("final"), "cell_type_fine")) is not None:
        fine, coarse = a.obs["cell_type_fine"].astype(str).to_numpy(), a.obs["cell_type"].astype(str).to_numpy()
        pal = label_palette(cfg, list(np.unique(fine)) + list(np.unique(coarse)), config_dir)
        log.info("colour table: %s", write_color_table(cfg, paths, pal, coarse=coarse, fine=fine))
        drawn.append(draw_fine_map(cfg, paths, a.obsm["spatial"], fine, pal, log))
    if not drawn:
        raise StageInputMissing(f"no checkpoint to draw from in {paths.checkpoints} — run the stages first "
                                f"(the final map needs stage 'merge')")
    return drawn


def export_palette(cfg: PipelineConfig, out: Path, config_dir=None, h5ad: Path | None = None, keys=(),
                   write_h5ad: bool = False, map_stem: Path | None = None, map_key: str | None = None) -> dict[str, str]:
    """Write the colour scheme as a CSV (label, hex, tier, legend_group) and optionally apply
    it to a dataset annotated elsewhere: the labels of ``keys`` in ``h5ad`` are added to the
    scheme (known labels keep their colour, new ones are placed by the same rules),
    ``write_h5ad`` stores ``uns[key + '_colors']`` in that file, and ``map_stem`` draws the
    hierarchical-legend tissue map of ``map_key``."""
    import anndata as ad

    cats: dict[str, list[str]] = {}
    xy = map_labels = None
    if h5ad is not None:
        a = ad.read_h5ad(h5ad, backed="r")
        try:
            for k in keys:
                if k not in a.obs:
                    raise KeyError(f"{h5ad}: obs has no column '{k}'")
                col = a.obs[k]
                cats[k] = ([str(c) for c in col.cat.categories] if isinstance(col.dtype, pd.CategoricalDtype)
                           else sorted(col.astype(str).unique()))
            if map_stem is not None:
                map_key = map_key or (keys[0] if keys else None)
                if map_key is None or "spatial" not in a.obsm:
                    raise ValueError("--map needs a label column (--keys / --map-key) and obsm['spatial']")
                map_labels = a.obs[map_key].astype(str).to_numpy()
                xy = np.asarray(a.obsm["spatial"])
        finally:
            a.file.close()
    pal = label_palette(cfg, [l for ls in cats.values() for l in ls], config_dir)
    table = P.scheme_table(pal, cfg.plot.legend_groups, cfg.plot.qc_group, cfg.plot.generic_labels)
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(out, index=False)
    if write_h5ad and h5ad is not None:
        import h5py
        try:
            from anndata.io import write_elem
        except ImportError:                                  # anndata < 0.11
            from anndata.experimental import write_elem
        with h5py.File(h5ad, "r+") as f:
            g = f.require_group("uns")
            for k, labels in cats.items():
                name = f"{k}_colors"
                if name in g:
                    del g[name]
                write_elem(g, name, np.array([pal[l] for l in labels], dtype=object))
    if map_stem is not None and xy is not None:
        P.save_fig(P.build_fine_map(xy, map_labels, pal, cfg.plot.legend_groups, title=str(map_key),
                                    qc_group=cfg.plot.qc_group), map_stem)
    return pal


RUNNERS = {"qc": run_qc, "lineage": run_lineage, "caf": run_caf, "csc": run_csc, "immune": run_immune_stage,
           "merge": run_merge, "spatial": run_spatial}


class StaleUpstream(RuntimeError):
    """A stage was asked to run on top of an upstream stage that is out of date."""


def preflight(cfg: PipelineConfig, stages=None, config_dir=None, check_inputs: bool = True) -> None:
    """Fail early, before any stage runs, with every problem that concerns ``stages``:
    structural config errors, lineage names that no lineage set defines (they would select
    zero cells), unreadable marker files, and (``check_inputs``) a missing Xenium folder.
    Raises :class:`ConfigError`."""
    stages = list(stages) if stages else list(STAGES)
    problems = list(validate_config(cfg))
    lineages = set(cfg.markers.lineage) if isinstance(cfg.markers.lineage, dict) else set()

    def known(names, where):
        for n in names:
            if n not in lineages:
                problems.append(f"{where}: '{n}' is not a lineage defined in markers.lineage "
                                f"(known: {sorted(lineages)})")

    if "qc" in stages and check_inputs:
        xd, sj_in = cfg.project.xenium_dir, cfg.project.spatioloji
        if xd and sj_in:
            problems.append("project.xenium_dir and project.spatioloji are both set; only one input can be used")
        elif sj_in:
            f = input_path(cfg, config_dir)
            if not f.exists():
                problems.append(f"project.spatioloji: {f} not found")
            else:
                try:
                    import spatioloji_s  # noqa: F401
                except ImportError:
                    problems.append("project.spatioloji: the spatioloji_s package is not installed (pip install spatioloji-s)")
        elif not xd:
            problems.append("project.xenium_dir / project.spatioloji: no input is set (the qc stage needs a Xenium "
                            "output folder or a saved spatioloji object)")
        else:
            for f in ("cell_feature_matrix.h5", "cells.parquet"):
                if not (Path(xd) / f).exists():
                    problems.append(f"project.xenium_dir: {Path(xd) / f} not found")
    if "lineage" in stages:
        if not lineages:
            problems.append("markers.lineage: no lineage sets are opted in, so the lineage stage has nothing to score. "
                            "Write `lineage: default` for the default library, or list your own sets")
        elif cfg.markers.proliferation:
            known([cfg.lineage.cycling_lineage], "lineage.cycling_lineage")
    # refinement stages without opted-in sets are skipped, so their settings are only checked when they will run
    if "caf" in stages:
        try:
            if resolve_caf_sets(cfg.caf, cfg.markers.caf, base=config_dir):
                known(cfg.caf.include_lineages, "caf.include_lineages")
        except (ValueError, OSError) as e:
            problems.append(f"caf: cannot read the CAF marker sets ({e})")
    if "csc" in stages:
        try:
            if resolve_csc_markers(cfg.csc, cfg.markers, base=config_dir)[0]:
                known(cfg.csc.tumor_lineages, "csc.tumor_lineages")
        except (ValueError, OSError) as e:
            problems.append(f"csc: cannot read the CSC marker set ({e})")
    if cfg.plot.palette_file:
        try:
            pinned_colors(cfg, config_dir)
        except (ValueError, OSError) as e:
            problems.append(f"plot.palette_file: cannot read the colour table ({e})")
    if "immune" in stages:
        for name in _selected_compartments(cfg):
            known(cfg.markers.compartments[name].lineages, f"markers.compartments.{name}.lineages")
    if problems:
        raise ConfigError("configuration problems:\n  - " + "\n  - ".join(dict.fromkeys(problems)))


def resolved_config(cfg: PipelineConfig, config_dir=None) -> tuple[PipelineConfig, dict]:
    """A copy of ``cfg`` that no longer depends on anything outside itself: paths are
    absolute and marker files are inlined into ``markers``. Returns ``(config, sources)``."""
    import copy
    import os

    out = copy.deepcopy(cfg)
    sources: dict = {}
    for key in ("outdir", "checkpoint_dir", "xenium_dir"):
        v = getattr(out.project, key)
        if v:
            setattr(out.project, key, os.path.abspath(v))
    if cfg.project.spatioloji:
        out.project.spatioloji = os.path.abspath(input_path(cfg, config_dir))
    try:
        csc, core, src = resolve_csc_markers(cfg.csc, cfg.markers, base=config_dir)
        out.markers.csc, out.markers.csc_core = list(csc), list(core)
        out.csc.marker_file = out.csc.core_marker_file = None
        sources["csc"] = src
    except (ValueError, OSError):
        pass
    try:
        if cfg.plot.palette_file:
            out.plot.colors = pinned_colors(cfg, config_dir)
            out.plot.palette_file = None
            sources["palette_file"] = cfg.plot.palette_file
    except (ValueError, OSError):
        pass
    try:
        out.markers.caf = resolve_caf_sets(cfg.caf, cfg.markers.caf, base=config_dir)
        if cfg.caf.marker_files:
            sources["caf_marker_files"] = {k: str(Path(v) if Path(v).is_absolute() or config_dir is None
                                                  else Path(config_dir) / v) for k, v in cfg.caf.marker_files.items()}
        out.caf.marker_files = {}
    except (ValueError, OSError):
        pass
    return out, sources


def write_config_used(cfg: PipelineConfig, paths: RunPaths, config_dir=None) -> Path:
    """``<outdir>/config_used.yaml``: the resolved configuration of the latest invocation.
    Re-running with this file reproduces the run without the original marker files. (What each
    individual stage actually used is in ``stage_state.json``.)"""
    import yaml

    resolved, sources = resolved_config(cfg, config_dir)
    header = [f"# Resolved configuration written by annotation_st {versions().get('annotation_st', '')} on {now()}.",
              "# Paths are absolute and marker files are inlined under 'markers'."]
    for k, v in sources.items():
        header.append(f"# source {k}: {v}")
    body = yaml.safe_dump(config_to_dict(resolved), sort_keys=False, default_flow_style=None, width=110)
    out = paths.outdir / "config_used.yaml"
    out.write_text("\n".join(header) + "\n" + body)
    return out


def status_table(cfg: PipelineConfig, paths: RunPaths | None = None, config_dir=None) -> list[tuple[str, str, str]]:
    """``[(stage, status, reason)]`` for every stage, without running anything."""
    paths = paths or RunPaths(cfg, mkdir=False)
    state = StageState(paths)
    return [(s, *stage_status(s, cfg, paths, state, config_dir)) for s in STAGES]


def run_stage(name: str, cfg: PipelineConfig, paths: RunPaths | None = None, force: bool = False,
              config_dir=None, log: logging.Logger | None = None, state: StageState | None = None,
              allow_stale: bool = False) -> bool:
    """Run one stage. Returns False when it was skipped because it is up to date.
    Raises :class:`StaleUpstream` when a stage it depends on is out of date (unless
    ``allow_stale``)."""
    import scanpy as sc

    if name not in RUNNERS:
        raise KeyError(f"unknown stage '{name}' (choose from {STAGES})")
    paths = paths or RunPaths(cfg)
    state = state or StageState(paths)
    log = log or setup_logging(paths, name)
    status, reason = stage_status(name, cfg, paths, state, config_dir)
    if not force and status == "up_to_date":
        log.info("stage %s: up to date, skipping (--force to redo)%s", name, f" [{reason}]" if reason else "")
        return False
    for dep in DEPENDS[name]:
        ds, dr = stage_status(dep, cfg, paths, state, config_dir)
        if ds == "stale":
            msg = (f"stage '{name}' depends on '{dep}', which is out of date ({dr}). Re-run it first, e.g. "
                   f"`annot-st run --stages {dep} {name}`, or pass --allow-stale to use the old result.")
            if not allow_stale:
                raise StaleUpstream(msg)
            log.warning("ALLOW-STALE: %s", msg)
    sc.settings.n_jobs = cfg.project.n_jobs
    sc.set_figure_params(dpi=110, dpi_save=130, frameon=False)
    log.info("stage %s: start (%s)", name, "forced" if status == "up_to_date" else reason)
    result = RUNNERS[name](cfg, paths, log, config_dir=config_dir, state=state)
    state.record(name, cfg, outputs=result["outputs"], status=result["status"], note=result["note"],
                 config_dir=config_dir)
    behind = [s for s in STAGES if s != name and state.get(s) is not None
              and stage_status(s, cfg, paths, state, config_dir)[0] == "stale"]
    log.info("stage %s: %s%s", name, result["status"],
             f"; now out of date: {', '.join(behind)}" if behind else "")
    return True


def run(cfg: PipelineConfig, stages=None, force: bool = False, config_dir=None,
        allow_stale: bool = False) -> dict[str, bool]:
    """Run the requested stages in pipeline order. A stage runs when it has never run, when
    its inputs changed, or when an upstream stage was re-run; otherwise it is skipped.
    Returns ``{stage: ran}``."""
    unknown = set(stages or []) - set(STAGES)
    if unknown:
        raise KeyError(f"unknown stage(s) {sorted(unknown)} (choose from {STAGES})")
    wanted = [s for s in STAGES if s in (stages or STAGES)]
    paths = RunPaths(cfg, mkdir=False)
    state = StageState(paths)
    qc_will_run = "qc" in wanted and (force or stage_status("qc", cfg, paths, state, config_dir)[0] != "up_to_date")
    preflight(cfg, wanted, config_dir=config_dir, check_inputs=qc_will_run)
    paths = RunPaths(cfg)
    write_config_used(cfg, paths, config_dir)
    write_json(paths.results / "marker_sources.json", dict(sorted(cfg.markers.sources.items())))
    out = {}
    for s in wanted:
        out[s] = run_stage(s, cfg, paths, force=force, config_dir=config_dir, state=state, allow_stale=allow_stale)
    return out
