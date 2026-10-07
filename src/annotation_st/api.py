"""Annotate an AnnData in memory, without checkpoints or a run folder.

``annotate_adata`` chains the same stage functions as the command-line pipeline (QC,
lineage, CAF, CSC, immune, merge and, on request, the spatial statistics) and gives the
same labels; it is the entry point for notebooks and for other packages' data objects
(see :mod:`annotation_st.spatioloji_bridge`).
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace

import pandas as pd

from .config import PipelineConfig


@dataclass
class AnnotationResult:
    adata: object                                   # all (QC-passed) cells with every label column
    palette: dict[str, str]                         # label -> hex colour (rough and fine types)
    stages: dict[str, str]                          # stage -> "completed" / "skipped: <reason>"
    tables: dict[str, pd.DataFrame] = field(default_factory=dict)
    lineage: object | None = None                   # LineageResult
    caf: object | None = None                       # CafResult
    csc: object | None = None                       # CscResult
    immune: dict = field(default_factory=dict)      # compartment -> ImmuneResult
    merge: object | None = None                     # MergeResult


def annotate_adata(adata, cfg: PipelineConfig, qc: bool = False, spatial: bool = False, normalize: bool = True,
                   umap: bool = False, copy: bool = True, config_dir=None) -> AnnotationResult:
    """Run the annotation on ``adata`` and return the labelled object.

    Args:
        adata: cells × genes with raw counts in ``X`` (or log-normalised values when
            ``normalize=False``, with raw counts in ``layers['counts']``) and centroids in
            ``obsm['spatial']`` (needed only for ``spatial=True``).
        cfg: the configuration; only marker sets it opts in to are used.
        qc: compute the QC metrics first (``obs`` gets gene_counts, control_frac, ...). No cell or
            gene is removed unless ``cfg.qc.apply_filters`` is True.
        spatial: also build the spatial graph and compute neighbourhood enrichment and the
            neighbourhood composition of ``cfg.spatial.composition_groups``.
        normalize: normalise and log-transform ``X`` (set False if it already is).
        umap: compute the lineage UMAP (only useful for plotting).
        copy: work on a copy and leave ``adata`` untouched.
        config_dir: folder against which marker / palette files of ``cfg`` resolve.

    Returns:
        AnnotationResult. Its ``adata.obs`` has ``leiden``, ``lineage``, ``cell_type``,
        ``cell_type_fine`` and the CAF / CSC / immune columns; stages without opted-in
        marker sets are listed as skipped and leave the rough label in place.
    """
    from .pipeline import colour_columns, label_palette, preflight
    from .plotting import color_table
    from .stages.caf import annotate_caf, resolve_caf_sets
    from .stages.csc import call_csc, resolve_csc_markers
    from .stages.immune import run_immune
    from .stages.lineage import annotate_lineage
    from .stages.merge import merge_labels
    from .stages.qc import compute_qc_metrics, qc_filter
    from .state import selected_compartments

    import scanpy as sc

    preflight(cfg, ["lineage", "caf", "csc", "immune", "merge"], config_dir=config_dir, check_inputs=False)
    sc.settings.n_jobs = cfg.project.n_jobs               # as the command-line pipeline does
    seed = cfg.project.seed
    stages: dict[str, str] = {}
    tables: dict[str, pd.DataFrame] = {}
    a = adata.copy() if copy else adata
    if qc:
        compute_qc_metrics(a)
        q = qc_filter(a, cfg.qc)
        a = q.adata
        tables["qc_filter_summary"] = q.flag_summary
        stages["qc"] = "completed"
    if "counts" not in a.layers:
        if not normalize:
            raise ValueError("normalize=False needs the raw counts in adata.layers['counts'] (the CSC rule counts detected markers)")
        a.layers["counts"] = a.X.copy()

    lin = annotate_lineage(a, replace(cfg.lineage, skip_umap=cfg.lineage.skip_umap or not umap), cfg.markers,
                           seed=seed, normalize=normalize)
    stages["lineage"] = "completed"
    tables["leiden_lineage_scores"], tables["cell_type_counts"] = lin.means, lin.counts

    caf = csc = None
    caf_sets = resolve_caf_sets(cfg.caf, cfg.markers.caf, base=config_dir)
    if not caf_sets:
        stages["caf"] = "skipped: no CAF marker sets opted in"
    else:
        fb = a[a.obs["lineage"].isin(list(cfg.caf.include_lineages)).values].copy()
        caf = annotate_caf(fb, cfg.caf, caf_sets, seed=seed, umap=False)
        stages["caf"] = "completed" if caf is not None else f"skipped: fewer than {cfg.caf.min_cells} cells"
        if caf is not None:
            tables["caf_subtype_counts"] = caf.counts

    csc_genes, core_genes, _ = resolve_csc_markers(cfg.csc, cfg.markers, base=config_dir)
    if not csc_genes:
        stages["csc"] = "skipped: no CSC marker set opted in"
    else:
        tum = a[a.obs["lineage"].isin(list(cfg.csc.tumor_lineages)).values].copy()
        csc = call_csc(tum, cfg.csc, seed=seed, umap=False, embed_cells=False, csc_genes=csc_genes, core_genes=core_genes)
        stages["csc"] = "completed" if csc is not None else f"skipped: fewer than {cfg.csc.min_cells} cells"
        if csc is not None:
            tables["csc_counts"] = csc.counts

    immune = {}
    if not selected_compartments(cfg):
        stages["immune"] = "skipped: no immune compartment opted in"
    else:
        immune = run_immune(a, cfg.immune, cfg.markers, seed=seed, umap=False)
        stages["immune"] = "completed" if immune else f"skipped: no compartment had {cfg.immune.min_cells} cells"
        for name, r in immune.items():
            tables[f"immune_{name}_subtype_counts"] = r.counts

    merged = merge_labels(a, caf.adata.obs if caf is not None else None, csc.adata.obs if csc is not None else None,
                          {name: r.adata.obs for name, r in immune.items()}, cfg.markers.compartments)
    stages["merge"] = "completed"
    tables["final_cell_type_counts"], tables["final_cell_type_coarse_counts"] = merged.fine_counts, merged.coarse_counts
    pal = label_palette(cfg, list(a.obs["cell_type_fine"].cat.categories) + list(a.obs["cell_type"].cat.categories), config_dir)
    colour_columns(cfg, a, ["cell_type", "cell_type_fine", "caf_subtype_cell", "csc_status", "immune_subtype"], config_dir)
    tables["cell_type_colors"] = color_table(pal, coarse=a.obs["cell_type"], fine=a.obs["cell_type_fine"],
                                             groups=cfg.plot.legend_groups, qc_group=cfg.plot.qc_group,
                                             generic=cfg.plot.generic_labels)

    if spatial:
        from .stages.spatial import build_spatial_graph, neighbourhood_composition, nhood_enrichment

        build_spatial_graph(a, cfg.spatial)
        z, count = nhood_enrichment(a, cfg.spatial.fine_key, cfg.spatial, seed=seed)
        tables["nhood_enrichment_zscore"], tables["nhood_enrichment_count"] = z, count
        tables["csc_neighbourhood_composition"] = neighbourhood_composition(
            a, cfg.spatial.fine_key, list(cfg.spatial.composition_groups))
        stages["spatial"] = "completed"
    return AnnotationResult(a, pal, stages, tables, lin, caf, csc, immune, merged)
