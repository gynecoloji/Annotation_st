"""Step 1 — normalise, embed, cluster and assign coarse lineages per Leiden cluster.

Rule: every lineage marker set is scored per cell (``score_genes``, z across all cells);
each Leiden cluster takes the set with the highest mean z if that mean > ``min_z``, else
``unassigned``. Cells of ``cycling_lineage`` whose proliferation z > ``cycling_z`` get
``cycling_label`` in ``obs['cell_type']`` (``obs['proliferating']`` holds the flag for all).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import pandas as pd

from ..config import LineageConfig, MarkerConfig
from ..markers import filter_to_panel
from ..scoring import cluster_assign, embed, score_sets, top_markers_table


@dataclass
class LineageResult:
    adata: object
    means: pd.DataFrame                       # leiden × set mean z + assigned, n_cells, mean_prolif_z
    mapping: dict                             # leiden cluster -> lineage
    counts: pd.DataFrame                      # cell_type counts
    sets_used: dict[str, list[str]]
    sets_dropped: dict[str, list[str]]
    prolif_used: list[str]
    top_markers: pd.DataFrame | None = None
    scored: list[str] = field(default_factory=list)


def annotate_lineage(adata, cfg: LineageConfig, markers: MarkerConfig, seed: int = 0, normalize: bool = True,
                     compute_top_markers: bool = True) -> LineageResult:
    """Operates in place on ``adata`` (counts in ``X`` when ``normalize``; log-normalised
    otherwise) and returns the tables the run folder needs."""
    import scanpy as sc

    if normalize:
        sc.pp.normalize_total(adata)          # median library size
        sc.pp.log1p(adata)
        adata.raw = adata
    embed(adata, cfg.n_pcs, cfg.n_neighbors, cfg.resolution, seed, "leiden",
          umap=not cfg.skip_umap, umap_cells=cfg.umap_cells)

    sets, dropped = filter_to_panel(markers.lineage, adata.var_names)
    prolif, _ = filter_to_panel({"Proliferating": markers.proliferation}, adata.var_names)
    scored = score_sets(adata, sets, prefix="score_", seed=seed)
    if not scored:
        raise ValueError("no lineage marker set has any gene on this panel")
    has_prolif = bool(score_sets(adata, prolif, prefix="score_", seed=seed))

    means, mapping = cluster_assign(adata, "leiden", scored, prefix="score_", min_z=cfg.min_z,
                                    unassigned=cfg.unassigned)
    means["assigned"] = pd.Series(mapping)
    means["n_cells"] = adata.obs["leiden"].value_counts()
    if has_prolif:
        means["mean_prolif_z"] = adata.obs.groupby("leiden", observed=True)["score_Proliferating_z"].mean()
        adata.obs["proliferating"] = adata.obs["score_Proliferating_z"] > cfg.cycling_z
    else:
        means["mean_prolif_z"] = float("nan")
        adata.obs["proliferating"] = False

    adata.obs["lineage"] = adata.obs["leiden"].map(mapping).astype("category")
    cell_type = adata.obs["lineage"].astype(str)
    cell_type[adata.obs["proliferating"] & (adata.obs["lineage"] == cfg.cycling_lineage)] = cfg.cycling_label
    adata.obs["cell_type"] = cell_type.astype("category")
    counts = adata.obs["cell_type"].value_counts().to_frame("n_cells")
    counts["fraction"] = counts["n_cells"] / adata.n_obs

    top = None
    if compute_top_markers:
        top = top_markers_table(adata, "leiden", mapping, use_raw=adata.raw is not None, label_col="lineage")
    return LineageResult(adata, means, mapping, counts, sets, dropped, prolif.get("Proliferating", []), top, scored)
