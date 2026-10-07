"""Step 3 — cancer stem-like cells among tumour epithelial cells.

The OCSC set is scored per cell and z-scored within tumour cells. ``CSC_like`` requires
z > ``z_thresh`` **and** ≥ ``min_core`` core markers with raw count > 0; otherwise
``Tumor_nonCSC``. Marker positivity and a CSC-vs-non-CSC Wilcoxon table are returned for
review. A tumour-only embedding (optionally on a subset of cells) is added for figures.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
import scipy.sparse as sp

from ..config import CscConfig, MarkerConfig
from ..markers import filter_to_panel, read_marker_file
from ..scoring import embed, score_sets


def resolve_csc_markers(cfg: CscConfig, markers: MarkerConfig, base=None) -> tuple[list[str], list[str], dict]:
    """Return ``(csc_genes, core_genes, source)`` according to the precedence documented on
    :class:`CscConfig`. Relative file paths resolve against ``base``."""
    from pathlib import Path

    def _path(p):
        p = Path(p)
        return p if p.is_absolute() or base is None else Path(base) / p

    source = {}
    if cfg.marker_file:
        csc = read_marker_file(_path(cfg.marker_file))
        source["csc"] = str(_path(cfg.marker_file))
    else:
        csc = list(markers.csc)
        source["csc"] = "config markers.csc"
    if not csc:
        return [], [], {"csc": "none (not opted in)", "core": "none"}
    if cfg.core_marker_file:
        core = read_marker_file(_path(cfg.core_marker_file))
        source["core"] = str(_path(cfg.core_marker_file))
    else:
        core = [g for g in markers.csc_core if g in csc]
        source["core"] = "config markers.csc_core restricted to the set"
        if not core:
            core = list(csc)
            source["core"] = "all genes of the set (no core marker overlaps it)"
    missing = [g for g in core if g not in csc]
    if missing:
        source["core_not_in_set"] = missing
    return csc, core, source


@dataclass
class CscResult:
    adata: object
    counts: pd.DataFrame
    positivity: pd.DataFrame
    de: pd.DataFrame | None
    sets_used: dict[str, list[str]]
    sets_dropped: dict[str, list[str]]
    n_csc: int


def _counts_matrix(tum):
    if "counts" in tum.layers:
        return tum.layers["counts"]
    return tum.X


def call_csc(tum, cfg: CscConfig, markers: MarkerConfig | None = None, seed: int = 0, umap: bool = True,
             embed_cells: bool = True, csc_genes: list[str] | None = None,
             core_genes: list[str] | None = None) -> CscResult | None:
    """``tum`` = tumour-epithelial subset (log-normalised ``X``, raw counts in
    ``layers['counts']``). The marker set comes from ``csc_genes`` / ``core_genes`` when given,
    else from ``markers`` via :func:`resolve_csc_markers`. Returns None when
    ``tum.n_obs < cfg.min_cells``."""
    import scanpy as sc

    if tum.n_obs < cfg.min_cells:
        return None
    if csc_genes is None or core_genes is None:
        if markers is None:
            raise ValueError("pass either markers or both csc_genes and core_genes")
        r_csc, r_core, _ = resolve_csc_markers(cfg, markers)
        csc_genes = list(csc_genes) if csc_genes is not None else r_csc
        core_genes = list(core_genes) if core_genes is not None else r_core
    sets, dropped = filter_to_panel({"OCSC": list(csc_genes), "OCSC_core": list(core_genes)}, tum.var_names)
    if not sets["OCSC"]:
        raise ValueError("no OCSC marker on this panel")
    score_sets(tum, {"OCSC": sets["OCSC"]}, prefix="csc_", seed=seed)

    counts = _counts_matrix(tum)
    core_idx = tum.var_names.get_indexer(sets["OCSC_core"])
    n_core = np.asarray((counts[:, core_idx] > 0).sum(1)).ravel() if len(core_idx) else np.zeros(tum.n_obs, int)
    tum.obs["csc_n_core_detected"] = n_core.astype(int)
    is_csc = (tum.obs["csc_OCSC_z"].to_numpy() > cfg.z_thresh) & (n_core >= cfg.min_core)
    tum.obs["csc_status"] = pd.Categorical(np.where(is_csc, cfg.positive_label, cfg.negative_label),
                                           categories=[cfg.negative_label, cfg.positive_label])
    stat = tum.obs["csc_status"].value_counts().to_frame("n_cells")
    stat["fraction"] = stat["n_cells"] / tum.n_obs

    idx = tum.var_names.get_indexer(sets["OCSC"])
    posm = counts[:, idx] > 0
    pos = pd.DataFrame(np.asarray(posm.todense()) if sp.issparse(posm) else np.asarray(posm),
                       columns=sets["OCSC"], index=tum.obs_names)
    posit = pos.groupby(tum.obs["csc_status"].values, observed=True).mean().T
    posit.columns = [f"frac_positive_{c}" for c in posit.columns]
    posit["frac_positive_all_tumor"] = pos.mean()

    de = None
    if is_csc.sum() >= 20 and (~is_csc).sum() >= 20:
        sc.tl.rank_genes_groups(tum, "csc_status", groups=[cfg.positive_label], reference=cfg.negative_label,
                                method="wilcoxon", use_raw=False)
        de = sc.get.rank_genes_groups_df(tum, group=cfg.positive_label).head(100)

    if embed_cells:
        if cfg.umap_cells and cfg.umap_cells < tum.n_obs:
            rng = np.random.default_rng(seed)
            sub_idx = np.sort(rng.choice(tum.n_obs, cfg.umap_cells, replace=False))
            sub = tum[sub_idx].copy()
            embed(sub, cfg.n_pcs, cfg.n_neighbors, cfg.resolution, seed, "tumor_leiden", umap=umap)
            tum.obs["tumor_leiden"] = "NA"
            tum.obs.loc[sub.obs_names, "tumor_leiden"] = sub.obs["tumor_leiden"].astype(str).values
            if umap:
                um = np.full((tum.n_obs, 2), np.nan, dtype=np.float32)
                um[sub_idx] = sub.obsm["X_umap"]
                tum.obsm["X_umap"] = um
            tum.obsm["X_pca"] = np.full((tum.n_obs, sub.obsm["X_pca"].shape[1]), np.nan, dtype=np.float32)
            tum.obsm["X_pca"][sub_idx] = sub.obsm["X_pca"]
        else:
            embed(tum, cfg.n_pcs, cfg.n_neighbors, cfg.resolution, seed, "tumor_leiden", umap=umap)
            tum.obs["tumor_leiden"] = tum.obs["tumor_leiden"].astype(str)
        tum.obs["tumor_leiden"] = tum.obs["tumor_leiden"].astype("category")
    return CscResult(tum, stat, posit, de, sets, dropped, int(is_csc.sum()))
