"""Step 6 — spatial graph, neighbourhood enrichment and neighbourhood composition.

* Delaunay graph on all centroids (or a fixed-radius graph when ``radius > 0``); edges
  longer than the ``edge_cut_pct`` percentile are removed (tissue gaps).
* ``squidpy.gr.nhood_enrichment`` on a label key with label permutations
  (``n_jobs=1``: parallel runs crash on read-only memmaps in squidpy 1.8.1); non-finite z
  set to 0.
* Neighbourhood composition: with adjacency A and one-hot labels G, ``A·G`` gives per-cell
  neighbour counts per label; fractions are averaged over each group of interest and a
  log2 ratio (pseudo-count 1e-4) is added when exactly two groups are given.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from ..config import SpatialConfig


def build_spatial_graph(adata, cfg: SpatialConfig) -> float:
    """Build ``obsp['spatial_connectivities'/'spatial_distances']``. Returns the edge cut (µm)."""
    import squidpy as sq

    if cfg.radius > 0:
        sq.gr.spatial_neighbors(adata, coord_type="generic", radius=cfg.radius, delaunay=False)
    else:
        sq.gr.spatial_neighbors(adata, coord_type="generic", delaunay=True)
    D = adata.obsp["spatial_distances"].tocsr()
    C = adata.obsp["spatial_connectivities"].tocsr()
    cut = float(np.percentile(D.data, cfg.edge_cut_pct)) if D.nnz else float("nan")
    if np.isfinite(cut):
        bad = (D > cut).astype(np.float32)
        C = (C - C.multiply(bad)).tocsr(); C.eliminate_zeros()
        D = (D - D.multiply(bad)).tocsr(); D.eliminate_zeros()
    adata.obsp["spatial_connectivities"] = C
    adata.obsp["spatial_distances"] = D
    return cut


def n_edges(adata) -> int:
    return int(adata.obsp["spatial_connectivities"].nnz // 2)


def nhood_enrichment(adata, key: str, cfg: SpatialConfig, seed: int = 0) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Permutation neighbourhood enrichment on ``obs[key]``. Returns ``(zscore, count)`` tables."""
    import squidpy as sq

    if not isinstance(adata.obs[key].dtype, pd.CategoricalDtype):
        adata.obs[key] = adata.obs[key].astype("category")
    adata.obs[key] = adata.obs[key].cat.remove_unused_categories()
    sq.gr.nhood_enrichment(adata, cluster_key=key, n_perms=cfg.n_perms, seed=seed, n_jobs=1,
                           show_progress_bar=False)
    cats = list(adata.obs[key].cat.categories)
    res = adata.uns[f"{key}_nhood_enrichment"]
    zs = np.asarray(res["zscore"], dtype=float)
    if not np.isfinite(zs).all():
        zs = np.nan_to_num(zs, nan=0.0, posinf=0.0, neginf=0.0)
        res["zscore"] = zs
    z = pd.DataFrame(zs, index=cats, columns=cats)
    count = pd.DataFrame(np.asarray(res["count"]), index=cats, columns=cats)
    return z, count


def neighbourhood_composition(adata, key: str, groups) -> pd.DataFrame:
    """Mean fraction of neighbours of each label around cells of each group in ``groups``
    (rows: labels; columns: groups [+ ``log2_ratio_<g0>_vs_<g1>`` for two groups])."""
    lab_s = adata.obs[key]
    if not isinstance(lab_s.dtype, pd.CategoricalDtype):
        lab_s = lab_s.astype("category")
    cats = list(lab_s.cat.categories)
    A = adata.obsp["spatial_connectivities"].tocsr()
    onehot = np.zeros((adata.n_obs, len(cats)), dtype=np.float32)
    onehot[np.arange(adata.n_obs), lab_s.cat.codes.to_numpy()] = 1
    neigh = A @ onehot
    comp = pd.DataFrame(np.asarray(neigh), index=adata.obs_names, columns=cats)
    rows = {}
    for grp in groups:
        m = (lab_s == grp).to_numpy()
        if m.sum():
            tot = comp[m].sum(1).replace(0, np.nan)
            rows[grp] = comp[m].div(tot, axis=0).mean()
    out = pd.DataFrame(rows)
    if len(groups) == 2 and out.shape[1] == 2:
        g0, g1 = groups
        out[f"log2_ratio_{g0}_vs_{g1}"] = np.log2((out[g0] + 1e-4) / (out[g1] + 1e-4))
    return out
