"""Gene-set scoring, cluster/per-cell assignment rules and the shared embedding recipe.

These primitives are what every stage is built from:

* :func:`score_sets` — ``scanpy.tl.score_genes`` per set, plus a z-score across the
  cells of the object it is called on (so "z within compartment" is achieved simply by
  calling it on the compartment subset).
* :func:`cluster_assign` — per-cluster argmax of mean z with a floor (lineage, CAF
  sub-cluster label, immune level 2).
* :func:`per_cell_assign` — per-cell argmax with floor and margin (CAF subtype).
* :func:`embed` — PCA → kNN → Leiden (igraph flavour) → optional UMAP on a subset.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def zscore(v: np.ndarray) -> np.ndarray:
    v = np.asarray(v, dtype=float)
    sd = v.std()
    return (v - v.mean()) / (sd if sd > 0 else 1.0)


def score_sets(adata, sets: dict[str, list[str]], prefix: str, seed: int = 0, ctrl_size: int = 50) -> list[str]:
    """Score every non-empty set. Writes ``obs[prefix+name]`` (raw score) and
    ``obs[prefix+name+'_z']`` (z across the cells of ``adata``). Returns the names scored."""
    import scanpy as sc

    scored = []
    for name, genes in sets.items():
        if len(genes) == 0:
            continue
        col = f"{prefix}{name}"
        kw = dict(score_name=col, random_state=seed, ctrl_size=max(ctrl_size, len(genes)))
        try:
            sc.tl.score_genes(adata, genes, **kw)
        except RuntimeError as e:
            # Only on a tiny gene universe (tests, custom mini panels): the set fills its whole
            # expression bin so scanpy finds no control genes. Retry with coarser bins, then
            # fall back to the plain mean expression of the set.
            if "No control genes" not in str(e):
                raise
            try:
                sc.tl.score_genes(adata, genes, n_bins=5, **kw)
            except RuntimeError:
                X = adata[:, genes].X
                adata.obs[col] = np.asarray(X.mean(axis=1)).ravel()
        adata.obs[col + "_z"] = zscore(adata.obs[col].to_numpy(dtype=float))
        scored.append(name)
    return scored


def cluster_assign(adata, cluster_key: str, sets: list[str], prefix: str, min_z: float = 0.0,
                   unassigned: str = "Unassigned") -> tuple[pd.DataFrame, dict]:
    """Each cluster takes the set with the highest mean z if that mean > ``min_z``.
    Returns ``(means, mapping)``: cluster × set table of mean z, and ``{cluster: label}``."""
    cols = [f"{prefix}{s}_z" for s in sets]
    means = adata.obs.groupby(cluster_key, observed=True)[cols].mean()
    means.columns = list(sets)
    mapping = {}
    for cl, row in means.iterrows():
        best = row.idxmax()
        mapping[cl] = best if row[best] > min_z else unassigned
    return means, mapping


def per_cell_assign(Z: np.ndarray, names: list[str], min_z: float, margin: float,
                    unassigned: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Per-cell argmax over the columns of ``Z`` (cells × sets). A call needs best z > ``min_z``
    and best − second best ≥ ``margin``. Returns ``(call, best_z, margin_arr)``."""
    Z = np.asarray(Z, dtype=float)
    if Z.ndim != 2 or Z.shape[1] != len(names):
        raise ValueError("Z must be cells x len(names)")
    order = np.argsort(-Z, axis=1)
    rows = np.arange(len(Z))
    best = Z[rows, order[:, 0]]
    second = Z[rows, order[:, 1]] if Z.shape[1] > 1 else np.full(len(Z), -np.inf)
    marg = best - second
    call = np.array(names, dtype=object)[order[:, 0]]
    call[(best <= min_z) | (marg < margin)] = unassigned
    return call, best, marg


def embed(adata, n_pcs: int, n_neighbors: int, resolution: float, seed: int, key: str,
          umap: bool = True, umap_cells: int = 0, n_iterations: int = 2) -> None:
    """PCA(n_pcs) → kNN(n_neighbors) → Leiden(resolution, igraph) into ``obs[key]``;
    optional UMAP, computed on a random subset of ``umap_cells`` cells when that is smaller
    than the object (NaN elsewhere; subset index stored in ``uns[key+'_umap_subset_idx']``)."""
    import scanpy as sc

    n_pcs = int(min(n_pcs, adata.n_vars - 1, max(adata.n_obs - 1, 1)))
    sc.pp.pca(adata, n_comps=n_pcs, random_state=seed)
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, random_state=seed)
    sc.tl.leiden(adata, resolution=resolution, flavor="igraph", n_iterations=n_iterations, directed=False,
                 random_state=seed, key_added=key)
    if not umap:
        return
    if umap_cells and umap_cells < adata.n_obs:
        rng = np.random.default_rng(seed)
        idx = np.sort(rng.choice(adata.n_obs, umap_cells, replace=False))
        sub = adata[idx].copy()
        sc.tl.umap(sub, random_state=seed)
        um = np.full((adata.n_obs, 2), np.nan, dtype=np.float32)
        um[idx] = sub.obsm["X_umap"]
        adata.obsm["X_umap"] = um
        adata.uns[f"{key}_umap_subset_idx"] = idx
    else:
        sc.tl.umap(adata, random_state=seed)


def top_markers_table(adata, cluster_key: str, labels: dict, n_genes: int = 15, use_raw: bool | None = None,
                      label_col: str = "lineage") -> pd.DataFrame:
    """Wilcoxon top genes per cluster (for manual review of an assignment)."""
    import scanpy as sc

    sc.tl.rank_genes_groups(adata, cluster_key, method="wilcoxon", n_genes=n_genes, use_raw=use_raw)
    rows = []
    for cl in adata.obs[cluster_key].cat.categories:
        df = sc.get.rank_genes_groups_df(adata, group=cl).head(n_genes)
        rows.append({"cluster": cl, label_col: labels.get(cl, ""), "top_genes": ",".join(map(str, df["names"]))})
    return pd.DataFrame(rows)
