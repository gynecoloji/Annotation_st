"""Step 0 — per-cell QC metrics and filtering.

Filters (a cell failing any is dropped): transcripts ≥ ``min_transcripts``, genes detected
≥ ``min_genes``, control fraction ≤ ``max_control_frac``, cell area inside the
``area_pct`` percentile window, ≥ 1 nucleus when ``require_nucleus``. Genes detected in
fewer than ``min_cells_per_gene`` cells are removed. Raw counts are kept in
``layers['counts']``.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from ..config import QCConfig

QC_METRICS = ["gene_counts", "n_genes_by_counts", "control_frac", "cell_area", "nucleus_area",
              "nucleus_ratio", "counts_per_area", "nucleus_count"]
PCTS = [.01, .05, .25, .5, .75, .95, .99]


@dataclass
class QCResult:
    adata: object
    flag_summary: pd.DataFrame
    area_bounds: tuple[float, float]
    summary_pre: pd.DataFrame
    summary_post: pd.DataFrame


def compute_qc_metrics(adata) -> None:
    """Adds gene_counts, n_genes_by_counts, control_counts, control_frac, nucleus_ratio,
    counts_per_area to ``obs`` (``X`` must hold raw counts)."""
    import scanpy as sc

    sc.pp.calculate_qc_metrics(adata, percent_top=None, inplace=True)
    adata.obs["gene_counts"] = adata.obs["total_counts"].astype(np.int64)
    ctrl = np.zeros(adata.n_obs)
    for c in ("neg_probe_counts", "neg_codeword_counts", "genomic_ctrl_counts"):
        if c in adata.obs:
            ctrl = ctrl + adata.obs[c].to_numpy(dtype=float)
    adata.obs["control_counts"] = ctrl.astype(np.int64)
    denom = (adata.obs["gene_counts"] + adata.obs["control_counts"]).clip(lower=1)
    adata.obs["control_frac"] = adata.obs["control_counts"] / denom
    if "cell_area" in adata.obs:
        area = adata.obs["cell_area"].clip(lower=1e-6)
        if "nucleus_area" in adata.obs:
            adata.obs["nucleus_ratio"] = adata.obs["nucleus_area"] / area
        adata.obs["counts_per_area"] = adata.obs["gene_counts"] / area


def _summary(adata) -> pd.DataFrame:
    cols = [c for c in QC_METRICS if c in adata.obs]
    return adata.obs[cols].describe(percentiles=PCTS).T


def qc_filter(adata, cfg: QCConfig) -> QCResult:
    """Apply the filters of ``cfg``. Returns the filtered copy plus the flag summary
    (rows: one per filter, ``ANY`` and ``KEPT``; columns ``n_cells_flagged``, ``fraction``)."""
    import scanpy as sc

    if "gene_counts" not in adata.obs:
        compute_qc_metrics(adata)
    summary_pre = _summary(adata)

    flags = pd.DataFrame(index=adata.obs_names)
    flags["low_transcripts"] = adata.obs["gene_counts"] < cfg.min_transcripts
    flags["low_genes"] = adata.obs["n_genes_by_counts"] < cfg.min_genes
    flags["high_control_frac"] = adata.obs["control_frac"] > cfg.max_control_frac
    if "cell_area" in adata.obs:
        lo, hi = np.percentile(adata.obs["cell_area"], cfg.area_pct)
        flags["area_outlier"] = (adata.obs["cell_area"] < lo) | (adata.obs["cell_area"] > hi)
    else:
        lo = hi = float("nan")
        flags["area_outlier"] = False
    if cfg.require_nucleus and "nucleus_count" in adata.obs:
        flags["no_nucleus"] = adata.obs["nucleus_count"] < 1
    else:
        flags["no_nucleus"] = False
    fail_any = flags.any(axis=1)
    fsum = flags.sum().to_frame("n_cells_flagged").astype(int)
    fsum["fraction"] = fsum["n_cells_flagged"] / adata.n_obs
    fsum.loc["ANY"] = [int(fail_any.sum()), float(fail_any.mean())]
    fsum.loc["KEPT"] = [int((~fail_any).sum()), float((~fail_any).mean())]
    fsum["n_cells_flagged"] = fsum["n_cells_flagged"].astype(int)

    if cfg.apply_filters:
        out = adata[~fail_any.values].copy()
        sc.pp.filter_genes(out, min_cells=cfg.min_cells_per_gene)
    else:                                   # already filtered upstream: keep everything, the summary still reports
        out = adata.copy()
    out.layers["counts"] = out.X.copy()
    return QCResult(out, fsum, (float(lo), float(hi)), summary_pre, _summary(out))
