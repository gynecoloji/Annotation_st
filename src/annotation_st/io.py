"""Input readers, run-folder layout, logging and ``params.json`` bookkeeping."""
from __future__ import annotations

import datetime as _dt
import json
import logging
import sys
from pathlib import Path

import numpy as np

from .config import PipelineConfig

# ----------------------------------------------------------------------------
# Xenium input
# ----------------------------------------------------------------------------
def read_xenium_h5(h5_path: str | Path, cell_idx=None):
    """Read a 10x ``cell_feature_matrix.h5`` (features × cells, CSC by cell) into an
    AnnData of cells × *gene-expression* features. Control probes / codewords are summed
    per cell into obs columns. ``cell_idx`` (sorted ints) restricts the columns read."""
    import anndata as ad
    import h5py
    import pandas as pd
    import scipy.sparse as sp

    with h5py.File(h5_path, "r") as f:
        m = f["matrix"]
        n_feat, n_cells = m["shape"][:]
        names = m["features"]["name"][:].astype(str)
        ftype = m["features"]["feature_type"][:].astype(str)
        fid = m["features"]["id"][:].astype(str)
        barcodes = m["barcodes"][:].astype(str)
        indptr = m["indptr"][:]
        if cell_idx is None:
            data, indices, ip, bc = m["data"][:], m["indices"][:], indptr, barcodes
        else:
            cell_idx = np.asarray(sorted(cell_idx))
            lens = indptr[cell_idx + 1] - indptr[cell_idx]
            ip = np.concatenate([[0], np.cumsum(lens)])
            data = np.empty(ip[-1], dtype=m["data"].dtype)
            indices = np.empty(ip[-1], dtype=m["indices"].dtype)
            dset_d, dset_i = m["data"], m["indices"]
            step = 20000
            for k in range(0, len(cell_idx), step):
                sub = cell_idx[k:k + step]
                s0, e0 = indptr[sub[0]], indptr[sub[-1] + 1]
                blk_d, blk_i = dset_d[s0:e0], dset_i[s0:e0]
                for j, c in enumerate(sub):
                    a, b = indptr[c] - s0, indptr[c + 1] - s0
                    o = ip[k + j]
                    data[o:o + (b - a)] = blk_d[a:b]
                    indices[o:o + (b - a)] = blk_i[a:b]
            bc = barcodes[cell_idx]
    X = sp.csc_matrix((data, indices, ip), shape=(n_feat, len(bc))).T.tocsr()
    keep = ftype == "Gene Expression"

    def _sum(kind):
        return np.asarray(X[:, ftype == kind].sum(1)).ravel()
    ctrl = {
        "neg_probe_counts": _sum("Negative Control Probe"),
        "neg_codeword_counts": _sum("Negative Control Codeword"),
        "genomic_ctrl_counts": _sum("Genomic Control"),
        "unassigned_codeword_counts_h5": _sum("Unassigned Codeword"),
        "deprecated_codeword_counts_h5": _sum("Deprecated Codeword"),
    }
    adata = ad.AnnData(X=X[:, keep].astype(np.float32), obs=pd.DataFrame(ctrl, index=bc),
                       var=pd.DataFrame({"gene_id": fid[keep]}, index=names[keep]))
    adata.var_names_make_unique()
    return adata


def read_xenium(xenium_dir: str | Path, subsample: int = 0, seed: int = 0, panel: str | None = None):
    """``cell_feature_matrix.h5`` + ``cells.parquet`` → AnnData with raw counts in ``X``,
    every ``cells.parquet`` column in ``obs`` and centroids (µm) in ``obsm['spatial']``."""
    import pandas as pd

    xenium_dir = Path(xenium_dir)
    cells = pd.read_parquet(xenium_dir / "cells.parquet").set_index("cell_id")
    cells.index = cells.index.astype(str)
    n_total = len(cells)
    cell_idx = None
    if subsample and subsample < n_total:
        cell_idx = np.sort(np.random.default_rng(seed).choice(n_total, subsample, replace=False))
    adata = read_xenium_h5(xenium_dir / "cell_feature_matrix.h5", cell_idx)
    obs = cells.loc[adata.obs_names]
    for c in obs.columns:
        adata.obs[c] = obs[c].values
    adata.obsm["spatial"] = obs[["x_centroid", "y_centroid"]].to_numpy(dtype=np.float64)
    adata.uns["xenium"] = {"dir": str(xenium_dir), "n_cells_in_run": int(n_total),
                           "panel": panel or "unknown"}
    return adata


# ----------------------------------------------------------------------------
# Run folder
# ----------------------------------------------------------------------------
STAGE_NUMBER = {"qc": "01", "lineage": "02", "caf": "03", "csc": "04", "immune": "05", "final": "06"}


class RunPaths:
    """Folder layout of one run.

    ``outdir/{results,figures,logs,qc}``, ``outdir/params.json``, ``outdir/analysis.log`` and
    checkpoints ``<checkpoint_dir>/<prefix>_<NN>_<stage>.h5ad`` (checkpoint_dir defaults to
    ``outdir/checkpoints``)."""

    def __init__(self, cfg: PipelineConfig, mkdir: bool = True):
        self.outdir = Path(cfg.project.outdir)
        self.results = self.outdir / "results"
        self.figures = self.outdir / "figures"
        self.logs = self.outdir / "logs"
        self.qc = self.outdir / "qc"
        self.checkpoints = Path(cfg.project.checkpoint_dir) if cfg.project.checkpoint_dir else self.outdir / "checkpoints"
        self.params = self.outdir / "params.json"
        self.analysis_log = self.outdir / "analysis.log"
        self.prefix = cfg.project.prefix
        if mkdir:
            for d in (self.results, self.figures, self.logs, self.qc, self.checkpoints):
                d.mkdir(parents=True, exist_ok=True)

    def checkpoint(self, stage: str) -> Path:
        base = stage.split("_")[0]
        num = STAGE_NUMBER.get(base)
        if num is None:
            raise KeyError(f"no checkpoint number for stage '{stage}'")
        return self.checkpoints / f"{self.prefix}_{num}_{stage}.h5ad"

    def immune_checkpoints(self, compartments) -> dict[str, Path]:
        return {c: self.checkpoint(f"immune_{c}") for c in compartments}


def setup_logging(paths: RunPaths, stage: str) -> logging.Logger:
    log = logging.getLogger(f"annotation_st.{stage}")
    log.setLevel(logging.INFO)
    log.propagate = False
    for h in list(log.handlers):
        log.removeHandler(h)
        h.close()
    fmt = logging.Formatter("%(asctime)s %(levelname)s %(message)s")
    for h in (logging.StreamHandler(sys.stdout), logging.FileHandler(paths.logs / f"{stage}.log")):
        h.setFormatter(fmt)
        log.addHandler(h)
    return log


def update_params(paths: RunPaths, stage: str, params: dict) -> None:
    """Merge one stage's parameters into ``params.json`` (one key per stage)."""
    data = {}
    if paths.params.exists():
        try:
            data = json.loads(paths.params.read_text())
        except json.JSONDecodeError:
            data = {}
    data[stage] = params
    paths.params.write_text(json.dumps(data, indent=2, default=str) + "\n")


def append_log(paths: RunPaths, message: str) -> None:
    with open(paths.analysis_log, "a") as fh:
        fh.write(f"{_dt.date.today().isoformat()} | {message}\n")


def now() -> str:
    return _dt.datetime.now().isoformat(timespec="seconds")


def versions() -> dict[str, str]:
    from importlib.metadata import PackageNotFoundError, version

    out = {}
    for pkg in ("annotation_st", "scanpy", "squidpy", "anndata", "numpy", "pandas"):
        try:
            out[pkg] = version(pkg)
        except PackageNotFoundError:
            pass
    return out


def write_json(path: Path, obj) -> None:
    Path(path).write_text(json.dumps(obj, indent=2, default=str) + "\n")
