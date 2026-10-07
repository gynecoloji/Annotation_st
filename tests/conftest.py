"""Synthetic Xenium-like data with planted populations.

Every population expresses a handful of marker genes at a high Poisson rate on top of a
uniform low background, sits in its own spatial blob, and carries per-cell QC columns, so
each stage's decision rule can be checked on data whose truth is known.
"""
from __future__ import annotations

import os

os.environ.setdefault("MPLBACKEND", "Agg")

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from annotation_st import markers as M

# Synthetic CAF sets (subset of the user-supplied lists that are on the 5K panel)
CAF_SETS = {
    "iCAF": ["CXCL12", "IL6", "CXCL2", "FGF7", "IGF1", "HAS2", "PI16", "MFAP5", "DPT", "THBS1"],
    "myCAF": ["MMP11", "HOPX", "CTHRC1", "INHBA", "COL10A1", "POSTN", "GRP", "CST1"],
    "apCAF": ["PDPN", "XBP1", "STAT1", "IRF5", "NFE2L2", "BCAM", "F11R"],
}

T_BASE = ["CD3E", "CD3G", "CD2", "CD247"]
NAIVE = ["CCR7", "SELL", "TCF7", "LEF1"]

# (truth label, truth_coarse (=lineage), n cells, {gene: rate})
POPULATIONS = [
    ("Epithelial_tumor", "Epithelial_tumor", 500, {g: 4 for g in M.LINEAGE_MARKERS["Epithelial_tumor"]}),
    ("CSC", "Epithelial_tumor", 100, {**{g: 4 for g in M.LINEAGE_MARKERS["Epithelial_tumor"]},
                                       **{g: 6 for g in ["CD44", "PROM1", "POU5F1", "CXCR4", "SOX2"]}}),
    ("Tumor_cycling", "Epithelial_tumor", 50, {**{g: 4 for g in M.LINEAGE_MARKERS["Epithelial_tumor"]},
                                                **{g: 5 for g in M.PROLIF_MARKERS}}),
    ("iCAF", "Fibroblast", 250, {**{g: 4 for g in M.LINEAGE_MARKERS["Fibroblast"]}, **{g: 6 for g in CAF_SETS["iCAF"]}}),
    ("myCAF", "Fibroblast", 250, {**{g: 4 for g in M.LINEAGE_MARKERS["Fibroblast"]}, **{g: 6 for g in CAF_SETS["myCAF"]}}),
    ("CD8_T", "T_NK", 130, {**{g: 5 for g in T_BASE + ["CD8A", "CD8B"]}, **{g: 4 for g in NAIVE}}),
    ("CD8_cyto", "T_NK", 120, {**{g: 5 for g in T_BASE + ["CD8A", "CD8B"]},
                               **{g: 5 for g in ["GZMB", "GZMA", "GZMH", "GZMK", "PRF1", "KLRG1"]}}),
    ("CD4_T", "T_NK", 150, {**{g: 5 for g in T_BASE + ["CD4", "CD40LG"]}, **{g: 4 for g in NAIVE}}),
    ("Treg", "T_NK", 100, {**{g: 5 for g in T_BASE + ["CD4"]},
                           **{g: 5 for g in ["FOXP3", "IL2RA", "CTLA4", "IKZF2", "TNFRSF18", "CCR8"]}}),
    ("NK", "T_NK", 150, {g: 5 for g in ["CD2", "CD247", "KLRD1", "KLRF1", "NCR1", "NCAM1", "FCGR3A", "KLRC1",
                                         "KLRK1", "XCL2", "GZMB", "PRF1"]}),
    # B cells are too rare to form their own lineage cluster in the real data and end up in T_NK
    ("B_cell", "T_NK", 120, {g: 5 for g in M.B_LEVEL1["B_cell"] + ["CD27", "TNFRSF13B"]}),
    ("Plasma", "Plasma", 200, {g: 5 for g in M.B_LEVEL1["Plasma"]}),
    ("Macrophage", "Macrophage_Mono", 250, {g: 4 for g in M.LINEAGE_MARKERS["Macrophage_Mono"]}),
    ("TREM2_TAM", "Macrophage_Mono", 150, {**{g: 4 for g in M.LINEAGE_MARKERS["Macrophage_Mono"]},
                                            **{g: 6 for g in ["TREM2", "LIPA", "LGMN", "CTSL", "GPNMB"]}}),
    ("Endothelial", "Endothelial", 200, {g: 4 for g in M.LINEAGE_MARKERS["Endothelial"]}),
    ("Ovarian_stroma", "Ovarian_stroma", 150, {g: 4 for g in M.LINEAGE_MARKERS["Ovarian_stroma"]}),
]


def _all_marker_genes() -> list[str]:
    genes: list[str] = []

    def add(gs):
        for g in gs:
            if g not in genes:
                genes.append(g)
    for v in M.LINEAGE_MARKERS.values():
        add(v)
    add(M.PROLIF_MARKERS); add(M.CSC_MARKERS)
    for v in M.SPILLOVER_SETS.values():
        add(v)
    for comp in M.IMMUNE_COMPARTMENTS.values():
        for v in comp["level1"].values():
            add(v)
        for states in comp["level2"].values():
            for v in states.values():
                add(v)
    for v in CAF_SETS.values():
        add(v)
    for _, _, _, rates in POPULATIONS:
        add(rates)
    return genes


def synthetic_counts(seed: int = 0, background: float = 0.3, n_filler: int = 600) -> ad.AnnData:
    """Marker genes + ``n_filler`` neutral genes whose rates span the expression range (like a
    real panel), so score_genes always finds neutral control genes in every expression bin."""
    rng = np.random.default_rng(seed)
    markers = _all_marker_genes()
    genes = markers + [f"FILL{i:03d}" for i in range(n_filler)]
    gidx = {g: i for i, g in enumerate(genes)}
    n_cells = sum(n for _, _, n, _ in POPULATIONS)
    lam = np.full((n_cells, len(genes)), background, dtype=np.float32)
    filler_rates = np.clip(rng.lognormal(np.log(0.3), 1.3, n_filler), 0.02, 8.0).astype(np.float32)
    lam[:, len(markers):] = filler_rates
    truth, coarse, xy = [], [], np.zeros((n_cells, 2))
    grid = int(np.ceil(np.sqrt(len(POPULATIONS))))
    start = 0
    for k, (name, lin, n, rates) in enumerate(POPULATIONS):
        sl = slice(start, start + n)
        for g, r in rates.items():
            lam[sl, gidx[g]] = r
        truth += [name] * n
        coarse += [lin] * n
        cx, cy = (k % grid) * 600 + 300, (k // grid) * 600 + 300
        xy[sl] = rng.normal([cx, cy], 70, size=(n, 2))
        start += n
    size = rng.lognormal(0, 0.3, n_cells).astype(np.float32)[:, None]
    X = rng.poisson(lam * size).astype(np.float32)
    obs = pd.DataFrame({"truth": truth, "truth_coarse": coarse}, index=[f"cell_{i:05d}" for i in range(n_cells)])
    obs["cell_area"] = rng.lognormal(np.log(80), 0.35, n_cells)
    obs["nucleus_area"] = obs["cell_area"] * rng.uniform(0.2, 0.5, n_cells)
    obs["nucleus_count"] = 1
    obs["x_centroid"], obs["y_centroid"] = xy[:, 0], xy[:, 1]
    for c in ["neg_probe_counts", "neg_codeword_counts", "genomic_ctrl_counts"]:
        obs[c] = rng.poisson(0.02, n_cells)
    obs["unassigned_codeword_counts_h5"] = rng.poisson(0.2, n_cells)
    obs["deprecated_codeword_counts_h5"] = 0
    a = ad.AnnData(X=sp.csr_matrix(X), obs=obs, var=pd.DataFrame({"gene_id": genes}, index=genes))
    a.obsm["spatial"] = xy
    for c in ["truth", "truth_coarse"]:
        a.obs[c] = a.obs[c].astype("category")
    return a


def plant_qc_failures(a: ad.AnnData, seed: int = 1) -> ad.AnnData:
    """20 cells with < 30 transcripts, 10 cells without a nucleus, 5 cells with high control fraction."""
    rng = np.random.default_rng(seed)
    a = a.copy()
    X = a.X.tolil()
    low = rng.choice(a.n_obs, 20, replace=False)
    for i in low:
        X.rows[i], X.data[i] = [], []
        for j in rng.choice(a.n_vars, 5, replace=False):
            X[i, j] = 1
    a.X = X.tocsr()
    rest = np.setdiff1d(np.arange(a.n_obs), low)
    a.obs.iloc[rng.choice(rest, 10, replace=False), a.obs.columns.get_loc("nucleus_count")] = 0
    a.obs.iloc[rng.choice(rest, 5, replace=False), a.obs.columns.get_loc("neg_probe_counts")] = 200
    a.uns["planted"] = {"low": 20, "no_nucleus": 10, "control": 5}
    return a


def normalise(a: ad.AnnData) -> ad.AnnData:
    import scanpy as sc

    a = a.copy()
    a.layers["counts"] = a.X.copy()
    sc.pp.normalize_total(a)
    sc.pp.log1p(a)
    a.raw = a
    return a


def write_synthetic_xenium(xdir: Path, a: ad.AnnData) -> None:
    """Write ``cell_feature_matrix.h5`` (10x layout, features × cells CSC) and ``cells.parquet``."""
    import h5py

    xdir = Path(xdir)
    xdir.mkdir(parents=True, exist_ok=True)
    genes = list(a.var_names)
    n_ctrl = 3
    feat_names = genes + [f"NegControlProbe_{i}" for i in range(n_ctrl)]
    feat_types = ["Gene Expression"] * len(genes) + ["Negative Control Probe"] * n_ctrl
    rng = np.random.default_rng(0)
    ctrl = sp.csr_matrix(rng.poisson(0.05, (a.n_obs, n_ctrl)).astype(np.float32))
    full = sp.hstack([sp.csr_matrix(a.X), ctrl]).tocsr()          # cells × features
    m = full.T.tocsc()                                            # features × cells, CSC by cell
    with h5py.File(xdir / "cell_feature_matrix.h5", "w") as f:
        g = f.create_group("matrix")
        g.create_dataset("shape", data=np.array(m.shape, dtype=np.int32))
        g.create_dataset("data", data=m.data.astype(np.int32))
        g.create_dataset("indices", data=m.indices.astype(np.int64))
        g.create_dataset("indptr", data=m.indptr.astype(np.int64))
        g.create_dataset("barcodes", data=np.array(a.obs_names, dtype="S"))
        fg = g.create_group("features")
        fg.create_dataset("name", data=np.array(feat_names, dtype="S"))
        fg.create_dataset("id", data=np.array(feat_names, dtype="S"))
        fg.create_dataset("feature_type", data=np.array(feat_types, dtype="S"))
    cells = a.obs[["x_centroid", "y_centroid", "cell_area", "nucleus_area", "nucleus_count"]].copy()
    cells["transcript_counts"] = np.asarray(a.X.sum(1)).ravel().astype(int)
    cells["segmentation_method"] = "Segmented by nucleus expansion"
    cells.index.name = "cell_id"
    cells.reset_index().to_parquet(xdir / "cells.parquet", index=False)


def ov_config():
    """Default parameters plus the marker sets of the Xenium OV run, opted in explicitly. The
    synthetic populations above express exactly these sets."""
    from annotation_st.config import config_from_dict

    return config_from_dict({"markers": M.ov_run_sets()})


# ----------------------------------------------------------------------------
# fixtures
# ----------------------------------------------------------------------------
@pytest.fixture(scope="session")
def synth_counts():
    return synthetic_counts()


@pytest.fixture(scope="session")
def synth_raw(synth_counts):
    return plant_qc_failures(synth_counts)


@pytest.fixture(scope="session")
def _synth_norm(synth_counts):
    return normalise(synth_counts)


@pytest.fixture
def synth(_synth_norm):
    return _synth_norm.copy()


@pytest.fixture
def synth_lineage(_synth_norm):
    a = _synth_norm.copy()
    a.obs["lineage"] = a.obs["truth_coarse"].astype(str).astype("category")
    a.obs["cell_type"] = a.obs["lineage"].astype(str)
    a.obs.loc[a.obs["truth"] == "Tumor_cycling", "cell_type"] = "Epithelial_tumor_cycling"
    a.obs["cell_type"] = a.obs["cell_type"].astype("category")
    a.obs["leiden"] = a.obs["truth"].cat.codes.astype(str).astype("category")
    return a


@pytest.fixture
def synth_fb(synth_lineage):
    return synth_lineage[synth_lineage.obs["lineage"] == "Fibroblast"].copy()


@pytest.fixture
def synth_tum(synth_lineage):
    return synth_lineage[synth_lineage.obs["lineage"] == "Epithelial_tumor"].copy()


@pytest.fixture
def caf_sets():
    return {k: list(v) for k, v in CAF_SETS.items()}
