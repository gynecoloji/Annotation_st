"""Compatibility with the ``spatioloji_s`` package (https://github.com/gynecoloji/spatioloji_s).

A ``spatioloji`` object keeps expression, cell metadata, coordinates and polygons under one
cell index. This module lets the annotation run on such an object and puts the result back
where spatioloji's own spatial, communication and plotting functions look for it:

* :func:`annotate` runs the annotation in memory and writes the labels into ``sp.cell_meta``
  as categorical columns (``cell_type``, ``cell_type_fine``, ...).
* :func:`palette` gives the colours of a label column as an ordered ``{label: colour}``
  mapping that works with every way spatioloji plots take colours (``colors=``,
  ``color_map=`` by label; ``palette=`` by position).
* :func:`to_anndata` / :func:`read_spatioloji` turn an object (or a saved one) into the
  AnnData the pipeline expects; :func:`add_annotations` writes the labels of a finished
  command-line run into an object.

``spatioloji_s`` is imported only when a function here is called.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from . import plotting as P
from .config import PipelineConfig

# Xenium control features, by 10x feature type and by name prefix -> QC column they are summed into
_CONTROL_TYPES = {"Negative Control Probe": "neg_probe_counts", "Negative Control Codeword": "neg_codeword_counts",
                  "Genomic Control": "genomic_ctrl_counts", "Unassigned Codeword": "unassigned_codeword_counts_h5",
                  "Deprecated Codeword": "deprecated_codeword_counts_h5"}
_CONTROL_PREFIXES = {"NegControlProbe": "neg_probe_counts", "NegControlCodeword": "neg_codeword_counts",
                     "Intergenic": "genomic_ctrl_counts", "UnassignedCodeword": "unassigned_codeword_counts_h5",
                     "DeprecatedCodeword": "deprecated_codeword_counts_h5"}
LABEL_COLUMNS = ["lineage", "cell_type", "cell_type_fine", "caf_subtype_cell", "caf_subtype_cluster", "csc_status",
                 "immune_compartment", "immune_subtype"]
NUMERIC_COLUMNS = ["csc_OCSC_z", "csc_n_core_detected"]
OTHER_COLUMNS = ["leiden", "proliferating"]
NA_LABEL = "NA"                          # cells without a call (filtered out, or not in that compartment)
NA_COLOUR = "#e6e6e6"


def _spatioloji_class():
    try:
        from spatioloji_s import spatioloji
    except ImportError as err:
        raise ImportError("spatioloji_s is needed for this function. Install it with: pip install spatioloji-s") from err
    return spatioloji


def to_anndata(sp, log_layer: str | None = None):
    """AnnData for the annotation from a ``spatioloji`` object.

    * ``X``: raw counts of the gene-expression features (CSR, float32). Control features
      (negative control probes / codewords, genomic controls, unassigned and deprecated
      codewords), recognised by ``gene_meta['feature_type']`` or by their name prefix, are
      not genes: they are summed per cell into the QC columns ``neg_probe_counts``,
      ``neg_codeword_counts``, ``genomic_ctrl_counts``, ... instead.
    * ``obs``: a copy of ``sp.cell_meta`` indexed by the cell ids.
    * ``obsm['spatial']``: global (x, y) only. (``sp.to_anndata()`` stores four columns,
      which a spatial neighbour graph would treat as four dimensions.)
    * ``log_layer``: name of a log-normalised layer of ``sp`` to use as ``X``; the raw
      counts then go to ``layers['counts']``.
    """
    import anndata as ad
    import scipy.sparse as sps

    expr = sp.expression
    X = (expr.get_sparse() if expr.is_sparse else sps.csr_matrix(expr.get_dense())).tocsr().astype(np.float32)
    genes = pd.Index(sp.gene_index.astype(str))
    kind = pd.Series("", index=genes)
    gm = sp.gene_meta
    if "feature_type" in gm.columns:
        ft = gm["feature_type"].astype(str).to_numpy()
        for t, col in _CONTROL_TYPES.items():
            kind[ft == t] = col
    for prefix, col in _CONTROL_PREFIXES.items():
        kind[(kind == "").to_numpy() & genes.str.startswith(prefix)] = col
    is_gene = (kind == "").to_numpy()

    obs = sp.cell_meta.copy()
    obs.index = pd.Index(sp.cell_index.astype(str), name=None)
    for col in dict.fromkeys(kind[~is_gene]):
        obs[col] = np.asarray(X[:, (kind == col).to_numpy()].sum(1)).ravel()
    var = gm.loc[is_gene].copy()
    var.index = genes[is_gene]
    a = ad.AnnData(X=X[:, is_gene].copy(), obs=obs, var=var)
    a.var_names_make_unique()
    a.obsm["spatial"] = np.column_stack([np.asarray(sp.spatial.x_global, dtype=np.float64),
                                         np.asarray(sp.spatial.y_global, dtype=np.float64)])
    if log_layer is not None:
        if log_layer not in sp.layers:
            raise KeyError(f"spatioloji object has no layer '{log_layer}' (available: {sorted(sp.layers)})")
        L = sp.layers[log_layer]
        L = (L.tocsr() if sps.issparse(L) else sps.csr_matrix(L))[:, is_gene].astype(np.float32)
        a.layers["counts"] = a.X.copy()
        a.X = L
    a.uns["source"] = {"package": "spatioloji_s", "n_control_features": int((~is_gene).sum())}
    return a


def load(path):
    """A saved ``spatioloji`` object: a pickle (``sp.to_pickle``) or a components folder (``sp.save_components``)."""
    cls = _spatioloji_class()
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"{path} not found")
    return cls.from_components(str(path)) if path.is_dir() else cls.from_pickle(str(path))


def read_spatioloji(path, subsample: int = 0, seed: int = 0, log_layer: str | None = None):
    """Load a saved ``spatioloji`` object and return the AnnData of :func:`to_anndata`."""
    a = to_anndata(load(path), log_layer=log_layer)
    if subsample and subsample < a.n_obs:
        idx = np.sort(np.random.default_rng(seed).choice(a.n_obs, subsample, replace=False))
        a = a[idx].copy()
    a.uns["source"]["path"] = str(path)
    return a


def _hierarchy_order(labels, groups, qc_group: str) -> list[str]:
    present = [l for l in dict.fromkeys(map(str, labels)) if l != NA_LABEL]
    return [l for _, ls in P.legend_layout(present, groups, qc_group) for l in ls]


def add_annotations(sp, labels, columns=None, prefix: str = "", overwrite: bool = True,
                    groups: dict[str, list[str]] | None = None, qc_group: str = P.QC_GROUP) -> list[str]:
    """Write annotation columns into ``sp.cell_meta``, aligned to the object's cells.

    Args:
        sp: the ``spatioloji`` object (changed in place).
        labels: a per-cell table indexed by cell id: a DataFrame, an AnnData (its ``obs``),
            or the path of ``results/final_annotation_by_cell.tsv.gz`` / a final h5ad.
        columns: columns to copy (default: every annotation column present).
        prefix: prepended to every column name, e.g. ``"ast_"``.
        overwrite: replace existing columns of the same name; with False that is an error.
        groups: legend hierarchy used to order the categories (default: the package's).

    Returns:
        The names of the columns written. Label columns are categorical with the categories
        in legend-hierarchy order; cells absent from ``labels`` get ``"NA"`` (labels) or NaN
        (scores). The lineage Leiden clusters go to ``annot_leiden`` so that a ``leiden``
        column of your own is never touched.
    """
    if isinstance(labels, (str, Path)):
        from .compare import load_labels
        labels = load_labels(labels)
    elif not isinstance(labels, pd.DataFrame):
        labels = labels.obs
    labels = labels.copy()
    labels.index = labels.index.astype(str)
    wanted = list(columns) if columns is not None else [c for c in LABEL_COLUMNS + NUMERIC_COLUMNS + OTHER_COLUMNS
                                                         if c in labels.columns]
    missing = [c for c in wanted if c not in labels.columns]
    if missing:
        raise KeyError(f"column(s) not in the annotation table: {missing}")
    target = {c: prefix + ("annot_leiden" if c == "leiden" else c) for c in wanted}
    meta = sp.cell_meta
    clash = [t for t in target.values() if t in meta.columns]
    if clash and not overwrite:
        raise ValueError(f"cell_meta already has column(s) {clash}; pass overwrite=True or a prefix")
    cells = pd.Index(sp.cell_index.astype(str))
    for col, name in target.items():
        s = labels[col].reindex(cells)
        if col in NUMERIC_COLUMNS:
            meta[name] = pd.to_numeric(s, errors="coerce").to_numpy(dtype=float)
        elif col == "proliferating":
            meta[name] = s.map(lambda v: str(v).lower() == "true" if pd.notna(v) else False).to_numpy(dtype=bool)
        else:
            vals = s.astype(object).where(s.notna(), NA_LABEL).astype(str)
            if col == "leiden":
                present = [v for v in dict.fromkeys(vals) if v != NA_LABEL]
                order = sorted(present, key=lambda v: (0, int(v)) if v.isdigit() else (1, v))
            else:
                order = _hierarchy_order(vals, groups, qc_group)
            if (vals == NA_LABEL).any():
                order = order + [NA_LABEL]
            meta[name] = pd.Categorical(vals.to_numpy(), categories=order)
    return list(target.values())


def palette(sp, column: str, cfg: PipelineConfig | None = None, config_dir=None) -> dict[str, str]:
    """Colours of the label column ``column`` of ``sp.cell_meta`` as an ordered mapping.

    The mapping has exactly the categories of the column, in their order, so it is correct
    for every spatioloji plot: functions that look colours up by label (``colors=``,
    ``color_map=``, ``palette=`` in the spatial plots) and functions that take the values
    of ``palette=`` by position (the embedding plots). The colours are the annotation's
    scheme (same label, same colour as in its own figures); ``"NA"`` is light grey."""
    from .pipeline import label_palette

    s = sp.cell_meta[column]
    cats = [str(c) for c in (s.cat.categories if isinstance(s.dtype, pd.CategoricalDtype) else sorted(s.astype(str).unique()))]
    real = [c for c in cats if c != NA_LABEL]
    pal = label_palette(cfg, real, config_dir) if cfg is not None else P.fine_palette(real)
    return {c: (NA_COLOUR if c == NA_LABEL else pal[c]) for c in cats}


@dataclass
class SpatiolojiAnnotation:
    result: object                       # annotation_st.api.AnnotationResult
    columns: list[str]                   # columns written to sp.cell_meta
    palettes: dict[str, dict[str, str]]  # label column -> ordered {label: colour}


def annotate(sp, cfg: PipelineConfig, qc: bool = False, spatial: bool = False, log_layer: str | None = None,
             prefix: str = "", overwrite: bool = True, config_dir=None) -> SpatiolojiAnnotation:
    """Annotate a ``spatioloji`` object and write the labels into ``sp.cell_meta``.

    Args:
        sp: the object (its cell metadata is changed in place; no cell is removed).
        cfg: annotation config (``annotation_st.load_config`` / ``preset_config``).
        qc: compute the QC metrics first. Nothing is removed unless ``cfg.qc.apply_filters`` is
            True (then cells that fail get ``"NA"``); QC filtering is spatioloji's job.
        spatial: also compute neighbourhood enrichment (returned in ``result.tables``).
        log_layer: use this log-normalised layer of ``sp`` instead of normalising the raw counts.
        prefix, overwrite: see :func:`add_annotations`.

    Returns:
        SpatiolojiAnnotation with the full result, the written columns, and one ordered
        palette per label column, ready for spatioloji's plots::

            out = annotate(sp, cfg)
            sj.visualization.plot_umap(sp, color_by="cell_type_fine", colors=out.palettes["cell_type_fine"])
    """
    from .api import annotate_adata

    res = annotate_adata(to_anndata(sp, log_layer=log_layer), cfg, qc=qc, spatial=spatial,
                         normalize=log_layer is None, copy=False, config_dir=config_dir)
    cols = add_annotations(sp, res.adata.obs, prefix=prefix, overwrite=overwrite, groups=cfg.plot.legend_groups,
                           qc_group=cfg.plot.qc_group)
    label_cols = [prefix + c for c in LABEL_COLUMNS if prefix + c in cols and c != "immune_compartment"]
    return SpatiolojiAnnotation(res, cols, {c: palette(sp, c, cfg, config_dir) for c in label_cols})
