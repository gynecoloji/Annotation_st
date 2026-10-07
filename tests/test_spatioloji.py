"""Compatibility with the spatioloji_s package: annotate a spatioloji object, write the labels
into its cell metadata, hand colours to its plots, and use a saved object as pipeline input."""
import json

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp_sparse

sj = pytest.importorskip("spatioloji_s")

from annotation_st import spatioloji_bridge as bridge  # noqa: E402
from annotation_st.api import annotate_adata  # noqa: E402
from annotation_st.cli import main  # noqa: E402
from annotation_st.config import ConfigError, dump_config, load_config  # noqa: E402
from annotation_st.io import read_xenium  # noqa: E402
from annotation_st.pipeline import label_palette, preflight  # noqa: E402
from conftest import write_synthetic_xenium  # noqa: E402
from test_cli import _small_config  # noqa: E402

N_CTRL = 3


def _spatioloji(a, feature_type=True, extra_meta=None):
    """A spatioloji object as its Xenium loader builds it: every feature in the matrix (control
    probes included), cells.parquet columns as cell metadata, local == global coordinates."""
    rng = np.random.default_rng(3)
    ctrl = sp_sparse.csr_matrix(rng.poisson(0.05, (a.n_obs, N_CTRL)).astype(np.float32))
    X = sp_sparse.hstack([sp_sparse.csr_matrix(a.X), ctrl]).tocsr()
    names = list(a.var_names) + [f"NegControlProbe_{i:05d}" for i in range(N_CTRL)]
    gene_meta = pd.DataFrame(index=names)
    if feature_type:
        gene_meta["feature_type"] = ["Gene Expression"] * a.n_vars + ["Negative Control Probe"] * N_CTRL
    meta = a.obs[["cell_area", "nucleus_area", "nucleus_count"]].copy()
    meta["fov"] = "xenium"
    if extra_meta:
        for k, v in extra_meta.items():
            meta[k] = v
    xy = np.asarray(a.obsm["spatial"])
    obj = sj.spatioloji(expression=X, cell_ids=list(a.obs_names), gene_names=names, cell_metadata=meta,
                        gene_metadata=gene_meta,
                        spatial_coords={"x_local": xy[:, 0], "y_local": xy[:, 1], "x_global": xy[:, 0], "y_global": xy[:, 1]})
    return obj, np.asarray(ctrl.sum(1)).ravel()


# ---- spatioloji -> AnnData ---------------------------------------------------------------------
@pytest.mark.parametrize("feature_type", [True, False])
def test_to_anndata_keeps_genes_and_counts_the_control_features(synth_raw, feature_type):
    obj, ctrl_sum = _spatioloji(synth_raw, feature_type=feature_type)
    a = bridge.to_anndata(obj)
    assert list(a.var_names) == list(synth_raw.var_names)                  # control probes are not genes
    assert list(a.obs_names) == list(synth_raw.obs_names)
    assert (a.X != sp_sparse.csr_matrix(synth_raw.X)).nnz == 0             # raw counts, untouched
    assert np.array_equal(a.obs["neg_probe_counts"].to_numpy(), ctrl_sum)  # ... but they feed the QC
    assert a.obsm["spatial"].shape == (a.n_obs, 2)                         # (x, y) only: a 4-column array breaks spatial graphs
    assert np.allclose(a.obsm["spatial"], synth_raw.obsm["spatial"])
    assert {"cell_area", "nucleus_count"} <= set(a.obs.columns)


def test_to_anndata_can_reuse_a_log_normalised_layer(synth_raw):
    obj, _ = _spatioloji(synth_raw)
    sj.processing.normalize_total(obj, inplace=True)
    sj.processing.log_transform(obj, layer="normalized_counts", inplace=True)
    a = bridge.to_anndata(obj, log_layer="log_normalized")
    assert a.n_vars == synth_raw.n_vars and "counts" in a.layers
    assert (a.layers["counts"] != sp_sparse.csr_matrix(synth_raw.X)).nnz == 0
    assert a.X.max() < 20 and a.X.max() != a.layers["counts"].max()        # X holds the log values
    with pytest.raises(KeyError, match="no_such_layer"):
        bridge.to_anndata(obj, log_layer="no_such_layer")


# ---- annotate in place -------------------------------------------------------------------------
@pytest.fixture
def cfg(tmp_path, synth_raw):
    return load_config(_small_config(tmp_path, synth_raw))


def test_annotate_writes_labels_into_cell_meta(synth_raw, cfg):
    obj, _ = _spatioloji(synth_raw, extra_meta={"leiden": "mine"})        # the user's own clustering column
    out = bridge.annotate(obj, cfg, qc=True)
    meta = obj.cell_meta
    assert len(meta) == synth_raw.n_obs                                    # the object keeps all its cells
    for col in ("lineage", "cell_type", "cell_type_fine", "caf_subtype_cell", "csc_status", "immune_subtype"):
        assert col in out.columns and isinstance(meta[col].dtype, pd.CategoricalDtype), col
    assert (meta["leiden"] == "mine").all() and "annot_leiden" in meta     # theirs is not overwritten

    ref = annotate_adata(bridge.to_anndata(obj_copy := _spatioloji(synth_raw)[0]), cfg, qc=True).adata.obs
    kept = meta.index.isin(ref.index)
    assert (meta.loc[ref.index, "cell_type_fine"].astype(str) == ref["cell_type_fine"].astype(str)).all()
    assert (meta.loc[~kept, "cell_type_fine"] == "NA").all() and (~kept).sum() > 0        # cells that failed QC
    assert meta.loc[~kept, "csc_OCSC_z"].isna().all()
    assert obj_copy.n_cells == obj.n_cells


def test_palette_matches_every_way_spatioloji_plots_take_colours(synth_raw, cfg):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.colors as mc

    obj, _ = _spatioloji(synth_raw)
    out = bridge.annotate(obj, cfg, qc=True)
    pal = bridge.palette(obj, "cell_type_fine", cfg)
    cats = list(obj.cell_meta["cell_type_fine"].cat.categories)
    assert list(pal) == cats                                               # same labels, same order, nothing extra
    full = label_palette(cfg, cats)
    assert all(pal[c] == full[c] for c in cats if c != "NA") and pal == out.palettes["cell_type_fine"]
    # categories follow the legend hierarchy (tumour first), not the alphabet
    assert cats.index("Tumor_nonCSC") < cats.index("iCAF") < cats.index("NK")

    obj.embeddings["X_umap"] = np.asarray(bridge.to_anndata(obj).obsm["spatial"])
    for kw in ({"colors": pal}, {"palette": pal}):                         # by key, and by position
        fig = sj.visualization.plot_umap(obj, color_by="cell_type_fine", show=False, **kw)
        ax = fig.axes[0]
        drawn = {t.get_text(): mc.to_hex(h.get_facecolor()[0]) for h, t in
                 zip(ax.get_legend().legend_handles, ax.get_legend().get_texts())}
        assert drawn == {c: mc.to_hex(pal[c]) for c in cats}, kw
        matplotlib.pyplot.close(fig)


def test_add_annotations_from_a_finished_run(tmp_path, synth_raw):
    cfgp = _small_config(tmp_path, synth_raw)
    assert main(["run", "-c", str(cfgp), "--stages", "qc", "lineage", "csc", "merge"]) == 0
    table = tmp_path / "run" / "results" / "final_annotation_by_cell.tsv.gz"
    obj, _ = _spatioloji(synth_raw)
    added = bridge.add_annotations(obj, table, prefix="ast_")
    assert "ast_cell_type_fine" in added and "ast_cell_type" in added
    labels = pd.read_csv(table, sep="\t", index_col=0, dtype=str, keep_default_na=False)
    assert (obj.cell_meta.loc[labels.index, "ast_cell_type_fine"].astype(str) == labels["cell_type_fine"]).all()
    with pytest.raises(ValueError, match="ast_cell_type"):
        bridge.add_annotations(obj, table, prefix="ast_", overwrite=False)


# ---- a saved spatioloji object as pipeline input ------------------------------------------------
def test_pipeline_accepts_a_spatioloji_pickle(tmp_path, synth_raw):
    cfgp = _small_config(tmp_path, synth_raw)                              # reads the Xenium folder
    assert main(["qc", "-c", str(cfgp)]) == 0
    import anndata as ad
    from_xenium = ad.read_h5ad(tmp_path / "run" / "checkpoints" / "sample_01_qc.h5ad")

    obj, _ = _spatioloji(read_xenium(tmp_path / "xenium"))                 # the same data as a spatioloji object
    obj.to_pickle(str(tmp_path / "sample.pkl"))
    cfg = load_config(cfgp)
    cfg.project.xenium_dir = None
    cfg.project.spatioloji = "sample.pkl"                                  # relative to the config file
    cfg.project.outdir = str(tmp_path / "run_sj")
    dump_config(cfg, cfgp)
    assert main(["run", "-c", str(cfgp), "--stages", "qc", "lineage", "csc", "merge"]) == 0
    from_sj = ad.read_h5ad(tmp_path / "run_sj" / "checkpoints" / "sample_01_qc.h5ad")
    assert list(from_sj.obs_names) == list(from_xenium.obs_names) and list(from_sj.var_names) == list(from_xenium.var_names)
    assert (from_sj.X != from_xenium.X).nnz == 0
    params = json.loads((tmp_path / "run_sj" / "params.json").read_text())
    assert params["qc_filter"]["input"].endswith("sample.pkl")

    # ... and the labels go back into a spatioloji object
    out = tmp_path / "annotated.pkl"
    assert main(["to-spatioloji", "-c", str(cfgp), "-o", str(out)]) == 0
    back = sj.spatioloji.from_pickle(str(out))
    fine = pd.read_csv(tmp_path / "run_sj" / "results" / "final_cell_type_counts.tsv", sep="\t", index_col=0)["n_cells"]
    got = back.cell_meta["cell_type_fine"].astype(str).value_counts()
    assert got.drop("NA", errors="ignore").to_dict() == fine.to_dict()


def test_input_must_be_one_of_xenium_folder_or_spatioloji_object(tmp_path, synth_raw):
    cfg = load_config(_small_config(tmp_path, synth_raw))
    cfg.project.xenium_dir = None
    with pytest.raises(ConfigError, match="project.spatioloji"):
        preflight(cfg, ["qc"])
    cfg.project.spatioloji = str(tmp_path / "missing.pkl")
    with pytest.raises(ConfigError, match="missing.pkl"):
        preflight(cfg, ["qc"])
    cfg.project.xenium_dir = str(tmp_path / "xenium")
    with pytest.raises(ConfigError, match="only one"):
        preflight(cfg, ["qc"])


def test_qc_can_keep_every_cell_of_an_already_filtered_object(synth_raw):
    from annotation_st.config import QCConfig
    from annotation_st.stages.qc import qc_filter
    a = bridge.to_anndata(_spatioloji(synth_raw)[0])
    res = qc_filter(a, QCConfig(apply_filters=False))
    assert res.adata.n_obs == a.n_obs and res.adata.n_vars == a.n_vars
    assert res.flag_summary.loc["ANY", "n_cells_flagged"] > 0 and "counts" in res.adata.layers   # still reported
