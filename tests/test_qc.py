import numpy as np

from annotation_st.config import QCConfig, default_config
from annotation_st.io import RunPaths, read_xenium, update_params
from annotation_st.stages.qc import QC_METRICS, compute_qc_metrics, qc_filter
from conftest import write_synthetic_xenium, ov_config


def test_no_cell_or_gene_is_removed_by_default(synth_raw):
    """QC filtering belongs to the upstream tools (e.g. spatioloji_s): this package only reports."""
    assert QCConfig().apply_filters is False
    a = synth_raw.copy()
    res = qc_filter(a, QCConfig())
    assert res.adata.n_obs == a.n_obs and res.adata.n_vars == a.n_vars
    assert list(res.adata.obs_names) == list(a.obs_names) and list(res.adata.var_names) == list(a.var_names)
    assert res.flag_summary.loc["ANY", "n_cells_flagged"] >= 35                 # the flags are still reported ...
    assert res.flag_summary.loc["KEPT", "n_cells_flagged"] < a.n_obs            # ... as what filtering would do
    assert "counts" in res.adata.layers and {"gene_counts", "control_frac"} <= set(res.adata.obs.columns)


def test_qc_filter_removes_planted_failures_when_asked(synth_raw):
    a = synth_raw.copy()
    compute_qc_metrics(a)
    assert {"gene_counts", "control_frac", "nucleus_ratio", "counts_per_area"} <= set(a.obs.columns)
    res = qc_filter(a, QCConfig(area_pct=(0, 100), apply_filters=True))
    fs = res.flag_summary
    assert fs.loc["low_transcripts", "n_cells_flagged"] >= 20
    assert fs.loc["no_nucleus", "n_cells_flagged"] == 10
    assert fs.loc["high_control_frac", "n_cells_flagged"] >= 5
    assert fs.loc["area_outlier", "n_cells_flagged"] == 0
    assert fs.loc["KEPT", "n_cells_flagged"] == res.adata.n_obs == a.n_obs - fs.loc["ANY", "n_cells_flagged"]
    assert "counts" in res.adata.layers
    assert res.summary_pre.shape[0] == len(QC_METRICS) and res.summary_post.shape[0] == len(QC_METRICS)


def test_area_percentile_window(synth_counts):
    a = synth_counts.copy()
    res = qc_filter(a, QCConfig(area_pct=(5, 95), min_transcripts=0, min_genes=0, apply_filters=True))
    frac = res.flag_summary.loc["area_outlier", "fraction"]
    assert 0.08 < frac < 0.12
    lo, hi = res.area_bounds
    assert res.adata.obs["cell_area"].between(lo, hi).all()


def test_read_xenium_roundtrip(tmp_path, synth_counts):
    write_synthetic_xenium(tmp_path, synth_counts)
    a = read_xenium(tmp_path)
    assert a.shape == synth_counts.shape
    assert np.allclose(a.X.sum(), synth_counts.X.sum())
    assert list(a.var_names) == list(synth_counts.var_names)
    assert {"neg_probe_counts", "cell_area", "nucleus_count"} <= set(a.obs.columns)
    assert a.obsm["spatial"].shape == (a.n_obs, 2)
    b = read_xenium(tmp_path, subsample=100, seed=1)
    assert b.n_obs == 100 and np.allclose(b.X.sum(1).ravel(), synth_counts[b.obs_names].X.sum(1).ravel())


def test_runpaths_and_params(tmp_path):
    cfg = ov_config()
    cfg.project.outdir = str(tmp_path / "run")
    cfg.project.prefix = "ov"
    p = RunPaths(cfg)
    assert p.checkpoint("lineage").name == "ov_02_lineage.h5ad"
    assert p.checkpoint("immune_TNK").name == "ov_05_immune_TNK.h5ad"
    assert p.checkpoint("final").parent == tmp_path / "run" / "checkpoints"
    assert p.results.is_dir() and p.qc.is_dir()
    update_params(p, "a", {"x": 1})
    update_params(p, "b", {"y": 2})
    import json
    assert json.loads(p.params.read_text()) == {"a": {"x": 1}, "b": {"y": 2}}
    cfg.project.checkpoint_dir = str(tmp_path / "proc")
    assert RunPaths(cfg).checkpoint("qc").parent == tmp_path / "proc"
