"""In-memory annotation of an AnnData: same labels as the checkpointed command-line run."""
import pandas as pd
import pytest

from annotation_st.api import annotate_adata
from annotation_st.cli import main
from annotation_st.config import ConfigError, load_config
from annotation_st.io import read_xenium
from test_cli import _small_config


def test_in_memory_run_gives_the_same_labels_as_the_pipeline(tmp_path, synth_raw):
    cfgp = _small_config(tmp_path, synth_raw)
    assert main(["run", "-c", str(cfgp), "--stages", "qc", "lineage", "caf", "csc", "immune", "merge"]) == 0
    piped = pd.read_csv(tmp_path / "run" / "results" / "final_annotation_by_cell.tsv.gz", sep="\t", index_col=0, dtype=str,
                        keep_default_na=False)

    cfg = load_config(cfgp)
    raw = read_xenium(tmp_path / "xenium")                                # exactly what the pipeline read
    res = annotate_adata(raw, cfg, qc=True)
    obs = res.adata.obs
    assert list(obs.index) == list(piped.index)                           # same cells survive QC
    for col in ("lineage", "cell_type", "cell_type_fine", "caf_subtype_cell", "csc_status", "immune_subtype"):
        assert (obs[col].astype(str).to_numpy() == piped[col].to_numpy()).all(), col
    assert raw.n_obs > res.adata.n_obs and "lineage" not in raw.obs and raw.X.max() > 20   # the input is left alone
    cats = list(obs["cell_type_fine"].cat.categories)
    assert list(res.adata.uns["cell_type_fine_colors"]) == [res.palette[c] for c in cats]
    assert set(res.stages) == {"qc", "lineage", "caf", "csc", "immune", "merge"}


def test_stages_without_opted_in_markers_are_reported_as_skipped(synth_raw):
    from annotation_st.config import config_from_dict
    cfg = config_from_dict({"qc": {"area_pct": [0, 100]}, "lineage": {"n_pcs": 20},
                            "markers": {"lineage": "default", "proliferation": "default"}})
    res = annotate_adata(synth_raw, cfg, qc=True)
    assert res.stages["caf"].startswith("skipped") and res.stages["csc"].startswith("skipped")
    assert res.stages["immune"].startswith("skipped") and res.stages["lineage"] == "completed"
    obs = res.adata.obs
    assert (obs["cell_type_fine"].astype(str) == obs["cell_type"].astype(str)).all()
    with pytest.raises(ConfigError, match="markers.lineage"):
        annotate_adata(synth_raw, config_from_dict({}))


def test_spatial_statistics_on_request(tmp_path, synth_raw):
    cfg = load_config(_small_config(tmp_path, synth_raw))
    res = annotate_adata(synth_raw, cfg, qc=True, spatial=True)
    assert "spatial_connectivities" in res.adata.obsp
    z = res.tables["nhood_enrichment_zscore"]
    assert list(z.index) == list(res.adata.obs["cell_type_fine"].cat.categories) and z.shape[0] == z.shape[1]
    assert "csc_neighbourhood_composition" in res.tables
