import numpy as np
import pytest
import pandas as pd

from annotation_st.config import CscConfig, default_config
from annotation_st.stages.csc import call_csc
from conftest import ov_config


def test_csc_rule(synth_tum):
    res = call_csc(synth_tum, CscConfig(n_pcs=15, umap_cells=0), ov_config().markers, umap=False)
    obs = res.adata.obs
    rule = (obs["csc_OCSC_z"] > 1.5) & (obs["csc_n_core_detected"] >= 2)
    assert ((obs["csc_status"] == "CSC_like") == rule).all()
    ct = pd.crosstab(obs["truth"], obs["csc_status"])
    assert ct.loc["CSC", "CSC_like"] / ct.loc["CSC"].sum() > 0.5
    assert ct.loc["Epithelial_tumor", "CSC_like"] / ct.loc["Epithelial_tumor"].sum() < 0.1
    assert res.n_csc == int(rule.sum()) == res.counts.loc["CSC_like", "n_cells"]
    assert set(res.positivity.columns) == {"frac_positive_CSC_like", "frac_positive_Tumor_nonCSC", "frac_positive_all_tumor"}
    assert res.positivity.loc["CD44", "frac_positive_CSC_like"] > res.positivity.loc["CD44", "frac_positive_Tumor_nonCSC"]
    assert res.de is not None and "CD44" in set(res.de["names"].head(10))
    assert obs["tumor_leiden"].nunique() >= 2 and "X_umap" not in res.adata.obsm


def test_csc_umap_subset(synth_tum):
    res = call_csc(synth_tum, CscConfig(n_pcs=15, umap_cells=200), ov_config().markers, umap=True)
    assert np.isnan(res.adata.obsm["X_umap"][:, 0]).sum() == synth_tum.n_obs - 200
    assert (res.adata.obs["tumor_leiden"] == "NA").sum() == synth_tum.n_obs - 200


def test_csc_too_few_cells(synth_tum):
    assert call_csc(synth_tum[:5].copy(), CscConfig(), ov_config().markers) is None


def test_custom_marker_file_replaces_default(tmp_path, synth_tum):
    from annotation_st.stages.csc import resolve_csc_markers
    f = tmp_path / "my_csc.txt"
    f.write_text("# my stemness list\nCD44\nPROM1\nCXCR4\nSOX2\nNOTAGENE\n")
    cfg = CscConfig(n_pcs=15, umap_cells=0, marker_file=str(f))
    csc, core, src = resolve_csc_markers(cfg, ov_config().markers)
    assert csc == ["CD44", "PROM1", "CXCR4", "SOX2", "NOTAGENE"]
    assert core == ["CD44", "PROM1", "SOX2"]               # default core restricted to the custom set
    assert src["csc"] == str(f) and "restricted" in src["core"]
    res = call_csc(synth_tum, cfg, ov_config().markers, umap=False)
    assert res.sets_used["OCSC"] == ["CD44", "PROM1", "CXCR4", "SOX2"] and res.sets_dropped["OCSC"] == ["NOTAGENE"]
    assert res.sets_used["OCSC_core"] == ["CD44", "PROM1", "SOX2"]
    assert list(res.positivity.index) == ["CD44", "PROM1", "CXCR4", "SOX2"]


def test_custom_core_file_and_fallback(tmp_path):
    from annotation_st.stages.csc import resolve_csc_markers
    (tmp_path / "set.txt").write_text("GENEA\nGENEB\nGENEC\n")
    (tmp_path / "core.txt").write_text("GENEA\n")
    m = ov_config().markers
    # relative paths resolve against base
    csc, core, src = resolve_csc_markers(CscConfig(marker_file="set.txt", core_marker_file="core.txt"), m, base=tmp_path)
    assert csc == ["GENEA", "GENEB", "GENEC"] and core == ["GENEA"]
    # no core file and no overlap with the default core -> every gene of the set is core
    csc, core, src = resolve_csc_markers(CscConfig(marker_file="set.txt"), m, base=tmp_path)
    assert core == csc and src["core"].startswith("all genes")
    # inline YAML list without a file
    m2 = ov_config().markers; m2.csc = ["CD44", "KIT"]; m2.csc_core = ["KIT"]
    assert resolve_csc_markers(CscConfig(), m2)[:2] == (["CD44", "KIT"], ["KIT"])


def test_call_csc_explicit_gene_lists(synth_tum):
    res = call_csc(synth_tum, CscConfig(n_pcs=15, umap_cells=0), umap=False, embed_cells=False,
                   csc_genes=["CD44", "PROM1", "POU5F1"], core_genes=["CD44", "PROM1"])
    assert res.sets_used["OCSC"] == ["CD44", "PROM1", "POU5F1"] and res.n_csc > 0
    with pytest.raises(ValueError):
        call_csc(synth_tum, CscConfig(), csc_genes=["CD44"])
