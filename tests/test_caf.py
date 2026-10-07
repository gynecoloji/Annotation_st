import pandas as pd
import pytest

from annotation_st.config import CafConfig
from annotation_st.stages.caf import annotate_caf


def test_caf_per_cell_and_cluster_labels(synth_fb, caf_sets):
    res = annotate_caf(synth_fb, CafConfig(n_pcs=15), caf_sets, umap=False)
    obs = res.adata.obs
    ct = pd.crosstab(obs["truth"], obs["caf_subtype_cell"])
    assert ct.loc["iCAF", "iCAF"] > ct.loc["iCAF"].drop("iCAF").max()
    assert ct.loc["myCAF", "myCAF"] > ct.loc["myCAF"].drop("myCAF").max()
    assert list(obs["caf_subtype_cell"].cat.categories) == ["iCAF", "myCAF", "apCAF", "FB_unassigned"]
    # per-cell rule holds exactly
    called = obs["caf_subtype_cell"] != "FB_unassigned"
    assert ((obs["caf_best_z"] > 0) & (obs["caf_margin"] >= 0.25)).equals(called)
    assert set(res.counts.columns) == {"per_cell", "per_subcluster", "frac_per_cell"}
    assert {"assigned", "n_cells", "frac_cell_iCAF", "frac_cell_FB_unassigned"} <= set(res.means.columns)
    assert res.scored == ["iCAF", "myCAF", "apCAF"]
    assert "X_umap" not in res.adata.obsm


def test_caf_strict_margin_gives_more_unassigned(synth_fb, caf_sets):
    loose = annotate_caf(synth_fb.copy(), CafConfig(n_pcs=15, margin=0.0), caf_sets, umap=False)
    strict = annotate_caf(synth_fb.copy(), CafConfig(n_pcs=15, margin=2.0), caf_sets, umap=False)
    assert (strict.adata.obs["caf_subtype_cell"] == "FB_unassigned").sum() > \
           (loose.adata.obs["caf_subtype_cell"] == "FB_unassigned").sum()


def test_caf_too_few_cells_returns_none(synth_fb, caf_sets):
    assert annotate_caf(synth_fb[:10].copy(), CafConfig(), caf_sets) is None


def test_caf_no_genes_on_panel_raises(synth_fb):
    with pytest.raises(ValueError):
        annotate_caf(synth_fb, CafConfig(n_pcs=15), {"iCAF": ["NOPE1"], "myCAF": ["NOPE2"]}, umap=False)
