import numpy as np

from annotation_st.markers import LINEAGE_MARKERS, filter_to_panel
from annotation_st.scoring import cluster_assign, embed, per_cell_assign, score_sets, zscore


def test_zscore_constant_vector_is_zero():
    assert np.allclose(zscore(np.ones(5)), 0)


def test_per_cell_assign_margin_rule():
    Z = np.array([[1.0, 0.5], [1.0, 0.9], [-0.1, -0.5], [0.0, -1.0]])
    call, best, marg = per_cell_assign(Z, ["a", "b"], min_z=0.0, margin=0.25, unassigned="U")
    assert list(call) == ["a", "U", "U", "U"]          # margin fail, floor fail, floor (<=) fail
    assert np.allclose(best, [1.0, 1.0, -0.1, 0.0]) and np.allclose(marg, [0.5, 0.1, 0.4, 1.0])


def test_score_sets_writes_raw_and_z(synth):
    sets, _ = filter_to_panel({"Tumor": LINEAGE_MARKERS["Epithelial_tumor"], "Empty": ["NOPE"]}, synth.var_names)
    scored = score_sets(synth, sets, "s_")
    assert scored == ["Tumor"]
    assert abs(synth.obs["s_Tumor_z"].mean()) < 1e-6 and abs(synth.obs["s_Tumor_z"].std() - 1) < 1e-2
    tum = synth.obs["truth_coarse"] == "Epithelial_tumor"
    assert synth.obs.loc[tum, "s_Tumor_z"].mean() > 1 > synth.obs.loc[~tum, "s_Tumor_z"].mean()


def test_cluster_assign_threshold(synth):
    sets, _ = filter_to_panel(LINEAGE_MARKERS, synth.var_names)
    scored = score_sets(synth, sets, "score_")
    means, mapping = cluster_assign(synth, "truth_coarse", scored, "score_", min_z=0.5)
    for lin in ["Fibroblast", "Epithelial_tumor", "T_NK", "Macrophage_Mono", "Plasma", "Endothelial", "Ovarian_stroma"]:
        assert mapping[lin] == lin
    assert list(means.columns) == scored
    # a floor above every mean gives Unassigned
    _, m2 = cluster_assign(synth, "truth_coarse", scored, "score_", min_z=100)
    assert set(m2.values()) == {"Unassigned"}


def test_embed_umap_subset(synth):
    embed(synth, 20, 15, 0.5, 0, "leiden", umap_cells=500)
    assert "leiden" in synth.obs and synth.obs["leiden"].nunique() >= 5
    assert np.isnan(synth.obsm["X_umap"][:, 0]).sum() == synth.n_obs - 500
    assert len(synth.uns["leiden_umap_subset_idx"]) == 500


def test_embed_no_umap(synth):
    embed(synth, 10, 10, 0.5, 0, "cl", umap=False)
    assert "X_umap" not in synth.obsm and "cl" in synth.obs
