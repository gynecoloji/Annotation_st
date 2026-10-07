import pandas as pd

from annotation_st.config import LineageConfig, default_config
from annotation_st.stages.lineage import annotate_lineage
from conftest import ov_config


def test_lineage_recovers_planted_populations(synth_counts):
    a = synth_counts.copy()
    res = annotate_lineage(a, LineageConfig(n_pcs=20, skip_umap=True), ov_config().markers)
    obs = res.adata.obs
    assert {"leiden", "lineage", "cell_type", "proliferating", "score_Epithelial_tumor_z"} <= set(obs.columns)
    ct = pd.crosstab(obs["truth_coarse"], obs["lineage"])
    for lin in ["Epithelial_tumor", "Fibroblast", "T_NK", "Macrophage_Mono", "Plasma", "Endothelial", "Ovarian_stroma"]:
        assert ct.loc[lin, lin] / ct.loc[lin].sum() > 0.8, ct
    cyc = obs["cell_type"] == "Epithelial_tumor_cycling"
    assert cyc.sum() > 0 and (obs.loc[obs["truth"] == "Tumor_cycling", "cell_type"] == "Epithelial_tumor_cycling").mean() > 0.8
    assert set(res.means.columns) >= {"assigned", "n_cells", "mean_prolif_z"}
    assert set(res.mapping) == set(res.means.index)
    assert res.counts["n_cells"].sum() == a.n_obs
    assert res.top_markers is not None and {"cluster", "lineage", "top_genes"} <= set(res.top_markers.columns)
    assert res.adata.raw is not None and "Epithelial_tumor" in res.sets_used


def test_lineage_unassigned_floor(synth):
    res = annotate_lineage(synth, LineageConfig(n_pcs=10, skip_umap=True, min_z=50), ov_config().markers,
                           normalize=False, compute_top_markers=False)
    assert set(res.adata.obs["lineage"]) == {"Unassigned"}
    assert res.top_markers is None
