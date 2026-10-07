import numpy as np

from annotation_st.config import SpatialConfig
from annotation_st.stages.spatial import build_spatial_graph, n_edges, neighbourhood_composition, nhood_enrichment


def test_graph_edge_cut_and_composition(synth):
    cut = build_spatial_graph(synth, SpatialConfig())
    D = synth.obsp["spatial_distances"]
    assert D.data.max() <= cut and n_edges(synth) > synth.n_obs
    synth.obs["lab"] = synth.obs["truth_coarse"]
    comp = neighbourhood_composition(synth, "lab", ("Epithelial_tumor", "Fibroblast"))
    assert "log2_ratio_Epithelial_tumor_vs_Fibroblast" in comp.columns
    assert comp.loc["Epithelial_tumor", "Epithelial_tumor"] > comp.loc["Epithelial_tumor", "Fibroblast"]
    assert comp.loc["Epithelial_tumor", "log2_ratio_Epithelial_tumor_vs_Fibroblast"] > 0
    assert np.allclose(comp[["Epithelial_tumor", "Fibroblast"]].sum(), 1, atol=1e-5)


def test_nhood_enrichment(synth):
    build_spatial_graph(synth, SpatialConfig())
    synth.obs["lab"] = synth.obs["truth_coarse"].astype(str)
    z, c = nhood_enrichment(synth, "lab", SpatialConfig(n_perms=20), seed=0)
    assert np.isfinite(z.values).all() and z.shape == c.shape == (7, 7)
    assert z.loc["Epithelial_tumor", "Epithelial_tumor"] > 0 > z.loc["Epithelial_tumor", "Endothelial"]


def test_radius_graph(synth):
    build_spatial_graph(synth, SpatialConfig(radius=40.0, edge_cut_pct=100))
    assert synth.obsp["spatial_distances"].data.max() <= 40.0
