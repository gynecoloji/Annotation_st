import pandas as pd

from annotation_st.config import CompartmentConfig, ImmuneConfig, default_config
from annotation_st.stages.immune import assign_level1, compartment_order, run_immune
from conftest import ov_config


def test_assign_level1_rules():
    comp = CompartmentConfig(lineages=["T_NK"], level1={}, level2={}, core=["CD8_T"], unresolved="TNK_unresolved",
                             cross={"cross_Myeloid": "Macrophage_Mono"}, unspecified={})
    cfg = ImmuneConfig()
    means = pd.DataFrame({"CD8_T": [-0.2, -0.9, -0.9, -0.9, -0.4, 0.5],
                          "B_cell": [0.1, 0.8, 0.1, 0.1, 0.5, 0.0],
                          "cross_Myeloid": [0, 0, 0, 0, 0, 0.9],
                          "spillover_Tumor": [0.2, 1.2, 0.7, 0.1, 0.0, 0.0]},
                         index=["c0", "c1", "c2", "c3", "c4", "c5"])
    lv, flag = assign_level1(means, ["CD8_T", "B_cell", "cross_Myeloid"], ["spillover_Tumor"], comp, cfg, "TNK")
    assert lv == {"c0": "CD8_T",                  # core above -0.5 beats minority below 0.5
                  "c1": "B_cell",                 # minority > 0.5 wins over depleted core
                  "c2": "TNK_spillover_Tumor",    # no candidate, spillover > 0.5
                  "c3": "TNK_unresolved",         # nothing qualifies
                  "c4": "B_cell",                 # minority floor is inclusive (>= 0.5), as in the source script
                  "c5": "cross_Myeloid"}          # cross set > 0.5 and highest
    assert flag == {"c0": False, "c1": True, "c2": False, "c3": False, "c4": False, "c5": False}


def test_compartment_order_puts_handoff_first():
    comps = ov_config().markers.compartments
    assert compartment_order(comps)[0] == "TNK"
    assert compartment_order(comps, ["Myeloid", "TNK"]) == ["TNK", "Myeloid"]


def test_run_immune_handoff_and_states(synth_lineage):
    cfg = ImmuneConfig(n_pcs=15, umap_cells=0, resolution=1.0)
    res = run_immune(synth_lineage, cfg, ov_config().markers, umap=False)
    assert set(res) == {"TNK", "Bcell", "Myeloid"}
    tnk, b, my = res["TNK"], res["Bcell"], res["Myeloid"]

    # planted B cells sat in the T_NK lineage and are handed to the B compartment
    assert "Bcell" in tnk.handoff and len(tnk.handoff["Bcell"]) > 60
    moved = list(tnk.handoff["Bcell"])
    assert not tnk.adata.obs_names.isin(moved).any()
    assert set(b.adata.obs.loc[moved, "immune_level1"]) <= {"B_cell", "Plasma"}
    assert b.n_prior_kept == len(moved)
    assert (b.adata.obs.loc[b.adata.obs["truth"] == "B_cell", "immune_level1"] == "B_cell").mean() > 0.6
    assert (b.adata.obs.loc[b.adata.obs["truth"] == "Plasma", "immune_subtype"] == "Plasma").mean() > 0.8

    obs = tnk.adata.obs
    assert {"TNK_leiden", "immune_level1", "immune_spillover", "immune_subtype", "immune_state_cell"} <= set(obs.columns)
    assert (obs.loc[obs["truth"] == "Treg", "immune_subtype"] == "Treg").mean() > 0.5
    assert (obs.loc[obs["truth"] == "CD8_cyto", "immune_subtype"] == "CD8_cytotoxic").mean() > 0.5
    assert (obs.loc[obs["truth"] == "NK", "immune_subtype"] == "NK").mean() > 0.8
    assert set(tnk.level2_tables) >= {"CD8_T", "CD4_T"}
    assert {"assigned", "spillover_flag", "n_cells", "handed_off", "final_label_majority"} <= set(tnk.means1.columns)
    assert tnk.top_markers is not None and {"subcluster", "level1", "final", "top_genes"} <= set(tnk.top_markers.columns)

    mo = my.adata.obs
    assert (mo.loc[mo["truth"] == "TREM2_TAM", "immune_subtype"] == "TREM2_TAM").mean() > 0.5
    assert set(mo.loc[mo["truth"] == "Macrophage", "immune_subtype"]) <= {
        "Macrophage_unpolarized", "M1_macrophage", "M2_macrophage", "TREM2_TAM", "Hypoxic_TAM"}
    assert my.sets_used["level1"]["dropped_not_on_panel"]["Macrophage"] == []


def test_run_immune_skips_small_compartment(synth_lineage):
    res = run_immune(synth_lineage, ImmuneConfig(n_pcs=15, umap_cells=0, min_cells=100_000),
                     ov_config().markers, umap=False)
    assert res == {}
