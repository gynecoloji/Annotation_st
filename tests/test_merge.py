import anndata as ad
import numpy as np
import pandas as pd

from annotation_st.config import default_config
from annotation_st.stages.merge import merge_labels
from conftest import ov_config


def _tiny(cell_types):
    n = len(cell_types)
    a = ad.AnnData(X=np.zeros((n, 2), dtype=np.float32),
                   obs=pd.DataFrame({"cell_type": pd.Categorical(cell_types)}, index=[f"c{i}" for i in range(n)]))
    return a


def test_merge_precedence_and_coarse_update():
    full = _tiny(["Fibroblast", "Epithelial_tumor_cycling", "T_NK", "T_NK", "Plasma", "Endothelial"])
    fb = pd.DataFrame({"caf_subtype_cell": ["myCAF"], "caf_subtype_cluster": ["iCAF"]}, index=["c0"])
    tum = pd.DataFrame({"csc_status": ["CSC_like"], "csc_OCSC_z": [2.0], "csc_n_core_detected": [3]}, index=["c1"])
    imm = {"TNK": pd.DataFrame({"immune_subtype": ["Treg"], "immune_level1": ["CD4_T"]}, index=["c2"]),
           "Bcell": pd.DataFrame({"immune_subtype": ["Memory_B", "Plasma"], "immune_level1": ["B_cell", "Plasma"]},
                                 index=["c3", "c4"])}
    res = merge_labels(full, fb, tum, imm, ov_config().markers.compartments)
    obs = res.adata.obs
    assert list(obs["cell_type_fine"]) == ["myCAF", "CSC_like", "Treg", "Memory_B", "Plasma", "Endothelial"]
    assert list(obs["cell_type"]) == ["Fibroblast", "Epithelial_tumor_cycling", "T_NK", "B_cell", "Plasma", "Endothelial"]
    assert list(obs["immune_compartment"]) == ["NA", "NA", "TNK", "Bcell", "Bcell", "NA"]
    assert list(obs["caf_subtype_cluster"]) == ["iCAF", "NA", "NA", "NA", "NA", "NA"]
    assert obs["csc_OCSC_z"].tolist()[0] != obs["csc_OCSC_z"].tolist()[0]     # NaN outside tumour
    assert obs["csc_OCSC_z"].iloc[1] == 2.0
    assert res.fine_counts["n_cells"].sum() == 6 and res.coarse_counts.loc["B_cell", "n_cells"] == 1


def test_merge_tolerates_missing_pieces():
    full = _tiny(["Fibroblast", "T_NK"])
    res = merge_labels(full, None, None, {}, ov_config().markers.compartments)
    assert list(res.adata.obs["cell_type_fine"]) == ["Fibroblast", "T_NK"]
