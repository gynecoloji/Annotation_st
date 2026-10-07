"""Step 5 — merge every refinement into ``obs['cell_type_fine']`` and update ``obs['cell_type']``.

``cell_type_fine`` = CAF subtype for fibroblasts, CSC status for tumour cells (cycling
tumour cells fold into CSC_like / Tumor_nonCSC; the cycling flag stays in ``cell_type``),
immune subtype for the immune compartments, coarse label otherwise. ``cell_type`` changes
only for cells a compartment's ``coarse`` map moves between lineages (B / plasma / DC /
cross-compartment).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from ..config import CompartmentConfig


@dataclass
class MergeResult:
    adata: object
    fine_counts: pd.DataFrame
    coarse_counts: pd.DataFrame


def _common(idx, obs_names):
    return idx[idx.isin(obs_names)]


def merge_labels(full, fb_obs: pd.DataFrame | None, tum_obs: pd.DataFrame | None,
                 immune: dict[str, pd.DataFrame], compartments: dict[str, CompartmentConfig]) -> MergeResult:
    """``fb_obs`` / ``tum_obs`` / ``immune[name]`` are the ``obs`` tables of the stage outputs
    (indexed by cell). Operates in place on ``full``."""
    fine = full.obs["cell_type"].astype(str).copy()
    coarse = fine.copy()

    for col in ("caf_subtype_cell", "caf_subtype_cluster", "csc_status", "immune_subtype", "immune_compartment"):
        full.obs[col] = "NA"
    for col in ("csc_OCSC_z", "csc_n_core_detected"):
        full.obs[col] = np.nan

    if fb_obs is not None and len(fb_obs):
        ix = _common(fb_obs.index, full.obs_names)
        fine.loc[ix] = fb_obs.loc[ix, "caf_subtype_cell"].astype(str).values
        for col in ("caf_subtype_cell", "caf_subtype_cluster"):
            if col in fb_obs:
                full.obs.loc[ix, col] = fb_obs.loc[ix, col].astype(str).values
    if tum_obs is not None and len(tum_obs):
        ix = _common(tum_obs.index, full.obs_names)
        fine.loc[ix] = tum_obs.loc[ix, "csc_status"].astype(str).values
        full.obs.loc[ix, "csc_status"] = tum_obs.loc[ix, "csc_status"].astype(str).values
        for col in ("csc_OCSC_z", "csc_n_core_detected"):
            if col in tum_obs:
                full.obs.loc[ix, col] = tum_obs.loc[ix, col].to_numpy(dtype=float)
    for name, imm in immune.items():
        if imm is None or not len(imm):
            continue
        ix = _common(imm.index, full.obs_names)
        fine.loc[ix] = imm.loc[ix, "immune_subtype"].astype(str).values
        full.obs.loc[ix, "immune_subtype"] = imm.loc[ix, "immune_subtype"].astype(str).values
        full.obs.loc[ix, "immune_compartment"] = name
        cmap = compartments[name].coarse if name in compartments else {}
        if cmap and "immune_level1" in imm:
            lv1 = imm.loc[ix, "immune_level1"].astype(str)
            hit = lv1.isin(list(cmap))
            coarse.loc[ix[hit.to_numpy()]] = lv1[hit].map(cmap).values

    full.obs["cell_type"] = pd.Categorical(coarse)
    full.obs["cell_type_fine"] = pd.Categorical(fine)
    for col in ("caf_subtype_cell", "caf_subtype_cluster", "csc_status", "immune_subtype", "immune_compartment"):
        full.obs[col] = full.obs[col].astype("category")

    fine_counts = full.obs["cell_type_fine"].value_counts().to_frame("n_cells")
    fine_counts["fraction"] = fine_counts["n_cells"] / full.n_obs
    coarse_counts = full.obs["cell_type"].value_counts().to_frame("n_cells")
    coarse_counts["fraction"] = coarse_counts["n_cells"] / full.n_obs
    return MergeResult(full, fine_counts, coarse_counts)
