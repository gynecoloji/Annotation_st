"""Step 2 — CAF subtypes (iCAF / myCAF / apCAF) within the fibroblast compartment.

Fibroblasts are re-embedded; each CAF set is scored per cell and z-scored *within
fibroblasts*. Primary label ``caf_subtype_cell``: per-cell argmax with best z > ``min_z``
and best − second ≥ ``margin`` (else ``unassigned``). Secondary label
``caf_subtype_cluster``: per-sub-cluster argmax of mean z (> ``min_z``).
"""
from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from ..config import CafConfig
from ..markers import filter_to_panel
from ..scoring import cluster_assign, embed, per_cell_assign, score_sets


@dataclass
class CafResult:
    adata: object
    means: pd.DataFrame        # fb_leiden × set mean z + assigned, n_cells, frac_cell_<s>
    counts: pd.DataFrame       # per_cell, per_subcluster, frac_per_cell
    sets_used: dict[str, list[str]]
    sets_dropped: dict[str, list[str]]
    scored: list[str]


def annotate_caf(fb, cfg: CafConfig, caf_sets: dict[str, list[str]], seed: int = 0,
                 umap: bool = True) -> CafResult | None:
    """``fb`` = fibroblast subset with log-normalised ``X``. Returns None when
    ``fb.n_obs < cfg.min_cells``."""
    if fb.n_obs < cfg.min_cells:
        return None
    embed(fb, cfg.n_pcs, cfg.n_neighbors, cfg.resolution, seed, "fb_leiden", umap=umap)

    sets, dropped = filter_to_panel(caf_sets, fb.var_names)
    scored = score_sets(fb, sets, prefix="caf_", seed=seed)
    if not scored:
        raise ValueError("no CAF marker set has any gene on this panel")

    Z = fb.obs[[f"caf_{s}_z" for s in scored]].to_numpy()
    call, best, marg = per_cell_assign(Z, scored, cfg.min_z, cfg.margin, cfg.unassigned)
    fb.obs["caf_subtype_cell"] = pd.Categorical(call, categories=scored + [cfg.unassigned])
    fb.obs["caf_best_z"] = best
    fb.obs["caf_margin"] = marg

    means, mapping = cluster_assign(fb, "fb_leiden", scored, prefix="caf_", min_z=cfg.min_z,
                                    unassigned=cfg.unassigned)
    means["assigned"] = pd.Series(mapping).astype(str)
    means["n_cells"] = fb.obs["fb_leiden"].value_counts()
    for s in scored + [cfg.unassigned]:
        means[f"frac_cell_{s}"] = fb.obs.groupby("fb_leiden", observed=True)["caf_subtype_cell"].apply(
            lambda x, s=s: float((x == s).mean()))
    fb.obs["caf_subtype_cluster"] = fb.obs["fb_leiden"].map(mapping).astype("category")

    counts = pd.concat({"per_cell": fb.obs["caf_subtype_cell"].value_counts(),
                        "per_subcluster": fb.obs["caf_subtype_cluster"].value_counts()}, axis=1).fillna(0).astype(int)
    counts["frac_per_cell"] = counts["per_cell"] / fb.n_obs
    return CafResult(fb, means, counts, sets, dropped, scored)


def resolve_caf_sets(cfg: CafConfig, inline: dict[str, list[str]] | None, base=None) -> dict[str, list[str]]:
    """Sets of ``markers.caf`` when given, else the files in ``caf.marker_files`` (relative
    paths resolve against ``base``; ``FileNotFoundError`` for a missing file). An empty
    mapping means the CAF stage was not opted in."""
    from ..markers import load_marker_files

    if inline:
        return {k: list(v) for k, v in inline.items()}
    if cfg.marker_files:
        return load_marker_files(cfg.marker_files, base=base)
    return {}
