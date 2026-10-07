"""Step 4 — immune compartments: identity (level 1) → functional state (level 2).

Compartments are processed so that those with a ``handoff`` run first (T/NK hands B and
plasma sub-clusters to the B compartment, where they keep that identity).

Level 1 (per sub-cluster, mean z within the compartment): core identities qualify at
mean z ≥ ``low_identity``, minority and ``cross_*`` identities at > ``min_z``; the highest
qualifying set wins. Otherwise a spillover set > ``min_z`` gives
``<compartment>_spillover_<X>``; otherwise ``unresolved``. ``immune_spillover`` is set on
any sub-cluster whose best spillover z > ``spill_flag``.

Level 2 (per state sub-cluster within a parent identity, re-clustered at ``resolution2``,
z within the parent): highest mean z > ``min_z2`` else the parent's ``unspecified`` label.
``cross_*`` labels map to the other compartment's coarse label.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from ..config import CompartmentConfig, ImmuneConfig, MarkerConfig
from ..markers import filter_to_panel
from ..scoring import cluster_assign, embed, score_sets


@dataclass
class ImmuneResult:
    name: str
    adata: object
    means1: pd.DataFrame                              # sub-cluster × set mean z + assigned/spillover_flag/n_cells
    level2_tables: dict[str, pd.DataFrame]            # parent -> state sub-cluster table
    counts: pd.DataFrame                              # immune_subtype counts
    handoff: dict[str, dict[str, str]]                # target compartment -> {cell: level-1 label}
    sets_used: dict                                   # {"level1": {...}, "level2": {parent: {...}}}
    top_markers: pd.DataFrame | None = None
    n_prior_kept: int = 0
    scored1: list[str] = field(default_factory=list)


def assign_level1(means1: pd.DataFrame, identities: list[str], spills: list[str], comp: CompartmentConfig,
                  cfg: ImmuneConfig, name: str) -> tuple[dict, dict]:
    """Pure decision rule on the sub-cluster × set mean-z table. Returns
    ``(level1: {cluster: label}, flag: {cluster: bool})``."""
    level1, flag = {}, {}
    for cl, row in means1.iterrows():
        # both floors are inclusive, exactly as in the source script
        cand = {i: row[i] for i in identities
                if row[i] >= (cfg.low_identity if i in comp.core else cfg.min_z)}
        best_sp = row[spills].idxmax() if spills else None
        if cand:
            level1[cl] = max(cand, key=cand.get)
        elif best_sp is not None and row[best_sp] > cfg.min_z:
            level1[cl] = f"{name}_{best_sp}"
        else:
            level1[cl] = comp.unresolved
        flag[cl] = bool(best_sp is not None and row[best_sp] > cfg.spill_flag)
    return level1, flag


def annotate_compartment(sub, name: str, comp: CompartmentConfig, cfg: ImmuneConfig,
                         spillover_sets: dict[str, list[str]], seed: int = 0, prior: dict | None = None,
                         umap: bool = True, compute_top_markers: bool = True) -> ImmuneResult | None:
    """Annotate one compartment. ``sub`` must already contain the compartment's cells (plus any
    handed-off cells, whose upstream identity is given in ``prior`` as ``{cell: label}``).
    Returns None when ``sub.n_obs < cfg.min_cells``. Cells handed off *from* this compartment
    are removed from the returned object and listed in ``result.handoff``."""
    import scanpy as sc

    if sub.n_obs < cfg.min_cells:
        return None
    key = f"{name}_leiden"
    embed(sub, cfg.n_pcs, cfg.n_neighbors, cfg.resolution, seed, key, umap=umap, umap_cells=cfg.umap_cells)

    # ---- level 1 -----------------------------------------------------------
    l1_all = {**comp.level1, **spillover_sets}
    l1_sets, l1_dropped = filter_to_panel(l1_all, sub.var_names)
    pre1 = f"imm1_{name}_"
    scored1 = score_sets(sub, l1_sets, prefix=pre1, seed=seed)
    identities = [s for s in scored1 if not s.startswith("spillover_")]
    spills = [s for s in scored1 if s.startswith("spillover_")]
    if not identities:
        raise ValueError(f"[{name}] no identity marker set has genes on this panel")
    means1 = sub.obs.groupby(key, observed=True)[[f"{pre1}{s}_z" for s in scored1]].mean()
    means1.columns = scored1
    level1, flag = assign_level1(means1, identities, spills, comp, cfg, name)
    sub.obs["immune_level1"] = sub.obs[key].map(level1).astype(str)
    n_prior = 0
    if prior:
        pr = pd.Series(prior).reindex(sub.obs_names)
        has = pr.notna()
        n_prior = int(has.sum())
        sub.obs.loc[has, "immune_level1"] = pr[has].values
    sub.obs["immune_spillover"] = sub.obs[key].map(flag).astype(bool)
    means1["assigned"] = pd.Series(level1)
    means1["spillover_flag"] = pd.Series(flag)
    means1["n_cells"] = sub.obs[key].value_counts()

    # ---- hand-off ------------------------------------------------------------
    handoff: dict[str, dict[str, str]] = {}
    for lab, target in comp.handoff.items():
        ids = sub.obs_names[(sub.obs["immune_level1"] == lab).to_numpy()]
        if len(ids):
            handoff.setdefault(target, {}).update({i: lab for i in ids})
    means1["handed_off"] = means1["assigned"].isin(list(comp.handoff)) if comp.handoff else False
    if handoff:
        moved = np.concatenate([np.array(list(v)) for v in handoff.values()])
        sub = sub[~sub.obs_names.isin(moved)].copy()

    # ---- level 2 -------------------------------------------------------------
    final = sub.obs["immune_level1"].astype(str).copy()
    state_cell = pd.Series("NA", index=sub.obs_names, dtype=object)
    l2_used, l2_tables = {}, {}
    for parent, states in comp.level2.items():
        pm = (sub.obs["immune_level1"] == parent).to_numpy()
        if pm.sum() < cfg.min_cells_level2:
            continue
        par = sub[pm].copy()
        key2 = f"{name}_{parent}_leiden"
        embed(par, cfg.n_pcs, cfg.n_neighbors, cfg.resolution2, seed, key2, umap=False)
        sub.obs[key2] = "NA"
        sub.obs.loc[par.obs_names, key2] = par.obs[key2].astype(str).values
        st_sets, st_dropped = filter_to_panel(states, par.var_names)
        pre2 = f"imm2_{name}_{parent}_"
        scored2 = score_sets(par, st_sets, prefix=pre2, seed=seed)
        l2_used[parent] = {"used_on_panel": st_sets, "dropped_not_on_panel": st_dropped}
        if not scored2:
            continue
        default = comp.unspecified.get(parent, parent)
        m2, map2 = cluster_assign(par, key2, scored2, prefix=pre2, min_z=cfg.min_z2, unassigned=default)
        m2["assigned"] = pd.Series(map2)
        m2["n_cells"] = par.obs[key2].value_counts()
        m2.insert(0, "parent", parent)
        l2_tables[parent] = m2
        final.loc[par.obs_names] = par.obs[key2].map(map2).astype(str).values
        # per-cell state argmax (review only; the label is the sub-cluster call). As in the
        # source script the per-cell floor is cfg.min_z, not min_z2.
        Z = par.obs[[f"{pre2}{s}_z" for s in scored2]].to_numpy()
        best = Z.argmax(1)
        cc = np.array(scored2, dtype=object)[best]
        cc[Z[np.arange(len(Z)), best] <= cfg.min_z] = default
        state_cell.loc[par.obs_names] = cc
        for s in scored2:
            col = f"{pre2}{s}_z"
            sub.obs[col] = np.nan
            sub.obs.loc[par.obs_names, col] = par.obs[col].values

    for cross, lab in comp.cross.items():
        final[final == cross] = lab
    sub.obs["immune_subtype"] = pd.Categorical(final)
    sub.obs["immune_state_cell"] = pd.Categorical(state_cell)

    # ---- tables --------------------------------------------------------------
    means1["final_label_majority"] = sub.obs.groupby(key, observed=True)["immune_subtype"].agg(
        lambda x: x.value_counts().index[0])
    counts = sub.obs["immune_subtype"].value_counts().to_frame("n_cells")
    counts["fraction"] = counts["n_cells"] / sub.n_obs
    counts["level1"] = [sub.obs.loc[sub.obs["immune_subtype"] == i, "immune_level1"].mode().iat[0] for i in counts.index]
    counts["n_spillover_flagged"] = sub.obs.groupby("immune_subtype", observed=True)["immune_spillover"].sum()

    top = None
    if compute_top_markers and sub.obs[key].nunique() > 1:
        sc.tl.rank_genes_groups(sub, key, method="wilcoxon", n_genes=15, use_raw=False)
        top = pd.DataFrame([{"subcluster": cl, "level1": level1[cl], "final": means1.loc[cl, "final_label_majority"],
                             "top_genes": ",".join(map(str, sc.get.rank_genes_groups_df(sub, group=cl).head(15)["names"]))}
                            for cl in sub.obs[key].cat.categories])
    return ImmuneResult(name, sub, means1, l2_tables, counts, handoff,
                        {"level1": {"used_on_panel": l1_sets, "dropped_not_on_panel": l1_dropped}, "level2": l2_used},
                        top, n_prior, scored1)


def compartment_order(compartments: dict[str, CompartmentConfig], selected=None) -> list[str]:
    names = list(selected) if selected else list(compartments)
    return [c for c in names if compartments[c].handoff] + [c for c in names if not compartments[c].handoff]


def run_immune(full, cfg: ImmuneConfig, markers: MarkerConfig, seed: int = 0, compartments=None,
               umap: bool = True, compute_top_markers: bool = True, log=None) -> dict[str, ImmuneResult]:
    """Run every (selected) compartment on the lineage-annotated object ``full`` in hand-off order."""
    comps = markers.compartments
    selected = compartments or cfg.compartments or list(comps)
    results: dict[str, ImmuneResult] = {}
    pending: dict[str, dict[str, str]] = {}
    for name in compartment_order(comps, selected):
        comp = comps[name]
        recv = pending.pop(name, None)
        mask = full.obs["lineage"].isin(comp.lineages).values
        if recv:
            mask = mask | full.obs_names.isin(list(recv))
        sub = full[mask].copy()
        if log:
            log.info("[%s] %d cells from lineages %s (+%d handed off)", name, sub.n_obs, comp.lineages,
                     0 if not recv else len(recv))
        res = annotate_compartment(sub, name, comp, cfg, markers.spillover, seed=seed, prior=recv,
                                   umap=umap, compute_top_markers=compute_top_markers)
        if res is None:
            if log:
                log.warning("[%s] fewer than %d cells, skipped", name, cfg.min_cells)
            continue
        for target, ids in res.handoff.items():
            pending.setdefault(target, {}).update(ids)
            if log:
                log.info("[%s] handing off %d cells to %s", name, len(ids), target)
        results[name] = res
    return results
