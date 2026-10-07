"""Figure helpers. Every figure is written as PNG *and* PDF; tissue plots are drawn in
image orientation (y axis inverted, matching Xenium Explorer)."""
from __future__ import annotations

from pathlib import Path

import numpy as np


def save_fig(fig, stem, dpi: int = 150) -> None:
    """Save ``fig`` as ``<stem>.png`` and ``<stem>.pdf`` (scatter collections rasterised)."""
    import matplotlib.pyplot as plt

    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    for ax in fig.axes:
        for coll in ax.collections:
            coll.set_rasterized(True)
    fig.savefig(str(stem) + ".png", dpi=dpi, bbox_inches="tight")
    fig.savefig(str(stem) + ".pdf", dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def _invert_tissue_axes(fig) -> None:
    for ax in fig.axes:
        if ax.get_label() == "<colorbar>" or not ax.collections:
            continue
        if not ax.yaxis_inverted():
            ax.invert_yaxis()
        ax.set_aspect("equal")


def _legend_defaults(adata, color, kw):
    cols = [color] if isinstance(color, str) else list(color)
    n = max((adata.obs[c].nunique() for c in cols
             if c in adata.obs and adata.obs[c].dtype.name in ("category", "object")), default=0)
    kw.setdefault("legend_loc", "right margin")
    kw.setdefault("legend_fontsize", 9 if n <= 12 else 8 if n <= 24 else 7)
    return kw


def spatial_plot(adata, color, stem, size=1, title=None, figsize=(7.5, 12), **kw) -> None:
    """Tissue scatter of ``obsm['spatial']`` coloured by ``color`` → PNG + PDF."""
    import matplotlib.pyplot as plt
    import scanpy as sc

    kw = _legend_defaults(adata, color, kw)
    with plt.rc_context({"figure.figsize": figsize}):
        fig = sc.pl.embedding(adata, basis="spatial", color=color, size=size, show=False, return_fig=True,
                              title=title, **kw)
    _invert_tissue_axes(fig)
    save_fig(fig, stem)


def umap_plot(adata, color, stem, figsize=(6.5, 6.5), **kw) -> None:
    import matplotlib.pyplot as plt
    import scanpy as sc

    if "X_umap" not in adata.obsm:
        return
    ok = ~np.isnan(adata.obsm["X_umap"][:, 0])
    sub = adata[ok] if not ok.all() else adata
    kw = _legend_defaults(sub, color, kw)
    with plt.rc_context({"figure.figsize": figsize}):
        fig = sc.pl.umap(sub, color=color, show=False, return_fig=True, **kw)
    save_fig(fig, stem)


def grouped_bar(df, stem, ylabel="", title="", colors=None) -> None:
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(max(7, 0.55 * len(df) + 2), 5.5))
    df.plot.bar(ax=ax, width=0.8, **({"color": list(colors)} if colors else {}))
    ax.set_xticklabels(df.index, rotation=45, ha="right", rotation_mode="anchor", fontsize=9)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False)
    fig.tight_layout()
    save_fig(fig, stem)


def dotplot_save(adata, var_names, groupby, stem, **kw) -> None:
    import scanpy as sc

    var_names = {k: v for k, v in var_names.items() if v} if isinstance(var_names, dict) else var_names
    if not var_names:
        return
    dp = sc.pl.dotplot(adata, var_names, groupby=groupby, show=False, return_fig=True, **kw)
    dp.make_figure()
    save_fig(dp.fig, stem)


def score_heatmap(means, cols, stem, title="", vmin=-2, vmax=2) -> None:
    """Cluster × marker-set mean-z heatmap with the assigned label on the y axis."""
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(max(6, 0.45 * len(cols) + 3), max(4, 0.32 * len(means) + 1)))
    im = ax.imshow(means[cols].values, cmap="RdBu_r", vmin=vmin, vmax=vmax, aspect="auto")
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(cols, rotation=90)
    ax.set_yticks(range(len(means)))
    ax.set_yticklabels([f"{i} → {means.loc[i, 'assigned']} (n={int(means.loc[i, 'n_cells']):,})" for i in means.index])
    plt.colorbar(im, ax=ax, label="mean z-score")
    ax.set_title(title)
    fig.tight_layout()
    save_fig(fig, stem)


def violin_save(adata, keys, groupby, stem, **kw) -> None:
    import matplotlib.pyplot as plt
    import scanpy as sc

    sc.pl.violin(adata, keys, groupby=groupby, rotation=45, show=False, stripplot=False, **kw)
    save_fig(plt.gcf(), stem)


def qc_histograms(adata, stem, title, cols, colour="#4C72B0", ncols=3) -> None:
    import matplotlib.pyplot as plt

    n = len(cols)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows), squeeze=False)
    for ax, (col, logx) in zip(axes.ravel(), cols):
        if col not in adata.obs:
            ax.axis("off")
            continue
        v = adata.obs[col].to_numpy(dtype=float)
        if logx:
            v = np.log10(v + 1)
            ax.set_xlabel(f"log10({col}+1)")
        else:
            ax.set_xlabel(col)
        ax.hist(v, bins=100, color=colour)
        ax.set_ylabel("cells")
    for ax in axes.ravel()[n:]:
        ax.axis("off")
    fig.suptitle(title)
    fig.tight_layout()
    save_fig(fig, stem)


def qc_spatial_maps(adata, stem, cols=("gene_counts", "n_genes_by_counts", "cell_area"), max_cells=300_000, seed=0) -> None:
    import matplotlib.pyplot as plt

    sub = np.random.default_rng(seed).choice(adata.n_obs, min(max_cells, adata.n_obs), replace=False)
    xy = adata.obsm["spatial"][sub]
    fig, axes = plt.subplots(1, len(cols), figsize=(7 * len(cols), 7), squeeze=False)
    for ax, col in zip(axes.ravel(), cols):
        v = np.log10(adata.obs[col].to_numpy(dtype=float)[sub] + 1)
        s = ax.scatter(xy[:, 0], xy[:, 1], c=v, s=0.3, cmap="viridis", rasterized=True)
        ax.set_aspect("equal"); ax.invert_yaxis(); ax.set_title(f"log10 {col}"); ax.axis("off")
        plt.colorbar(s, ax=ax, shrink=0.6)
    fig.tight_layout()
    save_fig(fig, stem)


def histogram_with_threshold(values, thresh, stem, xlabel, ylabel, label) -> None:
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(np.asarray(values, dtype=float), bins=150, color="#4C72B0")
    ax.axvline(thresh, color="red", ls="--", label=label)
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel); ax.legend()
    fig.tight_layout()
    save_fig(fig, stem)


def nhood_heatmap(adata, key, stem, vmin=-50, vmax=50) -> None:
    import matplotlib.pyplot as plt
    import squidpy as sq

    sq.pl.nhood_enrichment(adata, cluster_key=key, method="average", cmap="RdBu_r", vmin=vmin, vmax=vmax,
                           figsize=(11, 10))
    save_fig(plt.gcf(), stem)


# ----------------------------------------------------------------------------
# Final annotation map: hierarchical legend and computed colours
# ----------------------------------------------------------------------------
# Legend groups of the final annotation map, in display order (group title -> labels in
# order). Colours are computed from this hierarchy: labels of one group are pushed as far
# apart as possible, so subtypes that sit together in the tissue can be told apart.
QC_GROUP = "Unresolved / QC"
OTHER_GROUP = "Other"
LEGEND_HIERARCHY: dict[str, list[str]] = {
    "Tumour epithelium": ["Tumor_nonCSC", "CSC_like", "Epithelial_tumor", "Epithelial_tumor_cycling"],
    "Fibroblast / CAF": ["iCAF", "myCAF", "apCAF", "FB_unassigned", "Fibroblast"],
    "Other stroma": ["Ovarian_stroma", "Endothelial", "Lymphatic_EC", "Pericyte_SMC", "Mesothelial", "Adipocyte", "Schwann"],
    "T cell": ["CD8_cytotoxic", "CD8_exhausted", "CD8_naive_memory", "CD8_T_unspecified", "CD8_T", "CD4_Th_effector",
               "Treg", "Tfh", "CD4_naive_memory", "CD4_T_unspecified", "CD4_T", "gdT_MAIT", "Proliferating_T", "T_NK"],
    "NK cell": ["NK"],
    "B lineage": ["Naive_B", "Memory_B", "GC_B", "B_cell", "Proliferating_B", "Plasmablast", "Plasma"],
    "Myeloid: macrophage / monocyte": ["M1_macrophage", "M2_macrophage", "TREM2_TAM", "Hypoxic_TAM",
                                       "Macrophage_unpolarized", "Macrophage", "Monocyte_classical",
                                       "Monocyte_nonclassical", "Monocyte_unspecified", "Monocyte",
                                       "Proliferating_myeloid", "Osteoclast_like_giant", "Macrophage_Mono"],
    "Myeloid: dendritic cell": ["mregDC", "cDC1", "cDC2", "cDC_unspecified", "cDC", "pDC", "Dendritic"],
    "Myeloid: granulocyte": ["Neutrophil", "Mast"],
    QC_GROUP: ["Unassigned", "Myeloid_spillover_Tumor", "Myeloid_spillover_Fibroblast", "Bcell_spillover_Fibroblast",
               "TNK_spillover_Tumor", "TNK_spillover_Fibroblast", "Bcell_spillover_Tumor", "Immune_mixed",
               "Myeloid_spillover_Endothelial", "TNK_spillover_Endothelial", "Bcell_spillover_Endothelial",
               "T_NK_unresolved", "B_lineage_unresolved", "Myeloid_unresolved"],
}
# Labels that name a lineage or identity without a specific fine type. They are drawn in
# deep tones; labels with a specific fine type are drawn bright; QC classes in greys.
GENERIC_LABELS: list[str] = [
    "Epithelial_tumor", "FB_unassigned", "Fibroblast",
    "CD8_T_unspecified", "CD8_T", "CD4_T_unspecified", "CD4_T", "T_NK", "B_cell",
    "Macrophage_unpolarized", "Macrophage", "Monocyte_unspecified", "Monocyte", "Macrophage_Mono",
    "cDC_unspecified", "cDC", "Dendritic",
]
GENERIC_SUFFIXES = ("_unspecified", "_unpolarized", "_unassigned")   # for labels the hierarchy does not list
# Colour-distance settings (OKLab distance x 100; about 15 is the floor for telling two
# colours apart at a glance). Inside a legend group colours are placed by farthest-point
# selection (each new colour is the candidate farthest from the group's colours so far).
# Candidates closer than a radius to a colour of another group are excluded first; the
# largest radius of CROSS_RADII is used that still leaves the group at least GROUP_MIN and
# KEEP_FRACTION of its best achievable separation.
GROUP_MIN = 15.0
QC_MIN = 3.0
KEEP_FRACTION = 0.85
CROSS_RADII = (12.0, 10.0, 8.0, 6.0, 4.0, 2.0)
ALSO_DISTINCT_FLOOR = 15.0               # between the members of an `also_distinct` set
QC_NEUTRAL = "#bdbdbd"
NO_REUSE_RADIUS = 0.5                    # minimum distance to any colour of another group
SEED_HINT = "#1f77b4"                    # the very first colour is the candidate nearest to this blue
# Three tiers (OKLCH lightness and chroma grids, hue step in degrees):
#   bright - labels with a specific fine type: light and saturated
#   deep   - lineage / identity labels without a specific fine type: dark, still clearly coloured
#   grey   - unassigned / spillover / unresolved classes: greys with next to no hue
TIER_POOLS = {
    "bright": ([0.58, 0.62, 0.66, 0.70, 0.74, 0.78, 0.82, 0.86, 0.90], [0.13, 0.16, 0.19, 0.22, 0.25, 0.28, 0.31], 10),
    "deep": ([0.30, 0.34, 0.38, 0.42, 0.46, 0.50], [0.07, 0.10, 0.13, 0.16, 0.19, 0.22], 10),
    "grey": ([0.55, 0.585, 0.62, 0.655, 0.69, 0.725, 0.76, 0.795, 0.83, 0.865, 0.90], [0.0, 0.008, 0.015], 60),
}

_M1 = np.array([[0.4122214708, 0.5363325363, 0.0514459929],
                [0.2119034982, 0.6806995451, 0.1073969566],
                [0.0883024619, 0.2817188376, 0.6299787005]])
_M2 = np.array([[0.2104542553, 0.7936177850, -0.0040720468],
                [1.9779984951, -2.4285922050, 0.4505937099],
                [0.0259040371, 0.7827717662, -0.8086757660]])


def hex_to_oklab(colors) -> np.ndarray:
    """OKLab coordinates (n × 3) of hex colours."""
    import matplotlib.colors as mc

    rgb = np.array([mc.to_rgb(c) for c in colors], dtype=float).reshape(-1, 3)
    lin = np.where(rgb <= 0.04045, rgb / 12.92, ((rgb + 0.055) / 1.055) ** 2.4)
    return np.cbrt(lin @ _M1.T) @ _M2.T


def _oklch_pool(lightness, chroma, hue_step) -> list[str]:
    """In-gamut sRGB hex colours on an OKLCH grid."""
    import matplotlib.colors as mc

    m1i, m2i = np.linalg.inv(_M1), np.linalg.inv(_M2)
    out = []
    for L in lightness:
        for C in chroma:
            for h in ([0] if C == 0 else range(0, 360, hue_step)):
                lab = np.array([L, C * np.cos(np.radians(h)), C * np.sin(np.radians(h))])
                lin = ((lab @ m2i.T) ** 3) @ m1i.T
                if (lin < -1e-4).any() or (lin > 1 + 1e-4).any():
                    continue
                lin = np.clip(lin, 0, 1)
                rgb = np.where(lin <= 0.0031308, 12.92 * lin, 1.055 * lin ** (1 / 2.4) - 0.055)
                out.append(mc.to_hex(np.clip(rgb, 0, 1)))
    return list(dict.fromkeys(out))


_POOLS: dict[str, tuple[list[str], np.ndarray]] = {}


def _pool(tier: str) -> tuple[list[str], np.ndarray]:
    if tier not in _POOLS:
        lightness, chroma, step = TIER_POOLS[tier]
        cols = _oklch_pool(lightness, chroma, step)
        _POOLS[tier] = (cols, hex_to_oklab(cols))
    return _POOLS[tier]


def is_qc_label(label: str) -> bool:
    """Labels that say 'not resolved' rather than name a cell type."""
    return label == "Unassigned" or "_spillover_" in label or label.endswith("_unresolved") or label == "Immune_mixed"


def label_tier(label: str, group: str, groups: dict[str, list[str]] | None = None, qc_group: str = QC_GROUP,
               generic=None) -> str:
    """``grey`` (QC class), ``deep`` (no specific fine type) or ``bright`` (specific fine type)."""
    groups = LEGEND_HIERARCHY if groups is None else groups
    if group == qc_group:
        return "grey"
    generic_set = set(GENERIC_LABELS if generic is None else generic)
    listed = any(label in ls for ls in groups.values())
    if label in generic_set or (not listed and label.endswith(GENERIC_SUFFIXES)):
        return "deep"
    return "bright"


def _group_assignment(labels, groups: dict[str, list[str]], qc_group: str) -> list[tuple[str, str]]:
    """``[(label, group)]``: every label of ``groups`` in hierarchy order, then the other
    given labels (sorted) under the QC group when they look like QC labels, else ``Other``."""
    order, seen = [], set()
    for g, members in groups.items():
        for lab in members:
            if lab not in seen:
                order.append((lab, g))
                seen.add(lab)
    for lab in sorted(set(map(str, labels)) - seen):
        order.append((lab, qc_group if (is_qc_label(lab) and qc_group in groups) else OTHER_GROUP))
    return order


def fine_palette(labels=(), groups: dict[str, list[str]] | None = None, pinned: dict[str, str] | None = None,
                 qc_group: str = QC_GROUP, also_distinct=(), generic=None) -> dict[str, str]:
    """Label -> hex colour for every label of ``groups`` plus ``labels``.

    Three tiers: labels of ``qc_group`` (unassigned, spillover, unresolved) are greys;
    ``generic`` labels (a lineage or identity without a specific fine type, default
    ``GENERIC_LABELS``) get deep tones; every other label, i.e. a specific fine type, gets
    a bright colour. Within a legend group all colours are placed mutually far apart in
    OKLab (at least ``GROUP_MIN`` whenever the group can be packed that way), and no colour
    of one group lies close to a colour of another. ``also_distinct`` lists further label
    sets whose members must be mutually distinct although they sit in different groups (the
    rough cell types that share one coarse figure). The result depends only on the
    arguments other than ``labels``' composition, so a label keeps its colour across figures
    and runs. ``pinned`` fixes the colour of given labels (e.g. a saved colour table); the
    other labels are placed around them."""
    import matplotlib.colors as mc

    groups = LEGEND_HIERARCHY if groups is None else groups
    extra_sets = [set(map(str, s_)) for s_ in also_distinct]
    order = _group_assignment(list(labels) + [l for s_ in extra_sets for l in sorted(s_)], groups, qc_group)
    grp = dict(order)
    pal: dict[str, str] = {lab: mc.to_hex(pinned[lab]) for lab, _ in order if pinned and lab in pinned}
    lab_cache: dict[str, np.ndarray] = {}

    def tier(lab: str) -> str:
        return label_tier(lab, grp[lab], groups, qc_group, generic)

    def coords(hexes) -> np.ndarray:
        for h in hexes:
            if h not in lab_cache:
                lab_cache[h] = hex_to_oklab([h])[0]
        return np.array([lab_cache[h] for h in hexes]).reshape(-1, 3)

    def dist(clab: np.ndarray, hexes) -> np.ndarray:
        if not len(hexes):
            return np.full(len(clab), np.inf)
        return np.linalg.norm(clab[:, None] - coords(hexes)[None], axis=2).min(1) * 100

    def spread(hexes) -> float:
        c = coords(hexes)
        if len(c) < 2:
            return float("inf")
        d = np.linalg.norm(c[:, None] - c[None], axis=2) * 100
        return float(d[np.triu_indices(len(c), 1)].min())

    def place(g: str, members: list[str], radius: float) -> tuple[dict[str, str] | None, int]:
        """Farthest-point placement of the unpinned members of group ``g``. Returns the colours
        and the number of labels whose ``also_distinct`` distance could not be honoured."""
        floor = QC_MIN if g == qc_group else GROUP_MIN
        elsewhere = [pal[x] for x in pal if grp[x] != g]
        placed: dict[str, str] = {}
        missed = 0
        for lab in members:
            if lab in pal:
                continue
            t = tier(lab)
            cands, clab = _pool(t)
            in_group = [pal[x] for x in pal if grp[x] == g] + list(placed.values())
            if t == "grey" and not in_group:
                placed[lab] = QC_NEUTRAL
                continue
            d_group, d_else = dist(clab, in_group), dist(clab, elsewhere)
            ok = d_else >= radius
            if not ok.any():
                if radius > NO_REUSE_RADIUS:
                    return None, 0                               # this radius leaves the label no colour at all
                ok = np.ones(len(cands), dtype=bool)
            partners = [pal[x] for x in pal if grp[x] != g and any(lab in s_ and x in s_ for s_ in extra_sets)]
            if partners:
                near = dist(clab, partners) >= ALSO_DISTINCT_FLOOR
                if (ok & near & (d_group >= floor)).any():       # never at the cost of the group's own separation
                    ok &= near
                else:
                    missed += 1
            if in_group:
                score = d_group
            elif elsewhere:
                score = d_else                                   # first colour of the group: far from the others
            else:
                score = -dist(clab, [SEED_HINT])                 # very first colour of all
            placed[lab] = cands[int(np.argmax(np.where(ok, score, -np.inf)))]
        return placed, missed

    by_group: dict[str, list[str]] = {}
    for lab, g in order:
        by_group.setdefault(g, []).append(lab)
    for g, members in by_group.items():
        fixed = [pal[x] for x in pal if grp[x] == g]
        best, missed0 = place(g, members, NO_REUSE_RADIUS)        # a colour is never reused in another group
        spread0 = spread(fixed + list(best.values()))
        floor = QC_MIN if g == qc_group else GROUP_MIN
        # (radius, separation the group must keep): first insist on the floor, then only on
        # (nearly) the best the group can do
        plan = [(r, max(floor, KEEP_FRACTION * min(spread0, 40.0))) for r in CROSS_RADII] + \
               [(r, 0.93 * min(spread0, 40.0)) for r in (4.0, 2.0)]
        if g == qc_group:
            plan = []        # greys are apart from every coloured tier anyway: spend the whole grey scale on them
        for radius, need in plan:
            trial, missed = place(g, members, radius)
            if trial is not None and missed <= missed0 and spread(fixed + list(trial.values())) >= need:
                best = trial
                break
        pal.update(best)
    return pal


def read_palette_file(path) -> dict[str, str]:
    """Colour table with columns ``label`` and ``hex`` (``.csv`` comma-separated, anything
    else tab-separated; other columns are ignored). A table written by this package
    (``results/cell_type_colors.csv``) can be read back to reproduce its colours."""
    import pandas as pd

    path = Path(path)
    df = pd.read_csv(path, sep="," if path.suffix.lower() == ".csv" else "\t", dtype=str)
    if not {"label", "hex"} <= set(df.columns):
        raise ValueError(f"{path}: a palette file needs the columns 'label' and 'hex'")
    df = df.dropna(subset=["label", "hex"])
    return dict(zip(df["label"], df["hex"]))


def apply_palette(adata, key: str, pal: dict[str, str]) -> None:
    """Store the colours of ``obs[key]``'s categories in ``uns[key + '_colors']`` so that
    scanpy plots (and anything else reading the h5ad) use them."""
    import pandas as pd

    if not isinstance(adata.obs[key].dtype, pd.CategoricalDtype):
        adata.obs[key] = adata.obs[key].astype("category")
    adata.uns[f"{key}_colors"] = [pal.get(str(c), "#000000") for c in adata.obs[key].cat.categories]


def scheme_table(pal: dict[str, str], groups=None, qc_group: str = QC_GROUP, generic=None):
    """The whole colour scheme, one row per label of ``pal``: ``label``, ``hex``, ``tier``
    (grey / deep / bright), ``legend_group``. Portable: save it and point ``plot.palette_file``
    of any other dataset at it."""
    import pandas as pd

    rows = [{"label": lab, "hex": pal[lab], "tier": label_tier(lab, group, groups, qc_group, generic),
             "legend_group": group}
            for group, labs in legend_layout(list(pal), groups, qc_group) for lab in labs]
    return pd.DataFrame(rows, columns=["label", "hex", "tier", "legend_group"])


def color_table(pal: dict[str, str], coarse=None, fine=None, groups=None, qc_group: str = QC_GROUP, generic=None):
    """One row per cell-type label with its colour, for reuse and reproduction.

    ``coarse`` / ``fine`` are per-cell label vectors (either may be None). Columns:
    ``label``, ``hex``, ``level`` (``fine``, ``coarse`` or ``fine+coarse``), ``tier`` (grey /
    deep / bright), ``legend_group``
    (the big group of the hierarchical legend), ``coarse_type`` (for a fine label: the rough
    type most of its cells carry), ``n_cells_fine``, ``n_cells_coarse``."""
    import pandas as pd

    fine_s = None if fine is None else pd.Series(np.asarray(fine).astype(str))
    coarse_s = None if coarse is None else pd.Series(np.asarray(coarse).astype(str))
    n_fine = fine_s.value_counts() if fine_s is not None else pd.Series(dtype=int)
    n_coarse = coarse_s.value_counts() if coarse_s is not None else pd.Series(dtype=int)
    parent = {}
    if fine_s is not None and coarse_s is not None and len(fine_s) == len(coarse_s):
        ct = pd.crosstab(fine_s, coarse_s)
        parent = ct.idxmax(axis=1).to_dict()
    labels = list(n_fine.index) + [l for l in n_coarse.index if l not in n_fine.index]
    layout = legend_layout(labels, groups, qc_group)
    rows = []
    for group, labs in layout:
        for lab in labs:
            is_f, is_c = lab in n_fine.index, lab in n_coarse.index
            rows.append({"label": lab, "hex": pal.get(lab, "#000000"),
                         "level": "fine+coarse" if is_f and is_c else "fine" if is_f else "coarse",
                         "tier": label_tier(lab, group, groups, qc_group, generic),
                         "legend_group": group,
                         "coarse_type": parent.get(lab, lab if is_c else ""),
                         "n_cells_fine": int(n_fine.get(lab, 0)), "n_cells_coarse": int(n_coarse.get(lab, 0))})
    return pd.DataFrame(rows, columns=["label", "hex", "level", "tier", "legend_group", "coarse_type",
                                       "n_cells_fine", "n_cells_coarse"])


def legend_layout(labels, groups: dict[str, list[str]] | None = None,
                  qc_group: str = QC_GROUP) -> list[tuple[str, list[str]]]:
    """``[(group title, [labels])]`` for the labels that are present: groups and labels in
    the order of ``groups``; unlisted QC-like labels join the QC group, anything else goes
    under ``Other`` (sorted)."""
    groups = LEGEND_HIERARCHY if groups is None else groups
    present = set(map(str, labels))
    by_group: dict[str, list[str]] = {}
    for lab, g in _group_assignment(present, groups, qc_group):
        if lab in present:
            by_group.setdefault(g, []).append(lab)
    titles = [g for g in groups if g in by_group] + ([OTHER_GROUP] if OTHER_GROUP in by_group and OTHER_GROUP not in groups else [])
    return [(g, by_group[g]) for g in titles]


def build_fine_map(xy, labels, pal: dict[str, str], groups: dict[str, list[str]] | None = None,
                   title: str = "Final annotation (tissue)", size: float = 0.6, figsize=(15, 12),
                   legend_fontsize: float = 9, marker_size: float = 9, qc_group: str = QC_GROUP):
    """Tissue map of the fine annotation in image orientation with a single-column
    hierarchical legend: bold group headers, each followed by its labels with cell counts.
    Abundant labels are drawn first so rare ones are not hidden. Returns the figure."""
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    xy = np.asarray(xy)
    labels = np.asarray(labels).astype(str)
    layout = legend_layout(np.unique(labels), groups, qc_group)
    counts = {l: int((labels == l).sum()) for _, ls in layout for l in ls}
    fig, ax = plt.subplots(figsize=figsize)
    for l in sorted(counts, key=lambda l: -counts[l]):
        m = labels == l
        ax.scatter(xy[m, 0], xy[m, 1], s=size, c=pal.get(l, "#000000"), rasterized=True, linewidths=0)
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.axis("off")
    ax.set_title(title, fontsize=13)
    handles, texts, is_head = [], [], []
    for head, ls in layout:
        handles.append(Patch(facecolor="none", edgecolor="none"))
        texts.append(head)
        is_head.append(True)
        for l in ls:
            handles.append(Line2D([0], [0], marker="o", color="none", markerfacecolor=pal.get(l, "#000000"),
                                  markeredgecolor="none", markersize=marker_size))
            texts.append(f"{l}  ({counts[l]:,})")
            is_head.append(False)
    leg = ax.legend(handles, texts, loc="center left", bbox_to_anchor=(1.0, 0.5), ncol=1, frameon=False,
                    fontsize=legend_fontsize, handlelength=1.2, handletextpad=0.6, labelspacing=0.25,
                    borderaxespad=0.5)
    for t, h in zip(leg.get_texts(), is_head):
        if h:
            t.set_fontweight("bold")
            t.set_fontsize(legend_fontsize + 1)
    return fig
