"""Final annotation map: hierarchical legend, and colours that are distinct inside a group."""
import itertools
import textwrap

import matplotlib.colors as mc
import numpy as np
import pandas as pd
import pytest

from annotation_st.cli import main
from annotation_st.config import ConfigError, default_config, dump_config, load_config
from annotation_st.pipeline import label_palette
from annotation_st.plotting import (LEGEND_HIERARCHY, QC_GROUP, apply_palette, build_fine_map, color_table,
                                    fine_palette, legend_layout, read_palette_file)
from test_cli import _small_config
from conftest import ov_config

# the 43 fine labels of the OV validation run
OV_LABELS = ["Tumor_nonCSC", "Unassigned", "Ovarian_stroma", "myCAF", "Endothelial", "apCAF", "FB_unassigned",
             "Myeloid_spillover_Tumor", "iCAF", "CSC_like", "Macrophage_unpolarized", "NK", "TREM2_TAM", "Plasma",
             "CD4_T_unspecified", "Proliferating_myeloid", "M1_macrophage", "M2_macrophage",
             "Myeloid_spillover_Fibroblast", "CD8_T_unspecified", "Hypoxic_TAM", "Macrophage_Mono", "CD8_naive_memory",
             "Treg", "Bcell_spillover_Fibroblast", "CD8_cytotoxic", "Plasmablast", "Monocyte_unspecified",
             "Proliferating_T", "CD4_Th_effector", "mregDC", "CD8_exhausted", "Proliferating_B", "Monocyte_classical",
             "Memory_B", "T_NK", "Naive_B", "B_cell", "GC_B", "Osteoclast_like_giant", "Mesothelial", "pDC", "cDC1"]


def _oklab(hex_colour):
    """OKLab from the published matrices (Ottosson 2020), written out independently of the package."""
    r, g, b = [(c / 12.92) if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4 for c in mc.to_rgb(hex_colour)]
    l = np.cbrt(0.4122214708 * r + 0.5363325363 * g + 0.0514459929 * b)
    m = np.cbrt(0.2119034982 * r + 0.6806995451 * g + 0.1073969566 * b)
    s = np.cbrt(0.0883024619 * r + 0.2817188376 * g + 0.6299787005 * b)
    return np.array([0.2104542553 * l + 0.7936177850 * m - 0.0040720468 * s,
                     1.9779984951 * l - 2.4285922050 * m + 0.4505937099 * s,
                     0.0259040371 * l + 0.7827717662 * m - 0.8086757660 * s])


def _delta_e(a, b):
    return float(np.linalg.norm(_oklab(a) - _oklab(b)) * 100)


def _closest_pair(pal, labels):
    return min(((_delta_e(pal[a], pal[b]), a, b) for a, b in itertools.combinations(labels, 2)), default=(np.inf, "", ""))


def _run_palette(labels=()):
    """The palette a run with the default config uses (rough and fine types together)."""
    return label_palette(ov_config(), labels)


def test_cell_types_of_one_group_are_clearly_distinguishable():
    """Colour distance >= 15 (OKLab x 100, the floor for telling two colours apart) inside every cell-type group."""
    pal = _run_palette()
    for group, labels in LEGEND_HIERARCHY.items():
        if group == QC_GROUP:
            continue
        d, a, b = _closest_pair(pal, labels)
        assert d >= 15, f"{group}: {a} {pal[a]} vs {b} {pal[b]} are only {d:.1f} apart"


def test_rough_types_are_clearly_distinguishable_from_each_other():
    cfg = ov_config()
    rough = list(cfg.markers.lineage) + [cfg.lineage.cycling_label]
    d, a, b = _closest_pair(_run_palette(), rough)
    assert d >= 14.99, f"{a} vs {b} are only {d:.1f} apart"


def _chroma(hex_colour):
    return float(np.hypot(*_oklab(hex_colour)[1:]))


def _lightness(hex_colour):
    return float(_oklab(hex_colour)[0])


def _tiers():
    cfg = ov_config()
    qc = LEGEND_HIERARCHY[QC_GROUP]
    generic = [l for l in cfg.plot.generic_labels]
    cell_types = [l for g, ls in LEGEND_HIERARCHY.items() if g != QC_GROUP for l in ls]
    specific = [l for l in cell_types if l not in generic]
    return qc, generic, specific


def test_unassigned_and_spillover_labels_are_greys():
    pal = _run_palette()
    qc, _, _ = _tiers()
    assert pal["Unassigned"] == "#bdbdbd"
    assert max(_chroma(pal[l]) for l in qc) <= 0.025, {l: pal[l] for l in qc}      # grey family: next to no hue
    assert _closest_pair(pal, qc)[0] >= 3                                           # still separate keys
    # the classes seen in practice come first in the list and get greys that are easy to tell apart
    assert _closest_pair(pal, qc[:4])[0] >= 10, {l: pal[l] for l in qc[:4]}


def test_specific_fine_types_are_bright():
    pal = _run_palette()
    _, _, specific = _tiers()
    dull = {l: pal[l] for l in specific if _lightness(pal[l]) < 0.57 or _chroma(pal[l]) < 0.125}
    assert not dull, dull
    for l in ("Treg", "CSC_like", "iCAF", "myCAF", "M1_macrophage", "TREM2_TAM", "NK", "Plasma", "Endothelial"):
        assert l in specific


def test_labels_without_a_specific_fine_type_are_coloured_but_not_bright():
    pal = _run_palette()
    _, generic, _ = _tiers()
    for l in ("CD8_T_unspecified", "CD4_T_unspecified", "Macrophage_unpolarized", "Monocyte_unspecified",
              "FB_unassigned", "T_NK", "Macrophage_Mono"):
        assert l in generic
    odd = {l: pal[l] for l in generic if _lightness(pal[l]) > 0.52 or _chroma(pal[l]) < 0.06}
    assert not odd, odd
    # unknown labels follow the same rule by their name
    extra = _run_palette(["NewType_unspecified", "NewType_activated"])
    assert _lightness(extra["NewType_unspecified"]) <= 0.52 and _lightness(extra["NewType_activated"]) >= 0.57


def test_no_cell_type_can_pass_for_a_qc_class():
    pal = _run_palette()
    qc, generic, specific = _tiers()
    closest = min((_delta_e(pal[a], pal[b]), a, b) for a in generic + specific for b in qc)
    assert closest[0] >= 6, closest


def test_no_colour_is_used_twice_anywhere():
    pal = _run_palette(OV_LABELS)
    assert len(set(pal.values())) == len(pal)


def test_colour_of_a_label_does_not_depend_on_which_other_labels_are_present():
    full = _run_palette(OV_LABELS)
    assert _run_palette(["Treg", "NK"])["Treg"] == full["Treg"]
    assert _run_palette()["TREM2_TAM"] == full["TREM2_TAM"]


def test_unknown_labels_get_their_own_distinct_colours():
    pal = _run_palette(["Mystery_A", "Mystery_B", "Mystery_C"])
    assert _closest_pair(pal, ["Mystery_A", "Mystery_B", "Mystery_C"])[0] >= 15


def test_colour_table_lists_fine_and_rough_types_with_their_hierarchy():
    fine = ["Treg", "Treg", "NK", "CSC_like", "Tumor_nonCSC", "Plasma", "Unassigned"]
    coarse = ["T_NK", "T_NK", "T_NK", "Epithelial_tumor", "Epithelial_tumor_cycling", "Plasma", "Unassigned"]
    pal = {"Treg": "#111111", "NK": "#222222", "CSC_like": "#333333", "Tumor_nonCSC": "#444444", "Plasma": "#555555",
           "Unassigned": "#bdbdbd", "T_NK": "#666666", "Epithelial_tumor": "#777777", "Epithelial_tumor_cycling": "#888888"}
    t = color_table(pal, coarse=coarse, fine=fine).set_index("label")
    assert list(t.columns) == ["hex", "level", "tier", "legend_group", "coarse_type", "n_cells_fine", "n_cells_coarse"]
    assert t.loc["Treg"].tolist() == ["#111111", "fine", "bright", "T cell", "T_NK", 2, 0]
    assert t.loc["T_NK"].tolist() == ["#666666", "coarse", "deep", "T cell", "T_NK", 0, 3]
    assert t.loc["Plasma"].tolist() == ["#555555", "fine+coarse", "bright", "B lineage", "Plasma", 1, 1]
    assert t.loc["Unassigned", "tier"] == "grey"
    assert t.loc["CSC_like", "coarse_type"] == "Epithelial_tumor"
    assert t.loc["Unassigned", "legend_group"] == "Unresolved / QC"
    assert len(t) == 9
    # rough types only (lineage stage, before any refinement)
    t2 = color_table(pal, coarse=coarse)
    assert set(t2["level"]) == {"coarse"} and set(t2["label"]) == set(coarse)


def test_saved_colour_table_reproduces_the_palette(tmp_path):
    pal = _run_palette(OV_LABELS)
    table = color_table(pal, fine=OV_LABELS)
    table.loc[table["label"] == "Treg", "hex"] = "#ff0000"               # a hand-edited colour is respected too
    table.to_csv(tmp_path / "cell_type_colors.csv", index=False)
    cfg = ov_config()
    cfg.plot.palette_file = str(tmp_path / "cell_type_colors.csv")
    again = label_palette(cfg, OV_LABELS)
    assert again["Treg"] == "#ff0000"
    assert {l: again[l] for l in OV_LABELS if l != "Treg"} == {l: pal[l] for l in OV_LABELS if l != "Treg"}
    # the project's tab-separated table (label, hex, family) is readable as well
    (tmp_path / "old.tsv").write_text("label\thex\tfamily\nCSC_like\t#ff7f0e\tTumour\n")
    assert read_palette_file(tmp_path / "old.tsv") == {"CSC_like": "#ff7f0e"}


def test_apply_palette_stores_colours_in_category_order():
    import anndata as ad
    a = ad.AnnData(X=np.zeros((3, 1), dtype=np.float32), obs=pd.DataFrame({"cell_type": ["T_NK", "Plasma", "T_NK"]}, index=list("abc")))
    apply_palette(a, "cell_type", {"T_NK": "#010101", "Plasma": "#020202"})
    assert list(a.obs["cell_type"].cat.categories) == ["Plasma", "T_NK"]
    assert list(a.uns["cell_type_colors"]) == ["#020202", "#010101"]


def test_pinned_colours_are_kept_and_neighbours_move_away():
    pal = fine_palette(pinned={"Treg": "#ff0000", "Tumor_nonCSC": "#1F77B4"})
    assert pal["Treg"] == "#ff0000" and pal["Tumor_nonCSC"] == "#1f77b4"
    d, a, b = _closest_pair(pal, LEGEND_HIERARCHY["T cell"])
    assert d >= 13, (d, a, b)                                            # still distinct around the pinned red


def test_legend_layout_groups_present_labels_in_hierarchy_order():
    layout = legend_layout(["Treg", "CSC_like", "Tumor_nonCSC", "NK", "Unassigned", "Mystery", "iCAF", "pDC",
                            "Stroma_spillover_Tumor"])
    assert layout == [("Tumour epithelium", ["Tumor_nonCSC", "CSC_like"]),
                      ("Fibroblast / CAF", ["iCAF"]),
                      ("T cell", ["Treg"]),
                      ("NK cell", ["NK"]),
                      ("Myeloid: dendritic cell", ["pDC"]),
                      ("Unresolved / QC", ["Unassigned", "Stroma_spillover_Tumor"]),   # unlisted QC-like label
                      ("Other", ["Mystery"])]


def test_every_ov_label_has_a_legend_group():
    layout = legend_layout(OV_LABELS)
    assert "Other" not in [g for g, _ in layout]
    assert sorted(itertools.chain.from_iterable(ls for _, ls in layout)) == sorted(OV_LABELS)


def test_drawn_legend_has_bold_group_headers_followed_by_their_labels_with_counts():
    labels = np.array(["Tumor_nonCSC"] * 5 + ["CSC_like"] * 2 + ["Treg"] * 3 + ["NK"])
    xy = np.random.default_rng(0).normal(size=(len(labels), 2))
    pal = {"Tumor_nonCSC": "#112233", "CSC_like": "#ff0000", "Treg": "#00aa00", "NK": "#0000ff"}
    fig = build_fine_map(xy, labels, pal)
    ax = fig.axes[0]
    texts = ax.get_legend().get_texts()
    assert [t.get_text() for t in texts] == ["Tumour epithelium", "Tumor_nonCSC  (5)", "CSC_like  (2)",
                                             "T cell", "Treg  (3)", "NK cell", "NK  (1)"]
    assert [t.get_fontweight() == "bold" for t in texts] == [True, False, False, True, False, True, False]
    assert ax.yaxis_inverted()                                           # image orientation
    drawn = {mc.to_hex(c.get_facecolor()[0][:3]) for c in ax.collections}
    assert drawn == set(pal.values())                                    # one scatter layer per label, in its colour
    keys = [mc.to_hex(h.get_markerfacecolor()) for h in ax.get_legend().legend_handles if hasattr(h, "get_markerfacecolor")]
    assert keys == ["#112233", "#ff0000", "#00aa00", "#0000ff"]          # legend keys match their labels


def test_generic_labels_are_configurable(tmp_path):
    p = tmp_path / "c.yaml"
    p.write_text("plot: {generic_labels: []}\n")                      # nothing is generic: every cell type bright
    cfg = load_config(p)
    pal = label_palette(cfg)
    assert _lightness(pal["CD8_T_unspecified"]) >= 0.57
    p.write_text("plot: {generic_labels: Treg}\n")
    with pytest.raises(ConfigError, match="plot.generic_labels"):
        load_config(p)


def test_plot_section_is_configurable_and_validated(tmp_path):
    p = tmp_path / "c.yaml"
    p.write_text(textwrap.dedent("""
        plot:
          legend_groups: {Immune: [Treg, NK], Tumour: [Tumor_nonCSC]}
          colors: {NK: "#123456"}
    """))
    cfg = load_config(p)
    assert legend_layout(["NK", "Tumor_nonCSC", "Treg"], cfg.plot.legend_groups) == \
        [("Immune", ["Treg", "NK"]), ("Tumour", ["Tumor_nonCSC"])]
    pal = fine_palette(["NK", "Treg"], cfg.plot.legend_groups, cfg.plot.colors)
    assert pal["NK"] == "#123456" and _delta_e(pal["Treg"], pal["NK"]) >= 15
    p.write_text("plot: {legend_groups: {Immune: Treg}}\n")
    with pytest.raises(ConfigError, match="plot.legend_groups"):
        load_config(p)
    p.write_text("plot: {colors: {NK: blueish}}\n")
    with pytest.raises(ConfigError, match="plot.colors"):
        load_config(p)


def test_figures_command_redraws_the_map_without_rerunning_stages(tmp_path, synth_raw):
    cfgp = _small_config(tmp_path, synth_raw)                           # figures off
    assert main(["run", "-c", str(cfgp), "--stages", "qc", "lineage", "csc", "merge"]) == 0
    run = tmp_path / "run"
    assert not (run / "figures" / "spatial_cell_type_fine.png").exists()
    before = (run / "stage_state.json").read_text()

    assert main(["figures", "-c", str(cfgp)]) == 0
    assert (run / "figures" / "spatial_cell_type_fine.png").exists() and (run / "figures" / "spatial_cell_type_fine.pdf").exists()
    assert (run / "stage_state.json").read_text() == before            # nothing was re-run
    colours = pd.read_csv(run / "results" / "cell_type_colors.csv")
    fine = pd.read_csv(run / "results" / "final_cell_type_counts.tsv", sep="\t", index_col=0)
    coarse = pd.read_csv(run / "results" / "final_cell_type_coarse_counts.tsv", sep="\t", index_col=0)
    assert list(colours.columns) == ["label", "hex", "level", "tier", "legend_group", "coarse_type", "n_cells_fine", "n_cells_coarse"]
    # the rough-type and CSC maps come back too, in the same colours; stages that never ran are skipped quietly
    assert (run / "figures" / "spatial_cell_type.png").exists() and (run / "figures" / "spatial_tumor_csc.png").exists()
    assert not (run / "figures" / "spatial_fb_caf_subtype.png").exists()
    assert set(colours.loc[colours["level"].str.contains("fine"), "label"]) == set(fine.index)
    assert set(colours.loc[colours["level"].str.contains("coarse"), "label"]) == set(coarse.index)
    row = colours.set_index("label").loc["CSC_like"]
    assert row["legend_group"] == "Tumour epithelium" and row["coarse_type"].startswith("Epithelial_tumor")
    assert row["n_cells_fine"] == fine.loc["CSC_like", "n_cells"]

    # an accepted colour table can be pinned; the redraw then uses it
    (tmp_path / "accepted.tsv").write_text("label\thex\tfamily\nCSC_like\t#ff7f0e\tTumour\n")
    cfg = load_config(cfgp)
    cfg.plot.palette_file = "accepted.tsv"                              # relative to the config file
    dump_config(cfg, cfgp)
    assert main(["figures", "-c", str(cfgp)]) == 0
    colours = pd.read_csv(run / "results" / "cell_type_colors.csv").set_index("label")
    assert colours.loc["CSC_like", "hex"] == "#ff7f0e"


def test_figures_command_needs_the_merged_checkpoint(tmp_path, synth_raw):
    from annotation_st.pipeline import StageInputMissing
    cfgp = _small_config(tmp_path, synth_raw)
    with pytest.raises(StageInputMissing, match="merge"):
        main(["figures", "-c", str(cfgp)])


# ---- the colour scheme on other datasets ---------------------------------------
def _other_dataset(path):
    """A dataset annotated elsewhere: some labels the scheme knows, some it has never seen."""
    import anndata as ad
    labels = ["Treg"] * 4 + ["Hepatocyte"] * 5 + ["Kupffer_cell"] * 3 + ["Unassigned"] * 2 + ["Hepatocyte_unspecified"] * 2
    rough = ["T_NK"] * 4 + ["Parenchyma"] * 5 + ["Macrophage_Mono"] * 3 + ["Unassigned"] * 2 + ["Parenchyma"] * 2
    a = ad.AnnData(X=np.zeros((len(labels), 2), dtype=np.float32),
                   obs=pd.DataFrame({"ct": labels, "rough": rough}, index=[f"c{i}" for i in range(len(labels))]))
    a.obsm["spatial"] = np.random.default_rng(0).normal(size=(len(labels), 2))
    a.write_h5ad(path)
    return path


def test_palette_command_exports_the_whole_scheme(tmp_path):
    out = tmp_path / "scheme.csv"
    assert main(["palette", "-o", str(out)]) == 0
    t = pd.read_csv(out).set_index("label")
    assert list(t.columns) == ["hex", "tier", "legend_group"]
    assert {"Treg", "Tumor_nonCSC", "Unassigned", "Macrophage_Mono", "Endothelial"} <= set(t.index)
    assert t.loc["Treg", "tier"] == "bright" and t.loc["Unassigned", "tier"] == "grey" and t.loc["T_NK", "tier"] == "deep"
    assert t.loc["Treg", "hex"] == _run_palette()["Treg"]


def test_scheme_applies_to_a_dataset_annotated_elsewhere(tmp_path):
    import anndata as ad
    h5 = _other_dataset(tmp_path / "liver.h5ad")
    out = tmp_path / "liver_colors.csv"
    assert main(["palette", "--h5ad", str(h5), "--keys", "ct", "rough", "--write-h5ad", "-o", str(out),
                 "--map", str(tmp_path / "liver_map"), "--map-key", "ct"]) == 0
    t = pd.read_csv(out).set_index("label")
    ours = _run_palette()
    assert t.loc["Treg", "hex"] == ours["Treg"]                          # a known label keeps its colour everywhere
    assert t.loc["Hepatocyte", "legend_group"] == "Other" and t.loc["Hepatocyte", "tier"] == "bright"
    assert t.loc["Hepatocyte_unspecified", "tier"] == "deep" and t.loc["Unassigned", "tier"] == "grey"
    assert _delta_e(t.loc["Hepatocyte", "hex"], t.loc["Kupffer_cell", "hex"]) >= 15
    a = ad.read_h5ad(h5)
    assert list(a.uns["ct_colors"]) == [t.loc[c, "hex"] for c in a.obs["ct"].cat.categories]
    assert list(a.uns["rough_colors"]) == [t.loc[c, "hex"] for c in a.obs["rough"].cat.categories]
    assert (tmp_path / "liver_map.png").exists() and (tmp_path / "liver_map.pdf").exists()


def test_saved_scheme_gives_another_dataset_identical_colours(tmp_path):
    """Project workflow: export once, point every dataset at the file."""
    scheme = tmp_path / "project_palette.csv"
    assert main(["palette", "-o", str(scheme)]) == 0
    edited = pd.read_csv(scheme)
    edited.loc[edited["label"] == "Treg", "hex"] = "#ff0000"
    edited.to_csv(scheme, index=False)
    h5 = _other_dataset(tmp_path / "liver.h5ad")
    out = tmp_path / "liver_colors.csv"
    assert main(["palette", "--palette-file", str(scheme), "--h5ad", str(h5), "--keys", "ct", "-o", str(out)]) == 0
    t = pd.read_csv(out).set_index("label")
    assert t.loc["Treg", "hex"] == "#ff0000" and t.loc["Unassigned", "hex"] == "#bdbdbd"
    import anndata as ad
    assert "ct_colors" not in ad.read_h5ad(h5).uns                        # the file is only changed on request
