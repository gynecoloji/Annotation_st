import json

import pandas as pd
import pytest

from annotation_st.cli import main
from annotation_st.config import dump_config, load_config
from annotation_st.pipeline import StageInputMissing, run
from conftest import CAF_SETS, ov_config, write_synthetic_xenium


def _small_config(tmp_path, synth_raw, figures=False):
    xdir = tmp_path / "xenium"
    write_synthetic_xenium(xdir, synth_raw)
    cfgp = tmp_path / "cfg.yaml"
    assert main(["init-config", "-o", str(cfgp), "--xenium-dir", str(xdir)]) == 0
    cfg = load_config(cfgp)
    assert cfg.project.xenium_dir == str(xdir)
    cfg.project.outdir = str(tmp_path / "run")
    cfg.project.figures = figures
    cfg.project.n_jobs = 2
    cfg.qc.area_pct = (0, 100)
    cfg.qc.apply_filters = True                       # like the validated OV run; the package default removes nothing
    cfg.lineage.n_pcs = 20; cfg.lineage.skip_umap = not figures; cfg.lineage.umap_cells = 400
    cfg.caf.n_pcs = 15
    cfg.csc.n_pcs = 15; cfg.csc.umap_cells = 300
    cfg.immune.n_pcs = 15; cfg.immune.umap_cells = 300; cfg.immune.resolution = 1.0
    cfg.spatial.n_perms = 10
    cfg.markers = ov_config().markers                 # the sets the synthetic populations express
    cfg.markers.caf = {k: list(v) for k, v in CAF_SETS.items()}
    dump_config(cfg, cfgp)
    return cfgp


def test_cli_full_run_on_synthetic_xenium(tmp_path, synth_raw):
    cfgp = _small_config(tmp_path, synth_raw)
    assert main(["run", "-c", str(cfgp)]) == 0
    run_dir = tmp_path / "run"
    ck = run_dir / "checkpoints"
    for f in ["sample_01_qc.h5ad", "sample_02_lineage.h5ad", "sample_03_caf.h5ad", "sample_04_csc.h5ad",
              "sample_05_immune_TNK.h5ad", "sample_05_immune_Bcell.h5ad", "sample_05_immune_Myeloid.h5ad",
              "sample_06_final.h5ad"]:
        assert (ck / f).exists(), f
    res = run_dir / "results"
    for f in ["marker_sets_used.json", "leiden_lineage_scores.tsv", "cell_type_counts.tsv", "leiden_top_markers.tsv",
              "caf_subtype_counts.tsv", "caf_subtype_by_cell.tsv.gz", "csc_counts.tsv", "csc_marker_positivity.tsv",
              "immune_TNK_subcluster_scores.tsv", "immune_subtype_by_cell.tsv.gz", "immune_marker_sets_used.json",
              "final_cell_type_counts.tsv", "final_cell_type_coarse_counts.tsv", "final_annotation_by_cell.tsv.gz",
              "cell_type_colors.csv",
              "nhood_enrichment_zscore.tsv", "csc_neighbourhood_composition.tsv"]:
        assert (res / f).exists(), f
    assert (run_dir / "qc" / "qc_filter_summary.tsv").exists()
    params = json.loads((run_dir / "params.json").read_text())
    assert set(params) >= {"qc_filter", "cluster_annotate", "fibroblast_caf", "cancer_csc", "immune_subcluster",
                           "merge", "spatial_squidpy"}
    assert params["qc_filter"]["n_cells_kept"] < params["qc_filter"]["n_cells_input"]
    fine = pd.read_csv(res / "final_cell_type_counts.tsv", sep="\t", index_col=0)
    assert {"iCAF", "myCAF", "CSC_like", "Tumor_nonCSC", "NK"} <= set(fine.index)
    assert (run_dir / "analysis.log").read_text().count("\n") >= 7
    assert not list((run_dir / "figures").glob("*.png"))
    # the colours travel with the data: recorded in the table and stored in the checkpoints
    import anndata as ad
    colours = pd.read_csv(res / "cell_type_colors.csv").set_index("label")["hex"]
    final = ad.read_h5ad(ck / "sample_06_final.h5ad", backed="r")
    assert list(final.uns["cell_type_fine_colors"]) == [colours[c] for c in final.obs["cell_type_fine"].cat.categories]
    assert list(final.uns["cell_type_colors"]) == [colours[c] for c in final.obs["cell_type"].cat.categories]
    final.file.close()
    lin = ad.read_h5ad(ck / "sample_02_lineage.h5ad", backed="r")
    assert list(lin.uns["cell_type_colors"]) == [colours[c] for c in lin.obs["cell_type"].cat.categories]
    lin.file.close()
    # every subtype column of every stage carries the same colours as the final map
    for h5, keys in (("sample_03_caf.h5ad", ["caf_subtype_cell", "caf_subtype_cluster"]),
                     ("sample_04_csc.h5ad", ["csc_status"]),
                     ("sample_05_immune_TNK.h5ad", ["immune_subtype"]),
                     ("sample_05_immune_Myeloid.h5ad", ["immune_subtype"])):
        a = ad.read_h5ad(ck / h5, backed="r")
        for key in keys:
            cats = list(a.obs[key].cat.categories)
            assert list(a.uns[f"{key}_colors"]) == [colours[c] for c in cats], (h5, key)
        a.file.close()
    tnk = ad.read_h5ad(ck / "sample_05_immune_TNK.h5ad", backed="r")
    lv1 = dict(zip(tnk.obs["immune_level1"].cat.categories, tnk.uns["immune_level1_colors"]))
    assert lv1["NK"] == colours["NK"]                                    # identities too
    tnk.file.close()

    # second run skips every stage
    cfg = load_config(cfgp)
    assert run(cfg) == {s: False for s in ["qc", "lineage", "caf", "csc", "immune", "merge", "spatial"]}
    # forcing one stage re-runs only it
    assert run(cfg, stages=["merge"], force=True) == {"merge": True}


def test_pipeline_keeps_every_cell_and_gene_by_default(tmp_path, synth_raw):
    cfgp = _small_config(tmp_path, synth_raw)
    cfg = load_config(cfgp)
    cfg.qc.apply_filters = False
    dump_config(cfg, cfgp)
    assert main(["run", "-c", str(cfgp), "--stages", "qc", "lineage", "merge"]) == 0
    params = json.loads((tmp_path / "run" / "params.json").read_text())
    assert params["qc_filter"]["n_cells_kept"] == params["qc_filter"]["n_cells_input"] == synth_raw.n_obs
    assert params["qc_filter"]["n_genes_kept"] == synth_raw.n_vars and params["qc_filter"]["apply_filters"] is False
    assert (tmp_path / "run" / "qc" / "qc_filter_summary.tsv").exists()
    fine = pd.read_csv(tmp_path / "run" / "results" / "final_cell_type_counts.tsv", sep="\t", index_col=0)
    assert fine["n_cells"].sum() == synth_raw.n_obs
    text = (tmp_path / "run" / "config_used.yaml").read_text()
    assert "apply_filters: false" in text
    (tmp_path / "fresh.yaml").unlink(missing_ok=True)
    assert main(["init-config", "-o", str(tmp_path / "fresh.yaml")]) == 0
    assert load_config(tmp_path / "fresh.yaml").qc.apply_filters is False


def test_single_stage_needs_its_input(tmp_path, synth_raw):
    cfgp = _small_config(tmp_path, synth_raw)
    with pytest.raises(StageInputMissing, match="run stage 'qc' first"):
        main(["lineage", "-c", str(cfgp)])
    assert main(["qc", "-c", str(cfgp), "--subsample", "800"]) == 0
    params = json.loads((tmp_path / "run" / "params.json").read_text())
    assert params["qc_filter"]["n_cells_input"] == 800


def test_init_config_refuses_overwrite(tmp_path):
    p = tmp_path / "c.yaml"
    assert main(["init-config", "-o", str(p)]) == 0
    assert main(["init-config", "-o", str(p)]) == 1
    assert main(["init-config", "-o", str(p), "--force"]) == 0


@pytest.mark.slow
def test_cli_with_figures(tmp_path, synth_raw):
    cfgp = _small_config(tmp_path, synth_raw, figures=True)
    assert main(["run", "-c", str(cfgp)]) == 0
    figs = tmp_path / "run" / "figures"
    for stem in ["umap_leiden", "spatial_cell_type", "heatmap_cluster_lineage_scores", "umap_fb_caf_subtype",
                 "spatial_tumor_csc", "umap_TNK_subtype", "heatmap_Myeloid_subcluster_scores",
                 "spatial_cell_type_fine", "nhood_enrichment_fine", "csc_neighbour_caf_fraction"]:
        assert (figs / f"{stem}.png").exists() and (figs / f"{stem}.pdf").exists(), stem
    assert (tmp_path / "run" / "qc" / "qc_hist_prefilter.png").exists()
    # all palette-dependent maps can be redrawn from the checkpoints
    for f in figs.glob("*"):
        f.unlink()
    assert main(["figures", "-c", str(cfgp)]) == 0
    for stem in ["umap_cell_type", "spatial_cell_type", "umap_fb_caf_subtype", "spatial_fb_caf_subtype",
                 "umap_tumor_csc", "spatial_tumor_csc", "umap_TNK_subtype", "spatial_TNK_subtype",
                 "spatial_Myeloid_subtype", "spatial_cell_type_fine"]:
        assert (figs / f"{stem}.png").exists() and (figs / f"{stem}.pdf").exists(), stem
    assert not (figs / "dotplot_lineage_markers.png").exists()          # needs expression data: not a redraw


def test_cli_csc_marker_file_relative_to_config(tmp_path, synth_raw):
    cfgp = _small_config(tmp_path, synth_raw)
    (tmp_path / "my_csc.txt").write_text("CD44\nPROM1\nPOU5F1\nCXCR4\n")
    cfg = load_config(cfgp)
    cfg.csc.marker_file = "my_csc.txt"        # relative to the config file's folder
    dump_config(cfg, cfgp)
    assert main(["run", "-c", str(cfgp), "--stages", "qc", "lineage", "csc"]) == 0
    used = json.loads((tmp_path / "run" / "results" / "csc_marker_sets_used.json").read_text())
    assert used["used_on_panel"]["OCSC"] == ["CD44", "PROM1", "POU5F1", "CXCR4"]
    assert used["used_on_panel"]["OCSC_core"] == ["CD44", "PROM1", "POU5F1"]
    assert used["source"]["csc"].endswith("my_csc.txt")
    params = json.loads((tmp_path / "run" / "params.json").read_text())
    assert params["cancer_csc"]["marker_source"]["csc"].endswith("my_csc.txt")
