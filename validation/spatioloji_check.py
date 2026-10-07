#!/usr/bin/env python
"""Real-data check of the spatioloji_s bridge on the Xenium 5K OV section.

1. Load the section with spatioloji's own loader and with annotation_st's reader; the bridge's
   AnnData must have the same cells, genes, counts and control-probe sums.
2. Write the validated labels (validation/ov_full) into the spatioloji object.
3. Draw the fine annotation with a spatioloji plotting function and the shared colours.
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import spatioloji_s as sj

from annotation_st import spatioloji_bridge as bridge
from annotation_st.config import load_config
from annotation_st.io import read_xenium

HERE = Path(__file__).resolve().parent
OUT = HERE / "spatioloji_check"
XEN = "/easley/scratch/projects/amitra/amitra2016502/Xenium_5k_Human_OV_FFPE_outs"
report = {"spatioloji_s": sj.__version__}

t = time.time()
sp = sj.spatioloji.from_xenium(XEN, load_boundaries=False)
report["from_xenium_s"] = round(time.time() - t, 1)
a = bridge.to_anndata(sp)
b = read_xenium(XEN)
report["n_cells"] = [int(a.n_obs), int(b.n_obs)]
report["n_genes"] = [int(a.n_vars), int(b.n_vars)]
report["n_control_features"] = int(a.uns["source"]["n_control_features"])
same_cells = set(a.obs_names) == set(b.obs_names)
same_genes = set(a.var_names) == set(b.var_names)
report["same_cells"], report["same_genes"] = bool(same_cells), bool(same_genes)
if same_cells and same_genes:
    b = b[a.obs_names, a.var_names]
    report["counts_identical"] = bool((a.X != b.X).nnz == 0)
    for col in ("neg_probe_counts", "neg_codeword_counts", "genomic_ctrl_counts",
                "unassigned_codeword_counts_h5", "deprecated_codeword_counts_h5"):
        report[f"{col}_identical"] = bool(col in a.obs and np.array_equal(a.obs[col].to_numpy(float), b.obs[col].to_numpy(float)))
    report["coordinates_identical"] = bool(np.allclose(a.obsm["spatial"], b.obsm["spatial"]))
del a, b

cfg = load_config(HERE.parent / "configs" / "xenium_ov_5k.yaml")
table = HERE / "ov_full" / "results" / "final_annotation_by_cell.tsv.gz"
cols = bridge.add_annotations(sp, table, groups=cfg.plot.legend_groups, qc_group=cfg.plot.qc_group)
labels = pd.read_csv(table, sep="\t", index_col=0, dtype=str, keep_default_na=False)
fine = sp.cell_meta["cell_type_fine"].astype(str)
report["columns_written"] = cols
report["n_labelled"] = int((fine != "NA").sum())
report["n_NA"] = int((fine == "NA").sum())
report["labels_identical"] = bool((fine.reindex(labels.index).to_numpy() == labels["cell_type_fine"].to_numpy()).all())
pal = bridge.palette(sp, "cell_type_fine", cfg)
saved = pd.read_csv(HERE / "ov_full" / "results" / "cell_type_colors.csv").set_index("label")["hex"]
report["palette_matches_run_colour_table"] = bool(all(pal[k] == saved[k] for k in pal if k != "NA"))
report["n_palette_entries"] = len(pal)

try:
    import matplotlib
    matplotlib.use("Agg")
    sj.visualization.xenium_plot_spatial(sp, "cell_type_fine", color_dict=pal, dot_size=0.6, figsize=(12, 14),
                                         save_dir=str(OUT), filename="spatioloji_xenium_plot_spatial_cell_type_fine.png",
                                         dpi=150, show=False)
    report["spatioloji_plot"] = "ok"
except Exception as e:                      # the data checks above stand on their own
    report["spatioloji_plot"] = f"failed: {type(e).__name__}: {e}"

(OUT / "report.json").write_text(json.dumps(report, indent=1) + "\n")
print(json.dumps(report, indent=1))
ok = all(report.get(k) for k in ("same_cells", "same_genes", "counts_identical", "neg_probe_counts_identical",
                                 "coordinates_identical", "labels_identical", "palette_matches_run_colour_table"))
sys.exit(0 if ok else 1)
