#!/usr/bin/env python
"""Full in-memory annotation of the Xenium 5K OV section through the spatioloji_s bridge,
compared with the validated command-line run (validation/ov_full).

    sp = sj.spatioloji.from_xenium(...);  bridge.annotate(sp, cfg, qc=True)

Thread counts mirror the validated run because Leiden partitions can shift with them: QC and
lineage with 32 threads, every later stage with 16. The pipeline did that with two processes;
here the switch happens right after the lineage stage returns.
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import spatioloji_s as sj

import annotation_st.stages.lineage as lineage_stage
from annotation_st import spatioloji_bridge as bridge
from annotation_st.compare import compare_labels, format_summary, load_labels, write_report
from annotation_st.config import load_config

HERE = Path(__file__).resolve().parent
OUT = HERE / "spatioloji_full"
XEN = "/easley/scratch/projects/amitra/amitra2016502/Xenium_5k_Human_OV_FFPE_outs"
REF = HERE / "ov_full" / "results" / "final_annotation_by_cell.tsv.gz"

cfg = load_config(HERE.parent / "configs" / "xenium_ov_5k.yaml")
cfg.project.n_jobs = 32

_lineage = lineage_stage.annotate_lineage


def _lineage_then_16_threads(*args, **kwargs):
    res = _lineage(*args, **kwargs)
    import numba
    from threadpoolctl import threadpool_limits

    numba.set_num_threads(16)
    threadpool_limits(16)
    sc.settings.n_jobs = 16
    print(f"[{time.strftime('%X')}] lineage done; switched to 16 threads", flush=True)
    return res


lineage_stage.annotate_lineage = _lineage_then_16_threads

t0 = time.time()
sp = sj.spatioloji.from_xenium(XEN, load_boundaries=False)
print(f"[{time.strftime('%X')}] loaded {sp.n_cells} cells", flush=True)
out = bridge.annotate(sp, cfg, qc=True)
elapsed = time.time() - t0
print(f"[{time.strftime('%X')}] annotated in {elapsed / 60:.1f} min; stages: {out.result.stages}", flush=True)

cols = ["annot_leiden", "lineage", "cell_type", "cell_type_fine", "caf_subtype_cell", "immune_compartment",
        "immune_subtype", "csc_status", "csc_OCSC_z"]
new = sp.cell_meta[cols].rename(columns={"annot_leiden": "leiden"}).copy()
new.index = sp.cell_index.astype(str)
new = new[new["cell_type_fine"].astype(str) != "NA"]            # the cells that passed QC
new.to_csv(OUT / "annotation_by_cell.tsv.gz", sep="\t")
res = compare_labels(load_labels(OUT / "annotation_by_cell.tsv.gz"), load_labels(REF))
write_report(res, OUT / "compare")
print(format_summary(res), flush=True)

report = {"spatioloji_s": sj.__version__, "minutes": round(elapsed / 60, 1), "stages": out.result.stages,
          "n_cells_object": int(sp.n_cells), "n_cells_annotated": int(len(new)),
          "n_shared": res.n_shared, "n_only_new": res.n_only_new, "n_only_ref": res.n_only_ref,
          "agreement": {c: float(res.summary.loc[c, "agreement"]) for c in res.summary.index},
          "ari": {c: float(res.summary.loc[c, "ari"]) for c in res.summary.index},
          "n_disagree": {c: int(res.summary.loc[c, "n_disagree"]) for c in res.summary.index},
          "numeric": res.numeric.to_dict(orient="index")}
try:
    import matplotlib
    matplotlib.use("Agg")
    sj.visualization.xenium_plot_spatial(sp, "cell_type_fine", color_dict=out.palettes["cell_type_fine"], dot_size=0.6,
                                         figsize=(12, 14), save_dir=str(OUT), filename="spatioloji_cell_type_fine.png",
                                         dpi=150, show=False)
    report["plot"] = "ok"
except Exception as e:
    report["plot"] = f"failed: {type(e).__name__}: {e}"
(OUT / "report.json").write_text(json.dumps(report, indent=1, default=str) + "\n")
print(json.dumps({k: report[k] for k in ("minutes", "n_cells_annotated", "agreement", "n_disagree")}, indent=1))
sys.exit(0)
