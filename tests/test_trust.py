"""End-to-end checks that a run folder can be trusted: changed inputs re-run the right
stages, stale upstream is refused, skipped stages leave no stale labels behind, and the
resolved config is kept with the results."""
import json

import pandas as pd
import pytest

from annotation_st.cli import main
from annotation_st.config import ConfigError, dump_config, load_config
from annotation_st.pipeline import StaleUpstream, run
from test_cli import _small_config

ALL = ["qc", "lineage", "caf", "csc", "immune", "merge", "spatial"]


def _fine_counts(tmp_path):
    return pd.read_csv(tmp_path / "run" / "results" / "final_cell_type_counts.tsv", sep="\t", index_col=0)["n_cells"]


def test_changed_parameter_reruns_that_stage_and_everything_downstream(tmp_path, synth_raw):
    cfgp = _small_config(tmp_path, synth_raw)
    assert run(load_config(cfgp)) == {s: True for s in ALL}
    assert _fine_counts(tmp_path).get("CSC_like", 0) > 0

    cfg = load_config(cfgp)
    cfg.csc.z_thresh = 100.0                         # nothing can pass -> no CSC_like may survive in the merged labels
    dump_config(cfg, cfgp)
    assert run(load_config(cfgp)) == {"qc": False, "lineage": False, "caf": False, "csc": True, "immune": False,
                                      "merge": True, "spatial": True}
    assert _fine_counts(tmp_path).get("CSC_like", 0) == 0
    assert run(load_config(cfgp)) == {s: False for s in ALL}


def test_skipped_stage_leaves_no_stale_labels_and_is_not_retried(tmp_path, synth_raw):
    cfgp = _small_config(tmp_path, synth_raw)
    stages = ["qc", "lineage", "caf", "merge"]
    assert run(load_config(cfgp), stages=stages) == {s: True for s in stages}
    assert _fine_counts(tmp_path).get("iCAF", 0) > 0

    cfg = load_config(cfgp)
    cfg.caf.min_cells = 10_000_000                   # CAF stage now skips itself; its old checkpoint is still on disk
    dump_config(cfg, cfgp)
    assert run(load_config(cfgp), stages=stages) == {"qc": False, "lineage": False, "caf": True, "merge": True}
    fine = _fine_counts(tmp_path)
    assert fine.get("iCAF", 0) == 0 and fine.get("Fibroblast", 0) > 0
    assert run(load_config(cfgp), stages=stages) == {s: False for s in stages}


def test_single_stage_refuses_stale_upstream_unless_allowed(tmp_path, synth_raw):
    cfgp = _small_config(tmp_path, synth_raw)
    assert main(["run", "-c", str(cfgp), "--stages", "qc", "lineage", "csc", "merge"]) == 0
    cfg = load_config(cfgp)
    cfg.csc.min_core = 3
    dump_config(cfg, cfgp)
    with pytest.raises(StaleUpstream, match="csc"):
        main(["merge", "-c", str(cfgp), "--force"])
    assert main(["merge", "-c", str(cfgp), "--force", "--allow-stale"]) == 0


def test_forced_upstream_invalidates_downstream(tmp_path, synth_raw, capsys):
    cfgp = _small_config(tmp_path, synth_raw)
    cfg = load_config(cfgp)
    assert run(cfg, stages=["qc", "lineage", "csc"]) == {"qc": True, "lineage": True, "csc": True}
    assert run(cfg, stages=["lineage"], force=True) == {"lineage": True}
    capsys.readouterr()
    assert main(["status", "-c", str(cfgp)]) == 0
    out = capsys.readouterr().out
    row = {line.split()[0]: line for line in out.splitlines() if line.split() and line.split()[0] in ALL}
    assert "up to date" in row["lineage"] and "stale" in row["csc"] and "lineage" in row["csc"]
    assert "not run" in row["caf"]
    assert run(cfg, stages=["csc"]) == {"csc": True}


def test_run_folder_keeps_a_self_contained_config(tmp_path, synth_raw):
    cfgp = _small_config(tmp_path, synth_raw)
    (tmp_path / "my_csc.txt").write_text("CD44\nPROM1\nPOU5F1\nCXCR4\n")
    cfg = load_config(cfgp)
    cfg.csc.marker_file = "my_csc.txt"               # relative to the config file
    dump_config(cfg, cfgp)
    stages = ["qc", "lineage", "csc"]
    assert main(["run", "-c", str(cfgp), "--stages", *stages]) == 0

    used = tmp_path / "run" / "config_used.yaml"
    kept = load_config(used)
    assert kept.markers.csc == ["CD44", "PROM1", "POU5F1", "CXCR4"] and kept.csc.marker_file is None
    assert kept.project.outdir == str(tmp_path / "run")
    (tmp_path / "my_csc.txt").unlink()               # the saved config no longer needs the marker file ...
    assert run(kept, stages=stages) == {s: False for s in stages}        # ... and describes the same run

    state = json.loads((tmp_path / "run" / "stage_state.json").read_text())
    assert state["csc"]["inputs"]["csc_genes"] == ["CD44", "PROM1", "POU5F1", "CXCR4"]
    assert state["csc"]["upstream"]["lineage"] == state["lineage"]["run_id"]


def test_bad_config_fails_before_any_stage_runs(tmp_path, synth_raw):
    cfgp = _small_config(tmp_path, synth_raw)
    cfg = load_config(cfgp)
    cfg.csc.marker_file = "does_not_exist.txt"
    dump_config(cfg, cfgp)
    with pytest.raises(ConfigError, match="does_not_exist.txt"):
        main(["run", "-c", str(cfgp)])
    assert not (tmp_path / "run" / "checkpoints" / "sample_01_qc.h5ad").exists()
    assert main(["validate", "-c", str(cfgp), "--stages", "qc", "lineage"]) == 0
    assert main(["validate", "-c", str(cfgp)]) == 1
