"""Staleness tracking: a stage is up to date only if its outputs exist, the inputs that
determine its result are unchanged, and no upstream stage has been re-run since."""
import json

import pytest

from annotation_st.config import default_config
from annotation_st.io import RunPaths
from annotation_st.state import StageState, stage_inputs, stage_status
from conftest import ov_config


@pytest.fixture
def env(tmp_path):
    cfg = ov_config()
    cfg.project.outdir = str(tmp_path / "run")
    cfg.project.xenium_dir = str(tmp_path / "xenium")
    cfg.markers.caf = {"iCAF": ["IL6"], "myCAF": ["POSTN"]}
    paths = RunPaths(cfg)
    return cfg, paths, StageState(paths)


def _finish(stage, cfg, paths, state, status="completed"):
    """Pretend ``stage`` just ran: create its output file and record it."""
    out = []
    if status == "completed":
        f = paths.checkpoint("final" if stage in ("merge", "spatial") else stage if stage != "immune" else "immune_TNK")
        f.write_text("x")
        out = [str(f)]
    state.record(stage, cfg, outputs=out, status=status)


def test_never_run_stage_is_not_run(env):
    cfg, paths, state = env
    assert stage_status("qc", cfg, paths, state)[0] == "not_run"


def test_recorded_stage_with_outputs_is_up_to_date_and_persists(env):
    cfg, paths, state = env
    _finish("qc", cfg, paths, state)
    assert stage_status("qc", cfg, paths, state)[0] == "up_to_date"
    assert stage_status("qc", cfg, paths, StageState(paths))[0] == "up_to_date"      # reloaded from disk
    assert "qc" in json.loads((paths.outdir / "stage_state.json").read_text())


def test_deleted_output_means_not_run(env):
    cfg, paths, state = env
    _finish("qc", cfg, paths, state)
    paths.checkpoint("qc").unlink()
    assert stage_status("qc", cfg, paths, state)[0] == "not_run"


def test_changed_parameter_makes_stage_stale_and_names_the_key(env):
    cfg, paths, state = env
    for s in ("qc", "lineage", "csc"):
        _finish(s, cfg, paths, state)
    cfg.csc.z_thresh = 2.0
    status, reason = stage_status("csc", cfg, paths, state)
    assert status == "stale" and "z_thresh" in reason
    assert stage_status("lineage", cfg, paths, state)[0] == "up_to_date"             # other stages unaffected


def test_changed_marker_list_makes_only_its_stage_stale(env):
    cfg, paths, state = env
    for s in ("qc", "lineage", "caf", "csc"):
        _finish(s, cfg, paths, state)
    cfg.markers.csc = ["CD44", "PROM1"]
    assert stage_status("csc", cfg, paths, state)[0] == "stale"
    assert stage_status("caf", cfg, paths, state)[0] == "up_to_date"
    assert stage_status("lineage", cfg, paths, state)[0] == "up_to_date"


def test_edited_marker_file_with_same_path_makes_stage_stale(env, tmp_path):
    cfg, paths, state = env
    f = tmp_path / "csc.txt"
    f.write_text("CD44\nPROM1\n")
    cfg.csc.marker_file = str(f)
    for s in ("qc", "lineage", "csc"):
        _finish(s, cfg, paths, state)
    assert stage_status("csc", cfg, paths, state)[0] == "up_to_date"
    f.write_text("CD44\nPROM1\nKIT\n")
    assert stage_status("csc", cfg, paths, state)[0] == "stale"


def test_same_genes_from_file_or_inline_are_the_same_input(env, tmp_path):
    cfg, paths, state = env
    f = tmp_path / "csc.txt"
    f.write_text("CD44\nPROM1\n")
    cfg.csc.marker_file = str(f)
    a = stage_inputs("csc", cfg)
    cfg.csc.marker_file = None
    cfg.markers.csc = ["CD44", "PROM1"]
    cfg.markers.csc_core = ["CD44", "PROM1"]
    assert stage_inputs("csc", cfg) == a


@pytest.mark.parametrize("mutate", [
    lambda c: setattr(c.project, "figures", False),
    lambda c: setattr(c.project, "n_jobs", 3),
    lambda c: setattr(c.lineage, "umap_cells", 10),
    lambda c: setattr(c.lineage, "skip_umap", True),
    lambda c: setattr(c.csc, "umap_cells", 10),
    lambda c: setattr(c.immune, "umap_cells", 10),
])
def test_figure_and_thread_settings_do_not_invalidate(env, mutate):
    cfg, paths, state = env
    for s in ("qc", "lineage", "caf", "csc", "immune", "merge", "spatial"):
        _finish(s, cfg, paths, state)
    mutate(cfg)
    assert {s: stage_status(s, cfg, paths, state)[0] for s in ("qc", "lineage", "csc", "immune", "spatial")} == \
           {s: "up_to_date" for s in ("qc", "lineage", "csc", "immune", "spatial")}


def test_rerun_upstream_makes_downstream_stale(env):
    cfg, paths, state = env
    for s in ("qc", "lineage", "csc", "merge", "spatial"):
        _finish(s, cfg, paths, state)
    _finish("lineage", cfg, paths, state)                         # lineage re-run (e.g. --force), same config
    status, reason = stage_status("csc", cfg, paths, state)
    assert status == "stale" and "lineage" in reason
    assert stage_status("merge", cfg, paths, state)[0] == "stale"                    # transitively
    assert stage_status("spatial", cfg, paths, state)[0] == "stale"
    assert stage_status("qc", cfg, paths, state)[0] == "up_to_date"


def test_stale_upstream_propagates_before_it_is_rerun(env):
    cfg, paths, state = env
    for s in ("qc", "lineage", "csc", "merge"):
        _finish(s, cfg, paths, state)
    cfg.lineage.min_z = 0.9                                       # lineage out of date but not yet re-run
    assert stage_status("lineage", cfg, paths, state)[0] == "stale"
    status, reason = stage_status("csc", cfg, paths, state)
    assert status == "stale" and "lineage" in reason
    assert stage_status("merge", cfg, paths, state)[0] == "stale"


def test_stage_that_never_ran_does_not_invalidate_merge(env):
    """merge tolerates a missing CAF stage; that absence is not staleness."""
    cfg, paths, state = env
    for s in ("qc", "lineage", "csc", "merge"):
        _finish(s, cfg, paths, state)
    assert stage_status("merge", cfg, paths, state)[0] == "up_to_date"
    _finish("caf", cfg, paths, state)                             # CAF added later -> merge must be redone
    assert stage_status("merge", cfg, paths, state)[0] == "stale"


def test_skipped_stage_counts_as_done_without_outputs(env):
    cfg, paths, state = env
    for s in ("qc", "lineage"):
        _finish(s, cfg, paths, state)
    _finish("caf", cfg, paths, state, status="skipped")
    assert stage_status("caf", cfg, paths, state)[0] == "up_to_date"
    assert state.completed("caf") is False and state.completed("lineage") is True
