"""Per-stage run records and staleness checks.

``<outdir>/stage_state.json`` keeps, for every stage that ran, the inputs that determine
its result (parameters and the *resolved* gene lists, so editing a marker file counts),
a unique ``run_id``, the ``run_id`` of each upstream stage it consumed, and the files it
wrote. A stage is **up to date** only when

* its recorded outputs still exist,
* its current inputs equal the recorded ones, and
* every upstream stage is itself up to date and has not been re-run since.

Settings that only affect figures or speed (``project.figures``, ``project.n_jobs``,
``*.umap_cells``, ``lineage.skip_umap``) are deliberately not inputs. Note that Leiden
partitions can shift with the thread count; ``n_jobs`` is recorded for reference.
"""
from __future__ import annotations

import dataclasses
import datetime as _dt
import json
import os
import uuid
from pathlib import Path

from .config import PipelineConfig
from .io import RunPaths, versions

DEPENDS: dict[str, list[str]] = {
    "qc": [],
    "lineage": ["qc"],
    "caf": ["lineage"],
    "csc": ["lineage"],
    "immune": ["lineage"],
    "merge": ["lineage", "caf", "csc", "immune"],
    "spatial": ["merge"],
}
ORDER = list(DEPENDS)

_NOT_INPUTS = {
    "lineage": {"skip_umap", "umap_cells"},
    "caf": {"marker_files"},
    "csc": {"umap_cells", "marker_file", "core_marker_file"},
    "immune": {"umap_cells", "compartments"},
}


def _section(cfg: PipelineConfig, name: str) -> dict:
    d = dataclasses.asdict(getattr(cfg, name))
    return {k: v for k, v in d.items() if k not in _NOT_INPUTS.get(name, ())}


def selected_compartments(cfg: PipelineConfig) -> list[str]:
    from .stages.immune import compartment_order
    return compartment_order(cfg.markers.compartments, cfg.immune.compartments)


def stage_inputs(name: str, cfg: PipelineConfig, config_dir=None) -> dict:
    """Everything that determines the result of stage ``name``, as plain JSON data."""
    from .stages.caf import resolve_caf_sets
    from .stages.csc import resolve_csc_markers

    seed = cfg.project.seed
    if name == "qc":
        xd = cfg.project.xenium_dir
        qc = _section(cfg, "qc")
        if qc.get("apply_filters", True):                    # keys added after 0.1.0 are omitted at their
            qc.pop("apply_filters", None)                    # default so that earlier runs stay up to date
        d = {"xenium_dir": os.path.abspath(xd) if xd else None, "seed": seed, "qc": qc}
        if cfg.project.spatioloji:
            f = Path(cfg.project.spatioloji)
            d["spatioloji"] = str(f if f.is_absolute() or config_dir is None else Path(config_dir) / f)
    elif name == "lineage":
        d = {"seed": seed, "lineage": _section(cfg, "lineage"), "lineage_sets": cfg.markers.lineage,
             "proliferation": cfg.markers.proliferation}
    elif name == "caf":
        try:
            sets = resolve_caf_sets(cfg.caf, cfg.markers.caf, base=config_dir)
        except (ValueError, OSError) as e:
            sets = {"unresolved": str(e)}
        d = {"seed": seed, "caf": _section(cfg, "caf"), "caf_sets": sets}
    elif name == "csc":
        try:
            csc, core, _ = resolve_csc_markers(cfg.csc, cfg.markers, base=config_dir)
        except (ValueError, OSError) as e:
            csc, core = {"unresolved": str(e)}, None
        d = {"seed": seed, "csc": _section(cfg, "csc"), "csc_genes": csc, "core_genes": core}
    elif name == "immune":
        comps = {c: dataclasses.asdict(cfg.markers.compartments[c]) for c in selected_compartments(cfg)}
        d = {"seed": seed, "immune": _section(cfg, "immune"), "spillover": cfg.markers.spillover, "compartments": comps}
    elif name == "merge":
        d = {"coarse": {c: cfg.markers.compartments[c].coarse for c in selected_compartments(cfg)}}
    elif name == "spatial":
        d = {"seed": seed, "spatial": _section(cfg, "spatial")}
    else:
        raise KeyError(f"unknown stage '{name}'")
    return json.loads(json.dumps(d, sort_keys=True, default=list))      # tuples -> lists, canonical


def _flatten(d, prefix=""):
    out = {}
    if isinstance(d, dict):
        for k, v in d.items():
            out.update(_flatten(v, f"{prefix}{k}."))
    else:
        out[prefix[:-1]] = d
    return out


def changed_keys(old: dict, new: dict) -> list[str]:
    a, b = _flatten(old), _flatten(new)
    return sorted(k for k in set(a) | set(b) if a.get(k, "<absent>") != b.get(k, "<absent>"))


class StageState:
    """Reads and writes ``<outdir>/stage_state.json``."""

    def __init__(self, paths: RunPaths):
        self.path = Path(paths.outdir) / "stage_state.json"
        self.data: dict[str, dict] = {}
        if self.path.exists():
            try:
                self.data = json.loads(self.path.read_text())
            except json.JSONDecodeError:
                self.data = {}

    def get(self, stage: str) -> dict | None:
        return self.data.get(stage)

    def completed(self, stage: str) -> bool:
        """True when the stage ran and produced outputs (False if never run or skipped)."""
        rec = self.data.get(stage)
        return bool(rec) and rec.get("status") == "completed"

    def outputs(self, stage: str) -> list[str]:
        rec = self.data.get(stage)
        return list(rec.get("outputs", [])) if rec else []

    def record(self, stage: str, cfg: PipelineConfig, outputs, status: str = "completed", note: str = "",
               config_dir=None) -> dict:
        rec = {
            "run_id": f"{_dt.datetime.now().strftime('%Y%m%dT%H%M%S')}-{uuid.uuid4().hex[:8]}",
            "finished": _dt.datetime.now().isoformat(timespec="seconds"),
            "status": status,
            "note": note,
            "inputs": stage_inputs(stage, cfg, config_dir),
            "upstream": {d: self.data[d]["run_id"] for d in DEPENDS[stage] if d in self.data},
            "outputs": [str(o) for o in outputs],
            "n_jobs": cfg.project.n_jobs,
            "versions": versions(),
        }
        self.data[stage] = rec
        self.path.parent.mkdir(parents=True, exist_ok=True)
        tmp = self.path.with_suffix(".json.tmp")
        tmp.write_text(json.dumps(self.data, indent=2, default=str) + "\n")
        tmp.replace(self.path)
        return rec


def _own_status(name: str, cfg: PipelineConfig, state: StageState, config_dir) -> tuple[str, str]:
    rec = state.get(name)
    if rec is None:
        return "not_run", "never run in this folder"
    if not all(Path(o).exists() for o in rec.get("outputs", [])):
        return "not_run", "recorded outputs are missing"
    cur = stage_inputs(name, cfg, config_dir)
    if cur != rec.get("inputs"):
        keys = changed_keys(rec.get("inputs", {}), cur)
        shown = ", ".join(keys[:5]) + (f" (+{len(keys) - 5} more)" if len(keys) > 5 else "")
        return "stale", f"inputs changed: {shown}"
    for d in DEPENDS[name]:
        dep = state.get(d)
        if rec.get("upstream", {}).get(d) != (dep["run_id"] if dep else None):
            return "stale", f"upstream stage '{d}' was run after this stage"
    return "up_to_date", "skipped itself: " + rec["note"] if rec.get("status") == "skipped" else ""


def stage_status(name: str, cfg: PipelineConfig, paths: RunPaths, state: StageState | None = None,
                 config_dir=None) -> tuple[str, str]:
    """``(status, reason)`` with status in ``up_to_date`` / ``stale`` / ``not_run``. Staleness is
    transitive: a stage whose upstream is out of date is itself out of date."""
    state = state or StageState(paths)
    status, reason = _own_status(name, cfg, state, config_dir)
    if status != "up_to_date":
        return status, reason
    for d in DEPENDS[name]:
        ds, dr = stage_status(d, cfg, paths, state, config_dir)
        if ds == "stale":
            return "stale", f"upstream stage '{d}' is out of date ({dr})"
    return status, reason
