"""``annot-st`` command line.

    annot-st init-config [-o config.yaml] [--markers default|explicit|none]
    annot-st markers  [--set NAME] [--references refs.tsv]   the default marker library and its references
    annot-st validate -c config.yaml [--stages ...]         check the config before running
    annot-st run      -c config.yaml [--stages qc lineage ...] [--force] [--no-figures]
    annot-st qc|lineage|caf|csc|immune|merge|spatial -c config.yaml [...]
    annot-st status   -c config.yaml                         which stages are up to date / stale
    annot-st figures  -c config.yaml                         redraw the label maps from checkpoints (no recomputation)
    annot-st palette  [-c config.yaml] [-o scheme.csv]       export the colour scheme / apply it to another h5ad
    annot-st to-spatioloji -c config.yaml -o annotated.pkl   write a run's labels into a spatioloji_s object
    annot-st compare  --ref reference.tsv.gz (-c config.yaml | --new labels.tsv.gz)
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

from . import __version__
from .config import ConfigError, config_to_dict, default_config, dump_config, load_config, preset_config
from .pipeline import FIGURES, STAGES, export_palette, preflight, redraw_figures, run, status_table

os.environ.setdefault("MPLBACKEND", "Agg")

_STATUS_LABEL = {"up_to_date": "up to date", "stale": "stale", "not_run": "not run"}


def _add_overrides(p: argparse.ArgumentParser) -> None:
    p.add_argument("-c", "--config", required=True, type=Path, help="YAML config (see `annot-st init-config`)")
    p.add_argument("--outdir", help="override project.outdir")
    p.add_argument("--checkpoint-dir", help="override project.checkpoint_dir")
    p.add_argument("--xenium-dir", help="override project.xenium_dir")
    p.add_argument("--prefix", help="override project.prefix (checkpoint file prefix)")
    p.add_argument("--seed", type=int, help="override project.seed")
    p.add_argument("--subsample", type=int, help="override qc.subsample (keep N random cells; smoke tests)")


def _add_run_options(p: argparse.ArgumentParser) -> None:
    _add_overrides(p)
    p.add_argument("--n-jobs", type=int, help="override project.n_jobs (default: SLURM_CPUS_PER_TASK or config)")
    p.add_argument("--force", action="store_true", help="re-run even when the stage is up to date")
    p.add_argument("--allow-stale", action="store_true",
                   help="run on top of an out-of-date upstream stage instead of refusing")
    p.add_argument("--no-figures", action="store_true", help="skip every figure")


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="annot-st", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--version", action="version", version=f"annot-st {__version__}")
    sub = ap.add_subparsers(dest="cmd", required=True)

    ic = sub.add_parser("init-config", help="write the default (Xenium OV 5K) config for editing")
    ic.add_argument("-o", "--output", type=Path, default=Path("annot_config.yaml"))
    ic.add_argument("--xenium-dir", help="set project.xenium_dir in the written file")
    ic.add_argument("--force", action="store_true", help="overwrite an existing file")
    ic.add_argument("--markers", choices=["default", "explicit", "none"], default="default",
                    help="marker section to write: 'default' opts in to the default library by name (one editable "
                         "line per section), 'explicit' writes its gene lists out, 'none' leaves every set to you")

    mk = sub.add_parser("markers", help="show the default marker library and the reference of every gene")
    mk.add_argument("--set", dest="set_name", help="show the genes of one set with their references (e.g. lineage/Fibroblast)")
    mk.add_argument("--references", type=Path, help="write the full reference table (TSV)")
    mk.add_argument("--dropped", type=Path, help="write the candidate genes left out for lack of a reference (TSV)")

    va = sub.add_parser("validate", help="check the config (values, cross-references, marker files, input folder)")
    _add_overrides(va)
    va.add_argument("--stages", nargs="+", choices=STAGES, default=None, help="only check what these stages need")

    rp = sub.add_parser("run", help="run all (or selected) stages in order; up-to-date stages are skipped")
    rp.add_argument("--stages", nargs="+", choices=STAGES, default=None)
    _add_run_options(rp)
    for s in STAGES:
        _add_run_options(sub.add_parser(s, help=f"run only the '{s}' stage"))

    st = sub.add_parser("status", help="show which stages are up to date, stale or not run")
    _add_overrides(st)

    fg = sub.add_parser("figures", help="redraw figures from existing checkpoints without re-running stages")
    _add_overrides(fg)
    fg.add_argument("--which", nargs="+", choices=FIGURES, default=None, help="figures to redraw (default: all)")

    ts = sub.add_parser("to-spatioloji", help="write the labels of a finished run into a spatioloji_s object")
    _add_overrides(ts)
    ts.add_argument("--input", type=Path, help="saved spatioloji object to annotate (default: project.spatioloji)")
    ts.add_argument("-o", "--output", type=Path, required=True, help="pickle to write")
    ts.add_argument("--label-prefix", default="", help="prefix for the columns added to cell_meta")

    pa = sub.add_parser("palette", help="export the colour scheme, or apply it to a dataset annotated elsewhere")
    pa.add_argument("-c", "--config", type=Path, help="config whose plot section defines the scheme (default: built-in)")
    pa.add_argument("--palette-file", type=Path, help="saved colour table (label, hex) whose colours are kept")
    pa.add_argument("-o", "--output", type=Path, default=Path("cell_type_palette.csv"), help="CSV to write")
    pa.add_argument("--h5ad", type=Path, help="another dataset: its labels are added to the scheme")
    pa.add_argument("--keys", nargs="+", default=[], help="label columns of --h5ad (obs) to colour")
    pa.add_argument("--write-h5ad", action="store_true", help="store uns['<key>_colors'] in --h5ad")
    pa.add_argument("--map", type=Path, help="also draw the hierarchical-legend tissue map of --h5ad to this path stem")
    pa.add_argument("--map-key", help="label column for --map (default: the first of --keys)")

    cp = sub.add_parser("compare", help="compare per-cell labels of a run against a reference table")
    cp.add_argument("--ref", required=True, type=Path, help="reference per-cell table (.tsv[.gz], .csv, .parquet, .h5ad)")
    cp.add_argument("--new", type=Path, help="labels to check (default: <outdir>/results/final_annotation_by_cell.tsv.gz of -c)")
    cp.add_argument("-c", "--config", type=Path, help="config of the run to check (used when --new is not given)")
    cp.add_argument("--columns", nargs="+", help="columns to compare (default: every column both tables have)")
    cp.add_argument("-o", "--outdir", type=Path, help="where to write the report (default: <run outdir>/results/compare or ./compare)")
    cp.add_argument("--min-agreement", type=float, default=None,
                    help="exit with status 1 if any compared label column agrees on a smaller fraction of cells")
    return ap


def _apply_overrides(cfg, a: argparse.Namespace) -> None:
    for attr in ("outdir", "checkpoint_dir", "xenium_dir", "prefix", "seed"):
        v = getattr(a, attr, None)
        if v is not None:
            setattr(cfg.project, attr, v)
    if getattr(a, "n_jobs", None) is not None:
        cfg.project.n_jobs = a.n_jobs
    elif os.environ.get("SLURM_CPUS_PER_TASK"):
        cfg.project.n_jobs = int(os.environ["SLURM_CPUS_PER_TASK"])
    if getattr(a, "subsample", None) is not None:
        cfg.qc.subsample = a.subsample
    if getattr(a, "no_figures", False):
        cfg.project.figures = False


def _cmd_compare(a: argparse.Namespace) -> int:
    from .compare import compare_labels, format_summary, load_labels, write_report

    new_path, outdir = a.new, a.outdir
    if new_path is None:
        if a.config is None:
            print("compare: give --new or -c/--config", file=sys.stderr)
            return 2
        run_out = Path(load_config(a.config).project.outdir)
        new_path = run_out / "results" / "final_annotation_by_cell.tsv.gz"
        outdir = outdir or run_out / "results" / "compare"
    outdir = outdir or Path("compare")
    res = compare_labels(load_labels(new_path), load_labels(a.ref), columns=a.columns)
    write_report(res, outdir)
    print(f"new: {new_path}\nreference: {a.ref}")
    print(format_summary(res))
    print(f"report written to {outdir}")
    if a.min_agreement is not None and len(res.summary):
        low = res.summary.index[res.summary["agreement"] < a.min_agreement].tolist()
        if low:
            print(f"agreement below {a.min_agreement}: {', '.join(low)}", file=sys.stderr)
            return 1
    return 0


_MARKERS_DEFAULT = """\
# Marker sets are opt-in: only what is listed here is scored.
#   default            the package's default set(s); every gene has a literature reference
#                      (`annot-st markers` lists them, `annot-st markers --set <name>` shows the references)
#   [GENE1, GENE2]     your own gene list
#   path/to/file.txt   your own list, one gene per line (relative to this file)
# A section can be opted in as a whole (as below) or set by set, e.g.
#   lineage: {Epithelial_tumor: default, Fibroblast: default, Hepatocyte: [ALB, APOA1, TTR]}
#   compartments: {TNK: default}      # fine types of one immune compartment only
# Remove a line to switch that part off: without `csc` the CSC stage is skipped, and so on.
markers:
  lineage: default
  proliferation: default
  csc: default
  csc_core: default
  spillover: default
  caf: default
  compartments: default
"""
_MARKERS_NONE = """\
# Marker sets are opt-in: only what is listed here is scored. Nothing is listed yet.
# Use the word `default` for a set of the package's referenced default library
# (`annot-st markers` lists them), a gene list, or a file with one gene per line, e.g.
#   markers:
#     lineage: {Epithelial_tumor: default, Hepatocyte: [ALB, APOA1, TTR], Other: markers/other.txt}
#     proliferation: default
markers: {}
"""


def _cmd_init_config(a: argparse.Namespace) -> int:
    import yaml

    if a.output.exists() and not a.force:
        print(f"{a.output} exists (use --force to overwrite)", file=sys.stderr)
        return 1
    cfg = preset_config() if a.markers == "explicit" else default_config()
    if a.xenium_dir:
        cfg.project.xenium_dir = a.xenium_dir
    if a.markers == "explicit":
        dump_config(cfg, a.output)
    else:
        d = config_to_dict(cfg)
        d.pop("markers")
        body = yaml.safe_dump(d, sort_keys=False, default_flow_style=None, width=110)
        a.output.write_text(body + (_MARKERS_DEFAULT if a.markers == "default" else _MARKERS_NONE))
    print(f"wrote {a.output}")
    return 0


def _cmd_markers(a: argparse.Namespace) -> int:
    from . import library

    if a.references:
        library.references().to_csv(a.references, sep="\t", index=False)
        print(f"wrote {a.references}")
    if a.dropped:
        library.dropped().to_csv(a.dropped, sep="\t", index=False)
        print(f"wrote {a.dropped}")
    if a.set_name:
        ref = library.references()
        ref = ref[ref["set"] == a.set_name]
        if ref.empty:
            print(f"no set '{a.set_name}' in the default library (run `annot-st markers` for the list)", file=sys.stderr)
            return 1
        for gene, g in ref.groupby("gene", sort=False):
            print(gene)
            for _, r in g.iterrows():
                print(f"    PMID {r['pmid']}  {r['first_author']} {r['year']}, {r['journal']}  [{r['evidence']}]")
                print(f"        {r['detail'][:200]}")
        return 0
    if not (a.references or a.dropped):
        ov = library.overview()
        print(f"default marker library: {len(ov)} sets, {int(ov['n_genes'].sum())} set-gene pairs, each with a reference")
        print(f"{'set':<46}{'genes':>6}{'refs':>6}")
        for _, r in ov.iterrows():
            print(f"{r['set']:<46}{r['n_genes']:>6}{r['n_references']:>6}")
    return 0


def main(argv=None) -> int:
    a = build_parser().parse_args(argv)
    if a.cmd == "init-config":
        return _cmd_init_config(a)
    if a.cmd == "markers":
        return _cmd_markers(a)
    if a.cmd == "compare":
        return _cmd_compare(a)
    if a.cmd == "palette":
        cfg = load_config(a.config) if a.config else default_config()
        if a.palette_file:
            cfg.plot.palette_file = str(a.palette_file.resolve())
        if (a.keys or a.write_h5ad or a.map) and not a.h5ad:
            print("palette: --keys / --write-h5ad / --map need --h5ad", file=sys.stderr)
            return 2
        pal = export_palette(cfg, a.output, config_dir=a.config.resolve().parent if a.config else None,
                             h5ad=a.h5ad, keys=a.keys, write_h5ad=a.write_h5ad, map_stem=a.map, map_key=a.map_key)
        print(f"wrote {a.output} ({len(pal)} labels)")
        return 0

    config_dir = a.config.resolve().parent
    if a.cmd == "validate":
        try:
            cfg = load_config(a.config)
            _apply_overrides(cfg, a)
            preflight(cfg, a.stages, config_dir=config_dir)
        except ConfigError as e:
            print(e, file=sys.stderr)
            return 1
        print(f"{a.config}: configuration OK for stages {a.stages or STAGES}")
        return 0

    cfg = load_config(a.config)
    _apply_overrides(cfg, a)
    if a.cmd == "status":
        print(f"{'stage':<9}{'status':<12}reason")
        for stage, status, reason in status_table(cfg, config_dir=config_dir):
            print(f"{stage:<9}{_STATUS_LABEL[status]:<12}{reason}")
        return 0
    if a.cmd == "to-spatioloji":
        from . import spatioloji_bridge as bridge
        from .io import RunPaths
        from .pipeline import StageInputMissing, input_path

        table = RunPaths(cfg, mkdir=False).results / "final_annotation_by_cell.tsv.gz"
        if not table.exists():
            raise StageInputMissing(f"{table} not found — run stage 'merge' first")
        src = a.input or (input_path(cfg, config_dir) if cfg.project.spatioloji else None)
        if src is None:
            print("to-spatioloji: give --input, or set project.spatioloji in the config", file=sys.stderr)
            return 2
        obj = bridge.load(src)
        cols = bridge.add_annotations(obj, table, prefix=a.label_prefix, groups=cfg.plot.legend_groups,
                                      qc_group=cfg.plot.qc_group)
        obj.to_pickle(str(a.output))
        print(f"wrote {a.output} with cell_meta columns: {', '.join(cols)}")
        return 0
    if a.cmd == "figures":
        for stem in redraw_figures(cfg, a.which, config_dir=config_dir):
            print(f"wrote {stem}.png and {stem}.pdf")
        return 0
    stages = a.stages if a.cmd == "run" else [a.cmd]
    run(cfg, stages=stages, force=a.force, config_dir=config_dir, allow_stale=a.allow_stale)
    return 0


if __name__ == "__main__":
    sys.exit(main())
