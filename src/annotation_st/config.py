"""Configuration dataclasses and YAML loading.

Defaults reproduce the parameters of the documented run (analysis ``01_annot_Xenium_OV``,
run ``2026-09-08``, final job 1180895). A YAML config may be partial: only the keys
present override the defaults. Unknown keys raise ``ValueError`` so typos are caught.
"""
from __future__ import annotations

import copy
import dataclasses
import re
from dataclasses import dataclass, field, fields, is_dataclass
from pathlib import Path
from typing import Any

import yaml

from . import markers as _markers
from . import plotting as _plotting


class ConfigError(ValueError):
    """The configuration is invalid. The message lists every problem found."""


@dataclass
class ProjectConfig:
    name: str = "sample"
    xenium_dir: str | None = None            # input: a 10x Xenium output folder ...
    spatioloji: str | None = None            # ... or a saved spatioloji_s object (pickle or components folder)
    outdir: str = "annot_run"
    checkpoint_dir: str | None = None       # default: <outdir>/checkpoints
    prefix: str = "sample"                  # checkpoint file prefix
    seed: int = 0
    n_jobs: int = 8
    figures: bool = True


@dataclass
class QCConfig:
    subsample: int = 0                      # randomly keep N cells (0 = all); for tests
    min_transcripts: int = 30
    min_genes: int = 10
    max_control_frac: float = 0.02
    area_pct: tuple[float, float] = (1.0, 99.0)
    require_nucleus: bool = True
    min_cells_per_gene: int = 10
    # This package does not do QC filtering: the qc stage computes the metrics and reports what the
    # thresholds above would remove, but keeps every cell and gene. Set True only to reproduce a run
    # that did filter (the validated OV run); otherwise filter upstream (e.g. spatioloji_s QC).
    apply_filters: bool = False


@dataclass
class LineageConfig:
    n_pcs: int = 50
    n_neighbors: int = 15
    resolution: float = 1.0
    min_z: float = 0.5                      # cluster mean z needed to assign a lineage
    cycling_z: float = 1.0                  # proliferation z above which a cell is 'cycling'
    umap_cells: int = 300_000               # UMAP on a random subset (0 = all cells)
    skip_umap: bool = False
    cycling_lineage: str = "Epithelial_tumor"
    cycling_label: str = "Epithelial_tumor_cycling"
    unassigned: str = "Unassigned"


@dataclass
class CafConfig:
    include_lineages: tuple[str, ...] = ("Fibroblast",)
    n_pcs: int = 30
    n_neighbors: int = 15
    resolution: float = 0.5
    min_z: float = 0.0
    margin: float = 0.25                    # best z - second best z for a per-cell call
    min_cells: int = 50
    unassigned: str = "FB_unassigned"
    marker_files: dict[str, str] = field(default_factory=dict)   # {iCAF: path.txt, ...}


@dataclass
class CscConfig:
    tumor_lineages: tuple[str, ...] = ("Epithelial_tumor",)
    z_thresh: float = 1.5
    min_core: int = 2
    n_pcs: int = 30
    n_neighbors: int = 15
    resolution: float = 0.4
    umap_cells: int = 250_000
    min_cells: int = 50
    positive_label: str = "CSC_like"
    negative_label: str = "Tumor_nonCSC"
    # Replace the built-in OCSC set with your own list (one gene per line, '#' comments).
    # Precedence: csc.marker_file > markers.csc (inline list) > built-in default.
    marker_file: str | None = None
    # Core markers for the ">= min_core detected" rule. Precedence: csc.core_marker_file >
    # markers.csc_core restricted to the set in use; if none of those genes is in the set
    # (custom list), every gene of the set counts as core.
    core_marker_file: str | None = None


@dataclass
class ImmuneConfig:
    compartments: tuple[str, ...] | None = None   # None = all compartments in markers.compartments
    n_pcs: int = 30
    n_neighbors: int = 15
    resolution: float = 1.5
    min_z: float = 0.5                      # minority / cross / spillover threshold (mean z)
    low_identity: float = -0.5              # core identities qualify at mean z >= this
    resolution2: float = 1.0
    min_z2: float = 0.25                    # level-2 state threshold
    spill_flag: float = 1.0                 # immune_spillover flag threshold
    umap_cells: int = 200_000
    min_cells: int = 100
    min_cells_level2: int = 50


@dataclass
class SpatialConfig:
    n_perms: int = 100
    radius: float = 0.0                     # >0: fixed-radius graph instead of Delaunay
    edge_cut_pct: float = 99.0              # drop Delaunay edges longer than this percentile
    fine_key: str = "cell_type_fine"
    composition_groups: tuple[str, ...] = ("CSC_like", "Tumor_nonCSC")


@dataclass
class PlotConfig:
    """Figure-only settings for the final annotation map (they never invalidate a stage)."""
    # legend of spatial_cell_type_fine: group title -> labels, both in display order; colours
    # are computed so that the labels of one group are as distinct as possible
    legend_groups: dict[str, list[str]] = field(default_factory=lambda: copy.deepcopy(_plotting.LEGEND_HIERARCHY))
    qc_group: str = _plotting.QC_GROUP       # the group drawn in greys
    # labels without a specific fine type (deep tones); every other cell-type label is bright
    generic_labels: list[str] = field(default_factory=lambda: list(_plotting.GENERIC_LABELS))
    # fixed colours: a table with columns label, hex (e.g. an accepted fine_label_colors.tsv) ...
    palette_file: str | None = None
    # ... and/or inline {label: "#rrggbb"} (wins over the file); other labels are placed around them
    colors: dict[str, str] = field(default_factory=dict)


@dataclass
class CompartmentConfig:
    lineages: list[str]
    level1: dict[str, list[str]]
    level2: dict[str, dict[str, list[str]]]
    core: list[str]
    unresolved: str
    cross: dict[str, str] = field(default_factory=dict)
    unspecified: dict[str, str] = field(default_factory=dict)
    handoff: dict[str, str] = field(default_factory=dict)
    coarse: dict[str, str] = field(default_factory=dict)


@dataclass
class MarkerConfig:
    """Marker sets are opt-in: every field is empty unless the config lists it. A set can be the
    word ``default`` (the package's referenced default library), an inline gene list, or a file."""
    lineage: dict[str, list[str]] = field(default_factory=dict)
    proliferation: list[str] = field(default_factory=list)
    csc: list[str] = field(default_factory=list)
    csc_core: list[str] = field(default_factory=list)
    spillover: dict[str, list[str]] = field(default_factory=dict)
    compartments: dict[str, CompartmentConfig] = field(default_factory=dict)
    caf: dict[str, list[str]] = field(default_factory=dict)   # CAF sets (or caf.marker_files)
    # where each set came from ("default library", "inline", a file path); informational only
    sources: dict[str, str] = field(default_factory=dict, compare=False)


@dataclass
class PipelineConfig:
    project: ProjectConfig = field(default_factory=ProjectConfig)
    qc: QCConfig = field(default_factory=QCConfig)
    lineage: LineageConfig = field(default_factory=LineageConfig)
    caf: CafConfig = field(default_factory=CafConfig)
    csc: CscConfig = field(default_factory=CscConfig)
    immune: ImmuneConfig = field(default_factory=ImmuneConfig)
    spatial: SpatialConfig = field(default_factory=SpatialConfig)
    plot: PlotConfig = field(default_factory=PlotConfig)
    markers: MarkerConfig = field(default_factory=MarkerConfig)


# ----------------------------------------------------------------------------
# dict <-> dataclass
# ----------------------------------------------------------------------------
_TUPLE_FIELDS = {"area_pct", "include_lineages", "tumor_lineages", "compartments", "composition_groups"}


def _from_dict(cls, d: dict[str, Any] | None, section: str, base=None):
    """Build/merge dataclass ``cls`` from dict ``d`` on top of ``base`` (or the class defaults)."""
    obj = base if base is not None else cls()
    if not d:
        return obj
    names = {f.name for f in fields(cls)}
    for k, v in d.items():
        if k not in names:
            raise ConfigError(f"unknown key '{k}' in section '{section}' (allowed: {sorted(names)})")
        if k in _TUPLE_FIELDS and isinstance(v, list):
            v = tuple(v)
        setattr(obj, k, v)
    return obj


DEFAULT = "default"                      # the word that selects a set of the default library
_MARKER_KEYS = ("lineage", "proliferation", "csc", "csc_core", "spillover", "compartments", "caf")


class _MarkerResolver:
    """Turns the ``markers`` section as written (``default`` / gene lists / file paths) into gene lists."""

    def __init__(self, base_dir=None):
        self.base_dir = Path(base_dir) if base_dir is not None else None
        self.sources: dict[str, str] = {}
        self._lib = None

    @property
    def lib(self) -> dict:
        if self._lib is None:
            from . import library
            self._lib = library.default_library()
        return self._lib

    def genes(self, value, where: str, default_value):
        """One gene set: ``default`` -> ``default_value`` (None when the library has none), list, or file."""
        if isinstance(value, str):
            if value == DEFAULT:
                if default_value is None:
                    raise ConfigError(f"markers.{where}: the default library has no such set")
                self.sources[where] = "default library"
                return list(default_value)
            f = Path(value)
            if not f.is_absolute() and self.base_dir is not None:
                f = self.base_dir / f
            if not f.is_file():
                raise ConfigError(f"markers.{where}: '{value}' is neither the word 'default' nor a readable "
                                  f"marker file ({f} not found)")
            self.sources[where] = str(f)
            return _markers.read_marker_file(f)
        self.sources[where] = "inline"
        return value

    def set_dict(self, value, section: str, lib_sets: dict):
        """A mapping of set name -> gene set, or ``default`` for every library set of the section."""
        if value == DEFAULT:
            for name in lib_sets:
                self.sources[f"{section}.{name}"] = "default library"
            return copy.deepcopy(lib_sets)
        if not isinstance(value, dict):
            return value                                       # reported by validate_config
        out = {}
        for name, v in value.items():
            if v == DEFAULT and name not in lib_sets:
                raise ConfigError(f"markers.{section}.{name}: the default library has no set '{name}' "
                                  f"(it offers: {sorted(lib_sets)})")
            out[name] = self.genes(v, f"{section}.{name}", lib_sets.get(name))
        return out

    def compartment(self, name: str, value):
        lib_comps = self.lib["compartments"]
        if value == DEFAULT:
            if name not in lib_comps:
                raise ConfigError(f"markers.compartments.{name}: the default library has no compartment '{name}' "
                                  f"(it offers: {sorted(lib_comps)})")
            self.sources[f"compartments.{name}"] = "default library"
            return copy.deepcopy(lib_comps[name])
        if not isinstance(value, dict):
            raise ConfigError(f"markers.compartments.{name}: must be 'default' or a mapping")
        c = dict(value)
        # identity and state sets may come from any library compartment (a custom compartment can reuse them)
        l1_lib: dict[str, list[str]] = {}
        l2_lib: dict[str, dict[str, list[str]]] = {}
        for cname in [name] + [n for n in lib_comps if n != name]:
            for k, v in lib_comps.get(cname, {}).get("level1", {}).items():
                l1_lib.setdefault(k, v)
            for parent, states in lib_comps.get(cname, {}).get("level2", {}).items():
                l2_lib.setdefault(parent, states)
        w = f"compartments.{name}"
        if "level1" in c:
            c["level1"] = self.set_dict(c["level1"], f"{w}.level1", l1_lib)
        if isinstance(c.get("level2"), dict):
            level2 = {}
            for parent, states in c["level2"].items():
                if states == DEFAULT and parent not in l2_lib:
                    raise ConfigError(f"markers.{w}.level2.{parent}: the default library has no states for '{parent}'")
                level2[parent] = self.set_dict(states, f"{w}.level2.{parent}", l2_lib.get(parent, {}))
            c["level2"] = level2
        return c


def _marker_config(d, base: MarkerConfig | None = None, base_dir=None) -> MarkerConfig:
    """Resolve the ``markers`` section onto ``base``. Only what is listed is used (opt-in)."""
    if d == DEFAULT:
        d = {"preset": DEFAULT}
    if not isinstance(d, dict):
        raise ConfigError("markers: must be a mapping, or the word 'default' for the whole default library")
    allowed = set(_MARKER_KEYS) | {"preset"}
    for k in d:
        if k not in allowed:
            raise ConfigError(f"unknown key '{k}' in section 'markers' (allowed: {sorted(allowed)})")
    r = _MarkerResolver(base_dir)
    merged = dataclasses.asdict(base) if base is not None else dataclasses.asdict(MarkerConfig())
    r.sources.update(merged.pop("sources", {}))
    preset = d.get("preset")
    if preset is not None:
        if preset != DEFAULT:
            raise ConfigError(f"markers.preset: unknown preset '{preset}' (available: '{DEFAULT}')")
        d = {**{k: DEFAULT for k in _MARKER_KEYS}, **{k: v for k, v in d.items() if k != "preset"}}
    lib = r.lib if any(v == DEFAULT or isinstance(v, dict) for v in d.values()) else {}
    for key in ("lineage", "spillover", "caf"):
        if key in d:
            merged[key] = r.set_dict(d[key] if d[key] is not None else {}, key, lib.get(key, {}))
    for key in ("proliferation", "csc", "csc_core"):
        if key in d:
            merged[key] = r.genes(d[key] if d[key] is not None else [], key, lib.get(key) or None)
    if "compartments" in d:
        v = d["compartments"]
        if v == DEFAULT:
            v = {name: DEFAULT for name in r.lib["compartments"]}
        if not isinstance(v, dict) and v is not None:
            raise ConfigError("markers.compartments: must be 'default' or a mapping of compartment name -> definition")
        merged["compartments"] = {name: r.compartment(name, c) for name, c in (v or {}).items()}
    comps = {}
    allowed_c = {f.name for f in fields(CompartmentConfig)}
    for name, c in merged["compartments"].items():
        for k in c:
            if k not in allowed_c:
                raise ConfigError(f"unknown key '{k}' in markers.compartments.{name} (allowed: {sorted(allowed_c)})")
        missing = [f.name for f in fields(CompartmentConfig)
                   if f.default is dataclasses.MISSING and f.default_factory is dataclasses.MISSING and f.name not in c]
        if missing:
            raise ConfigError(f"markers.compartments.{name}: missing key(s) {missing}")
        comps[name] = CompartmentConfig(**c)
    merged["compartments"] = comps
    return MarkerConfig(**merged, sources=r.sources)


# ----------------------------------------------------------------------------
# validation
# ----------------------------------------------------------------------------
def _is_int(v) -> bool:
    return isinstance(v, int) and not isinstance(v, bool)


def _is_num(v) -> bool:
    return isinstance(v, (int, float)) and not isinstance(v, bool)


def _is_strs(v) -> bool:
    return isinstance(v, (list, tuple)) and all(isinstance(x, str) for x in v)


# annotation (as written on the dataclass) -> (checker, human description)
_KINDS = {
    "int": (_is_int, "an integer"),
    "float": (_is_num, "a number"),
    "bool": (lambda v: isinstance(v, bool), "true or false"),
    "str": (lambda v: isinstance(v, str), "a string"),
    "str | None": (lambda v: v is None or isinstance(v, str), "a string or null"),
    "tuple[float, float]": (lambda v: isinstance(v, (list, tuple)) and len(v) == 2 and all(_is_num(x) for x in v),
                            "two numbers"),
    "tuple[str, ...]": (_is_strs, "a list of strings"),
    "tuple[str, ...] | None": (lambda v: v is None or _is_strs(v), "a list of strings or null"),
    "dict[str, str]": (lambda v: isinstance(v, dict) and all(isinstance(k, str) and isinstance(x, str) for k, x in v.items()),
                       "a mapping of strings"),
}

# (section, field) -> (predicate on a correctly-typed value, requirement text)
_RANGES = {
    ("project", "n_jobs"): (lambda v: v >= 1, "must be >= 1"),
    ("project", "prefix"): (lambda v: bool(v) and "/" not in v and "\\" not in v, "must be a non-empty name without path separators"),
    ("project", "outdir"): (lambda v: bool(v), "must not be empty"),
    ("qc", "subsample"): (lambda v: v >= 0, "must be >= 0"),
    ("qc", "min_transcripts"): (lambda v: v >= 0, "must be >= 0"),
    ("qc", "min_genes"): (lambda v: v >= 0, "must be >= 0"),
    ("qc", "min_cells_per_gene"): (lambda v: v >= 0, "must be >= 0"),
    ("qc", "max_control_frac"): (lambda v: 0 <= v <= 1, "must be between 0 and 1"),
    ("qc", "area_pct"): (lambda v: 0 <= v[0] < v[1] <= 100, "must be [low, high] with 0 <= low < high <= 100"),
    ("caf", "margin"): (lambda v: v >= 0, "must be >= 0"),
    ("csc", "min_core"): (lambda v: v >= 0, "must be >= 0"),
    ("immune", "resolution2"): (lambda v: v > 0, "must be > 0"),
    ("immune", "min_cells_level2"): (lambda v: v >= 1, "must be >= 1"),
    ("spatial", "n_perms"): (lambda v: v >= 1, "must be >= 1"),
    ("spatial", "radius"): (lambda v: v >= 0, "must be >= 0"),
    ("spatial", "edge_cut_pct"): (lambda v: 0 < v <= 100, "must be in (0, 100]"),
}
for _sec in ("lineage", "caf", "csc", "immune"):
    _RANGES[(_sec, "n_pcs")] = (lambda v: v >= 1, "must be >= 1")
    _RANGES[(_sec, "n_neighbors")] = (lambda v: v >= 1, "must be >= 1")
    _RANGES[(_sec, "resolution")] = (lambda v: v > 0, "must be > 0")
    _RANGES[(_sec, "umap_cells")] = (lambda v: v >= 0, "must be >= 0")
    _RANGES[(_sec, "min_cells")] = (lambda v: v >= 1, "must be >= 1")


def _check_gene_sets(sets, where: str, problems: list[str]) -> None:
    if not isinstance(sets, dict):
        problems.append(f"{where}: must be a mapping of set name -> gene list")
        return
    for name, genes in sets.items():
        if not isinstance(genes, (list, tuple)) or not genes or not all(isinstance(g, str) and g for g in genes):
            problems.append(f"{where}.{name}: must be a non-empty list of gene symbols (strings)")


def _check_gene_list(genes, where: str, problems: list[str], allow_empty: bool = False) -> None:
    if not isinstance(genes, (list, tuple)) or not all(isinstance(g, str) and g for g in genes) \
            or (not genes and not allow_empty):
        problems.append(f"{where}: must be a {'' if allow_empty else 'non-empty '}list of gene symbols (strings)")


def validate_config(cfg: PipelineConfig) -> list[str]:
    """Structural checks that do not depend on which stages will run: value types and
    ranges, gene-set shapes, and the cross-references inside the immune compartments.
    Returns a list of problems (empty when the config is valid). Stage-specific checks
    (lineage names, marker files, input folder) are done by ``pipeline.preflight``."""
    problems: list[str] = []
    for sec in _SECTIONS:
        obj = getattr(cfg, sec)
        for f in fields(obj):
            v = getattr(obj, f.name)
            check, what = _KINDS.get(str(f.type), (None, None))
            if check is not None and not check(v):
                problems.append(f"{sec}.{f.name}: expected {what}, got {v!r}")
                continue
            rng = _RANGES.get((sec, f.name))
            if rng is not None and not rng[0](v):
                problems.append(f"{sec}.{f.name}: {rng[1]} (got {v!r})")
    if isinstance(cfg.csc.positive_label, str) and cfg.csc.positive_label == cfg.csc.negative_label:
        problems.append("csc.positive_label: must differ from csc.negative_label")

    lg = cfg.plot.legend_groups
    if not isinstance(lg, dict) or not all(isinstance(k, str) and isinstance(ls, (list, tuple)) and ls
                                           and all(isinstance(x, str) and x for x in ls) for k, ls in lg.items()):
        problems.append("plot.legend_groups: must be a mapping of group title -> non-empty list of labels")
    gl = cfg.plot.generic_labels
    if not isinstance(gl, (list, tuple)) or not all(isinstance(x, str) and x for x in gl):
        problems.append("plot.generic_labels: must be a list of labels")
    pc = cfg.plot.colors
    if not isinstance(pc, dict) or not all(isinstance(k, str) and isinstance(v, str)
                                           and re.fullmatch(r"#[0-9a-fA-F]{6}", v) for k, v in pc.items()):
        problems.append("plot.colors: must be a mapping of label -> '#rrggbb' colour")

    m = cfg.markers
    _check_gene_sets(m.lineage, "markers.lineage", problems)
    _check_gene_list(m.proliferation, "markers.proliferation", problems, allow_empty=True)
    _check_gene_list(m.csc, "markers.csc", problems, allow_empty=True)
    _check_gene_list(m.csc_core, "markers.csc_core", problems, allow_empty=True)
    _check_gene_sets(m.caf, "markers.caf", problems)
    _check_gene_sets(m.spillover, "markers.spillover", problems)
    if isinstance(m.spillover, dict):
        for name in m.spillover:
            if not str(name).startswith("spillover_"):
                problems.append(f"markers.spillover.{name}: spillover set names must start with 'spillover_'")

    comps = m.compartments
    str_map = _KINDS["dict[str, str]"][0]
    malformed = set()
    for name, c in comps.items():
        w = f"markers.compartments.{name}"
        shape_ok = True
        if not isinstance(c.core, (list, tuple)) or not _is_strs(c.core):
            problems.append(f"{w}.core: must be a list of identity names")
            shape_ok = False
        for key in ("cross", "unspecified", "handoff", "coarse"):
            if not str_map(getattr(c, key)):
                problems.append(f"{w}.{key}: must be a mapping of label -> name (strings); got {getattr(c, key)!r}")
                shape_ok = False
        if not shape_ok:                      # cross-references below assume well-shaped values
            malformed.add(name)
            continue
        if not _is_strs(c.lineages) or not c.lineages:
            problems.append(f"{w}.lineages: must be a non-empty list of lineage names")
        _check_gene_sets(c.level1, f"{w}.level1", problems)
        l1 = set(c.level1) if isinstance(c.level1, dict) else set()
        if not l1:
            problems.append(f"{w}.level1: at least one identity set is required")
        for ident in l1:
            if str(ident).startswith("spillover_"):
                problems.append(f"{w}.level1.{ident}: identity names must not start with 'spillover_' (reserved for spillover sets)")
        for ident in c.core:
            if ident not in l1:
                problems.append(f"{w}.core: '{ident}' has no level1 marker set")
        if not isinstance(c.unresolved, str) or not c.unresolved:
            problems.append(f"{w}.unresolved: must be a non-empty label")
        if isinstance(c.level2, dict):
            for parent, states in c.level2.items():
                if parent not in l1:
                    problems.append(f"{w}.level2: parent '{parent}' has no level1 marker set")
                _check_gene_sets(states, f"{w}.level2.{parent}", problems)
        else:
            problems.append(f"{w}.level2: must be a mapping of parent identity -> state sets")
        for key in c.unspecified:
            if key not in l1:
                problems.append(f"{w}.unspecified: '{key}' has no level1 marker set")
        for key in c.cross:
            if key not in l1:
                problems.append(f"{w}.cross: '{key}' has no level1 marker set")
        for key in c.coarse:
            if key not in l1:
                problems.append(f"{w}.coarse: '{key}' has no level1 marker set")
        for label, target in c.handoff.items():
            if label not in l1:
                problems.append(f"{w}.handoff: label '{label}' has no level1 marker set")
            if target == name:
                problems.append(f"{w}.handoff: '{label}' is handed off to the compartment itself")
            elif target not in comps:
                problems.append(f"{w}.handoff: target compartment '{target}' is not defined")
    for name, c in comps.items():
        if name in malformed:
            continue
        senders = [s for s, sc in comps.items() if s not in malformed and s != name and name in sc.handoff.values()]
        if senders and c.handoff:
            problems.append(f"markers.compartments.{name}: receives hand-offs (from {', '.join(senders)}) and also "
                            f"hands cells off; only one hop is supported")
    if cfg.immune.compartments is not None and _is_strs(cfg.immune.compartments):
        for name in cfg.immune.compartments:
            if name not in comps:
                problems.append(f"immune.compartments: '{name}' is not defined in markers.compartments")
    return problems


_SECTIONS = {"project": ProjectConfig, "qc": QCConfig, "lineage": LineageConfig, "caf": CafConfig,
             "csc": CscConfig, "immune": ImmuneConfig, "spatial": SpatialConfig, "plot": PlotConfig}


def config_from_dict(d: dict[str, Any] | None, base: PipelineConfig | None = None, base_dir=None) -> PipelineConfig:
    """Merge a (possibly partial) config dict onto ``base`` (default: :func:`default_config`).
    Marker files named in the ``markers`` section resolve against ``base_dir``."""
    cfg = base if base is not None else PipelineConfig()
    d = d or {}
    for k in d:
        if k not in _SECTIONS and k != "markers":
            raise ConfigError(f"unknown top-level section '{k}' (allowed: {sorted([*_SECTIONS, 'markers'])})")
    for sec, cls in _SECTIONS.items():
        setattr(cfg, sec, _from_dict(cls, d.get(sec), sec, base=getattr(cfg, sec)))
    if "markers" in d:
        cfg.markers = _marker_config(d["markers"] if d["markers"] is not None else {}, base=cfg.markers, base_dir=base_dir)
    problems = validate_config(cfg)
    if problems:
        raise ConfigError("invalid configuration:\n  - " + "\n  - ".join(problems))
    return cfg


def default_config() -> PipelineConfig:
    """Default parameters and NO marker sets: markers are opt-in (``preset_config`` opts in to the library)."""
    return PipelineConfig()


def preset_config(preset: str = DEFAULT) -> PipelineConfig:
    """Default parameters with every set of the default marker library opted in."""
    return config_from_dict({"markers": {"preset": preset}})


def load_config(path: str | Path) -> PipelineConfig:
    with open(path) as fh:
        d = yaml.safe_load(fh)
    if d is not None and not isinstance(d, dict):
        raise ConfigError(f"{path}: top level must be a mapping")
    return config_from_dict(d, base_dir=Path(path).resolve().parent)


def config_to_dict(cfg: PipelineConfig) -> dict[str, Any]:
    d = dataclasses.asdict(cfg)
    d["markers"].pop("sources", None)                       # informational, not part of the config

    def _tuples(o):
        if isinstance(o, tuple):
            return list(o)
        if isinstance(o, dict):
            return {k: _tuples(v) for k, v in o.items()}
        if isinstance(o, list):
            return [_tuples(v) for v in o]
        return o
    return _tuples(d)


def dump_config(cfg: PipelineConfig, path: str | Path) -> None:
    Path(path).write_text(yaml.safe_dump(config_to_dict(cfg), sort_keys=False, default_flow_style=None, width=110))
