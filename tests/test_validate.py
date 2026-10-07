"""Config validation: structural problems are reported at load time, all at once;
stage-specific problems (lineage names, marker files, input folder) at preflight."""
import textwrap

import pytest

from annotation_st.config import ConfigError, default_config, load_config, validate_config
from conftest import ov_config


def _load(tmp_path, text):
    p = tmp_path / "c.yaml"
    p.write_text(textwrap.dedent(text))
    return load_config(p)


def test_default_config_is_valid():
    assert validate_config(ov_config()) == []


@pytest.mark.parametrize("yaml_text, expected", [
    ("qc: {min_transcripts: '30'}", "qc.min_transcripts"),                       # string where a number is needed
    ("qc: {require_nucleus: 1}", "qc.require_nucleus"),                          # int where a bool is needed
    ("qc: {area_pct: [99, 1]}", "qc.area_pct"),                                  # window upside down
    ("qc: {area_pct: [1, 50, 99]}", "qc.area_pct"),                              # wrong length
    ("qc: {max_control_frac: 1.5}", "qc.max_control_frac"),
    ("lineage: {resolution: 0}", "lineage.resolution"),
    ("lineage: {n_pcs: 0}", "lineage.n_pcs"),
    ("immune: {resolution2: -1}", "immune.resolution2"),
    ("spatial: {edge_cut_pct: 150}", "spatial.edge_cut_pct"),
    ("spatial: {n_perms: 0}", "spatial.n_perms"),
    ("csc: {positive_label: X, negative_label: X}", "csc.positive_label"),
    ("project: {n_jobs: 0}", "project.n_jobs"),
    ("project: {prefix: 'a/b'}", "project.prefix"),
    ("markers: {lineage: {Foo: []}}", "markers.lineage.Foo"),                    # empty gene list
    ("markers: {lineage: {Foo: [A, 3]}}", "markers.lineage.Foo"),                # non-string gene
    ("markers: {spillover: {Tumor: [EPCAM]}}", "markers.spillover.Tumor"),       # must be named spillover_*
    ("immune: {compartments: [TNK, Nope]}", "immune.compartments"),
])
def test_structural_problem_is_reported_with_its_key(tmp_path, yaml_text, expected):
    with pytest.raises(ConfigError) as e:
        _load(tmp_path, yaml_text)
    assert expected in str(e.value)


def _comp(**over):
    base = dict(lineages=["T_NK"], level1={"A": ["G1"], "B": ["G2"], "cross_X": ["G3"]}, level2={"A": {"A1": ["G4"]}},
                core=["A"], unresolved="C_unresolved", cross={"cross_X": "Macrophage_Mono"}, unspecified={"A": "A_unspecified"},
                handoff={}, coarse={})
    base.update(over)
    return base


@pytest.mark.parametrize("over, expected", [
    (dict(core=["Z"]), "core"),                                   # core identity without a level-1 set
    (dict(level2={"Z": {"Z1": ["G"]}}), "level2"),                # state parent without a level-1 set
    (dict(cross={"cross_Q": "T_NK"}), "cross"),
    (dict(unspecified={"Z": "Z_unspecified"}), "unspecified"),
    (dict(handoff={"Z": "Other"}), "handoff"),                    # label that is not a level-1 identity
    (dict(handoff={"B": "Nowhere"}), "Nowhere"),                  # target compartment does not exist
    (dict(handoff={"B": "C"}), "itself"),                         # hands off to itself
    (dict(coarse={"Z": "B_cell"}), "coarse"),
    (dict(level1={"spillover_A": ["G1"]}, core=[], level2={}, cross={}, unspecified={}), "spillover_"),
    (dict(lineages=[]), "lineages"),
    (dict(level1={}), "level1"),
])
def test_compartment_cross_references(over, expected):
    import yaml
    from annotation_st.config import config_from_dict
    other = _comp(lineages=["Plasma"])
    with pytest.raises(ConfigError) as e:
        config_from_dict({"markers": {"compartments": {"C": _comp(**over), "Other": other}}})
    assert "markers.compartments.C" in str(e.value) and expected in str(e.value)


def test_handoff_chain_is_rejected():
    """Hand-off ordering supports one hop: a receiving compartment must not hand off itself."""
    from annotation_st.config import config_from_dict
    a = _comp(handoff={"B": "Other"})
    b = _comp(lineages=["Plasma"], handoff={"B": "Third"})
    c = _comp(lineages=["Macrophage_Mono"])
    with pytest.raises(ConfigError, match="receives hand-offs"):
        config_from_dict({"markers": {"compartments": {"C": a, "Other": b, "Third": c}}})


def test_all_problems_are_listed_together(tmp_path):
    with pytest.raises(ConfigError) as e:
        _load(tmp_path, """
            qc: {min_genes: -1}
            lineage: {resolution: 0}
        """)
    assert "qc.min_genes" in str(e.value) and "lineage.resolution" in str(e.value)


def test_replacing_lineage_sets_alone_is_a_valid_config(tmp_path):
    """Another tissue: only the lineage stage is wanted, so dangling CAF/CSC lineage names are not a load error."""
    cfg = _load(tmp_path, "markers: {lineage: {Foo: [A, B], Bar: [C]}}")
    assert list(cfg.markers.lineage) == ["Foo", "Bar"]


# ---- preflight (stage-specific) ------------------------------------------------
def test_preflight_catches_lineage_name_typo_only_for_selected_stages():
    from annotation_st.pipeline import preflight
    cfg = ov_config()
    cfg.markers.caf = {"iCAF": ["IL6"]}
    cfg.caf.include_lineages = ("Fibroblasts",)                   # typo: would silently select 0 cells
    assert preflight(cfg, ["lineage", "csc"], check_inputs=False) is None
    with pytest.raises(ConfigError, match="Fibroblasts"):
        preflight(cfg, ["caf"], check_inputs=False)


def test_preflight_missing_marker_file_and_missing_caf_sets(tmp_path):
    from annotation_st.pipeline import preflight
    cfg = ov_config()
    cfg.csc.marker_file = "nope.txt"
    with pytest.raises(ConfigError, match="nope.txt"):
        preflight(cfg, ["csc"], config_dir=tmp_path, check_inputs=False)
    cfg = ov_config()                                         # no CAF sets opted in: the stage will skip itself,
    cfg.caf.include_lineages = ("Typo",)                      # so neither that nor its settings are a problem
    assert preflight(cfg, ["caf"], check_inputs=False) is None
    cfg.caf.marker_files = {"iCAF": "missing_icaf.txt"}
    with pytest.raises(ConfigError, match="missing_icaf.txt"):
        preflight(cfg, ["caf"], config_dir=tmp_path, check_inputs=False)


def test_preflight_checks_xenium_folder(tmp_path):
    from annotation_st.pipeline import preflight
    cfg = ov_config()
    cfg.project.xenium_dir = str(tmp_path / "missing")
    with pytest.raises(ConfigError, match="cell_feature_matrix.h5"):
        preflight(cfg, ["qc"])
    cfg.project.xenium_dir = None
    with pytest.raises(ConfigError, match="xenium_dir"):
        preflight(cfg, ["qc"])


@pytest.mark.parametrize("over, expected", [
    (dict(handoff={"B": ["Other", "B"]}), "handoff"),        # the reference scripts' (target, identity) form is not supported
    (dict(cross={"cross_X": ["T_NK"]}), "cross"),
    (dict(coarse=["B"]), "coarse"),
    (dict(core="A"), "core"),
    (dict(unspecified={"A": 3}), "unspecified"),
])
def test_wrongly_shaped_compartment_values_are_config_errors_not_crashes(over, expected):
    from annotation_st.config import config_from_dict
    with pytest.raises(ConfigError) as e:
        config_from_dict({"markers": {"compartments": {"C": _comp(**over), "Other": _comp(lineages=["Plasma"])}}})
    assert "markers.compartments.C" in str(e.value) and expected in str(e.value)
