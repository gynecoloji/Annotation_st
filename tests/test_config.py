import textwrap

import pytest

from annotation_st.config import default_config, dump_config, load_config
from conftest import ov_config


def test_defaults_match_documented_run():
    cfg = ov_config()
    assert cfg.qc.min_transcripts == 30 and cfg.qc.max_control_frac == 0.02
    assert cfg.qc.area_pct == (1.0, 99.0) and cfg.qc.require_nucleus
    assert cfg.lineage.resolution == 1.0 and cfg.lineage.min_z == 0.5 and cfg.lineage.n_pcs == 50
    assert cfg.caf.margin == 0.25 and cfg.caf.min_z == 0.0 and cfg.caf.resolution == 0.5
    assert cfg.csc.z_thresh == 1.5 and cfg.csc.min_core == 2
    assert cfg.immune.low_identity == -0.5 and cfg.immune.min_z2 == 0.25 and cfg.immune.resolution == 1.5
    assert cfg.spatial.n_perms == 100 and cfg.spatial.edge_cut_pct == 99.0
    assert set(cfg.markers.compartments) == {"TNK", "Bcell", "Myeloid"}
    assert cfg.markers.compartments["TNK"].handoff == {"B_cell": "Bcell", "Plasma": "Bcell"}
    assert cfg.markers.compartments["Bcell"].coarse["pDC"] == "Dendritic"
    assert len(cfg.markers.lineage) == 16


def test_partial_yaml_overrides_only_given_keys(tmp_path):
    p = tmp_path / "c.yaml"
    p.write_text(textwrap.dedent("""
        project: {name: demo, outdir: out}
        qc: {min_transcripts: 300}
        markers:
          lineage: {Foo: [A, B], Bar: [C]}
    """))
    cfg = load_config(p)
    assert cfg.project.name == "demo" and cfg.project.seed == 0
    assert cfg.qc.min_transcripts == 300 and cfg.qc.min_genes == 10
    assert cfg.markers.lineage == {"Foo": ["A", "B"], "Bar": ["C"]}
    # marker sets are opt-in: what the file does not list is not used
    assert cfg.markers.proliferation == [] and cfg.markers.compartments == {} and cfg.markers.csc == []


def test_compartment_override(tmp_path):
    p = tmp_path / "c.yaml"
    p.write_text(textwrap.dedent("""
        markers:
          compartments:
            Only:
              lineages: [X]
              level1: {A: [G1]}
              level2: {}
              core: [A]
              unresolved: Only_unresolved
    """))
    cfg = load_config(p)
    assert list(cfg.markers.compartments) == ["Only"]
    assert cfg.markers.compartments["Only"].handoff == {}


def test_unknown_key_raises(tmp_path):
    p = tmp_path / "c.yaml"
    p.write_text("qc: {min_trnscripts: 3}\n")
    with pytest.raises(ValueError, match="unknown key"):
        load_config(p)
    p.write_text("qcc: {min_transcripts: 3}\n")
    with pytest.raises(ValueError, match="unknown top-level"):
        load_config(p)


def test_roundtrip(tmp_path):
    cfg = ov_config()
    p = tmp_path / "c.yaml"
    dump_config(cfg, p)
    assert load_config(p) == cfg
