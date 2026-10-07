"""Marker sets are opt-in: nothing is scored unless the config lists it. The package ships a
default library whose every gene carries a literature reference; a config selects from it with
the word ``default``, or supplies its own lists or files."""
import json
import re
import textwrap

import pandas as pd
import pytest

from annotation_st import library
from annotation_st.cli import main
from annotation_st.config import ConfigError, config_from_dict, default_config, dump_config, load_config, validate_config
from annotation_st.pipeline import StageState, preflight, run
from annotation_st.io import RunPaths
from conftest import write_synthetic_xenium

LIB = library.default_library()


def _load(tmp_path, text):
    p = tmp_path / "c.yaml"
    p.write_text(textwrap.dedent(text))
    return load_config(p)


# ---- nothing is implicit ---------------------------------------------------------------------
def test_a_config_without_markers_scores_nothing():
    m = default_config().markers
    assert m.lineage == {} and m.proliferation == [] and m.csc == [] and m.csc_core == []
    assert m.spillover == {} and m.compartments == {} and m.caf == {}


def test_only_the_listed_sets_are_used(tmp_path):
    cfg = _load(tmp_path, """
        markers:
          lineage:
            Epithelial_tumor: default
            Hepatocyte: [ALB, APOA1, TTR]
    """)
    assert list(cfg.markers.lineage) == ["Epithelial_tumor", "Hepatocyte"]
    assert cfg.markers.lineage["Epithelial_tumor"] == LIB["lineage"]["Epithelial_tumor"]
    assert cfg.markers.lineage["Hepatocyte"] == ["ALB", "APOA1", "TTR"]
    assert cfg.markers.csc == [] and cfg.markers.compartments == {} and cfg.markers.proliferation == []
    assert cfg.markers.sources["lineage.Epithelial_tumor"].startswith("default library")
    assert cfg.markers.sources["lineage.Hepatocyte"] == "inline"


def test_a_whole_section_or_the_whole_library_can_be_selected(tmp_path):
    cfg = _load(tmp_path, "markers: {lineage: default, csc: default, csc_core: default}")
    assert cfg.markers.lineage == LIB["lineage"] and cfg.markers.csc == LIB["csc"] and cfg.markers.compartments == {}
    everything = _load(tmp_path, "markers: default")
    assert everything.markers.lineage == LIB["lineage"] and set(everything.markers.compartments) == set(LIB["compartments"])
    assert everything.markers.caf == LIB["caf"] and everything.markers.spillover == LIB["spillover"]
    assert validate_config(everything) == []
    # preset plus a replacement: the given key wins, the rest stays default
    mixed = _load(tmp_path, "markers: {preset: default, csc: [CD44, KIT], csc_core: [KIT]}")
    assert mixed.markers.csc == ["CD44", "KIT"] and mixed.markers.lineage == LIB["lineage"]


def test_own_sets_from_files_resolve_next_to_the_config(tmp_path):
    (tmp_path / "markers").mkdir()
    (tmp_path / "markers" / "chol.txt").write_text("KRT19\nSOX9\n# comment\n")
    (tmp_path / "prolif.txt").write_text("MKI67\nTOP2A\n")
    cfg = _load(tmp_path, """
        markers:
          lineage: {Cholangiocyte: markers/chol.txt}
          proliferation: prolif.txt
    """)
    assert cfg.markers.lineage == {"Cholangiocyte": ["KRT19", "SOX9"]}
    assert cfg.markers.proliferation == ["MKI67", "TOP2A"]
    assert cfg.markers.sources["lineage.Cholangiocyte"].endswith("chol.txt")


@pytest.mark.parametrize("text, expected", [
    ("markers: {lineage: {Hepatocyte: default}}", "Hepatocyte"),             # the library has no such set
    ("markers: {preset: liver}", "preset"),
    ("markers: {lineage: {X: nope.txt}}", "nope.txt"),
    ("markers: {compartments: {Nope: default}}", "Nope"),
    ("markers: {compartments: {TNK: {lineages: [T_NK], level1: {Mystery: default}, level2: {}, core: [], unresolved: U}}}", "Mystery"),
])
def test_unknown_defaults_and_missing_files_are_config_errors(tmp_path, text, expected):
    with pytest.raises(ConfigError, match=expected):
        _load(tmp_path, text)


def test_fine_types_are_selected_one_by_one(tmp_path):
    cfg = _load(tmp_path, """
        markers:
          compartments:
            Bcell: default
            Myeloid:
              lineages: [Macrophage_Mono]
              core: [Macrophage]
              unresolved: Myeloid_unresolved
              unspecified: {Macrophage: Macrophage_unpolarized}
              level1: {Macrophage: default, Monocyte: default, MyCell: [GENE1, GENE2]}
              level2:
                Macrophage: {TREM2_TAM: default, M1_macrophage: default, Mine: [GENE3]}
    """)
    my = cfg.markers.compartments["Myeloid"]
    lib_my = LIB["compartments"]["Myeloid"]
    assert list(my.level1) == ["Macrophage", "Monocyte", "MyCell"] and my.level1["Macrophage"] == lib_my["level1"]["Macrophage"]
    assert list(my.level2["Macrophage"]) == ["TREM2_TAM", "M1_macrophage", "Mine"]              # M2 was not opted in
    assert my.level2["Macrophage"]["TREM2_TAM"] == lib_my["level2"]["Macrophage"]["TREM2_TAM"]
    assert cfg.markers.compartments["Bcell"].level1 == LIB["compartments"]["Bcell"]["level1"]
    assert "TNK" not in cfg.markers.compartments
    whole_parent = _load(tmp_path, """
        markers:
          compartments:
            TNK: {lineages: [T_NK], core: [CD8_T], unresolved: U, level1: {CD8_T: default}, level2: {CD8_T: default}}
    """)
    assert whole_parent.markers.compartments["TNK"].level2["CD8_T"] == LIB["compartments"]["TNK"]["level2"]["CD8_T"]


def test_explicit_config_round_trips(tmp_path):
    cfg = _load(tmp_path, "markers: default")
    dump_config(cfg, tmp_path / "explicit.yaml")
    assert load_config(tmp_path / "explicit.yaml") == cfg


# ---- the default library: every gene has a reference -------------------------------------------
def _library_pairs():
    out = []
    for name, genes in LIB["lineage"].items(): out += [(f"lineage/{name}", g) for g in genes]
    out += [("proliferation/Proliferating", g) for g in LIB["proliferation"]]
    out += [("csc/CSC", g) for g in LIB["csc"]]
    for name, genes in LIB["spillover"].items(): out += [(f"spillover/{name}", g) for g in genes]
    for name, genes in LIB["caf"].items(): out += [(f"caf/{name}", g) for g in genes]
    for comp, c in LIB["compartments"].items():
        for name, genes in c["level1"].items(): out += [(f"{comp}.level1/{name}", g) for g in genes]
        for parent, states in c["level2"].items():
            for name, genes in states.items(): out += [(f"{comp}.level2.{parent}/{name}", g) for g in genes]
    return out


def test_every_gene_of_the_default_library_has_a_reference():
    ref = library.references()
    assert list(ref.columns) == ["set", "gene", "evidence", "pmid", "first_author", "year", "journal", "title", "doi", "detail"]
    have = set(zip(ref["set"], ref["gene"]))
    missing = [p for p in _library_pairs() if p not in have]
    assert not missing, missing[:20]
    assert len(_library_pairs()) > 400


def test_references_are_well_formed():
    ref = library.references()
    assert ref["pmid"].str.fullmatch(r"\d{6,9}").all()
    assert set(ref["evidence"]) <= {"curated", "primary_text", "literature_text"}
    assert (ref["title"].str.len() > 10).all() and (ref["year"].str.fullmatch(r"(19|20)\d\d")).all()
    assert (ref["detail"].str.len() > 10).all()
    text = ref[ref["evidence"] != "curated"]
    # a text reference quotes a sentence that actually names the gene (or its protein)
    assert all(len(d) > 40 for d in text["detail"])


def test_genes_without_a_reference_are_not_in_the_library():
    dropped = library.dropped()
    pairs = set(_library_pairs())
    assert not [(s, g) for s, g in zip(dropped["set"], dropped["gene"]) if (s, g) in pairs]
    assert len(dropped) > 0 and (dropped["reason"].str.len() > 10).all()


def test_library_is_structurally_sound():
    assert set(LIB["csc_core"]) <= set(LIB["csc"]) and len(LIB["csc_core"]) >= 2
    assert all(genes for genes in LIB["lineage"].values())
    for comp, c in LIB["compartments"].items():
        assert all(c["level1"].values()), comp
        assert set(c["core"]) <= set(c["level1"])
        assert set(c["level2"]) <= set(c["level1"])
        for states in c["level2"].values():
            assert states and all(states.values())
    assert set(LIB["caf"]) == {"iCAF", "myCAF", "apCAF"} and all(len(v) >= 5 for v in LIB["caf"].values())


# ---- pipeline behaviour ------------------------------------------------------------------------
def _run_config(tmp_path, synth_raw, markers_yaml):
    xdir = tmp_path / "xenium"
    write_synthetic_xenium(xdir, synth_raw)
    p = tmp_path / "cfg.yaml"
    p.write_text(textwrap.dedent(f"""
        project: {{xenium_dir: {xdir}, outdir: {tmp_path / 'run'}, figures: false, n_jobs: 2}}
        qc: {{area_pct: [0, 100]}}
        lineage: {{n_pcs: 20, skip_umap: true}}
        csc: {{n_pcs: 15, umap_cells: 0}}
        immune: {{n_pcs: 15, umap_cells: 0, resolution: 1.0}}
        spatial: {{n_perms: 10}}
    """) + textwrap.dedent(markers_yaml))
    return p


def test_lineage_stage_needs_lineage_sets(tmp_path, synth_raw):
    p = _run_config(tmp_path, synth_raw, "")
    with pytest.raises(ConfigError, match="markers.lineage"):
        main(["run", "-c", str(p)])
    assert main(["run", "-c", str(p), "--stages", "qc"]) == 0             # QC needs no markers


def test_refinements_that_were_not_opted_in_are_skipped_and_recorded(tmp_path, synth_raw):
    p = _run_config(tmp_path, synth_raw, "markers: {lineage: default, proliferation: default}\n")
    cfg = load_config(p)
    assert run(cfg) == {s: True for s in ["qc", "lineage", "caf", "csc", "immune", "merge", "spatial"]}
    state = StageState(RunPaths(cfg))
    for stage in ("caf", "csc", "immune"):
        rec = state.get(stage)
        assert rec["status"] == "skipped" and "opted in" in rec["note"], stage
    fine = pd.read_csv(tmp_path / "run" / "results" / "final_cell_type_counts.tsv", sep="\t", index_col=0)
    coarse = pd.read_csv(tmp_path / "run" / "results" / "final_cell_type_coarse_counts.tsv", sep="\t", index_col=0)
    assert set(fine.index) == set(coarse.index)                            # no refinement: fine labels are the lineages
    assert run(load_config(p)) == {s: False for s in ["qc", "lineage", "caf", "csc", "immune", "merge", "spatial"]}
    sources = json.loads((tmp_path / "run" / "results" / "marker_sources.json").read_text())
    assert sources["lineage.Fibroblast"].startswith("default library")


def test_default_library_runs_end_to_end(tmp_path, synth_raw):
    p = _run_config(tmp_path, synth_raw, "markers: default\n")
    assert main(["run", "-c", str(p), "--stages", "qc", "lineage", "csc", "immune", "merge"]) == 0
    fine = pd.read_csv(tmp_path / "run" / "results" / "final_cell_type_counts.tsv", sep="\t", index_col=0)
    assert {"CSC_like", "Tumor_nonCSC", "NK"} <= set(fine.index)
    used = json.loads((tmp_path / "run" / "results" / "marker_sets_used.json").read_text())
    assert used["lineage"]["Fibroblast"] == [g for g in LIB["lineage"]["Fibroblast"] if g in synth_raw.var_names]


# ---- command line ------------------------------------------------------------------------------
def test_init_config_opts_in_to_the_defaults_visibly(tmp_path):
    p = tmp_path / "c.yaml"
    assert main(["init-config", "-o", str(p)]) == 0
    text = p.read_text()
    assert re.search(r"^  lineage: default$", text, re.M) and re.search(r"^  compartments: default$", text, re.M)
    assert load_config(p).markers.lineage == LIB["lineage"]
    assert main(["init-config", "-o", str(p), "--force", "--markers", "none"]) == 0
    assert load_config(p).markers.lineage == {}
    assert main(["init-config", "-o", str(p), "--force", "--markers", "explicit"]) == 0
    assert "EPCAM" in p.read_text() and load_config(p).markers.lineage == LIB["lineage"]


def test_markers_command_lists_sets_and_exports_references(tmp_path, capsys):
    assert main(["markers"]) == 0
    out = capsys.readouterr().out
    assert "lineage/Fibroblast" in out and "TNK.level2.CD4_T/Treg" in out
    assert main(["markers", "--set", "TNK.level2.CD4_T/Treg"]) == 0
    out = capsys.readouterr().out
    assert "FOXP3" in out and "PMID" in out
    assert main(["markers", "--references", str(tmp_path / "refs.tsv"), "--dropped", str(tmp_path / "dropped.tsv")]) == 0
    refs = pd.read_csv(tmp_path / "refs.tsv", sep="\t", dtype=str)
    assert {"set", "gene", "pmid", "evidence"} <= set(refs.columns) and len(refs) > 400
    assert len(pd.read_csv(tmp_path / "dropped.tsv", sep="\t")) > 0
