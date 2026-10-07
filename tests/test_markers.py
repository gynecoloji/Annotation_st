from annotation_st.markers import (LINEAGE_MARKERS, PROLIF_MARKERS, ov_run_sets,
                                   filter_to_panel, load_marker_files, read_marker_file)


def test_read_marker_file_dedups_and_skips_comments(tmp_path):
    f = tmp_path / "x.txt"
    f.write_text("A\n# c\nB\n\nA\n")
    assert read_marker_file(f) == ["A", "B"]


def test_load_marker_files_relative_to_base(tmp_path):
    (tmp_path / "iCAF.txt").write_text("IL6\nCXCL12\n")
    assert load_marker_files({"iCAF": "iCAF.txt"}, base=tmp_path) == {"iCAF": ["IL6", "CXCL12"]}


def test_filter_to_panel():
    kept, dropped = filter_to_panel({"S": ["A", "Z"], "T": ["Q"]}, ["A", "B"])
    assert kept == {"S": ["A"], "T": []} and dropped == {"S": ["Z"], "T": ["Q"]}


def test_defaults_have_16_lineages_and_prolif():
    assert len(LINEAGE_MARKERS) == 16 and PROLIF_MARKERS[0] == "MKI67"
    assert set(ov_run_sets()) == {"lineage", "proliferation", "csc", "csc_core", "spillover", "compartments", "caf"}


def test_defaults_returns_copy():
    d = ov_run_sets()
    d["lineage"]["Epithelial_tumor"].append("FAKE")
    assert "FAKE" not in LINEAGE_MARKERS["Epithelial_tumor"]
