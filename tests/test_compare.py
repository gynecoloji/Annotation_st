"""Label comparison between a new run and a reference per-cell table."""
import anndata as ad
import numpy as np
import pandas as pd
import pytest

from annotation_st.cli import main
from annotation_st.compare import adjusted_rand_index, compare_labels, load_labels, write_report


def _df(index, **cols):
    return pd.DataFrame(cols, index=pd.Index(index, name="cell_id"))


def test_identical_tables_agree_fully():
    a = _df(list("abcd"), cell_type=["A", "A", "B", "B"])
    res = compare_labels(a, a.copy())
    assert res.summary.loc["cell_type", "agreement"] == 1.0 and res.summary.loc["cell_type", "ari"] == 1.0
    assert res.n_shared == 4 and res.n_only_new == 0 and res.n_only_ref == 0


def test_hand_computed_agreement_per_label_and_ari():
    ref = _df(list("abcd"), cell_type=["A", "A", "B", "B"])
    new = _df(list("abcd"), cell_type=["A", "B", "B", "B"])
    res = compare_labels(new, ref)
    s = res.summary.loc["cell_type"]
    assert s["agreement"] == 0.75 and s["n_labels_ref"] == 2 and s["n_labels_new"] == 2
    assert s["ari"] == pytest.approx(0.0)                                     # contingency [[1,1],[0,2]]
    pl = res.per_label["cell_type"]
    assert pl.loc["A", ["n_ref", "n_new", "n_both"]].tolist() == [2, 1, 1]
    assert pl.loc["A", "recall"] == 0.5 and pl.loc["A", "precision"] == 1.0
    assert pl.loc["B", "recall"] == 1.0 and pl.loc["B", "precision"] == pytest.approx(2 / 3)
    assert res.confusion["cell_type"].loc["A", "B"] == 1                      # rows = reference, columns = new


def test_renumbered_clusters_have_low_agreement_but_perfect_ari():
    ref = _df(list("abcdef"), leiden=["0", "0", "1", "1", "2", "2"])
    new = _df(list("abcdef"), leiden=["2", "2", "0", "0", "1", "1"])
    s = compare_labels(new, ref).summary.loc["leiden"]
    assert s["agreement"] == 0.0 and s["ari"] == 1.0


def test_cells_are_matched_by_id_not_position_and_unshared_cells_are_counted():
    ref = _df(list("abcd"), cell_type=["A", "A", "B", "B"])
    new = _df(list("dcbx"), cell_type=["B", "B", "A", "Q"])
    res = compare_labels(new, ref)
    assert res.n_shared == 3 and res.n_only_new == 1 and res.n_only_ref == 1
    assert res.summary.loc["cell_type", "agreement"] == 1.0


def test_label_present_in_one_table_only():
    ref = _df(list("abc"), cell_type=["A", "A", "B"])
    new = _df(list("abc"), cell_type=["A", "A", "C"])
    pl = compare_labels(new, ref).per_label["cell_type"]
    assert pl.loc["B", ["n_ref", "n_new", "n_both"]].tolist() == [1, 0, 0] and np.isnan(pl.loc["B", "precision"])
    assert pl.loc["C", ["n_ref", "n_new", "n_both"]].tolist() == [0, 1, 0] and np.isnan(pl.loc["C", "recall"])


def test_numeric_columns_get_correlation_and_max_difference():
    ref = _df(list("abcd"), score=[0.0, 1.0, 2.0, 3.0], cell_type=list("AABB"))
    new = _df(list("abcd"), score=[0.0, 1.0, 2.0, 3.5], cell_type=list("AABB"))
    res = compare_labels(new, ref)
    assert "score" not in res.summary.index
    assert res.numeric.loc["score", "max_abs_diff"] == 0.5 and res.numeric.loc["score", "pearson_r"] > 0.99


def test_missing_values_are_compared_as_a_label():
    ref = _df(list("abc"), csc_status=["NA", "CSC_like", None])
    new = _df(list("abc"), csc_status=["NA", "CSC_like", "NA"])
    assert compare_labels(new, ref).summary.loc["csc_status", "agreement"] == 1.0


def test_only_requested_or_shared_columns_are_compared():
    ref = _df(list("ab"), x=["A", "B"], y=["A", "B"])
    new = _df(list("ab"), x=["A", "B"], z=["A", "B"])
    assert list(compare_labels(new, ref).summary.index) == ["x"]
    with pytest.raises(KeyError, match="y"):
        compare_labels(new, ref, columns=["y"])


def test_ari_matches_definition_on_a_known_case():
    # two clusterings of 6 items; value from the Hubert-Arabie formula worked by hand
    x = ["a", "a", "a", "b", "b", "b"]
    y = ["p", "p", "q", "q", "r", "r"]
    # n_ij: a:[2,1,0] b:[0,1,2] -> sum C(nij,2)=2; rows C(3,2)*2=6; cols 3*C(2,2)=3; C(6,2)=15
    # expected=6*3/15=1.2; max=(6+3)/2=4.5; ARI=(2-1.2)/(4.5-1.2)
    assert adjusted_rand_index(pd.Series(x), pd.Series(y)) == pytest.approx(0.8 / 3.3)


def test_loaders_read_tsv_gz_and_h5ad_obs(tmp_path):
    t = _df(list("abc"), cell_type=["A", "B", "B"], leiden=["0", "1", "1"])
    t.to_csv(tmp_path / "t.tsv.gz", sep="\t")
    got = load_labels(tmp_path / "t.tsv.gz")
    assert list(got.index) == list("abc") and list(got["leiden"]) == ["0", "1", "1"]   # cluster ids stay strings
    a = ad.AnnData(X=np.zeros((3, 1), dtype=np.float32), obs=t.copy())
    a.write_h5ad(tmp_path / "t.h5ad")
    assert list(load_labels(tmp_path / "t.h5ad")["cell_type"]) == ["A", "B", "B"]


def test_report_files_and_cli_exit_code(tmp_path):
    ref = _df(list("abcd"), cell_type=["A", "A", "B", "B"], lineage=["A", "A", "B", "B"])
    new = _df(list("abcd"), cell_type=["A", "B", "B", "B"], lineage=["A", "A", "B", "B"])
    ref.to_csv(tmp_path / "ref.tsv.gz", sep="\t")
    new.to_csv(tmp_path / "new.tsv.gz", sep="\t")
    out = tmp_path / "cmp"
    assert main(["compare", "--new", str(tmp_path / "new.tsv.gz"), "--ref", str(tmp_path / "ref.tsv.gz"), "-o", str(out)]) == 0
    summ = pd.read_csv(out / "compare_summary.tsv", sep="\t", index_col=0)
    assert summ.loc["cell_type", "agreement"] == 0.75 and summ.loc["lineage", "agreement"] == 1.0
    assert (out / "compare_per_label_cell_type.tsv").exists() and (out / "compare_confusion_cell_type.tsv").exists()
    # a threshold turns disagreement into a failing exit code
    assert main(["compare", "--new", str(tmp_path / "new.tsv.gz"), "--ref", str(tmp_path / "ref.tsv.gz"),
                 "-o", str(out), "--min-agreement", "0.9"]) == 1
    assert main(["compare", "--new", str(tmp_path / "new.tsv.gz"), "--ref", str(tmp_path / "ref.tsv.gz"),
                 "-o", str(out), "--min-agreement", "0.9", "--columns", "lineage"]) == 0
