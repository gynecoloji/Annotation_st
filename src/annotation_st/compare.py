"""Compare the per-cell labels of a run against a reference table.

Used to validate a re-run (or a refactor, or a parameter change) against an accepted
annotation: for every label column both tables share it reports the agreement on shared
cells, the adjusted Rand index (robust to renumbered clusters), per-label recall and
precision, and the confusion table; numeric columns get a correlation and the largest
difference. Tables may be ``.tsv[.gz]`` / ``.csv[.gz]`` / ``.parquet`` files indexed by
cell id, or ``.h5ad`` files (``obs`` is used).
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


@dataclass
class CompareResult:
    summary: pd.DataFrame                    # index: label column
    per_label: dict[str, pd.DataFrame]       # column -> label table
    confusion: dict[str, pd.DataFrame]       # column -> reference (rows) × new (columns)
    numeric: pd.DataFrame                    # index: numeric column
    n_shared: int
    n_only_new: int
    n_only_ref: int


def _looks_float(s: pd.Series) -> bool:
    vals = s[(s != "") & (s.str.lower() != "nan") & (s != "NA")]
    if vals.empty or not vals.str.contains(r"[.eE]", regex=True).any():
        return False                          # integer-like columns (cluster ids) stay labels
    return pd.to_numeric(vals, errors="coerce").notna().all()


def load_labels(path: str | Path) -> pd.DataFrame:
    """Per-cell table indexed by cell id. Text columns stay strings (cluster ids are not
    turned into numbers, the literal ``NA`` is kept); columns of decimals become floats."""
    path = Path(path)
    name = path.name.lower()
    if name.endswith(".h5ad"):
        import anndata as ad

        df = ad.read_h5ad(path, backed="r").obs.copy()
    elif name.endswith(".parquet"):
        df = pd.read_parquet(path)
        if isinstance(df.index, pd.RangeIndex):
            df = df.set_index(df.columns[0])
    else:
        sep = "," if ".csv" in name else "\t"
        df = pd.read_csv(path, sep=sep, index_col=0, dtype=str, keep_default_na=False)
        for c in df.columns:
            if _looks_float(df[c]):
                df[c] = pd.to_numeric(df[c].replace({"": np.nan, "NA": np.nan}), errors="coerce")
    df.index = df.index.astype(str)
    return df


def _is_numeric(s: pd.Series) -> bool:
    return pd.api.types.is_float_dtype(s.dtype)


def _labels(s: pd.Series) -> pd.Series:
    s = s.astype(object).where(s.notna(), "NA")
    return s.astype(str)


def adjusted_rand_index(a: pd.Series, b: pd.Series) -> float:
    """Hubert–Arabie adjusted Rand index from the contingency table."""
    ct = pd.crosstab(np.asarray(a), np.asarray(b)).to_numpy(dtype=np.float64)
    n = ct.sum()
    if n < 2:
        return float("nan")

    def comb2(x):
        return x * (x - 1) / 2.0
    sum_ij = comb2(ct).sum()
    sum_a, sum_b = comb2(ct.sum(1)).sum(), comb2(ct.sum(0)).sum()
    expected = sum_a * sum_b / comb2(n)
    maximum = 0.5 * (sum_a + sum_b)
    if maximum == expected:
        return 1.0
    return float((sum_ij - expected) / (maximum - expected))


def compare_labels(new: pd.DataFrame, ref: pd.DataFrame, columns=None) -> CompareResult:
    """Compare ``new`` against ``ref`` on the cells they share (matched by index)."""
    new, ref = new.copy(), ref.copy()
    new.index, ref.index = new.index.astype(str), ref.index.astype(str)
    shared = ref.index.intersection(new.index)
    if columns is None:
        columns = [c for c in ref.columns if c in new.columns]
    else:
        missing = [c for c in columns if c not in ref.columns or c not in new.columns]
        if missing:
            raise KeyError(f"column(s) not present in both tables: {missing}")
    r, n = ref.loc[shared], new.loc[shared]

    rows, per_label, confusion, num_rows = {}, {}, {}, {}
    for c in columns:
        if _is_numeric(r[c]) and _is_numeric(n[c]):
            ok = r[c].notna() & n[c].notna()
            x, y = r.loc[ok, c].to_numpy(dtype=float), n.loc[ok, c].to_numpy(dtype=float)
            diff = np.abs(x - y)
            num_rows[c] = {"n_compared": int(ok.sum()),
                           "pearson_r": float(np.corrcoef(x, y)[0, 1]) if len(x) > 1 and x.std() > 0 and y.std() > 0 else float("nan"),
                           "max_abs_diff": float(diff.max()) if len(diff) else float("nan"),
                           "mean_abs_diff": float(diff.mean()) if len(diff) else float("nan"),
                           "n_missing_in_one": int((r[c].isna() != n[c].isna()).sum())}
            continue
        a, b = _labels(r[c]), _labels(n[c])
        same = (a.to_numpy() == b.to_numpy())
        ct = pd.crosstab(a.rename("reference"), b.rename("new"))
        labels = sorted(set(ct.index) | set(ct.columns))
        n_ref = ct.sum(1).reindex(labels, fill_value=0)
        n_new = ct.sum(0).reindex(labels, fill_value=0)
        n_both = pd.Series({l: int(ct.loc[l, l]) if l in ct.index and l in ct.columns else 0 for l in labels})
        pl = pd.DataFrame({"n_ref": n_ref.astype(int), "n_new": n_new.astype(int), "n_both": n_both.astype(int)})
        pl["recall"] = pl["n_both"] / pl["n_ref"].where(pl["n_ref"] > 0)
        pl["precision"] = pl["n_both"] / pl["n_new"].where(pl["n_new"] > 0)
        union = pl["n_ref"] + pl["n_new"] - pl["n_both"]
        pl["jaccard"] = pl["n_both"] / union.where(union > 0)
        pl.index.name = c
        per_label[c], confusion[c] = pl, ct
        rows[c] = {"n_shared": int(len(shared)), "n_only_new": int(len(new) - len(shared)),
                   "n_only_ref": int(len(ref) - len(shared)),
                   "agreement": float(same.mean()) if len(same) else float("nan"),
                   "ari": adjusted_rand_index(a, b),
                   "n_labels_ref": int(a.nunique()), "n_labels_new": int(b.nunique()),
                   "n_disagree": int((~same).sum())}
    summary = pd.DataFrame.from_dict(rows, orient="index")
    summary.index.name = "column"
    numeric = pd.DataFrame.from_dict(num_rows, orient="index")
    numeric.index.name = "column"
    return CompareResult(summary, per_label, confusion, numeric, int(len(shared)),
                         int(len(new) - len(shared)), int(len(ref) - len(shared)))


def write_report(res: CompareResult, outdir: str | Path) -> Path:
    """``compare_summary.tsv``, ``compare_numeric.tsv``, and per column
    ``compare_per_label_<col>.tsv`` / ``compare_confusion_<col>.tsv``."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    res.summary.to_csv(outdir / "compare_summary.tsv", sep="\t")
    if len(res.numeric):
        res.numeric.to_csv(outdir / "compare_numeric.tsv", sep="\t")
    for c, pl in res.per_label.items():
        pl.to_csv(outdir / f"compare_per_label_{c}.tsv", sep="\t")
        res.confusion[c].to_csv(outdir / f"compare_confusion_{c}.tsv", sep="\t")
    return outdir


def format_summary(res: CompareResult) -> str:
    lines = [f"cells: {res.n_shared:,} shared, {res.n_only_new:,} only in new, {res.n_only_ref:,} only in reference"]
    if len(res.summary):
        lines.append(res.summary[["agreement", "ari", "n_disagree", "n_labels_ref", "n_labels_new"]]
                     .round(4).to_string())
    if len(res.numeric):
        lines.append(res.numeric.round(6).to_string())
    return "\n".join(lines)
