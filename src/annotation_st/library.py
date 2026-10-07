"""The default marker library: gene sets a config can opt in to with the word ``default``.

Every gene of every set carries at least one literature reference (``references()``):

* ``curated``         a CellMarker 2.0 record (Hu et al., Nucleic Acids Res 2023, PMID 36300619)
                      for that cell type and gene; the record's own PMID is the reference.
* ``primary_text``    the gene named with the cell type in a sentence of one of the set's
                      primary references (the sentence is kept in ``detail``).
* ``literature_text`` the same kind of sentence from another open-access paper, read and
                      accepted by hand (``tools/reviewed_literature_evidence.tsv``).

Candidate genes for which no such reference was found are not in the library; they are listed
by ``dropped()``. The tables are built by ``tools/build_marker_references.py``.
"""
from __future__ import annotations

import copy
import json
from functools import lru_cache
from pathlib import Path

from . import markers as _m

LIBRARY_NAME = "default"
LIB_DIR = Path(__file__).resolve().parent / "data" / "marker_library"


@lru_cache(maxsize=1)
def _flat() -> dict[str, list[str]]:
    return json.loads((LIB_DIR / "default_sets.json").read_text())


@lru_cache(maxsize=1)
def _library() -> dict:
    flat = _flat()

    def section(prefix: str) -> dict[str, list[str]]:
        return {k.split("/", 1)[1]: list(v) for k, v in flat.items() if k.startswith(prefix + "/") and v}

    csc = list(flat.get("csc/CSC", []))
    comps = {}
    for name, c in _m.IMMUNE_COMPARTMENTS.items():
        level1 = section(f"{name}.level1")
        level2 = {}
        for parent in c["level2"]:
            states = section(f"{name}.level2.{parent}")
            if states and parent in level1:
                level2[parent] = states
        comps[name] = {
            "lineages": list(c["lineages"]), "level1": level1, "level2": level2,
            "core": [i for i in c["core"] if i in level1], "unresolved": c["unresolved"],
            "cross": {k: v for k, v in c["cross"].items() if k in level1},
            "unspecified": {k: v for k, v in c["unspecified"].items() if k in level1},
            "handoff": {k: v for k, v in c["handoff"].items() if k in level1},
            "coarse": {k: v for k, v in c["coarse"].items() if k in level1},
        }
    return {"lineage": section("lineage"), "proliferation": list(flat.get("proliferation/Proliferating", [])),
            "csc": csc, "csc_core": [g for g in _m.CSC_CORE if g in csc], "spillover": section("spillover"),
            "compartments": comps, "caf": section("caf")}


def default_library() -> dict:
    """The library in the shape of the ``markers`` config section (a fresh copy)."""
    return copy.deepcopy(_library())


def _table(name: str):
    import pandas as pd

    return pd.read_csv(LIB_DIR / name, sep="\t", dtype=str, keep_default_na=False)


def references():
    """One row per (set, gene, reference): set, gene, evidence, pmid, first_author, year, journal, title, doi, detail."""
    ref = _table("default_references.tsv")
    keep = {(k, g) for k, genes in _flat().items() for g in genes}
    return ref[[(s, g) in keep for s, g in zip(ref["set"], ref["gene"])]].reset_index(drop=True)


def dropped():
    """Candidate genes left out of the library because no reference was found: set, gene, reason."""
    return _table("default_dropped.tsv")


def overview():
    """One row per library set: set, n_genes, n_references, genes."""
    import pandas as pd

    ref = references()
    n_ref = ref.groupby("set")["pmid"].nunique()
    rows = [{"set": k, "n_genes": len(v), "n_references": int(n_ref.get(k, 0)), "genes": ", ".join(v)}
            for k, v in _flat().items() if v]
    return pd.DataFrame(rows)
