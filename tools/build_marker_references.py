#!/usr/bin/env python
"""Build the gene-level reference table of the default marker library.

For every (set, gene) of the candidate sets this script looks for a reference in two ways:

  curated       a CellMarker 2.0 record (Hu et al., Nucleic Acids Res 2023, PMID 36300619) whose
                cell name matches the set and whose gene symbol is the gene; the record's own PMID
                is the reference.
  primary_text  the gene and the cell type named in the same sentence of the full text (PubMed
                Central) of one of the set's primary references; the sentence is kept as evidence.
  literature_text  a sentence of the same kind from another open-access paper that a person has
                read and accepted: the rows of tools/reviewed_literature_evidence.tsv (candidates
                come from tools/find_literature_evidence.py).

A gene with neither is NOT included in the default library; it is listed in the "dropped" table.

Inputs : CellMarker 2.0 human table (downloaded on demand), PubMed Central full texts (NCBI E-utilities).
Outputs: <out>/default_references.tsv, <out>/default_dropped.tsv, <out>/default_sets.json
Usage  : python tools/build_marker_references.py --out src/annotation_st/data/marker_library --cache <dir>
"""
from __future__ import annotations

import argparse
import html
import json
import re
import time
import urllib.request
from pathlib import Path

import pandas as pd

from annotation_st import markers as M

CELLMARKER_URL = "http://117.50.127.228/CellMarker/CellMarker_download_files/file/Cell_marker_Human.xlsx"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
IDCONV = "https://pmc.ncbi.nlm.nih.gov/tools/idconv/api/v1/articles/"

# ---------------------------------------------------------------------------------------------
# What counts as "the same cell type": regex on the CellMarker cell name / on a sentence.
# (pattern, must_also_match or None). Case-insensitive. EXCLUDE lists names that must not count.
# ---------------------------------------------------------------------------------------------
T_NAME = r"\bT cell|T\(|T-cell|\bT helper|T lymphocyte"
PAT: dict[str, tuple[str, str | None]] = {
    # lineages
    "Epithelial_tumor": (r"^epithelial cell$|^cancer cell$|malignant|tumou?r cell|carcinoma cell|ovarian.*epitheli", None),
    "Fibroblast": (r"fibroblast", None),
    "Pericyte_SMC": (r"pericyte|smooth muscle|mural", None),
    "Endothelial": (r"endotheli", None),
    "Lymphatic_EC": (r"lymphatic endotheli", None),
    "T_NK": (T_NAME + r"|natural killer|\bNK\b", None),
    "B_cell": (r"\bB cell|B\(|B-cell|B lymphocyte", None),
    "Plasma": (r"plasma cell|plasmablast|plasmacyte", None),
    "Macrophage_Mono": (r"macrophage|monocyte", None),
    "Dendritic": (r"dendritic", None),
    "Neutrophil": (r"neutrophil", None),
    "Mast": (r"mast cell", None),
    "Mesothelial": (r"mesotheli", None),
    "Ovarian_stroma": (r"granulosa|theca|ovarian strom|luteal", None),
    "Adipocyte": (r"adipocyte", None),
    "Schwann": (r"schwann", None),
    "Proliferating": (r"prolifer|cycling|mitotic|dividing", None),
    "CSC": (r"cancer stem|stem-like|tumou?r.initiating|side population", None),
    # immune identities
    "CD8_T": (r"CD8", None),
    "CD4_T": (r"CD4|helper|regulatory T|treg", None),
    "NK": (r"natural killer|\bNK\b", None),
    "gdT_MAIT": (r"gamma.?delta|γδ|MAIT|mucosal.associated|invariant", None),
    "mregDC": (r"dendritic", r"mature|LAMP3|migratory|mreg|activated|regulatory|CCR7"),
    "Macrophage": (r"macrophage", None),
    "Monocyte": (r"monocyte", None),
    "cDC": (r"dendritic", None),
    "pDC": (r"plasmacytoid", None),
    "Osteoclast_like_giant": (r"osteoclast|giant cell", None),
    # T / B / myeloid states
    "CD8_cytotoxic": (r"cytotoxic|effector|TEMRA", r"CD8|\bT cell|T\(|T lymphocyte"),
    "CD8_exhausted": (r"exhaust|dysfunction", None),
    "CD8_naive_memory": (r"na[iï]ve|memory", T_NAME),          # CD4-named records excluded below
    "CD4_naive_memory": (r"na[iï]ve|memory", T_NAME),
    "Treg": (r"regulatory T|treg", None),
    "Tfh": (r"follicular helper|Tfh", None),
    "CD4_Th_effector": (r"helper|Th1|Th2|Th17", None),
    "Naive_B": (r"na[iï]ve B", None),
    "Memory_B": (r"memory B|atypical B|age.associated B", None),
    "GC_B": (r"germinal", None),
    "Plasmablast": (r"plasmablast", None),
    "M1_macrophage": (r"\bM1\b|pro-?inflammatory macro|inflammatory macro|classically activated", None),
    "M2_macrophage": (r"\bM2\b|alternatively activated|anti-?inflammatory macro", None),
    "TREM2_TAM": (r"TREM2|lipid.associated|tumou?r.associated macrophage|\bTAM", None),
    "Hypoxic_TAM": (r"hypoxi", None),
    "Monocyte_classical": (r"classical monocyte|CD14\+ monocyte", None),
    "Monocyte_nonclassical": (r"non-?classical|CD16\+ monocyte|patrolling", None),
    "cDC1": (r"cDC1|dendritic cell 1|CD141|CLEC9A|XCR1|type 1 conventional", None),
    "cDC2": (r"cDC2|dendritic cell 2|CD1C|type 2 conventional", None),
    # CAF subtypes
    "iCAF": (r"inflammatory.{0,40}fibroblast|iCAF|inflammatory CAF", None),
    "myCAF": (r"myofibroblastic cancer|myofibroblastic CAF|myCAF", None),
    "apCAF": (r"antigen.present.{0,40}(fibroblast|CAF)|apCAF", None),
}
EXCLUDE: dict[str, str] = {
    "CD8_cytotoxic": r"CD4|natural killer|\bNK\b|regulatory",
    "CD8_naive_memory": r"CD4|regulatory|\bB cell",
    "CD4_naive_memory": r"CD8|regulatory|\bB cell",
    "CD8_T": r"CD4",
    "CD4_T": r"CD8",
    "Proliferating": r"germ",
    "Plasma": r"plasmacytoid",
    "Monocyte_classical": r"non-?classical",
}
ALIAS_SET = {"Proliferating_T": "Proliferating", "Proliferating_B": "Proliferating",
             "Proliferating_myeloid": "Proliferating", "cross_Myeloid": "Macrophage_Mono", "cross_T": "T_NK",
             "spillover_Tumor": "Epithelial_tumor", "spillover_Fibroblast": "Fibroblast",
             "spillover_Endothelial": "Endothelial"}
# In running text a state is also named by these words (used for primary_text only).
TEXT_EXTRA = {
    "Hypoxic_TAM": r"hypoxi|HIF", "TREM2_TAM": r"TREM2|lipid|TAM", "M1_macrophage": r"M1|IFN|LPS|inflammatory",
    "M2_macrophage": r"M2|IL-?4|IL-?10|alternative", "Mesothelial": r"mesotheli", "Epithelial_tumor": r"epitheli|tumou?r|cancer|malignant",
    "Ovarian_stroma": r"granulosa|theca|strom", "Proliferating": r"prolifer|cycl|mitos|G2/M|cell.cycle",
    "Osteoclast_like_giant": r"osteoclast", "Plasma": r"plasma", "Plasmablast": r"plasmablast|plasma",
    "Naive_B": r"na[iï]ve", "Memory_B": r"memory|atypical", "CD4_Th_effector": r"Th1|Th2|Th17|helper|effector",
    "iCAF": r"iCAF|inflammatory", "myCAF": r"myCAF|myofibroblast", "apCAF": r"apCAF|antigen",
    "Fibroblast": r"fibroblast|CAF", "mregDC": r"mreg|LAMP3|mature|migratory",
}
# Primary references of each set (PMIDs; from the project's annotation evidence document).
PRIMARY: dict[str, list[str]] = {
    "Epithelial_tumor": ["32572264", "34238352", "33961783", "36517593"],
    "Fibroblast": ["33981032", "32769974", "33961783", "36517593", "31197017", "32434947"],
    "Pericyte_SMC": ["32769974"], "Endothelial": ["32059779"], "Lymphatic_EC": ["32059779"],
    "T_NK": ["34914499", "30413361", "35549406"], "B_cell": ["33579751", "35549406"], "Plasma": ["33579751", "35549406"],
    "Macrophage_Mono": ["33545035", "35690521", "32783918"], "Dendritic": ["28428369", "32269339"],
    "Neutrophil": ["33545035"], "Mast": ["27135604"], "Mesothelial": ["21984916", "32572264"],
    "Ovarian_stroma": ["31320652", "32123174"], "Adipocyte": ["33981032"], "Schwann": ["35549406"],
    "Proliferating": ["27124452", "34914499"],
    "CSC": ["18519691", "18836486", "19816957", "21498635", "16849428", "31810474", "34645505"],
    "CD8_T": ["34914499"], "CD4_T": ["34914499"], "NK": ["30413361"], "gdT_MAIT": ["34914499", "35549406"],
    "mregDC": ["32269339"], "Macrophage": ["33545035", "35690521"], "Monocyte": ["28428369"], "cDC": ["28428369", "32269339"],
    "pDC": ["28428369"], "Osteoclast_like_giant": ["12479813"],
    "CD8_cytotoxic": ["34914499", "32024970"], "CD8_exhausted": ["34914499", "32024970", "27124452"],
    "CD8_naive_memory": ["34914499", "32024970"], "CD4_naive_memory": ["34914499"],
    "Treg": ["27851913", "34914499"], "Tfh": ["31117010", "34914499"], "CD4_Th_effector": ["34914499"],
    "Naive_B": ["33579751"], "Memory_B": ["33579751", "35549406"], "GC_B": ["35113731", "33579751"], "Plasmablast": ["33579751"],
    "M1_macrophage": ["25035950", "33545035", "35690521"], "M2_macrophage": ["25035950", "33545035", "35690521"],
    "TREM2_TAM": ["35690521", "32783918"], "Hypoxic_TAM": ["27482883", "35690521"],
    "Monocyte_classical": ["28428369"], "Monocyte_nonclassical": ["28428369"],
    "cDC1": ["28428369", "32269339"], "cDC2": ["28428369"],
    "iCAF": ["28232471", "31197017", "32434947"], "myCAF": ["28232471", "31197017", "32434947"], "apCAF": ["31197017"],
}
SOURCE_RANK = {"Experiment": 0, "Single-cell sequencing": 1, "Review": 2}


def candidate_sets(caf_dir: Path | None) -> dict[str, list[str]]:
    """``category/set`` -> candidate genes: the sets of the Xenium OV run plus the CAF lists."""
    out: dict[str, list[str]] = {}
    for k, v in M.LINEAGE_MARKERS.items():
        out[f"lineage/{k}"] = list(v)
    out["proliferation/Proliferating"] = list(M.PROLIF_MARKERS)
    out["csc/CSC"] = list(M.CSC_MARKERS)
    for k, v in M.SPILLOVER_SETS.items():
        out[f"spillover/{k}"] = list(v)
    for comp, cfg in M.IMMUNE_COMPARTMENTS.items():
        for k, v in cfg["level1"].items():
            out[f"{comp}.level1/{k}"] = list(v)
        for parent, states in cfg["level2"].items():
            for k, v in states.items():
                out[f"{comp}.level2.{parent}/{k}"] = list(v)
    if caf_dir is not None:
        for k in ("iCAF", "myCAF", "apCAF"):
            f = caf_dir / f"{k}.txt"
            if f.exists():
                out[f"caf/{k}"] = M.read_marker_file(f)
    return out


def set_name(key: str) -> str:
    name = key.split("/")[1]
    return ALIAS_SET.get(name, name)


def _get(url: str, tries: int = 4, pause: float = 0.4) -> bytes:
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=90) as r:
                data = r.read()
            time.sleep(pause)
            return data
        except Exception:
            if i == tries - 1:
                raise
            time.sleep(2 + 3 * i)
    return b""


def load_cellmarker(cache: Path) -> pd.DataFrame:
    f = cache / "Cell_marker_Human.xlsx"
    if not f.exists():
        f.write_bytes(_get(CELLMARKER_URL))
    d = pd.read_excel(f)
    d = d[d["PMID"].notna() & d["Symbol"].notna()].copy()
    d["PMID"] = d["PMID"].astype(int).astype(str)
    for c in ("cell_name", "Symbol", "marker_source", "Title", "journal", "tissue_type", "cancer_type"):
        d[c] = d[c].astype(str)
    return d


def matches(name: str, text: str) -> bool:
    pat, also = PAT[name]
    if name in EXCLUDE and re.search(EXCLUDE[name], text, re.I):
        return False
    return bool(re.search(pat, text, re.I)) and (also is None or bool(re.search(also, text, re.I)))


def pmc_ids(pmids: list[str], cache: Path) -> dict[str, str]:
    f = cache / "pmid_to_pmcid.json"
    known = json.loads(f.read_text()) if f.exists() else {}
    todo = [p for p in pmids if p not in known]
    for i in range(0, len(todo), 100):
        chunk = todo[i:i + 100]
        rec = json.loads(_get(f"{IDCONV}?ids={','.join(chunk)}&format=json&tool=annotation_st&email=mitrabiolab@gmail.com"))
        got = {str(r.get("pmid")): r.get("pmcid", "") for r in rec.get("records", [])}
        for p in chunk:
            known[p] = got.get(p, "") or ""
    f.write_text(json.dumps(known, indent=1))
    return known


def pmc_sentences(pmcid: str, cache: Path) -> list[str]:
    """Sentences of the body (and figure / table captions) of an open PMC article; [] if no full text."""
    f = cache / f"{pmcid}.xml"
    if not f.exists():
        f.write_bytes(_get(f"{EUTILS}/efetch.fcgi?db=pmc&id={pmcid.replace('PMC', '')}"))
    xml = f.read_text(errors="ignore")
    m = re.search(r"<body.*?</body>", xml, re.S)
    if not m:
        return []
    body = re.sub(r"<ref-list.*?</ref-list>", " ", m.group(0), flags=re.S)
    body = re.sub(r"<xref[^>]*>.*?</xref>", " ", body, flags=re.S)
    text = html.unescape(re.sub(r"<[^>]+>", " ", body))
    text = re.sub(r"\s+", " ", text)
    return [s.strip() for s in re.split(r"(?<=[.;])\s+(?=[A-Z(])", text) if 20 < len(s) < 900]


def pubmed_summaries(pmids: list[str], cache: Path) -> dict[str, dict]:
    f = cache / "pubmed_summaries.json"
    known = json.loads(f.read_text()) if f.exists() else {}
    todo = [p for p in dict.fromkeys(pmids) if p not in known]
    for i in range(0, len(todo), 150):
        chunk = todo[i:i + 150]
        res = json.loads(_get(f"{EUTILS}/esummary.fcgi?db=pubmed&id={','.join(chunk)}&retmode=json"))["result"]
        for p in chunk:
            r = res.get(p, {})
            doi = next((a["value"] for a in r.get("articleids", []) if a.get("idtype") == "doi"), "")
            known[p] = {"title": r.get("title", ""), "journal": r.get("source", ""), "year": str(r.get("pubdate", ""))[:4],
                        "first_author": (r.get("authors") or [{}])[0].get("name", ""), "doi": doi,
                        "found": bool(r) and "error" not in r}
    f.write_text(json.dumps(known, indent=1))
    return known


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--cache", type=Path, required=True, help="download cache (CellMarker table, PMC XML)")
    ap.add_argument("--caf-dir", type=Path, default=None, help="folder with iCAF.txt / myCAF.txt / apCAF.txt candidates")
    ap.add_argument("--max-refs", type=int, default=3, help="references kept per (set, gene)")
    ap.add_argument("--reviewed", type=Path, default=Path(__file__).resolve().parent / "reviewed_literature_evidence.tsv",
                    help="hand-reviewed sentences: columns name (set name), gene, pmid, sentence")
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    a.cache.mkdir(parents=True, exist_ok=True)

    cm = load_cellmarker(a.cache)
    cands = candidate_sets(a.caf_dir)
    all_primary = sorted({p for ps in PRIMARY.values() for p in ps})
    pmc = pmc_ids(all_primary, a.cache)
    texts = {p: pmc_sentences(pmc[p], a.cache) if pmc.get(p) else [] for p in all_primary}

    reviewed = (pd.read_csv(a.reviewed, sep="\t", dtype=str) if a.reviewed.exists()
                else pd.DataFrame(columns=["name", "gene", "pmid", "sentence"]))
    rows, dropped = [], []
    for key, genes in cands.items():
        name = set_name(key)
        cm_set = cm[cm["cell_name"].map(lambda n: matches(name, n))]
        text_pat = re.compile(PAT[name][0] + ("|" + TEXT_EXTRA[name] if name in TEXT_EXTRA else ""), re.I)
        text_not = re.compile(EXCLUDE[name], re.I) if name in EXCLUDE and name.startswith(("CD8", "CD4")) else None
        for gene in genes:
            found = []
            rec = cm_set[cm_set["Symbol"] == gene].copy()
            if len(rec):
                rec["rank"] = rec["marker_source"].map(SOURCE_RANK).fillna(3)
                rec = rec.sort_values(["rank", "year"], ascending=[True, False]).drop_duplicates("PMID")
                for _, r in rec.head(a.max_refs).iterrows():
                    found.append({"evidence": "curated", "pmid": r["PMID"], "detail": f"CellMarker 2.0: {r['cell_name']}"
                                  f" [{r['tissue_type']}; {r['marker_source']}]"})
            if len(found) < a.max_refs:
                gpat = re.compile(rf"(?<![A-Za-z0-9-]){re.escape(gene)}(?![A-Za-z0-9])")
                for p in PRIMARY.get(name, []):
                    if any(f["pmid"] == p for f in found):
                        continue
                    hit = next((s for s in texts.get(p, []) if gpat.search(s) and text_pat.search(s)
                                and not (text_not and text_not.search(s))), None)
                    if hit:
                        found.append({"evidence": "primary_text", "pmid": p, "detail": hit[:400]})
                    if len(found) >= a.max_refs:
                        break
            if len(found) < a.max_refs:
                for _, r in reviewed[(reviewed["name"] == name) & (reviewed["gene"] == gene)].iterrows():
                    if not any(f["pmid"] == r["pmid"] for f in found):
                        found.append({"evidence": "literature_text", "pmid": r["pmid"], "detail": r["sentence"][:400]})
            if found:
                rows += [{"set": key, "gene": gene, **f} for f in found[: a.max_refs]]
            else:
                dropped.append({"set": key, "gene": gene,
                                "reason": "no CellMarker 2.0 record for this cell type, not named with it in a primary "
                                          "reference, and no reviewed sentence from other literature"})
    ref = pd.DataFrame(rows)
    summ = pubmed_summaries(ref["pmid"].tolist(), a.cache)
    for col in ("first_author", "year", "journal", "title", "doi"):
        ref[col] = ref["pmid"].map(lambda p: summ[p][col])
    bad = sorted({p for p in ref["pmid"] if not summ[p]["found"]})
    if bad:
        raise SystemExit(f"PMIDs not found in PubMed: {bad}")
    ref = ref[["set", "gene", "evidence", "pmid", "first_author", "year", "journal", "title", "doi", "detail"]]
    ref.to_csv(a.out / "default_references.tsv", sep="\t", index=False)
    pd.DataFrame(dropped, columns=["set", "gene", "reason"]).to_csv(a.out / "default_dropped.tsv", sep="\t", index=False)
    kept = {k: [g for g in genes if ((ref["set"] == k) & (ref["gene"] == g)).any()] for k, genes in cands.items()}
    (a.out / "default_sets.json").write_text(json.dumps(kept, indent=1) + "\n")
    n_pairs = sum(len(v) for v in cands.values())
    print(f"{len(ref)} references for {sum(len(v) for v in kept.values())}/{n_pairs} set-gene pairs; "
          f"{len(dropped)} dropped; {ref['pmid'].nunique()} distinct PMIDs; full text available for "
          f"{sum(bool(t) for t in texts.values())}/{len(texts)} primary references")


if __name__ == "__main__":
    main()
