#!/usr/bin/env python
"""Candidate sentences for set-gene pairs that have neither a curated record nor a primary-text hit.

For each pair in default_dropped.tsv this script searches Europe PMC (open-access full texts) for
papers that mention the gene (or a listed protein alias) together with the cell type, fetches the
PubMed Central text and collects sentences in which gene and cell type co-occur with a marker cue.
The output is a CANDIDATE list for human review; only sentences copied into
tools/reviewed_literature_evidence.tsv become references of the default library.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import urllib.parse
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_marker_references import _get, pmc_sentences, set_name  # noqa: E402

EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
# set -> (query phrases, sentence pattern)
CONTEXT = {
    "Epithelial_tumor": (["ovarian cancer", "high-grade serous"], r"ovarian|serous|epitheli|tumou?r cell|cancer cell"),
    "Fibroblast": (["cancer-associated fibroblasts", "fibroblasts"], r"fibroblast|CAF"),
    "Mesothelial": (["mesothelial cells"], r"mesotheli"),
    "Ovarian_stroma": (["granulosa cells", "theca cells", "ovarian stroma"], r"granulosa|theca|ovarian strom"),
    "Adipocyte": (["adipocytes"], r"adipocyte"),
    "Proliferating": (["proliferating cells", "cycling cells", "cell cycle"], r"prolifer|cycling|G2/M|cell.cycle|mitos"),
    "mregDC": (["mregDC", "LAMP3+ dendritic cells", "mature dendritic cells"], r"mreg|LAMP3|mature|migratory|dendritic"),
    "Macrophage_Mono": (["macrophages"], r"macrophage|monocyte"),
    "Macrophage": (["macrophages"], r"macrophage"),
    "CD4_Th_effector": (["T helper cells", "Th1"], r"helper|Th1|Th2|Th17|effector"),
    "Plasma": (["plasma cells"], r"plasma cell|plasmablast"),
    "Plasmablast": (["plasmablasts"], r"plasmablast"),
    "Naive_B": (["naive B cells"], r"na[iï]ve B"),
    "Memory_B": (["memory B cells"], r"memory B"),
    "Osteoclast_like_giant": (["osteoclasts", "giant cells"], r"osteoclast|giant cell"),
    "M1_macrophage": (["M1 macrophages", "interferon macrophages"], r"\bM1\b|IFN|interferon|inflammatory macro|classically"),
    "M2_macrophage": (["M2 macrophages", "alternatively activated macrophages", "resident macrophages"], r"\bM2\b|alternativ|resident|perivascular|anti-?inflam"),
    "TREM2_TAM": (["TREM2 macrophages", "lipid-associated macrophages"], r"TREM2|lipid.associated|\bLAM|TAM"),
    "Hypoxic_TAM": (["hypoxic macrophages", "tumor-associated macrophages hypoxia"], r"hypoxi"),
    "Monocyte_classical": (["classical monocytes"], r"classical monocyte|CD14\+"),
    "Monocyte_nonclassical": (["non-classical monocytes", "nonclassical monocytes"], r"non-?classical"),
    "myCAF": (["myCAF", "myofibroblastic CAFs"], r"myCAF|myofibroblast"),
    "apCAF": (["apCAF", "antigen-presenting CAFs"], r"apCAF|antigen.present"),
}
ALIASES = {"MSLN": ["mesothelin"], "CALB2": ["calretinin"], "NR5A1": ["SF-1", "SF1", "steroidogenic factor 1"],
           "INHA": ["inhibin"], "PPARG": ["PPARγ", "PPAR-γ"], "CCNB1": ["cyclin B1"], "AURKB": ["aurora B", "Aurora kinase B"],
           "SLC2A1": ["GLUT1", "GLUT-1"], "SLC2A3": ["GLUT3"], "HK2": ["hexokinase 2"], "VEGFA": ["VEGF"],
           "ADM": ["adrenomedullin"], "NFATC1": ["NFATc1"], "TENT5C": ["FAM46C"], "SELENOP": ["SEPP1"],
           "SLC40A1": ["ferroportin"], "LGALS3": ["galectin-3", "Gal-3"], "LGMN": ["legumain"], "CTSL": ["cathepsin L"],
           "IDO1": ["IDO"], "SELL": ["CD62L", "L-selectin"], "IL2": ["IL-2"], "F13A1": ["FXIIIA", "factor XIIIa"],
           "CDH11": ["cadherin-11", "cadherin 11"], "CD80": ["B7-1"], "LYVE1": ["LYVE-1"]}
CUE = re.compile(r"marker|express|signature|characteri[sz]|defin|identif|high|upregulat|enrich|positive|\+|specific", re.I)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dropped", type=Path, required=True)
    ap.add_argument("--cache", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--max-papers", type=int, default=8)
    ap.add_argument("--skip-sets", nargs="*", default=["caf/iCAF"])
    a = ap.parse_args()
    dr = pd.read_csv(a.dropped, sep="\t", dtype=str)
    dr = dr[~dr["set"].isin(a.skip_sets)]
    dr["name"] = dr["set"].map(set_name)
    pairs = dr.drop_duplicates(["name", "gene"])
    scache = a.cache / "epmc_search.json"
    searches = json.loads(scache.read_text()) if scache.exists() else {}
    rows = []
    for _, r in pairs.iterrows():
        name, gene = r["name"], r["gene"]
        if name not in CONTEXT:
            continue
        phrases, spat = CONTEXT[name]
        terms = [gene] + ALIASES.get(gene, [])
        q = "(" + " OR ".join(f'"{t}"' for t in terms) + ") AND (" + " OR ".join(f'"{p}"' for p in phrases) + \
            ") AND OPEN_ACCESS:y AND IN_EPMC:y AND SRC:MED"
        if q not in searches:
            url = f"{EPMC}?query={urllib.parse.quote(q)}&format=json&pageSize=25&sort=CITED%20desc&resultType=lite"
            res = json.loads(_get(url)).get("resultList", {}).get("result", [])
            searches[q] = [{"pmid": x.get("pmid", ""), "pmcid": x.get("pmcid", ""), "cited": x.get("citedByCount", 0),
                            "title": x.get("title", ""), "year": x.get("pubYear", ""), "journal": x.get("journalTitle", "")}
                           for x in res if x.get("pmcid") and x.get("pmid")]
            scache.write_text(json.dumps(searches, indent=1))
        gpat = re.compile(r"(?<![A-Za-z0-9-])(" + "|".join(re.escape(t) for t in terms) + r")(?![A-Za-z0-9])")
        cpat = re.compile(spat, re.I)
        n = 0
        for hit in searches[q][: a.max_papers * 2]:
            try:
                sents = pmc_sentences(hit["pmcid"], a.cache)
            except Exception:
                continue
            good = [s for s in sents if gpat.search(s) and cpat.search(s) and CUE.search(s) and len(s) < 600]
            if good:
                best = min(good, key=len) if len(good) > 3 else good[0]
                rows.append({"name": name, "gene": gene, "pmid": hit["pmid"], "pmcid": hit["pmcid"], "cited": hit["cited"],
                             "year": hit["year"], "journal": hit["journal"], "title": hit["title"], "sentence": best})
                n += 1
            if n >= 3:
                break
    out = pd.DataFrame(rows)
    out.to_csv(a.out, sep="\t", index=False)
    got = out.drop_duplicates(["name", "gene"]).shape[0] if len(out) else 0
    print(f"{len(pairs)} pairs searched; candidates for {got}; {len(out)} sentences -> {a.out}")


if __name__ == "__main__":
    main()
