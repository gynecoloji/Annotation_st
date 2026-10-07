"""Default marker sets and marker-file helpers.

The gene lists reproduce ``01_xenium_common.py`` of analysis ``01_annot_Xenium_OV``
(Xenium Human 5K Pan Tissue & Pathways panel, human ovarian adenocarcinoma). Every set
is intersected with the panel at run time (:func:`filter_to_panel`), so genes that are
not on a panel are simply dropped and recorded.

All of these are *defaults*: a YAML config can replace any of them (see
:mod:`annotation_st.config`).
"""
from __future__ import annotations

import copy
from pathlib import Path

# ----------------------------------------------------------------------------
# Lineage level (coarse annotation, per Leiden cluster)
# ----------------------------------------------------------------------------
LINEAGE_MARKERS: dict[str, list[str]] = {
    "Epithelial_tumor": ["EPCAM", "PAX8", "WT1", "MUC16", "SOX17", "CLDN7", "CDH1",
                          "MUC1", "FOLR1", "ERBB2", "MSLN"],
    "Fibroblast": ["PDGFRA", "FAP", "THY1", "COL11A1", "COL5A1", "COL5A2", "FBLN2",
                   "FBLN5", "POSTN", "CTHRC1", "SFRP4", "PI16", "MFAP5", "DPT",
                   "LRRC15", "MMP11", "PRRX1", "LOX", "ASPN", "TWIST2", "CDH11"],
    "Pericyte_SMC": ["RGS5", "MCAM", "NOTCH3", "CSPG4", "KCNJ8", "ABCC9", "PDGFRB"],
    "Endothelial": ["PECAM1", "CDH5", "CLDN5", "KDR", "PLVAP", "ENG", "TEK", "ERG",
                    "SOX18"],
    "Lymphatic_EC": ["PROX1", "FLT4", "PDPN", "PECAM1"],
    "T_NK": ["CD3E", "CD3G", "CD2", "CD247", "CD4", "CD8A", "CD8B", "GZMB", "GZMK",
             "GZMA", "PRF1", "KLRD1", "KLRB1", "FOXP3", "IL2RA", "CTLA4", "PDCD1",
             "TIGIT", "CCR7", "SELL"],
    "B_cell": ["MS4A1", "CD79A", "CD79B", "CD19", "BANK1", "PAX5", "CD22"],
    "Plasma": ["MZB1", "XBP1", "PRDM1", "TNFRSF17", "DERL3", "FCRL5", "CD38"],
    "Macrophage_Mono": ["CD68", "CD163", "CD14", "FCGR3A", "MRC1", "MSR1", "TREM2",
                        "CSF1R", "ITGAM", "MARCO", "FOLR2", "SIGLEC1"],
    "Dendritic": ["CD1C", "CLEC10A", "FCER1A", "CLEC9A", "XCR1", "LAMP3", "LILRA4",
                  "CLEC4C", "IL3RA"],
    "Neutrophil": ["CSF3R", "FCGR3B", "CD177"],
    "Mast": ["MS4A2", "HDC", "HPGDS", "GATA2", "KIT"],
    "Mesothelial": ["UPK3B", "CALB2", "PRG4", "PTGIS", "MSLN"],
    "Ovarian_stroma": ["FOXL2", "NR5A1", "STAR", "CYP17A1", "INHA", "AMH", "ARX",
                       "TCF21", "DLK1", "GATA4"],
    "Adipocyte": ["PLIN1", "ADIPOQ", "LEP", "PPARG"],
    "Schwann": ["MPZ", "SOX10"],
}
PROLIF_MARKERS = ["MKI67", "TOP2A", "CDK1", "BIRC5", "CCNB1"]

# ----------------------------------------------------------------------------
# Ovarian cancer stem-like cells
# ----------------------------------------------------------------------------
CSC_MARKERS = ["CD44", "PROM1", "KIT", "LGR5", "SOX2", "POU5F1", "NANOG", "ABCG2",
               "NES", "LIN28A", "LIN28B", "BMI1", "CXCR4"]
CSC_CORE = ["CD44", "PROM1", "KIT", "LGR5", "SOX2", "POU5F1", "NANOG", "ABCG2"]

# ----------------------------------------------------------------------------
# Immune compartments: two-level scheme + spillover / cross-compartment sets
# ----------------------------------------------------------------------------
SPILLOVER_SETS: dict[str, list[str]] = {
    "spillover_Tumor": ["EPCAM", "PAX8", "WT1", "MUC16", "SOX17", "CLDN7", "CDH1", "MUC1", "FOLR1", "MSLN", "BCAM"],
    "spillover_Fibroblast": ["PDGFRA", "FAP", "THY1", "COL11A1", "COL5A1", "COL5A2", "POSTN", "CTHRC1", "SFRP4",
                             "MMP11", "AEBP1", "LTBP2", "THBS2", "CDH11"],
    "spillover_Endothelial": ["PECAM1", "CDH5", "CLDN5", "KDR", "PLVAP", "ENG", "TEK", "ERG", "SOX18", "CD34"],
}
TNK_LEVEL1: dict[str, list[str]] = {
    "CD8_T": ["CD8A", "CD8B"],
    "CD4_T": ["CD4", "CD40LG", "FOXP3", "IL2RA", "CXCR5", "ICOS"],
    "NK": ["KLRD1", "KLRF1", "NCR1", "NCAM1", "FCGR3A", "KLRC1", "KLRK1", "XCL2", "XCL1"],
    "gdT_MAIT": ["TRDC", "TRGC1", "TRGC2", "KLRB1", "IL18R1", "SLC4A10"],
    "B_cell": ["MS4A1", "CD79A", "CD79B", "CD19", "BANK1", "PAX5", "CD22"],
    "Plasma": ["MZB1", "XBP1", "PRDM1", "TNFRSF17", "DERL3", "CD38"],
    "mregDC": ["LAMP3", "CCL19", "FSCN1", "CD83", "WDFY4", "CD274", "CCR7"],
    "Proliferating_T": ["MKI67", "TOP2A", "AURKB", "RRM2", "FOXM1"],
    "cross_Myeloid": ["CD14", "CD68", "CSF1R", "MSR1", "CD163", "SIGLEC1", "ITGAM", "MRC1", "SLCO2B1"],
}
TNK_LEVEL2: dict[str, dict[str, list[str]]] = {
    "CD8_T": {
        "CD8_cytotoxic": ["GZMB", "GZMA", "GZMH", "GZMK", "PRF1", "KLRG1", "CX3CR1", "FGFBP2", "IFNG", "NKG7", "GNLY"],
        "CD8_exhausted": ["PDCD1", "LAG3", "HAVCR2", "TIGIT", "CTLA4", "TOX", "CXCL13", "ENTPD1", "TNFRSF9", "LAYN"],
        "CD8_naive_memory": ["CCR7", "SELL", "TCF7", "LEF1", "IL7R"],
    },
    "CD4_T": {
        "Treg": ["FOXP3", "IL2RA", "CTLA4", "IKZF2", "TNFRSF18", "CCR8", "TNFRSF4"],
        "Tfh": ["CXCL13", "CXCR5", "BCL6", "PDCD1", "ICOS", "IL21"],
        "CD4_Th_effector": ["TBX21", "IFNG", "GATA3", "RORC", "IL17A", "IL2", "CD40LG"],
        "CD4_naive_memory": ["CCR7", "SELL", "TCF7", "LEF1", "IL7R"],
    },
}
MYELOID_LEVEL1: dict[str, list[str]] = {
    "Macrophage": ["CD68", "CD163", "MRC1", "MSR1", "CSF1R", "MARCO", "FOLR2", "TREM2", "SIGLEC1", "STAB1",
                   "SLCO2B1", "C1QA", "C1QB"],
    "Monocyte": ["FCN1", "CCR2", "SELL", "LILRB2", "CDKN1C", "CX3CR1", "S100A8", "S100A9", "VCAN"],
    "cDC": ["CD1C", "CLEC10A", "FCER1A", "CD1E", "CLEC9A", "XCR1", "BATF3", "LAMP3", "CCL19", "FSCN1", "CD83"],
    "pDC": ["LILRA4", "CLEC4C", "IL3RA", "TCF4", "IRF7"],
    "Neutrophil": ["CSF3R", "FCGR3B", "CD177", "CXCR2", "CXCR1", "S100A8"],
    "Mast": ["KIT", "MS4A2", "HDC", "HPGDS", "GATA2", "TPSAB1"],
    "Osteoclast_like_giant": ["ACP5", "CTSK", "NFATC1", "TNFRSF11A", "ATP6V0D2", "OCSTAMP", "DCSTAMP"],
    "Proliferating_myeloid": ["MKI67", "TOP2A", "AURKB", "RRM2", "FOXM1"],
    "cross_T": ["CD3E", "CD3G", "CD2", "CD247", "CD8A", "TRDC"],
}
MYELOID_LEVEL2: dict[str, dict[str, list[str]]] = {
    "Macrophage": {
        "M1_macrophage": ["CD86", "CD80", "IL1B", "TNF", "CXCL9", "CXCL10", "CXCL11", "IDO1", "STAT1", "NOS2",
                          "IL6", "SOCS3", "GBP1", "GBP5", "CCL5"],
        "M2_macrophage": ["CD163", "MRC1", "FOLR2", "CD209", "IL10", "TGFB1", "CCL18", "CCL22", "F13A1", "STAB1",
                          "SELENOP", "LYVE1", "SLC40A1", "GAS6"],
        "TREM2_TAM": ["TREM2", "LIPA", "GPNMB", "APOE", "APOC1", "SPP1", "LGALS3", "FABP5", "LGMN", "CTSL"],
        "Hypoxic_TAM": ["SLC2A1", "HK2", "ADM", "LDHA", "PGK1", "VEGFA", "PLIN2", "BNIP3", "SLC2A3"],
    },
    "Monocyte": {
        "Monocyte_classical": ["CD14", "FCN1", "CCR2", "SELL", "S100A8", "S100A9", "VCAN", "LYZ"],
        "Monocyte_nonclassical": ["FCGR3A", "CX3CR1", "LILRB2", "CDKN1C"],
    },
    "cDC": {
        "cDC1": ["CLEC9A", "XCR1", "BATF3", "IRF8", "CADM1"],
        "cDC2": ["CD1C", "CLEC10A", "FCER1A", "CD1E"],
        "mregDC": ["LAMP3", "CCR7", "CCL19", "FSCN1", "CD274", "CD83"],
    },
}
B_LEVEL1: dict[str, list[str]] = {
    "B_cell": ["MS4A1", "CD79A", "CD79B", "CD19", "BANK1", "PAX5", "CD22", "CIITA"],
    "Plasma": ["MZB1", "XBP1", "PRDM1", "TNFRSF17", "DERL3", "CD38", "IRF4", "TENT5C", "POU2AF1"],
    "Proliferating_B": ["MKI67", "TOP2A", "AURKB", "RRM2", "FOXM1"],
    "pDC": ["LILRA4", "CLEC4C", "IL3RA", "TCF4", "IRF7"],
    "cDC1": ["CLEC9A", "XCR1", "BATF3", "IRF8", "CADM1"],
    "cross_T": ["CD3E", "CD3G", "CD2", "CD247", "CD8A", "TRDC"],
    "cross_Myeloid": ["CD14", "CD68", "CSF1R", "MSR1", "CD163", "SIGLEC1", "ITGAM", "MRC1", "SLCO2B1"],
}
B_LEVEL2: dict[str, dict[str, list[str]]] = {
    "B_cell": {
        "Naive_B": ["TCL1A", "FCER2", "IL4R", "BACH2", "SELL", "CCR7", "IGHD", "IGHM"],
        "Memory_B": ["CD27", "TNFRSF13B", "CD80", "CD86", "ITGAX", "TBX21", "FCRL4", "FCRL5", "AIM2"],
        "GC_B": ["BCL6", "AICDA", "LMO2", "SUGCT", "SERPINA9", "CD38", "MKI67", "RGS13"],
    },
    "Plasma": {"Plasmablast": ["MKI67", "TOP2A", "AURKB", "RRM2"]},
}

# Compartment definitions. Processing order: compartments with a ``handoff`` run first so
# the receiving compartment gets the handed-off cells (they keep their identity there).
IMMUNE_COMPARTMENTS: dict[str, dict] = {
    "TNK": {
        "lineages": ["T_NK"], "level1": TNK_LEVEL1, "level2": TNK_LEVEL2, "core": ["CD8_T", "CD4_T", "NK"],
        "unresolved": "T_NK_unresolved", "cross": {"cross_Myeloid": "Macrophage_Mono"},
        "unspecified": {"CD8_T": "CD8_T_unspecified", "CD4_T": "CD4_T_unspecified"},
        "handoff": {"B_cell": "Bcell", "Plasma": "Bcell"},
        "coarse": {},
    },
    "Bcell": {
        "lineages": ["Plasma"], "level1": B_LEVEL1, "level2": B_LEVEL2, "core": ["Plasma"],
        "unresolved": "B_lineage_unresolved", "cross": {"cross_T": "T_NK", "cross_Myeloid": "Macrophage_Mono"},
        "unspecified": {"B_cell": "B_cell", "Plasma": "Plasma"},
        "handoff": {},
        "coarse": {"B_cell": "B_cell", "Proliferating_B": "B_cell", "Plasma": "Plasma",
                   "pDC": "Dendritic", "cDC1": "Dendritic", "cross_T": "T_NK", "cross_Myeloid": "Macrophage_Mono"},
    },
    "Myeloid": {
        "lineages": ["Macrophage_Mono"], "level1": MYELOID_LEVEL1, "level2": MYELOID_LEVEL2,
        "core": ["Macrophage"], "unresolved": "Myeloid_unresolved", "cross": {"cross_T": "T_NK"},
        "unspecified": {"Macrophage": "Macrophage_unpolarized", "Monocyte": "Monocyte_unspecified",
                        "cDC": "cDC_unspecified"},
        "handoff": {},
        "coarse": {},
    },
}


def ov_run_sets() -> dict:
    """The gene sets of the Xenium OV run (2026-09-08) in the shape of the ``markers`` config
    section (a deep copy). These are the *candidates* of the default library, not the library:
    ``annotation_st.library.default_library()`` keeps only the genes that carry a reference."""
    return copy.deepcopy({
        "lineage": LINEAGE_MARKERS,
        "proliferation": PROLIF_MARKERS,
        "csc": CSC_MARKERS,
        "csc_core": CSC_CORE,
        "spillover": SPILLOVER_SETS,
        "compartments": IMMUNE_COMPARTMENTS,
        "caf": {},
    })



# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------
def read_marker_file(path: str | Path) -> list[str]:
    """One gene per line; blank lines and ``#`` comments ignored; order kept, duplicates dropped."""
    genes: list[str] = []
    for line in Path(path).read_text().splitlines():
        g = line.strip()
        if g and not g.startswith("#") and g not in genes:
            genes.append(g)
    return genes


def load_marker_files(mapping: dict[str, str | Path], base: str | Path | None = None) -> dict[str, list[str]]:
    """``{set_name: path}`` -> ``{set_name: genes}``. Relative paths resolve against ``base``."""
    out = {}
    for name, p in mapping.items():
        p = Path(p)
        if base is not None and not p.is_absolute():
            p = Path(base) / p
        out[name] = read_marker_file(p)
    return out


def filter_to_panel(sets: dict[str, list[str]], var_names) -> tuple[dict[str, list[str]], dict[str, list[str]]]:
    """Intersect each set with ``var_names``. Returns ``(kept, dropped)`` with the same keys."""
    present = set(map(str, var_names))
    kept = {k: [g for g in v if g in present] for k, v in sets.items()}
    dropped = {k: [g for g in v if g not in present] for k, v in sets.items()}
    return kept, dropped
