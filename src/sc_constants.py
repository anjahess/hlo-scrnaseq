"""

Constants

@author: Anja Hess
@date: 2022-OCT-01

"""

import os
import seaborn as sns
import scanpy as sc

#########################################################################
# 1. PATHS ETC
#########################################################################
SCRIPT_PATH = str(os.path.dirname
                  (os.path.abspath(__file__)))
SCRIPT_PATH = SCRIPT_PATH.rsplit('/src')[0]
cell_cycle_file = f"{SCRIPT_PATH}/config/other/" \
                  f"regev_lab_cell_cycle_genes.txt"
EXTERNAL_SAVE_PATH = "/where/to/save"
MATRIXPATH = f"{EXTERNAL_SAVE_PATH}data/raw/datasets/"
DATA_PATH = f"{EXTERNAL_SAVE_PATH}data/"
TENX_scheme = ["barcodes.tsv.gz", "features.tsv.gz",
               "matrix.mtx.gz", 'barcodes.tsv',
               'features.tsv', 'matrix.mtx']
script_titles = ["start.py", "core_functions.py",
                 "lazy_functions.py", "sc_constants.py",
                 "data_parsing.py", "plots.py",
                 "cell_qc.py", "annotate.py"]
SEARCHTERMS = ["WikiPathway_2021_Human"]
REPLOT = True
RECLUSTER = False
MIN_GENES = 500
MIN_CELLS = 10
COUNT_LOWER, COUNT_UPPER = 500, 30000
MITO_CUTOFF = 0.2
RIBO_UPPER, RIBO_LOWER = 0.4, 1
DOUBLET_SCORE = "off"
STATSTEST = "wilcoxon"
GENES_BY_COUNT_FRAME = [0, 4000]  # Deprecated, just naming
RES = [0.1]  # leiden resolutions
keys_to_plot = [f"leiden_0.1"]
H5AD_DATA_PATH = f"{os.getcwd()}/data/"

#########################################################################
# 2. SCANPY PARAMS
#########################################################################
n_highest_expressed_to_show = 20
v_min, v_max = -4, 4
sc.settings.file_format_figs = 'pdf'
sc.settings.dpi_save = 300
sc.settings.autosave = True
sc.settings.n_jobs = 10
sc.settings.set_figure_params(dpi=400,
                              fontsize=20)

#########################################################################
# 3. ANNOTATION GENE LISTS
#########################################################################
CELLMARKER_PATH = f"{SCRIPT_PATH}/config/cell_marker/"
GO_PATH = f"{SCRIPT_PATH}/config/GO_terms/GO"
annoset2path = {
    "sctype": f"{CELLMARKER_PATH}"
              f"sctype.csv",
    "literature": f"{CELLMARKER_PATH}"
                  f"literature.csv",
    "wesley_pairwise": f"{CELLMARKER_PATH}"
                       f"wesley_pairwise.csv",
}
NMARKERS = [200]
MAIN_ANNO = "leiden_0.12sctype_200"

# Hepatocyte subclustering
TARGET_HEPATO = "leiden_0.42wesley_pairwise_200"
RES_HEPATO = "leiden_0.4"
ANNO_DICT_HEPATO = {"wesley_pairwise": f"{CELLMARKER_PATH}"
                                       f"wesley_pairwise.csv"}

# Cholangiocyte subclustering
TARGET_CHOL = "leiden_0.32sctype_200_simple"
RES_CHOL = "leiden_0.3"
ANNO_DICT_CHOL = {"sctype": f"{CELLMARKER_PATH}sctype.csv"}

#########################################################################
# 3. GO TERMS
#########################################################################
geneset2path = {
    "98-gene_signature": f"{SCRIPT_PATH}"
                         f"/config/signatures/"
                         f"98-gene_signature.txt",
    "26-gene_signature": f"{SCRIPT_PATH}"
                         f"/config/signatures/"
                         f"26-gene_signature.txt",
    "Fibroblast_act": f"{GO_PATH}0072537_"
                      f"fibroblast-activation",
    "Fibroblast_prol": f"{GO_PATH}0048144_"
                       f"fibroblast-proliferation",
    "Pos_ecm_orga": f"{GO_PATH}1903055_positive"
                    f"-regulation-of-extracellular"
                    f"-matrix-organization",
    "Pos_hep_stell_act": f"{GO_PATH}2000491_"
                         f"positive-regulation-"
                         f"of-hepatic-stellate-"
                         f"cell-activation",
    "Positive-regulation-of-collagen-fibril"
    "-organization": f"{GO_PATH}1904028_positive-"
                     f"regulation-of-collagen"
                     f"-fibril-organization",
    "Smooth_muscle_"
    "contractile_fiber": f"{GO_PATH}0030485_"
                         f"smooth-muscle-"
                         f"contractile-fiber",
    "Fibroblast-migration": f"{GO_PATH}0010761_"
                            f"fibroblast-migration",
    "Connective-tissue"
    "-response-to-"
    "inflammation": f"{GO_PATH}0002248_connective-"
                    f"tissue-response-to-"
                    f"inflammation",
    "Activation_of_innate"
    "_immune_response": f"{GO_PATH}0002218_activation"
                        f"-of-innate-immune-response",
    "Acute-"
    "inflammatory-response": f"{GO_PATH}0002526_acute"
                             f"-inflammatory-response",
    "Chronic-inflammatory-response": f"{GO_PATH}0002544_"
                                     f"chronic-inflammatory"
                                     f"-response",
    "Chemokine-activity": f"{GO_PATH}0008009_"
                          f"chemokine-activity",
    "Chemokine-production": f"{GO_PATH}0032602_"
                            f"chemokine-production",
    "Hepatic-immune-response": f"{GO_PATH}0002384"
                               f"_hepatic"
                               f"-immune-response",
    "Lps-immune-"
    "receptor-activity": f"{GO_PATH}0001875"
                         f"_lps-immune-"
                         f"receptor-activity",
}
#########################################################################
# 3.1 COLLECTION FIBROSIS
#########################################################################
geneset2path_fibrosis = {
    "Fibroblast-migration": f"{GO_PATH}0010761_fibroblast-migration",
    "Connective-tissue-response-to-inflammation":
        f"{GO_PATH}0002248_connective-tissue-response-to-inflammation",
    "Positive-regulation-of-collagen-fibril-organization":
        f"{GO_PATH}1904028_positive-regulation-of-collagen-"
        f"fibril-organization",
    "Pos_ecm_orga": f"{GO_PATH}1903055_positive-regulation-of-"
                    f"extracellular-matrix-organization",
    "Fibroblast_act": f"{GO_PATH}0072537_fibroblast-activation",
    "Fibroblast_prol": f"{GO_PATH}0048144_fibroblast-proliferation",
    "Pos_hep_stell_act": f"{GO_PATH}2000491_positive"
                         f"-regulation-of-hepatic-stellate"
                         f"-cell-activation",
    "Smooth_muscle_contractile_fiber": f"{GO_PATH}"
                                       f"0030485_smooth-muscle"
                                       f"-contractile-fiber"
}
#########################################################################
# 3.1 COLLECTION INFLAMMATION
#########################################################################
geneset2path_inflammation = {
    "Chemokine-activity": f"{GO_PATH}0008009_chemokine-activity",
    "Hepatic-immune-response": f"{GO_PATH}0002384_hepatic"
                               f"-immune-response",
    "Lps-immune-receptor-activity": f"{GO_PATH}0001875_lps"
                                    f"-immune-receptor-activity",
    "Acute-inflammatory-response": f"{GO_PATH}0002526_acute"
                                   f"-inflammatory-response",
    "Chronic-inflammatory-response": f"{GO_PATH}0002544_"
                                     f"chronic-inflammatory"
                                     f"-response",
    "Chemokine-production": f"{GO_PATH}0032602_chemokine-production",
    "Activation_of_innate_immune_response": f"{GO_PATH}"
                                            f"0002218_activation"
                                            f"-of-innate-immune"
                                            f"-response",
    "Activation_of_immune_response": f"{GO_PATH}0002253_"
                                     f"activation-"
                                     f"of-immune-response",
}

#########################################################################
# 4. COLOR SCHEMES
#########################################################################
cond_palette = ["#bfcfcd", "#fcd2a1", "#477b80", "#d56763",
                "#2d435b", "#f1e8d7", "#eacdcb", "lightblue",
                "#fbc27b", "cadetblue", "#fbc27b", 'lightslategrey',
                "#85ada3", "#d56763", "#fcd2a1", "#477b80",
                "#eacdcb", "#bfcfcd", "#2d435b", "#986960",
                "#f1e8d7", "#d56763", "#fcd2a1", "#477b80", "#bfcfcd",
                "#d56763", "#fcd2a1", "#477b80", "#2d435b", "#477b80",
                "#2d435b",
                "#986960", "#f1e8d7", "#d56763", "#fcd2a1", "#477b80",
                "lightblue", "#fbc27b", "cadetblue", "#fbc27b",
                "lightslategrey", "#85ada3", "#d56763", "#fcd2a1",
                "#477b80", "#eacdcb", "#bfcfcd","#2d435b", "#986960",
                "#f1e8d7", "#d56763", "#fcd2a1",
                "#477b80", "black", "brown", "blue",
                "orange", "yellow", "pink"]

palettediff2 = ["#d56763", "#fcd2a1", "#bfcfcd", "#477b80",
                "#2d435b","#f1e8d7", "#eacdcb", 'grey',
                "lightblue", "#fbc27b", "cadetblue",
                "#fbc27b", 'lightslategrey', "#85ada3",
                "#d56763", "#fcd2a1", "#eacdcb", "#bfcfcd",
                "#986960", "#d56763"]

#########################################################################
# 5. DICTS
#########################################################################
CELL_HYPERCATS = {
    'Liver progenitor cell': '0_Hepatocytes',
    'progenitor': '1_Hepatocytes',
    'Hepatocytes': '1_Hepatocytes',
    'HB1': '11_Hepatocytes',
    'HB2': '2_Hepatocytes',
    'FH1': "3_Hepatocytes",
    'FH2': '4_Hepatocytes',
    'cAH': '4_Hepatocytes',
    'AH': '5_Hepatocytes',
    'Cholangiocytes': '6_Cholangiocytes',
    'CHOLs': '6_Cholangiocytes',
    'Ductal cells': "7_Cholangiocytes",
    'DCs': "7_Cholangiocytes",
    'Hepatic stellate cells': '8_Fibroblasts',
    'HSCs': '8_Fibroblasts',
    'Activated hepatic stellate cells': '9_Fibroblasts',
    'HSCs-ACT': '9_Fibroblasts',
    'Fibroblasts': '99_Fibroblasts',
    'FIBs': '99_Fibroblasts',
    'Myofibroblasts': '999_Fibroblasts',
    'Smooth': '9999_Fibroblasts',
    'Smooth muscle cells': '9999_Fibroblasts',
    'SMCs': '9999_Fibroblasts',
}

cell2abbr = {"HB1": "HB1",
             "Hepatic stellate cells": "HSCs",
             "HB2": "HB2",
             "Fibroblasts": "FIBs",
             "Cholangiocytes": "CHOLs",
             "FH1": "FH1",
             "Smooth muscle cells": "SMCs",
             "Activated hepatic stellate cells": "HSCs-ACT",
             "FH2": "cAH",
             "Ductal cells": "DCs",
             "AH": "AH"}

CELL_COLS = {"KRT19": ["#d56763",
                       "#d25f5b", "#dd8481", "#e7aaa8",
                       "#f2d0ce"],
             "epato": ["#2d435b", "#436488", "#537ca9",
                       "#bcd4e6", "#bccae6", "#c3bce6",
                       "#bccae6", "#e4e1f4",
                       "#9184d1", "#ffb6c1",
                       "#d6dfea", "#ebb5b3"
                       ],
             "AH": ["#2d435b", "#436488", "#537ca9",
                    "#bcd4e6", "#bccae6", "#c3bce6",
                    "#bccae6", "#e4e1f4",
                    "#9184d1", "#ffb6c1",
                    "#d6dfea", "#ebb5b3"
                    ],

             "HB1": ["#ffb6c1",
                     "#d6dfea", "#ebb5b3"],
             "HB2": ["#bcd4e6"],
             "FH": ["#6b9ac8", "#96a2ae", "#3e74a8"],
             "Enteroendocrine": ["#76a3a3"],
             "yoepithelial": ["#477b80", "darkgreen", "#8b4513"],
             "epatoblast": ["#ffb6c1", "#d6dfea",
                            "#ebb5b3"],
             "holangio": ["#d56763", "#a02f2b", "#c73a35",
                          "#d25f5b", "#dd8481", "#e7aaa8",
                          "#f2d0ce"],
             "CHOLs": ["#a02f2b"],
             "ibroblast": ["#fcd2a1", "#965304", "#c56d06",
                           "#f58707", "#f89e35", "#fab564",
                           "#fbcb94", "#fde2c3"],
             "stellate": ["#fcd2a1", "#965304", "#c56d06",
                          "#f58707", "#f89e35", "#fab564",
                          "#fbcb94", "#fde2c3"],
             "HSCs-ACT": ["orange", "yellow"],
             "HSCs": ["#fcd2a1", "#965304", "#c56d06",
                      "#f58707", "#f89e35", "#fab564",
                      "#fbcb94", "#fde2c3"],
             "FIBs": ["black"],
             "ancreatic": ["#b5c7da", "#d6dfea", "#5d9b9b",
                           "#82b4b4"],
             "stem": ["#bfcfcd", "#a7c6c9", "#edd5d3", "#dfb3b0",
                      "#d1928d"],
             "uscle": ["#85ada3", "#477b80", "#fbc27b",
                       "#663a03", "#965504"],
             "SMCs": ["#85ada3", "#477b80", "#fbc27b",
                      "#663a03", "#965504"],
             "emature": ["#b5c7da", "#d6dfea", "#5d9b9b",
                         "#82b4b4"],
             "progenitor": ["#b5c7da", "#d6dfea", "#5d9b9b",
                            "#82b4b4"],
             "ycling": ["grey", "lightgrey", "#d2e5e5"],
             "Ductal": ["#d56763", "#a02f2b", "#c73a35",
                        "#d25f5b", "#dd8481", "#e7aaa8",
                        "#f2d0ce"],
             "DCs": ["#d56763", "#a02f2b", "#c73a35", ],
             }

GROUP_COLS = {"CTRL": "lightgrey",
              "CTRL-OA": "lightgrey",
              "CTRLOA": "lightgrey",
              "CTRL-PA": "lightgrey",
              "CTRL-TGFB1": "lightgrey",
              "PA": "#477b80",
              "TGFB1": "#d56763",
              "OA500": "#fcd2a1",
              }

cmap = sns.diverging_palette(220, 20, as_cmap=True)

#########################################################################
# 6. MARKER GENES
#########################################################################
GENERAL_MARKERS = ["AFP", "ALB", "HNF4A",
                   "CEBPA", "HP", "SAA1",
                   "TAT", "F9", "CYP2C9",
                   "AZGP1", "APCS", "HPD",
                   "GPC3", "CES1", "CYP8B1",
                   "CYP2C92", "KRT19", "KRT7",
                   "VIM", "SPARC", "COL3A1",
                   "COL1A1"]

FIBROSIS_MARKERS = ["COL1A1", "COL3A1", "TAGLN2",
                    "TGFBI", "PDGFA", "PDGFB",
                    "ACTG1", "COL4A1", "COL5A1",
                    "COL6A1", "S100A11"]

INFLAMMATION_MARKERS = ["NFKBIA", "CX3CL1", "CXCL1",
                        "CD24", "IL1A", "IL32",
                        "CCL20", "COX7A2", "TNFAIP3",
                        "NFKB2", "S100A4"]

#########################################################################
# 7. CONDITIONS
#########################################################################
DEFAULT_CONTROL = "CTRL"
sample2condition = {
                    "50": "CTRL",
                    "51": "CTRL",
                    "45": "Day21",
                    "46": "Day21",
                    "47": "Day21",
                    "48": "Day34",
                    "49": "Day34",
                    "52": "TGFB1",
                    "53": "TGFB1",
                    "62": "CTRL",
                    "63": "CTRL",
                    "64": "PA",
                    "65": "PA",
                    "SM-KBV62": "CTRL",
                    "SM-KBV63": "CTRL",
                    "SM-KBV64": "PA",
                    "SM-KBV65": "PA",
                    "SM-L3XWE": "CTRL",
                    "SM-L3XWF": "CTRL",
                    "SM-L3XWG": "CTRLOA",
                    "SM-L3XWH": "CTRLOA",
                    "SM-L3XWI": "OA500",
                    "SM-L3XWJ": "OA500",
                    }

sample2controlstatus = {
    "45": "CTRL",
    "46": "CTRL",
    "47": "CTRL",
    "48": "CTRL",
    "49": "CTRL",
    "50": "CTRL",
    "51": "CTRL",
    "52": "TGFB1",
    "53": "TGFB1",
    "62": "CTRL",
    "63": "CTRL",
    "64": "PA",
    "65": "PA",
    "SM-KBV62": "CTRL",
    "SM-KBV63": "CTRL",
    "SM-KBV64": "PA",
    "SM-KBV65": "PA",
    "SM-L3XWE": "CTRL",
    "SM-L3XWF": "CTRL",
    "SM-L3XWG": "CTRL",
    "SM-L3XWH": "CTRL",
    "SM-L3XWI": "OA500",
    "SM-L3XWJ": "OA500",
}

sample2culture = {"50": "OS",
                  "51": "OS",
                  "45": "ULA",
                  "46": "ULA",
                  "47": "ULA",
                  "48": "ULA",
                  "49": "ULA",
                  "52": "OS",
                  "53": "OS",
                  "62": "OS",
                  "63": "OS",
                  "64": "OS",
                  "65": "OS",
                  "SM-KBV62": "OS",
                  "SM-KBV63": "OS",
                  "SM-KBV64": "OS",
                  "SM-KBV65": "OS",
                  "SM-L3XWE": "OS",
                  "SM-L3XWF": "OS",
                  "SM-L3XWG": "OS",
                  "SM-L3XWH": "OS",
                  "SM-L3XWI": "OS",
                  "SM-L3XWJ": "OS",
                  }

hep2loc = {
    "Periportal_C5": 0,
    "Periportal_C14": 1,
    "Periportal_C6": 2,
    "Pericentral_C15": 3,
    "Pericentral_C1": 4,
    "Pericentral_C3": 5
}

sample2detailedcondition = {
    "50": "CTRL-TGFB1",
    "51": "CTRL-TGFB1",
    "45": "Day21",
    "46": "Day21",
    "47": "Day21",
    "48": "Day34",
    "49": "Day34",
    "52": "TGFB1",
    "53": "TGFB1",
    "62": "CTRL-PA",
    "63": "CTRL-PA",
    "64": "PA",
    "65": "PA",
    "SM-L3XWE": "CTRL",
    "SM-L3XWF": "CTRL",
    "SM-L3XWG": "CTRL-OA",
    "SM-L3XWH": "CTRL-OA",
    "SM-L3XWI": "OA500",
    "SM-L3XWJ": "OA500",
}

sample2diffcond = {
}

diff2controlcondition = {
    "PA": "CTRL-PA",
    "OA500": "CTRL-OA",
    "TGFB1": "CTRL-TGFB1"
}

sample2controlcondition = {
    "45": "CTRL",
    "46": "CTRL",
    "47": "CTRL",
    "48": "CTRL",
    "49": "CTRL",
    "50": "CTRL",
    "51": "CTRL",
    "52": "CTRL-TGFB1",
    "53": "CTRL-TGFB1",
    "62": "CTRL",
    "63": "CTRL",
    "64": "CTRL-PA",
    "65": "CTRL-PA",
    "SM-L3XWE": "CTRL",
    "SM-L3XWF": "CTRL",
    "SM-L3XWG": "CTRL",
    "SM-L3XWH": "CTRL",
    "SM-L3XWI": "CTRL-OA",
    "SM-L3XWJ": "CTRL-OA",
}

#########################################################################
# QC, clustering
#########################################################################
GENECLASSES = {"MT-": "mito", "LNC": "lnc",
               ("RPS", "RPL"): "ribo",
               "MIR": "mir"}

umap_vars = ["culture", "condition", "sample",
             "phase", "total_counts",
             "n_genes_by_counts", "doublet_score",
             "pct_ribo", "pct_mito", "pct_mir",
             ]

genesets_lineage = {
    "HSC": ["COL1A2", "COL3A1",
            "COL6A1", "SPARC", "DES",
            "COL1A1"],
    "Adult_Hepatocyte": ["CES1", "CYP2C9"],
    "Portal_Hepatocyte:": ["HSD17B13"],
    "Premature_hepatocyte": ["AFP", "APOA1", "CEBPA"],
    "Proliferating": ["CDK1", "CENPF", "TYMS", "TOP2A"],
    "CHOL": ["KRT7", "KRT19", "CLDN4", ],
    "HB /FH": ["APOA1", "AFP", "HNF4A", "CEBPA", "APOM"],
}

genesets_plots = {
    "Embryonic_DE": ["NANOG", "POU5F1",
                      "UTF1", "SOX2", "SOX17"],
    "FH": ["IER2", "BTG2", "LIPC",
           "GC"],
    "Cholangiocyte": ["KRT19", "CLDN4", "KRT7", "EPCAM",
                      "TACSTD2", "CD24", "KRT8", "EPCAM", "CLDN4"],
    "Premature_hepatocyte": ["AFP", "APOA1", "HNF4A", "CEBPA", "APOA1", "APOM",
                             "CYP2C8", "CYP2C9", "CYP2C18", "CYP2A6"
                             ],
    "Adult_Hepatocyte": ["CYP2E1", "ADH1C", "FGL1", "CES1",
                         "APCS", "CFI", "APOF"
                         ],
    "Proliferating": ["CDK1", "CENPF", "TYMS", "TOP2A"],
    "(Myo)fibroblast": ["TAGLN", "MYL9", "TGFBI",
                        "ACTA2", "LOX"],
    "HSCs": ["SPARC", 'COL1A1', "COL3A1", 'COL1A2',
             'COL6A1', "DES"],
}

genesets_plots_big = {
    "Embryonic_DE": ["NANOG", "POU5F1",
                     "UTF1", "SOX2", "SOX17"],
    "HSC": ["SPARC", "PPARG", 'COL1A1',
            "COL3A1", 'COL6A1', "IGFBP3",
            "MYL9", "PDGFRB",
            "TGFBI", "ACTA2", "DES", ],
    "Adult_Hepatocyte": ["ADH1C", "FGL1", "SERPING1", "CLU",
                         "CFI", "NFIA",
                         "GSTM3", "CHST9", "GAL3ST1", "UGT2B15", "F5",
                         "SULT1E1", ],
    "Proliferating": ["CDK1", "CENPF", "TYMS",
                      "TOP2A"],
    "Cholangiocyte": ["KRT7", "KRT19", "CLDN4", "S100A14", "EPCAM",
                      "TACSTD2", "PLPP2",
                      "CHST4", ],
    "Myo_HSCa ": ["ITGB1", "MMP14",
                  "IGFBP7", "ACKR3", "S100A11", "VIM",
                  "TAGLN2", "FN1", "ITGAV", "IGFBP2",
                  "ITGB1", "MMP7", "TNFRSF12A"],
    "Fetal_hepatocyte": ["AHSG", "AFP", "APOA1",
                         "APOC1", "HNF4A", "CEBPA",
                         "ALB"],

}

# CAVE: for whatever reason max ~12 genes, otherwise numpy error
PALANTIR_GENES = {
    "Progenitor": ['MKI67', 'NANOG',
                   'SOX2', 'POU5F1'],
    "HB": ['AFP', 'ALB', 'AMN',
           'APOA1', 'APOB', 'APOM',
           'CEBPA', 'CYP2C9', 'C3'],
    "FH": ['F2', 'HNF4A', 'BAAT', 'VTN',
           "CEBPA", "APCS", "CES1",
           "KRT19", "KRT7", "APOC1"],
    "CTRL_HEP": ["GC", 'CES1', 'APCS',
                 'GATA4', 'GPAT4', 'GHR'],
    "QHSC": ["CYGB", "OLFML3", "HGF",
             "COLEC11", "SPARC",
             "NGFR", "NES", "DES", "RGS5",
             "BAMBI", "PPARG",
             "ADFP", "ADIPOR1"],
    "Hepatomesenchyme": ["SNAI1", "SNAI2",
                         "CDH2", "ALCAM",
                         "ZEB2", "LHX2",
                         "ITGB1", "FLNA",
                         "VCL", "MYL4",
                         "MYO9A",
                         "DLK1", "FOXA3",
                         "SLBP", "WASF1",
                         "MFAP2", "MFAP4",
                         "CNN2"],
    "Inflammation": ["IL1R1", "CX3CL1",
                     'CXCL12', "CXCL1",
                     'NFKBIA', 'TNFAIP3',
                     "CCL20", "CXCL5",
                     "FAM3C", "IL1A",
                     "IL32", "IL8",
                     'NFKB2', "HLA-A",
                     "HLA-C", "TLR2",
                     "TNFSF2", "TNFA"],
    "Inflammation_module": ["NFKBIA",
                            "CX3CL1",
                            "IL1A",
                            "IL32",
                            "CCL20",
                            "TLR2",
                            "HLA-DPB1",
                            "HLA-A",
                            "IL1R1", "CXCL5",
                            "TNFAIP3", "NFKB2"],
    "Fibrosis_module": ["COL1A1", "COL3A1",
                        "TAGLN", "TAGLN2",
                        "TGFBI", "PDGFA",
                        "PDGFB", "ACTG1",
                        "COL4A1", "COL5A1",
                        "COL6A1", "S100A11"],
    "Fibrosis": ["ITGB1", "TAGLN",
                 "TAGLN2", "IGFBP3",
                 "PDGFRA", "COL4A1",
                 "ACTA2", "AEBP12",
                 "COL1A1", "COL3A1",
                 "COL5A1",
                 "ITGAV", "LOXL2",
                 "DES", "PDGFRB", "PDGFB",
                 "S100A11", "S100A16",
                 'TGFBI',
                 'TGFBRI',
                 'TIMP1'],
    "Hepato_OS": ["CLU", "CDC42",
                  "GNG5", "KRT17",
                  "ENO1", "LDHA",
                  "CKB"],
    "1": ["COL1A1", "COL1A2", "CXCL6",
          "ACTA2", "TNFAIP3",
          "COL11A1",
          "COL3A1",
          "COL21A1",
          "COL20A1",
          "COL4A5",
          "COL9A1",
          "COL8A2",
          "COL9A3",
          "COL9A2",
          "COL11A2",
          "KPNA2",
          "SMLR1",
          "GREM1"
          "MMP7", "MMP9", "MMP10", "MMP11", "MMP13", "MMP24", "MXRA5"
                                                              "SPARCL1",
          "LAMA1", "LAMA2", "LAMA4",
          "MUC5AC", "PLAT",
          "OTOG",
          "AMTN",
          "IDH2", "INMT",
          "FBN3",
          "ANXA1", "S100A7", "S100A8", "S100A9", "SPARCL1",
          "KPNA2"
          "PDGFD",
          "SPHK1"
          "SCARA3",
          "DDR2"
          "CPB2",
          "DES",
          ],
    "Cancer": ["WNT6", ],
    "0000OACHOL": ["NFKBIA",
                   "IL1A",
                   "CXCL5",
                   "TNFAIP3",
                   "C1QTNF3",
                   "EPHA2",
                   "IL32",
                   "CXCL14",
                   "ADORA1",
                   "JAG1",
                   "TAGLN", "ACTA2", "COL1A1",
                   "ITGAV", "PDGFB", "S100A11", "TGFBI", "TIMP1"
                   ],
    "000HEP": [
        "THBS1",
        "CD84",
        "CCL28",
        "ALOX15B",
        "SAA2",
        "SAA1",
        "REG3A", "CREB3L3",
        "UTG1A1",
        "KLKB1"
        "APOA4",
        "PLCG2",
        "C3",
        "BLNK",
        "TOMM40L", "APOH",
        "PLG",
        "TGM2",
        "HTRA1",
        "TNXB",
        "AGT", "ARG1", "SLPI",
        "TRIM15", "TRIM26", "TRIM31", "MX1",
        "ISG15", ],

    "1_HEP": ["PDGFRB",
              "LEP"
              "COL4A1",
              "COL4A2",
              "CYTOR",
              "GSN",
              "COL11A1",
              "COL21A1",
              "COL9A2",
              "COL11A2",
              "SMLR1",
              "OTOG",
              "S100A9", "SPARCL1",
              "S100A4",
              "TRBC2",
              "VIM",
              "KRT1",
              "JAG1",
              "TAGLN", "ACTA2", "COL1A1",
              "ITGAV", "PDGFB", "S100A11", "TGFBI", "TIMP1"],
    "HEP_INFL": ["KLRC3",
                 "F2R"
                 "UL1A"
                 "SERPINA1"
                 "KLRD1",
                 "SRC",
                 "FFAR2",
                 "ADORA1",
                 "IFI16",
                 "CXCL12",
                 "VCAM1",
                 "TNF", "CXCL13", "S100A8", "CCL5",
                 "TNFAIP3",
                 ],

    "OA500CHOL_INFL": [
        "NFKBIA",
        "CX3CL1",
        "IL1A",
        "CXCL5",
        "TNFAIP3",
        "C1QTNF3",
        "CD74", "CXCL6",
        "EPHA2",
        "IL32",
        "CXCL14",
        "ADORA1",

    ],
    "2_STELLATE": ["IL1RL1",
                   "CXCL13",
                   "IL33",
                   "TWIST1",
                   "CXCL14",
                   "ADORA1",
                   "TYMS",
                   "BST2",
                   "C1QTNF3",
                   "FOXJ1",
                   "ADCYAP1"
                   "ANXA1",
                   "CCL21",
                   "INS",
                   "IL18BP", ],

    "2": ["FOXJ1",
          "TYMS",
          "BST2",
          "C1QTNF3",
          "CD74",
          "LPB",
          "CXCL6",
          "ALOX15B",
          "EPHA2",
          "GSTP1",
          "IL1RL1",
          "C5",
          "ADCYAP1"
          "ANXA1",
          "CCL21",
          "CXCL13",
          "IL32",
          "IL33",
          "TWIST1",
          "CXCL14",
          "IL18BP",
          "INS",
          "ADORA1",
          "APOA4",
          "FFAR2",
          "IVS",
          "SNAP1",
          "NPY",
          "NNS",
          ]
}
# END OF SCRIPT
