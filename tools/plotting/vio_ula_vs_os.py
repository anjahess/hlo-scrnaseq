"""

Python script for violin plots comparing gene
expression in ULA and OS HLOs
@author: Anja Hess
@date: 2023-JUN-09

"""

import scanpy as sc
import seaborn as sns
import warnings
import pandas as pd

sc.settings.file_format_figs = 'svg'
sc.settings.dpi_save = 10
sc.settings.verbosity = 0
sc.settings.autosave = True
sc.settings.set_figure_params(dpi=300,
                              fontsize=14,
                              figsize=[3, 3])
cmap = sns.diverging_palette(220, 20,
                             as_cmap=True)
warnings.filterwarnings("ignore")


sample2controlstatus = {
    "45": "CTRL",
    "46": "CTRL",
    "47": "CTRL",
    "50": "CTRL",
    "51": "CTRL",
    "52": "TGFB1",
    "53": "TGFB1",
    "62": "CTRL",
    "63": "CTRL",
    "64": "PA",
    "65": "PA",
    "SM-L3XWE": "CTRL",
    "SM-L3XWF": "CTRL",
    "SM-L3XWG": "CTRLOA",
    "SM-L3XWH": "CTRLOA",
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

genes = ["SPARC",  "COL3A1", 'COL1A1',
         'IGFBP7','COL1A2', 'COL6A1',
         "IGFBP3", "DCN", "DES",
         "PDGFRB", "TGFBI", "ACTA2"]


def vioplots(path_to_mtx="", file_id="", cluster_identifier="composite",
             title=""):
    """

    :param path_to_mtx: str
    :param file_id: str
    :param cluster_identifier: obs slot
    :return: Statistics for abundance

    """
    #######################################################################
    # Define names and dirs
    #######################################################################
    file_name = file_id.rsplit(".h5ad")[0]

    ###################################################################
    # 1. Load the query data
    ###################################################################
    adata = sc.read_h5ad(path_to_mtx)
    print(adata)

    if "culture" not in list(adata.obs):
        adata.obs["culture"] = adata.obs[
            "sample"].map(sample2culture).astype('category')
        adata.write_h5ad(path_to_mtx)

    if "controlstatus" not in list(adata.obs):
        adata.obs["controlstatus"] = adata.obs[
            "sample"].map(sample2controlstatus).astype('category')
        adata.write_h5ad(path_to_mtx)

    ###################################################################
    # 2. Select controls Load the query data
    ###################################################################
    keep = (adata.obs["controlstatus"]
            == "CTRL")
    adata = adata[keep, :]

    ###################################################################
    # 3. Plot
    ###################################################################
    sc.pl.violin(adata,
                 genes,
                 groupby='culture',
                 rotation=90,
                 save=f"_{title}",
                 stripplot=False,
                 layer="scaled",
                 use_raw=False)

    sc.pl.stacked_violin(adata,
                         genes,
                         groupby='culture',
                         rotation=90,
                         save=f"_{title}",
                         layer="scaled",
                         use_raw=False
                         )
    # ENF OF FUNCTION


############################################################################
# START
############################################################################

# 1. OS-HLOs
dir = "./"
file_id = "all_scrub_clean_0.5_qc_dr.h5ad"
vioplots(path_to_mtx=f"{dir}/{file_id}",
         file_id=file_id, title=file_id.rsplit(".h5ad")[0])

# END OF SCRIPT
