"""

Python script for the testing of cell type abundance
https://github.com/theislab/scCODA

@author: Anja Hess
@date: 2023-MAY-01

"""

import scanpy as sc
import seaborn as sns
import warnings
import pandas as pd
import pickle as pkl
from sccoda import util as scu
import os
from sccoda.util import comp_ana as mod

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
    "SM-L3XWG": "CTRL",
    "SM-L3XWH": "CTRL",
    "SM-L3XWI": "OA500",
    "SM-L3XWJ": "OA500",
}


def run_scCODA(path_to_mtx="", file_id="", cluster_identifier="composite"):
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
    outdir = f"{file_id}/"
    outpath = f"{file_id}/"
    save_file = f"{outpath}/{file_id}_ciara.h5ad"
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(outpath, exist_ok=True)
    sc.settings.figdir = outpath

    print(f"--- Searching rare celltypes for {file_name}.")
    if not os.path.isfile(save_file):
        ###################################################################
        # 1. Load the query data
        ###################################################################
        adata = sc.read_h5ad(path_to_mtx)

        ###################################################################
        # 2. Make covariate DataFrame
        ###################################################################
        samples = sample2controlstatus.keys()
        idx = sample2controlstatus.values()
        cov_df = pd.DataFrame({"controlstatus": idx},
                              index=samples)
        print(cov_df)
        print(adata.obs["sample"].value_counts())

        ###################################################################
        # 3. Load from scanpy
        ###################################################################
        data_scanpy_1 = scu.cell_composition_data.from_scanpy(
            adata,
            cell_type_identifier=cluster_identifier,
            sample_identifier="sample",
            covariate_df=cov_df
        )
        model_salm = mod.CompositionalAnalysis(
            data_scanpy_1,
            formula="controlstatus",
            reference_cell_type="AH_ns_0")

        # Run MCMC
        sim_results = model_salm.sample_hmc()
        sim_results.summary()
        print(sim_results.credible_effects())
        sim_results.set_fdr(est_fdr=0.4)
        sim_results.summary()

        # SAVE
        path = "test"
        sim_results.save(path)

        # loading
        with open(path, "rb") as f:
            sim_results_2 = pkl.load(f)
        sim_results_2.summary()
    else:
        print("scCODA version exists.")
    # ENF OF FUNCTION


############################################################################
# START
############################################################################

# 1. OS-HLOs
dir = "$HOME/"
file_id = "all_scrub_clean_0.5_qc_dr_OS_qc_clust.h5ad"
run_scCODA(path_to_mtx=f"{dir}/{file_id}",
           file_id=file_id)

# END OF SCRIPT
