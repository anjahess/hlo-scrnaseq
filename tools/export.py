"""

Python script for exporting h5ad files to other formats.

@author: Anja Hess
@date: 2023-MAY-01

"""

from __future__ import print_function
import os
import scanpy as sc
import scipy
import pandas as pd


def export(path_to_h5ad, out="./"):
    """

    :param path_to_mtx: str
    :param file_id: str
    :return: tables with cluster robustness metrics
    """
    ####################################################################
    # Define names and dirs
    ####################################################################
    ad = sc.read_h5ad(path_to_h5ad)
    print(ad)

    destination = "./"
    pd.DataFrame(ad.obsm["X_umap"], index=ad.obs_names).to_csv(
        f"{out}umap.tsv", sep="\t")
    pd.DataFrame(ad.obsm["X_draw_graph_fa"], index=ad.obs_names).to_csv(
        f"{out}fa.tsv", sep="\t")
    pd.DataFrame(ad.var.index).to_csv(f"{out}genes.tsv", sep="\t")
    pd.DataFrame(ad.obs.index).to_csv(f"{out}barcodes.tsv", sep="\t")
    ad.obs.to_csv(os.path.join(destination, f"{out}metadata.tsv"),
                  sep="\t")
    scipy.io.mmwrite(f"{out}matrix.mtx", ad.X)

    print("--- Done.")
    # ENF OF FUNCTION


############################################################################
# START
############################################################################
# 1. OS-HLOs
exp_path = "./EXPORT/"
os.makedirs(exp_path, exist_ok=True)
file_ids = ["all_scrub_clean_0.5_qc_dr_Day21_qc_clust.h5ad",
            "all_scrub_clean_0.5_qc_dr_OS_qc_clust.h5ad"]
for file_id in file_ids:
    export(file_id, out=exp_path)
