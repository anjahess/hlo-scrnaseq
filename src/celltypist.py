"""

Script for the automated, reference-based annotation of
10X scRNAseq data.
See end of script for starting instructions.

@author: Anja Hess
@date: 2023-MAY-01

"""

import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import numpy as np
import os
import gzip
import shutil
import scanpy.external as sce
import pickle
import seaborn as sns
import celltypist
from celltypist import models

sc.settings.verbosity = 3
sc.logging.print_header()
sc.settings.file_format_figs = 'svg'
sc.settings.dpi_save = 10
sc.settings.verbosity = 0
sc.settings.autosave = True
sc.settings.set_figure_params(dpi=300,
                              fontsize=10,
                              figsize=[3, 3])


def run_ct(path_to_mtx="", file_id="",
           ref="", mode="best match",
           ref_clusteridentifier="Manual_Annotation",
           feature_selection=True, min_prop=0.1,
           top_genes=100):
    """

    CellTypist implementation

    :param path_to_mtx: str
    :param file_id: str
    :param ref: str
    :param mode: str
    :param ref_clusteridentifier: str
    :param feature_selection: bool
    :param min_prop: float
    :param top_genes: int
    :return:
    """
    ref_id = ref.rsplit(".h5ad")[0]
    if "/" in ref_id:
        ref_id = ref_id.rsplit("/", 1)[1]
    print(f"--- Annotating {file_id} to {ref_id}")
    outdir = f"{ref_id}/"
    outpath = f"{ref_id}/{file_id}_{mode}_" \
              f"{top_genes}" \
              f"_{ref_clusteridentifier}_" \
              f"{feature_selection}_{min_prop}/"
    save_file = f"{outpath}/{file_id}_" \
                f"to{ref_id}_{mode}" \
                f"_SCN.h5ad"
    model_dir = f"{ref_id}/{ref_id}_" \
                f"{feature_selection}_{top_genes}_" \
                f"{ref_clusteridentifier}_model.pkl"
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(outpath, exist_ok=True)
    sc.settings.figdir = outpath

    if not os.path.isfile(save_file):
        ####################################################################
        # 1. Load the reference (= training) data.
        ####################################################################
        if not os.path.isfile(model_dir):
            adata_2000 = sc.read_h5ad(ref)
            sc.pp.normalize_total(adata_2000,
                                  target_sum=1e4)
            sc.pp.log1p(adata_2000)

            # Train model
            new_model = celltypist.train(
                adata_2000,
                labels=ref_clusteridentifier,
                n_jobs=40,
                top_genes=top_genes,
                feature_selection=feature_selection)

            # Save the model.
            new_model.write(model_dir)
            del new_model

        ####################################################################
        # 2. Load the query data
        ####################################################################
        print(f"--- Loading query sample {file_id}")

        if ".h5ad" in file_id:
            qDat = sc.read_h5ad(path_to_mtx)
            if "OS" in file_id:
                keep = (qDat.obs["controlstatus"]
                        == "CTRL")
                qDat = qDat[keep, :]
        else:
            qDat = sc.read_10x_mtx(
                path_to_mtx +
                "/filtered_feature_bc_matrix/",
                var_names='gene_symbols',
                cache=True,
                gex_only=True)

        qDat.obs["cluster"] = "Undefined"
        qDat.obs["sample"] = file_id
        qDat.obs["cluster"].fillna("Undefined",
                                   inplace=True)

        ####################################################################
        # 3. Annotate
        ####################################################################
        sc.pp.normalize_total(qDat, target_sum=1e4)
        sc.pp.log1p(qDat)

        predictions = celltypist.annotate(
            qDat,
            model=model_dir,
            mode=mode,
            majority_voting=True,
            min_prop=min_prop,
            over_clustering=
            "leiden_0.12sctype_200",
            p_thres=0.1
        )
        adata = predictions.to_adata()

        sc.pl.umap(
            adata,
            color=['leiden_0.12sctype_200',
                   'majority_voting'],
            save=ref_id + f"_voting")

        try:
            sc.pl.draw_graph(
                adata,
                color=['leiden_0.12sctype_200',
                       'majority_voting'],
                save=ref_id + f"_voting")
        except:
            print("No graph.")
    # ENF OF FUNCTION


############################################################################
# START
############################################################################

root = "$HOME/"
files = ["all_scrub_clean_0.5_qc_dr_Day21_qc_clust.h5ad"]

for file_id in files:
    run_ct(path_to_mtx=f"{dir}/{file_id}",
           ref=f"{root}anno/wesley_full.h5ad",
           ref_clusteridentifier="cell",
           file_id=file_id,
           mode="prob match",
           top_genes=100)

# END OF SCRIPT