"""

Preprocessing functions for scrnaseq analysis.
@author: Anja Hess
@date: 2022-OCT-01

"""
import os
import shutil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
from sys import exit
from sc_constants import GENECLASSES, DOUBLET_SCORE, \
    COUNT_LOWER, COUNT_UPPER, MIN_GENES, MIN_CELLS, \
    MITO_CUTOFF, RIBO_LOWER, RIBO_UPPER, \
    GENES_BY_COUNT_FRAME, geneset2path
from tools.lazy_functions import adjust_list, remove_duplicates
from tools.plotting.plots import qc_plot


def easy_scrublet(dirpath, file_names=[], outputpath="",
                  sample_name=False):
    """

    Implementation of doublet detection method "Scrublet"

    # pip install scrublet
    # Notebook: https://github.com/swolock/scrublet/blob/
    master/examples/scrublet_basics.ipynb

    :param dirpath -> cellranger out
    :param file_names: list
    :param outputpath: str
    :return: creates a csv (barcode2scrublet_outs)

    """
    import scrublet as scr
    import scipy.io
    import matplotlib.pyplot as plt
    import numpy as np
    import gzip
    import shutil
    import os

    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = 'Arial'
    plt.rc('font', size=14)
    plt.rcParams['pdf.fonttype'] = 42

    ############################################################################
    # 1. Set the outpath
    ############################################################################
    save_file = f"{outputpath}/{sample_name}_doublet-scores.csv"
    matrix_loc = dirpath+"/filtered_feature_bc_matrix"

    if os.path.isfile(save_file):
        print(f"--- Scublet file for {sample_name} exists.")
        return "DONE"
    print(f"--- {sample_name}: Scrublet, saving to {save_file}")
    print(f"--- Searching 10X {matrix_loc}")
    print(f"--- Making outdir {outputpath}")
    os.makedirs(outputpath, exist_ok=True)

    ############################################################################
    # 2. Scrublet wants uncompressed matrices
    ############################################################################
    for dirpath, dirnames, file_names in os.walk(matrix_loc):
        for elem in file_names:
            if "barcodes" in elem:
                if not os.path.isfile(dirpath + "/"
                                      + elem.split(".gz")[0]):
                    print("--- BARCODE FILE. Decompressing")
                    with gzip.open(dirpath
                                   + "/" + elem, 'rb') as f_in:
                        with open(dirpath
                                  + "/" + elem.split(
                            ".gz")[0], 'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out)
            if "matrix.mtx" in elem:
                matrixfile = elem
                if not sample_name:
                    sample_name = elem.split("_matrix.mtx")[0]
                if not os.path.isfile(dirpath
                                      + "/"
                                      + elem.split(".gz")[0]):
                    print("--- MTX FILE. Decompressing")
                    with gzip.open(dirpath
                                   + "/" + elem, 'rb') as f_in:
                        with open(dirpath
                                  + "/"
                                  + elem.split(".gz")[0],
                                  'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out)
            if ("genes" in elem or "features" in elem) \
                    and ".gz" in elem:
                if not os.path.isfile(dirpath +
                                      "/" +
                                      elem.split(".gz")[0]):
                    print("--- GENE FILE. Decompressing")
                    with gzip.open(dirpath
                                   + "/" + elem, 'rb') as f_in:
                        with open(dirpath
                                  + "/"
                                  + elem.split(".gz")[0],
                                  'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out)
                genefile = elem.split(".gz")[0]
        break

    ############################################################################
    # 3. Start scrublet
    ############################################################################
    input_dir = outputpath
    counts_matrix = scipy.io.mmread(matrix_loc +
                                    '/matrix.mtx').T.tocsc()
    genes = np.array(scr.load_genes(matrix_loc +
                                    '/features.tsv',
                                    delimiter='\t',
                                    column=1))
    barcodes = pd.read_table(matrix_loc +
                             '/barcodes.tsv',
                             header=None)
    barcodes.columns = ["barcodes"]
    barcodes = list(barcodes["barcodes"])

    # Expected fraction that are doublets, typically 0.05-0.1.
    scrub = scr.Scrublet(counts_matrix,
                         expected_doublet_rate=0.06)
    # Default pipeline
    doublet_scores, predicted_doublets = \
        scrub.scrub_doublets(min_counts=2,
                             min_cells=3,
                             min_gene_variability_pctl=85,
                             n_prin_comps=30)
    scrub.plot_histogram()
    plt.savefig(f"{outputpath}/{sample_name}_histogram.pdf")

    print("--- Saving the results")
    # Save the output to a df
    df = pd.DataFrame({'barcode':barcodes,
                       'doublet_score':
                           scrub.doublet_scores_obs_,
                       'predicted_doublet':
                           scrub.predicted_doublets_})
    df.to_csv(save_file)

    # Get 2D embdedding (just for vis.)
    print('--- Running UMAP...')
    scrub.set_embedding('UMAP',
                        scr.get_umap(scrub.manifold_obs_, 10,
                                             min_dist=0.3))

    scrub.plot_embedding('UMAP', order_points=True)
    plt.savefig(f"{outputpath}/{sample_name}_umap.pdf")
    print(f'--- Scrublet done for sample {sample_name}.')
    # END OF FUNCTION


def easy_qc(adata, run_id="", min_genes=MIN_GENES, min_cells=MIN_CELLS,
            outputpath="", plot_data=False, n_highest_expr=20,
            filtering=True, group_by=None, samplewise_qc_criteria=False,
            scoring=False, basic_qc_criteria={}):
    """

    Main function for QC of scRNAseq data.

    :param adata: AnnData
    :param run_id: str
    :param min_genes: int
    :param min_cells: int
    :param outputpath: str
    :param plot_data: bool
    :param n_highest_expr: int
    :param filtering: bool
    :param group_by: str
    :param samplewise_qc_criteria: bool
    :param scoring: bool
    :param basic_qc_criteria: dict
    :return: QC-filtered AnnData object

    """
    # 0.1 Load the h5ad file
    adata_path = adata
    adata = sc.read_h5ad(adata)

    # 0.2 Path settings
    os.makedirs(outputpath, exist_ok=True)
    sc.settings.figdir = outputpath
    os.makedirs(f"{outputpath}QC/", exist_ok=True)
    logbook_path = f"{outputpath}QC/logbook.txt"
    print(f"--- Saving metrics here: {logbook_path}")
    logbook = open(logbook_path, "w")
    print(adata)

    # 1. Cell filtering
    print(f"--- Filter cells.")
    print(adata)
    sc.pp.filter_cells(adata, min_genes=MIN_GENES)

    # 2. Gene filtering
    print("--- Filter genes.")
    sc.pp.filter_genes(adata, min_cells=MIN_CELLS)

    # 3. Calculate metrics:
    sc.pp.calculate_qc_metrics(adata,
                               inplace=True,
                               log1p=False)
    metrics_custom(adata)

    # 4. Plot situation before filtering
    try:
        qc_plot(adata,
                sample_id=run_id,
                outputpath=outputpath,
                foldername="/QC/PRE_FILTER",
                group_by=group_by)
    except:
        qc_plot(adata,
                sample_id=run_id,
                outputpath=outputpath,
                foldername="/QC/PRE_FILTER",
                )
    # 5. Count filtering (post plot inspection)
    print("count filtering")
    keep = (adata.obs['total_counts'] < COUNT_UPPER)
    adata = adata[keep, :]
    print(adata)

    # 6. Mito filtering
    print("Mito filtering")
    keep = (adata.obs['pct_mito'] < MITO_CUTOFF)
    adata = adata[keep, :]
    print(adata)

    # 7. Ribo filtering - inspected violin plots for 1.5 interquartile
    # range (< 0.4)
    print("Ribo filtering")
    keep = (adata.obs['pct_ribo'] < 0.4)
    adata = adata[keep, :]
    print(adata)

    # 8. Doublet filtering
    print("Doublet filtering")
    keep = (adata.obs['doublet_score'] < 0.5)
    adata = adata[keep, :]
    print(f"--- After doublet clearing: {adata.shape}")

    # 9. Plot situation after filtering
    qc_plot(adata,
            sample_id=run_id,
            outputpath=outputpath,
            foldername="/QC/POST_FILTER",
            group_by=group_by)

    # 10. Close logboook
    logbook.write(str(adata.shape))
    logbook.close()

    # 11. Save the qc filtered results
    save_file = f'{outputpath}{run_id}_qc.h5ad'
    adata.write_h5ad(save_file)
    print(f"--- Preprocessing (QC) complete.")
    print(f"--- Saved to {save_file}")
    return save_file


def metrics_custom(adata):
    """
    if a certain gene class needs to be inspected

    Function to calculate QC metrics for
    custom gene classes
    :param adata: AnnData object
    :return: new slot in AnnData object

    """

    for geneclass in GENECLASSES:
        genes_of_interest = \
            adata.var_names.str.startswith(geneclass)
        adata.obs[f'pct_{GENECLASSES[geneclass]}'] \
            = np.sum(adata[:,
                     genes_of_interest].X, axis=1) / \
              np.sum(adata.X, axis=1)
    # END OF FUNCTION

def preprocessing(adata, run_id="", min_genes=300, min_cells=3, outputpath="",
                  plot_data=False, n_highest_expr=20, filtering=True,
                  group_by=None, samplewise_qc_criteria=False, scoring=False,
                  basic_qc_criteria={'n_genes_by_counts': False,
                                     'total_counts': False,
                                     'percent_mito': 0.4,
                                     'percent_ribo': (0.1, 0.45)}):

    """
    Preprocessing to add gene set scores.
    QC previously incorporated now moved to easy_pp.

    :param adata: AnnData
    :param run_id: str
    :param min_genes: int, deprecated
    :param min_cells: int, deprecated
    :param outputpath: str
    :param plot_data: bool
    :param n_highest_expr: int, deprecated
    :param filtering: bool
    :param group_by: str
    :param samplewise_qc_criteria: bool, deprecated
    :param scoring: bool
    :param basic_qc_criteria: dict, deprecated
    :return: Each cell gets a score for gene set of interest

    """
    # 0. Load the h5ad file
    adata_path = adata
    adata = sc.read_h5ad(adata)

    if plot_data:
        os.makedirs(outputpath, exist_ok=True)
        sc.settings.figdir = outputpath

    if scoring:
        for i, geneset in enumerate(geneset2path):
            # 1. Read the custom genelist from config
            try:
                df = pd.read_csv(geneset2path[geneset])
            except:
                df = pd.read_table(geneset2path[geneset])

            # 2. generate a boolean variable in the df
            geneset_list = df.iloc[:, 0].tolist()
            try:
                geneset_list = remove_duplicates(
                    adjust_list(geneset_list,
                                adata.var['gene_ids']))
            except:
                print("No removal.")

            adata.var[geneset] = adata.var_names.isin(geneset_list)
            # Calculate top 5 genes with highest counts and exclude
            targets_total_counts = adata[:,
                                   geneset_list].var[
                "total_counts"].sort_values(ascending=False)
            minus_top5 = [e for e in
                          geneset_list if e
                          not in list(
                    targets_total_counts.index[:4])]

            scenarios = { "all": geneset_list,
                          "minus_top5": minus_top5}

            for scenario in scenarios:
                # PART A: calculate percentage of counts
                adata.obs[f'pct_counts_{geneset}_{scenario}'] = \
                    np.sum(adata[:,
                           adata.var_names.isin(scenarios[
                                                    scenario])].X,
                           axis=1) / np.sum(adata.X, axis=1)
                adata.obs[f'rel_pct_counts_' \
                          f'{geneset}_{scenario}'] = \
                    adata.obs[
                        f'pct_counts_{geneset}_{scenario}'
                    ] * 1000

                # PART B: calculate the number of different
                # inflammatory genes / all genes per cell
                adata.obs[f'diversity_{geneset}_{scenario}'] = \
                    np.count_nonzero(adata[:,
                                     adata.var_names.isin(scenarios
                                                          [scenario
                                                          ])].X.toarray(),
                                     axis=1)
                adata.obs[f'pct_diversity_{geneset}_{scenario}'] = \
                    adata.obs[f'diversity_{geneset}_{scenario}'] / \
                    (len(scenarios[scenario])/100)
                adata.obs[f'score_{geneset}_{scenario}'] = \
                    adata.obs[f'rel_pct_counts_{geneset}_{scenario}'] \
                    + adata.obs[f'pct_diversity_{geneset}_{scenario}']

    # 4. Save the qc filtered results
    save_file = f'{outputpath}{run_id}_qc.h5ad'
    adata.write_h5ad(save_file)
    print(f"--- Preprocessing complete.")
    print(f"--- Saved to {save_file}")
    return save_file

    # END OF FUNCTION

# END OF FILE