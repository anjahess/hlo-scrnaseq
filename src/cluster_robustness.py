"""

Python script for the evaluation of cluster
robustness parameters of 10X scRNAseq data.

Tutorials used:
https://evafast.github.io/blog/2019/06/28/example_content/
https://github.com/evafast/Jupyter_notebooks/blob/master/
DBindex_Silscore_clustering_tutorial.ipynb

@author: Anja Hess
@date: 2023-MAY-01

"""

from __future__ import print_function
import os
import scanpy as sc
from time import sleep
import pandas as pd
import warnings
import matplotlib.pyplot as plt
from sklearn.metrics import silhouette_samples, silhouette_score
from sklearn.metrics import davies_bouldin_score
sc.settings.verbosity = 0
sc.logging.print_header()
sc.settings.dpi_save = 10
sc.settings.verbosity = 0
sc.settings.autosave = True
sc.settings.set_figure_params(dpi=300,
                              fontsize=14,
                              figsize=[3, 3])
warnings.filterwarnings('ignore')


def run_ana(path_to_mtx="", file_id=""):
    """
    Example for implementation can be found here:
    https://github.com/evafast/Jupyter_notebooks/
    blob/master/DBindex_Silscore_clustering_tutorial.ipynb
    :param path_to_mtx: str
    :param file_id: str
    :return: tables with cluster robustness metrics
    """
    ####################################################################
    # Define names and dirs
    ####################################################################
    file_name = file_id.rsplit(".h5ad")[0]
    outpath = f"{file_id}/"
    save_file = f"{outpath}/{file_id}_results.tsv"
    os.makedirs(outpath, exist_ok=True)
    sc.settings.figdir = outpath

    print(f"--- Cluster robustness: {file_name}.")

    if not os.path.isfile(save_file):
        ################################################################
        # 0. Load data
        ################################################################
        if ".h5ad" in file_id:
            adata = sc.read_h5ad(path_to_mtx)
        else:
            adata = sc.read_10x_mtx(path_to_mtx +
                                   "/filtered_feature"
                                   "_bc_matrix/",
                                   var_names=
                                   'gene_symbols',
                                   cache=True,
                                   gex_only=True)

        ################################################################
        # 1. Evaluate clustering by varying Leiden res
        ################################################################
        results_df = pd.DataFrame(columns=
                                  ['resolution',
                                   'number_of_clusters',
                                   'sil',
                                   'davie_bould'])

        for i in range(1, 20, 1):
            res = i/10
            print(f"--- Testing leiden res {res}")
            sc.tl.leiden(adata, resolution=res)
            silhouette_avg = \
                silhouette_score(adata.obsm['X_pca'],
                                 adata.obs['leiden'])
            print(f"     Silhouette: {silhouette_avg}")
            davies_bouldin_avg = davies_bouldin_score(
                adata.obsm['X_pca'],
                adata.obs['leiden'])

            print(f"     DB: {davies_bouldin_avg}")
            results_df = results_df.append(
                pd.DataFrame([[res, max(adata.obs['leiden']),
                               silhouette_avg,
                               davies_bouldin_avg]],
                             columns=['resolution',
                                      'number_of_clusters',
                                      'sil',
                                      'davie_bould']))
            sleep(0.01)
            # END OF LOOP

        ################################################################
        # 2. Save results
        ################################################################
        results_df['number_of_clusters'] = \
            results_df['number_of_clusters'
            ].astype(float) + 1
        print(results_df)
        results_df.to_csv(save_file,
                          sep="\t")
        # END OF IF STATEMENT

    if os.path.isfile(save_file):
        ################################################################
        # 3. Load results
        ################################################################
        results_df = pd.read_table(save_file)
        print(results_df)

        # make some new dataframe for plotting
        max_sil = \
            results_df[results_df.groupby('number_of_clusters')[
                           'sil'].transform('max') ==
                       results_df['sil']].sort_values(
                by=['number_of_clusters'])

        min_db = results_df[
            results_df.groupby('number_of_clusters')[
                'davie_bould'].transform('min') ==
            results_df['davie_bould']].sort_values(
            by=['number_of_clusters'])

        # plot the top value for each cluster
        plt.plot(min_db['number_of_clusters'],
                 min_db['davie_bould'])
        plt.plot(max_sil['number_of_clusters'],
                 min_db['sil'] * 10)
        plt.legend(['Davies–Bouldin index',
                    'Silhouette score * 10'])
        plt.xlabel('number of clusters')
        plt.ylabel('clustering score')
        plt.savefig(f"{outpath}table_sil.pdf")
        plt.close()
        print("--- Done.")
    # ENF OF FUNCTION

############################################################################
# START ANALYSIS
############################################################################

# 1. OS-HLOs
dir = "$HOME/"
file_id = "all_scrub_clean_0.5_qc_dr_OS_qc_clust.h5ad"
run_ana(path_to_mtx=f"{dir}/{file_id}",
        file_id=file_id)

# 2. ULA-HLOs
file_id = "all_scrub_clean_0.5_qc_dr_Day21_qc_clust.h5ad"
run_ana(path_to_mtx=f"{dir}/{file_id}",
        file_id=file_id)

# END OF SCRIPT