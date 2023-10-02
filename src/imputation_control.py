"""

Compare imputed and normalized gene expression side-by-side.

@author: Anja Hess
@date: 2023-JUN-06

"""

import os
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib

sns.set_style('ticks')
matplotlib.rcParams['figure.figsize'] = [4, 4]
matplotlib.rcParams['figure.dpi'] = 500
#cmap = sns.diverging_palette(220, 20, as_cmap=True)
matplotlib.rcParams['image.cmap'] = 'Spectral_r'
sc.set_figure_params(scanpy=True, fontsize=6)


def side_by_side(adata, gene_list, layers, title="", save=""):
    """

    To compare two AnnData layers on a projection side-by-side

    :param adata: AnnData object
    :param gene_list: list
    :param layers: list
    :param title: str
    :param output_path: str

    returns: Plots are generated

    """

    fig, ax = plt.subplots(20, 2, figsize=(10, 60))  # nrows, ncols
    print("Plotting imputed / non-imp.")
    for i, layer in enumerate(layers):
        for x, gene in enumerate(gene_list):
            if x < 11:
                trend = "up"
            else:
                trend = "down"

            try:
                subtitle = f"{trend}_{gene}_{layer.rsplit('_')[0]}"
            except:
                subtitle = f"{trend}_{gene}_{layer}"
            sc.pl.draw_graph(adata,
                             color=gene,
                             layer=layer,
                             title=subtitle,
                             ncols=1,
                             ax=ax[x, i],
                             save=False)
    fig.suptitle(title)
    plt.subplots_adjust(top=0.95)

    plt.savefig(f"{save}.pdf")
    plt.close()
    print("Saved.")
    # END OF FUNCTION


def side_by_side_dot(adata, gene_list, layers, title="", save=""):
    """

    To compare two AnnData layers on a projection side-by-side

    :param adata: AnnData object
    :param gene_list: list
    :param layers: list
    :param title: str
    :param save: str

    returns: Plots are generated

    """

    fig, ax = plt.subplots(2, figsize=(15, 6))  # nrows, ncols
    print("Plotting imputed / non-imp.")

    for i, layer in enumerate(layers):

        try:
            subtitle = f"{layer.rsplit('_')[0]}"
        except:
            subtitle = layer

        sc.pl.matrixplot(adata,
                         gene_list,
                         groupby="composite_simple",
                         save=False,
                         # cmap="Reds",
                         layer=layer,
                         ax=ax[i],
                         standard_scale="var",
                         colorbar_title=
                         f'Mean {subtitle} '
                         f'\n expression '
                         f'\n in group',
                         dendrogram=True,
                         )
        plt.title(subtitle)

    plt.tight_layout()
    fig.suptitle(title)
    plt.savefig(f"{save}.pdf")
    plt.close()
    print("Saved.")

    # END OF FUNCTION

def plot_imp_norm(adata_path, cluster_identifier="", nn=20,
                  ncomp=20, num_wp=1200, use_start=False,
                  pca_proj="X_draw_graph_fa", hvg=False,
                  lineage=False,
                  lineage_info={"Hepatic": ["FH", "HB", "AH",
                                            "Hepatocyte", "DC",
                                            "Cholangiocytes",
                                            "Ductal"],
                                "Stellate": ["ellate", "ibro", "uscle"]},
                  impute_data=False, top=10
                  ):
    """

    Palantir implementation, see for https://nbviewer.org/github/
    dpeerlab/Palantir/blob/master/notebooks/
    Palantir_sample_notebook.ipynb

    :param adata_path: str
    :param cluster_identifier: str
    :param nn: int
    :param ncomp: int
    :param num_wp: int
    :param use_start: bool
    :param pca_proj: str
    :param hvg: bool
    :param lineage: bool
    :param id: bool / str
    :param lineage_info: python dictionary
    ;param impute_data: bool
    :return: trajectory inference analysis

    """

    ############################################################################
    # 1. Directories
    ############################################################################
    print("--- Imputation control for", lineage)
    suffix = adata_path.rsplit(".h5ad", 1)[0]
    save_file = f"{suffix}_palantir_{hvg}_{pca_proj}_{ncomp}_" \
                f"{lineage}_lineage.h5ad"
    save_folder = f"{suffix.rsplit('/', 1)[0]}/PALANTIR/"
    root_dir = f"{suffix.rsplit('/', 1)[0]}/PALANTIR/{pca_proj}_{hvg}_{ncomp}" \
               f"_{nn}_{num_wp}_{lineage}/"
    os.makedirs(save_folder, exist_ok=True)
    magic_key = f"MAGIC_imputed_data_X_draw_graph_fa_{nn}_{lineage}"
    gene_file_key = "Supplementary_Table_7_Inflammation-fibrosis-palantir_"

    ############################################################################
    # 2.  Load adata with MAGIC
    ############################################################################
    if os.path.isfile(save_file):
        adata = sc.read_h5ad(save_file)
        print(adata)
    else:
        print("File missing")
        exit()

    ############################################################################
    # 3. Load the genes dataframe
    ############################################################################
    for term in ["inflammation", "fibrosis"]:
        print(f"----------------------{term}----------------------------------")
        df = pd.read_csv(f"{gene_file_key}{term}.csv")

        ########################################################################
        # Cut to the top / bottom 10
        ########################################################################
        #df_top = df.head(n=top)
        #df_bottom = df.tail(n=top)
        #df = df_top.append(df_bottom)

        ########################################################################
        # 3.0 Get sorted genes for each cell type
        ########################################################################
        for terminal_state in df.columns:
            if lineage == "Stellate" and ("ycling" in terminal_state or
                                          "Ductal" in terminal_state):
                continue
            if lineage == "Hepatic" and ("ellate" in terminal_state or
                                         "ibrobl" in terminal_state or
                                         "uscle" in terminal_state):
                continue
            print(f"--- Plotting {terminal_state}")
            genes = list(df[terminal_state])
            side_by_side(adata, gene_list=genes, layers=[None, magic_key],
                         title=f"{term} \n {lineage} lineage \n {terminal_state}",
                         save=f"{term}_{lineage}_{terminal_state}")
            side_by_side_dot(adata, gene_list=genes, layers=[None, magic_key],
                         title=f"Term: {term}. Lineage: {lineage}. "
                               f"Terminal state: {terminal_state}.",
                         save=f"DOT_{term}_{lineage}_{terminal_state}")
    # END OF FUNCTION
# END OF FILE
