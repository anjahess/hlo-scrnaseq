"""

Plotting functions

@author: Anja Hess
@date: 2022-OCT-01

"""
import os
import scanpy as sc
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import sys
import scipy.stats as scipy_stats
import decimal
import numpy as np
import scikit_posthocs as sp
import pandas as pd

path_delim = "/"
script_path = str(os.path.dirname(os.path.abspath(__file__)))
maindir = script_path.split(path_delim + "tools/plotting")[0]
src_path = maindir + path_delim + "src"
bin_path = maindir + path_delim + "bin"
sys.path.insert(0, script_path)
sys.path.insert(0, maindir)
sys.path.insert(0, src_path)
sys.path.insert(0, bin_path)
sys.path.insert(0, os.path.abspath('.'))
sys.path.insert(0, os.path.abspath('..' + path_delim))
sys.path.insert(0, os.path.abspath('..' + path_delim + 'src'))
cwd = os.getcwd()
absFilePath = os.path.abspath(__file__)
script_path = str(os.path.dirname(os.path.abspath(__file__)))

from src.sc_constants import \
    sample2condition, palettediff2, \
    GENERAL_MARKERS, CELL_HYPERCATS, diff2controlcondition, \
    geneset2path, GROUP_COLS, cmap, CELL_COLS, sample2controlstatus, \
    umap_vars, genesets_plots, genesets_lineage, genesets_plots_big, \
    cond_palette
from src.data_parsing import find_cell
from src.core_functions import remove_duplicates, \
    adjust_list, score_geneset


def show_values(axs, orient="v", space=.01):
    """

    https://www.statology.org/seaborn-barplot-show-values/

    :param axs:
    :param orient:
    :param space:
    :return:

    """

    def _single(ax):
        if orient == "v":
            for p in ax.patches:
                _x = p.get_x() + p.get_width() / 2
                _y = p.get_y() + p.get_height() + (p.get_height() * 0.01)
                value = '{:.1f}'.format(p.get_height())
                ax.text(_x, _y, value, ha="center")
        elif orient == "h":
            for p in ax.patches:
                _x = p.get_x() + p.get_width() + float(space)
                _y = p.get_y() + p.get_height() - (p.get_height() * 0.5)
                value = '{:.1f}'.format(p.get_width())
                ax.text(_x, _y, value, ha="left")

    if isinstance(axs, np.ndarray):
        for idx, ax in np.ndenumerate(axs):
            _single(ax)
    else:
        _single(axs)
    # END OF FUNCTION


def get_label_rotation(angle, offset):
    """

    https://www.python-graph-gallery.com/circular-barplot-with-groups

    :param angle: int
    :param offset: int

    :return:

    """
    # Rotation must be specified in degrees :(
    rotation = np.rad2deg(angle + offset)
    if angle <= np.pi:
        alignment = "right"
        rotation = rotation + 180
    else:
        alignment = "left"
    return rotation, alignment


def add_labels(angles, values, labels, offset, ax):
    """
    https://www.python-graph-gallery.com/circular-barplot-with-groups
    :param angles:
    :param values:
    :param labels:
    :param offset:
    :param ax:
    :return:
    """

    # This is the space between the end of the bar and the label
    padding = 1

    # Iterate over angles, values, and labels, to add all of them.
    for angle, value, label, in zip(angles,
                                    values,
                                    labels):
        angle = angle
        label = label.rsplit(" ", 1)[0]

        words = label.count(" ")
        n_split = 4
        if words > n_split:
            label = label.rsplit(" ", words - n_split)[0]

        # Obtain text rotation and alignment
        rotation, alignment = get_label_rotation(angle, offset)

        # And finally add the text
        ax.text(
            color="black",
            size=20,
            x=angle,
            y=value + padding,
            s=label,
            ha=alignment,
            va="center",
            rotation=rotation,
            rotation_mode="anchor"
        )
    # END OF FUNCTION


def circular_barplot(df, outfolder, value="", name="", group="", id="", top=10,
                     title="CELLTYPE"):
    """

    For pathway enrichment visualization.

     	name 	value 	group
    0 	item 1   31     A
    https://www.python-graph-gallery.com/
    circular-barplot-with-groups

    :param df: pandas df
    :param outfolder: str
    :param value: str
    :param name: str
    :param group: str
    :param top: int
    :param title: str

    :return: plot generated

    """

    ############################################################################
    # 0. Load data, prepare dataframe
    ############################################################################
    if type(df) == str:
        df = pd.read_csv(df, index_col=0)
    df = df.loc[df[group] != "CTRL"]
    df = df.loc[df[group] != "CTRLOA"]
    df["-log10pval"] = -np.log10(df[value].astype(float))

    ############################################################################
    # 1. Plot only top n in list
    ############################################################################
    GROUP_NAMES = df[group].value_counts().index.tolist()
    GROUP_NAMES = sorted(GROUP_NAMES)
    if top:
        top_df = pd.DataFrame()
        for group_name in GROUP_NAMES:
            subdf = df.loc[df[group] == group_name][:top]
            top_df = pd.concat([top_df, subdf])
        # Update the df
        df = top_df
        df.reset_index(inplace=True, drop=True)

    ############################################################################
    # 2. Get the values
    ############################################################################
    VALUES = df["-log10pval"].values
    LABELS = df[name].values
    GROUP = df[group].values
    GROUP2N = df[group].value_counts().to_dict()
    GROUPS_SIZE = list(GROUP2N.values())
    ############################################################################
    # Set plot parameters
    ############################################################################
    PAD = 0  # Space between the groups
    ANGLES_N = len(VALUES) + PAD * len(
        np.unique(GROUP))
    ANGLES = np.linspace(0, 2 * np.pi,
                         num=ANGLES_N,
                         endpoint=False)
    WIDTH = (2 * np.pi) / len(ANGLES)

    ############################################################################
    # Determines where to place the first bar.
    ############################################################################
    OFFSET = np.pi / 2
    offset = 0
    IDXS = []
    for size in GROUPS_SIZE:
        IDXS += list(range(offset + PAD,
                           offset + size + PAD))
        offset += size + PAD
    fig, ax = plt.subplots(figsize=(25, 25),
                           subplot_kw={
                               "projection": "polar"})
    ax.set_theta_offset(OFFSET)
    ax.set_frame_on(False)
    ax.xaxis.grid(False)
    ax.yaxis.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])

    # BAR COLORS
    COLORS = [GROUP_COLS[e] for i, e in
              enumerate(GROUP_NAMES)
              for _ in range(GROUP2N[e])]

    # LEGEND
    patches = []
    for i, group_name in enumerate(
            GROUP_NAMES):
        group_patch = mpatches.Patch(
            color=GROUP_COLS[group_name],
            label=group_name)
        patches.append(group_patch)
    ax.legend(handles=patches)
    ax.bar(
        ANGLES[IDXS],
        VALUES,
        width=WIDTH,
        color=COLORS,
        linewidth=2,
        edgecolor="black"
    )

    ############################################################################
    # Extra customization
    ############################################################################
    offset = 0
    for group, size in zip(df[group
                           ].value_counts().index.tolist(),
                           GROUPS_SIZE):
        # Add line below bars
        x1 = np.linspace(ANGLES[offset + PAD],
                         ANGLES[offset + size + PAD - 1],
                         num=50)
        # Add text to indicate group
        ax.text(
            np.mean(x1), -.5,
            group,
            color="#333333",
            fontsize=40,
            ha="center",
            va="center"
        )
        # Add reference lines at 20, 40, 60, and 80
        x2 = np.linspace(ANGLES[offset],
                         ANGLES[offset + PAD - 1],
                         num=50)
        ax.plot(x2, [0.1] * 50,
                color="#bebebe",
                lw=0.8)
        ax.plot(x2, [0.5] * 50,
                color="#bebebe",
                lw=0.8)
        ax.plot(x2, [1] * 50,
                color="#bebebe",
                lw=0.8)

        offset += size + PAD

    add_labels(ANGLES[IDXS],
               VALUES,
               LABELS,
               OFFSET,
               ax)

    plt.tight_layout()
    plt.title(title + f" top {top} ", fontsize=40)
    plt.savefig(f"{outfolder}/circularbp_{id}_{top}.pdf")
    plt.close()
    print(f"--- Saved to {outfolder}")

    # END OF FUNCTION


def enrichment_barplot(df, outputpath="",
                       title="", cutoff_col="Adjusted P-value"):
    """

    Function to plot enrichment analysis.

    :param df: pandas df
    :param outputpath: str
    :param title: str
    :param cutoff_col: str

    :return: plot is generated

    """

    sns.set(rc={'figure.figsize': (8, 6)})
    sns.set_context("paper")
    sns.set_style("ticks")

    ############################################################################
    # 0. Prepare df
    ############################################################################
    if type(df) == str:
        df = pd.read_csv(df, index_col=0)
    df["-log10Padj"] = -np.log10(df["Adjusted P-value"].astype(float))
    df["Term"] = df["Term"].str.rsplit(" ", 1, expand=True)[0]

    ############################################################################
    # 1. Limit to significant results
    ############################################################################
    df = df.loc[df[cutoff_col] < 0.05]
    opposite = {"Combined Score": "-log10Padj"}
    # Sort by combined score OR adj-Pval
    for sorter in ["Combined Score"]:
        df.sort_values(by=[sorter], inplace=True, ascending=False)

        # Take top 10
        subdf = df.iloc[:10]
        sns.scatterplot(y="Term", x=sorter,
                        size=opposite[sorter],
                        hue="Odds Ratio",
                        palette="rocket_r",
                        sizes=(100, 200),
                        data=subdf)
        plt.title(title)
        plt.tight_layout()
        plt.savefig(outputpath + f"/{title}_{sorter}.pdf")
        plt.close()
    sns.set(rc={})
    # END OF FUNCTION


def assign_me_a_color(cell, colors_used=[]):
    """

    Function to choose color for cell type.

    :param cell: str
    :param colors_used: list

    :return: color str

    """
    counter = 0

    for celltype in CELL_COLS:
        if celltype in cell:
            color = CELL_COLS[celltype][0]
            while color in colors_used:
                counter += 1
                color = CELL_COLS[celltype][counter]
            return color

    # In case unknown:
    for col in cond_palette:
        if col not in colors_used:
            return col
    # END OF FUNCTION


def umap_cells_custom_color(adata, cluster_identifier):
    """

    Color cell types consistently

    :param adata: AnnData object
    :param cluster_identifier: str
    :return: plots and changes in Anndata uns slot

    """

    ############################################################################
    # 0. Assign default colors
    ############################################################################
    sc.pl.umap(adata,
               color=cluster_identifier,
               title=cluster_identifier,
               save=False)

    sc.pl.draw_graph(adata,
                     color=cluster_identifier,
                     title=cluster_identifier,
                     save=False)
    cells = adata.obs[cluster_identifier].value_counts().index.tolist()
    colors = []

    ############################################################################
    # 1. Assign custom colors
    ############################################################################
    for cell in cells:
        col = assign_me_a_color(cell, colors)
        colors.append(col)
        adata.uns[f"{cluster_identifier}_colors"][
            sorted(adata.obs[
                       cluster_identifier
                   ].value_counts().index.tolist()).index(cell)
        ] = col

    ############################################################################
    # 2. Plot with new colors
    ############################################################################
    sc.pl.umap(adata, color=cluster_identifier, title=cluster_identifier,
               save=cluster_identifier)
    sc.pl.draw_graph(adata, color=cluster_identifier, title=cluster_identifier,
                     save=cluster_identifier)

    for e in ["condition", "controlstatus", "culture", "phase"]:
        sc.pl.draw_graph(adata,
                         color=e,
                         title=e,
                         save=f"{e}_all",
                         palette=cond_palette)
    # END OF FUNCTION


def plot_cellproportions_per_cluster(adata, cluster_identifier="",
                                     plot_path="", dataset_name=""):
    """

    Function to plot cell fractions

    :param adata: anndata object
    :param cluster_identifier: adata.obs
    :param plot_path: str
    :param dataset_name

    :return: plots are created

    """
    samples = adata.obs["sample"].value_counts().index.tolist()
    adata.obs[f"{cluster_identifier}_counts"] = 1
    table = pd.pivot_table(adata.obs,
                           index=cluster_identifier,
                           columns="sample",
                           values=f"{cluster_identifier}"
                                  f"_counts",
                           aggfunc=np.sum)
    table.sort_values(by=samples,
                      inplace=True,
                      ascending=False)
    plt.rcParams["axes.grid"] = False
    table.plot.bar(stacked=True,
                   figsize=(7, 7),
                   color=palettediff2,
                   width=0.9,
                   linewidth=0.3,
                   edgecolor="black")

    plt.xticks(rotation=70,
               fontsize=12,
               fontweight='normal')

    plt.savefig(plot_path
                + f"abs_{dataset_name}.pdf")

    plt.close()
    # END OF FUNCTION


def plot_diff_props(adata, cluster_identifier="",
                    plot_path=""):
    """

    Bar graphs indicating the proportions of cells per cluster

    :param adata:
    :param cluster_identifier
    :param plot_path
    :param dataset_name

    :return: plots are created

    """

    sns.set_style("whitegrid", {'axes.grid': False})

    ############################################################################
    # 0. Prepare df
    ############################################################################
    for x, y in [("culture", cluster_identifier),
                 ("sample", cluster_identifier),
                 (cluster_identifier, "sample"),
                 ]:
        print(f"... Diff for {x}, {y}")
        adata.obs[f"counts"] = 1
        df = adata.obs.pivot_table(
            values=f"counts",
            index=x,
            columns=y,
            aggfunc=np.sum)

        ########################################################################
        # Calculate % of all cells in a sample for each cell type
        ########################################################################
        df = df.div(df.sum(axis=1), axis=0)

        ########################################################################
        # Calculate % of all cells in a sample for each cell type
        ########################################################################
        df.columns = df.columns.astype(str)
        df = df.reset_index()

        if x == cluster_identifier:
            df[x] = df[x].str.split("_", expand=True)
        df_for_stacked = df.transpose()

        df_for_stacked.to_csv(plot_path + f"props_{y}_{x}.csv")
        df_for_stacked = pd.read_csv(f"{plot_path}props_{y}_{x}.csv",
                                     index_col=0, header=[0, 1])
        df_for_stacked = df_for_stacked.fillna(0)
        df_for_stacked = df_for_stacked.transpose()

        ########################################################################
        # Get the colors from anndata uns slot
        ########################################################################
        try:
            palette = adata.uns[f"{cluster_identifier}_colors"]
        except:
            palette = palettediff2

        if x == "sample" and y == cluster_identifier and \
                "sctype" in cluster_identifier \
                or "composite" in cluster_identifier:
            try:
                # Reorder so colors match
                colors = {}
                for i, column in enumerate(df_for_stacked.columns):
                    match_key = [e for e in CELL_HYPERCATS if e in column]
                    new_term = CELL_HYPERCATS[match_key[0]] + "_" + column
                    df_for_stacked.rename(columns={column: new_term},
                                          inplace=True)
                    colors.update({new_term: palette[i]})
                df_for_stacked = df_for_stacked.reindex(
                    sorted(df_for_stacked), axis=1)
                palette = [colors[e] for e in df_for_stacked.columns]
            except:
                print("")

        ########################################################################
        # Plot.
        ########################################################################
        df_for_stacked.plot.bar(stacked=True,
                                figsize=(7, 7),
                                width=1,
                                linewidth=0.01,
                                edgecolor="black",
                                rot=45,
                                color=palette)
        plt.savefig(plot_path + f"props_{y}_{x}.pdf")
        plt.close()

        df_for_stacked.plot.bar(stacked=True,
                                figsize=(10, 10),
                                linewidth=0.01,
                                width=1,
                                color=palette,
                                edgecolor="black",
                                rot=45,
                                subplots=True,
                                legend=None)
        plt.savefig(plot_path + f"path_{y}_{x}.pdf")
        print(f"-- Stacked barplot saved to {plot_path}.")
        # END OF FUNCTION


def reduce_to_genes_expressed(adata, dictionary=""):
    """

    To take out genes from a list if not present in adata.var

    :param adata: h5ad
    :param dictionary: dict
    :return: geneset_adjusted: dict

    """

    try:
        geneset_adjusted = \
            adjust_dict(dictionary, adata.var['gene_ids'])
    except:
        try:
            geneset_adjusted = \
                adjust_dict(dictionary, adata.var['gene_ids-1'])
        except:
            try:
                geneset_adjusted = \
                    adjust_dict(dictionary, adata.var['gene_ids-0'])
            except:
                print(f"--- None is expressed: {dictionary}.")
                exit()

    return geneset_adjusted
    # END OF FUNCTION


def plot_gene_umaps(adata, dictionary="", run_id="",
                    plot_path="", dataset_name="", scaling=""):
    """

    To display gene expression on projections

    :param adata: h5ad object
    :param dictionary: dict
    :param run_id: str

    :return:

    """

    ########################################################################
    # Prepare marker gene structures
    ########################################################################
    dict_marker = {}
    all_genes = []

    ########################################################################
    # Plot
    ########################################################################
    for coloring in [True, False]:
        plot_rows = round(len(dictionary) + 1)
        max_geneset_len = max([len(v) for v in dictionary.values()])
        fig, axs = plt.subplots(plot_rows, max_geneset_len + 1,
                                figsize=(25, 35))
        plt.subplots_adjust(left=0.125, bottom=0.1, right=0.9, top=0.9,
                            wspace=1, hspace=1)
        fig.suptitle(f'Marker genes for {run_id}', fontsize=16)
        for elem in range(0, max_geneset_len + 1):
            for i in range(0, plot_rows):
                try:
                    axs[i, elem].axis('off')
                except:
                    continue
        ######################################################################
        # Iterate through gene sets
        ######################################################################
        for i, geneset in enumerate(dictionary):
            try:
                genes = dictionary[geneset]
                valid_genes = [x for x in genes if
                               x in adata.var['gene_ids'] or x in list(
                                   adata.var['gene_ids'].index.values)]
            except:
                try:
                    genes = dictionary[geneset]
                    valid_genes = [x for x in genes if
                                   x in adata.var['gene_ids-0'] or x in list(
                                       adata.var['gene_ids-0'].index.values)]
                except:
                    genes = dictionary[geneset]
                    valid_genes = [x for x in genes if
                                   x in adata.var['gene_ids-1'] or x in list(
                                       adata.var['gene_ids-1'].index.values)]
            dict_marker.update({geneset: valid_genes})
            all_genes += valid_genes

            for num, gene in enumerate(valid_genes):
                if coloring:
                    sub = sc.pl.umap(adata,
                                     color=gene,
                                     save=False,
                                     show=False,
                                     size=40,
                                     vmin=-2,
                                     vmax=4,
                                     ax=axs[i][num],
                                     legend_fontsize="xx-small")
                else:
                    sub = sc.pl.umap(adata,
                                     color=gene,
                                     save=False,
                                     show=False,
                                     size=40,
                                     ax=axs[i][num],
                                     legend_fontsize="xx-small")
                axs[i, num].axis('on')
        ########################################################################
        #  Append the geneset name
        ########################################################################
        pad = 5  # in points
        rows = list(dictionary.keys())
        for ax, row in zip(axs[:, 0], rows):
            ax.annotate(row, xy=(0, 0.5), xytext=(-ax.yaxis.labelpad - pad, 0),
                        xycoords=ax.yaxis.label, textcoords='offset points',
                        size='large', ha='right', va='center')
        ########################################################################
        # Plot umap
        ########################################################################
        sc.pl.umap(adata, color=all_genes, save=f"all_markers",
                   layer='scaled', vmin=-2, vmax=3)
        fig.savefig(f"{plot_path}/umap_genesets_{run_id}_"
                    f"{coloring}_{dataset_name}_{scaling}.pdf",
                    format="pdf")
        plt.close()
        # END OF FUNCTION


def adjust_dict(dictionary, list_valid):
    """

    :param dictionary: pyhton dict
    :param list_valid: list

    """
    myDict = {key: [e for e in dictionary[key]
                    if e in list_valid] for key in dictionary
              }
    return myDict


def score_stats(adata, cluster_identifier="", outputpath="",
                condition_obs="", scores=[]):
    """

    Function for scoring statistics

    :param adata: AnnData object
    :param cluster_identifier: str
    :param outputpath: str
    :param condition_obs: list/str
    :param scores: list

    :return: stats folder

    """

    ########################################################################
    # 0. Prepare folders
    ########################################################################
    outputpath = outputpath + "/statistics/"
    os.makedirs(outputpath, exist_ok=True)

    outputpath = outputpath + f"/{condition_obs}/"
    os.makedirs(outputpath, exist_ok=True)

    ########################################################################
    # 1. Iterate over bulk and individual cell populations
    ########################################################################
    adata.obs[scores + [condition_obs] + [cluster_identifier]].to_csv(f"{outputpath}/"
                                             f"sourcefile.csv")
    for celltype in ["all"] + adata.obs[
        cluster_identifier].value_counts().index.tolist():
        if celltype == "all":
            cell_df = adata.obs
        else:
            cell_df = adata.obs.loc[adata.obs[cluster_identifier]
                                    == celltype]

        ####################################################################
        # Export the source file
        ####################################################################
        cell_df[scores + [condition_obs]].to_csv(f"{outputpath}/"
                         f"sourcefile_{celltype}.csv")

        ####################################################################
        # Calc score statistics for each condition
        ####################################################################
        for score in scores:
            val_dict = {}
            kruskal_groups = []
            for e in cell_df[condition_obs].unique():
                sliced = cell_df.loc[
                    cell_df[condition_obs] == e]
                val_dict.update({e: sliced[score
                ].values.tolist()})
                kruskal_groups.append(
                    sliced[score].values.tolist())
            ################################################################
            # Kruskal Wallis test
            ################################################################
            if len(kruskal_groups) < 2:
                print("--- Too little groups")
                continue
            stats, p_value = scipy_stats.kruskal(*kruskal_groups)

            ################################################################
            # If groups differ, run post hoc test
            ################################################################
            if p_value < 0.05:
                p_adjust_test = 'bonferroni'
                kruskal_groups_for_posthoc = np.array(kruskal_groups)
                if len(kruskal_groups) < 3:
                    results = sp.posthoc_conover([kruskal_groups[0],
                         kruskal_groups[1]], p_adjust=p_adjust_test)
                else:
                    results = sp.posthoc_conover(kruskal_groups_for_posthoc,
                                                 p_adjust=p_adjust_test)
            results.columns = list(val_dict.keys())
            results.index = list(val_dict.keys())

            ################################################################
            # Export results
            ################################################################
            results.to_csv(f"{outputpath}/kruskal_{score}{celltype}.csv")
    # END OF FUNCTION


def plot_scores(adata, cluster_identifier="", on_raw=True,
                order=["CTRL", "OA500", "PA", "TGFB1"],
                interest=["sc_score_26-gene_signature",
                          "sc_score_98-gene_signature"],
                splitters=["controlstatus","detailed_condition", ],
                replot_umaps=False):
    """

    Visualize scores on projection.

    :param adata: AnnData object
    :param cluster_identifier: str
    :param on_raw: bool
    :param order: list
    :param splitters: list
    :param replot_umaps: bool

    :return: plots are being generated

    """
    sns.set_style("whitegrid", {'axes.grid': False})
    root = sc.settings.figdir

    ########################################################################
    # 1. Generate the scores
    ########################################################################
    for status in ["all_genes"]: #, "hvg_only"]:
        print(f"--- Scoring {status} data.")
        plot_path = f"{root}/{status}/"
        os.makedirs(plot_path, exist_ok=True)
        tmp_file = f"{root}/tmp_file_{status}.h5ad"
        os.makedirs(plot_path, exist_ok=True)
        sc.settings.figdir = plot_path

        ####################################################################
        # 1.1. Scores on all available genes
        ####################################################################
        if status == "all_genes":
            if not os.path.isfile(tmp_file):
                print("--- Temporary file not found")
                # Temporarily go back to raw count data (retrieve all genes)
                tmp = adata.raw.to_adata()
                tmp.raw = tmp  # Store raw slot
                # Normalize and scale
                sc.pp.normalize_total(tmp, exclude_highly_expressed=True,
                                      max_fraction=0.1)
                sc.pp.log1p(tmp)
                print("--- Scaling")
                tmp.layers['scaled'] = sc.pp.scale(tmp, copy=True,
                                                   max_value=5).X
                tmp.write_h5ad(tmp_file)
                print(f"Wrote to {tmp_file}")

            if os.path.isfile(tmp_file):
                print("--- Found all genes file (tmp)")
                tmp = sc.read_h5ad(tmp_file)
                if len([e for e in list(tmp.obs) if "sc_score" in e]) == 0:
                    # Score
                    tmp = score_geneset(tmp, on_raw=on_raw, interest=["26-gene_signature",
                                                                      "98-gene_signature"])
                    tmp.write_h5ad(tmp_file)
        else:
            if not os.path.isfile(tmp_file):
                ##############################################################
                # 1.2. Scores on hvg only
                ##############################################################
                if len([e for e in list(adata.obs) if "sc_score" in e]) == 0:
                    tmp = score_geneset(adata, on_raw=on_raw)
                    tmp.write_h5ad(tmp_file)
            if os.path.isfile(tmp_file):
                tmp = sc.read_h5ad(tmp_file)
        ####################################################################
        # Include minus top 5% genes to make sure they don't solely
        # account for the effect
        ####################################################################
        for scenario in ["all"]: #, "minus_top5"]:
            interest_sub = [f"{e}_{scenario}_{on_raw}"
                            for e in interest]
            plot_path = f"{root}/{status}/{scenario}/"
            os.makedirs(plot_path, exist_ok=True)

            #################################################################
            # 1.3 Calculate score statistics (26, 98 signatures)
            #################################################################
            for splitter in splitters:
                score_stats(tmp, scores=interest_sub,
                            cluster_identifier=cluster_identifier,
                            condition_obs=splitter,
                            outputpath=plot_path)

                ##############################################################
                # 1.4 Plot overview UMAPs (all scores in one plot)
                ##############################################################
                plot_path = f"{root}/{status}/{scenario}/{splitter}/"
                if replot_umaps:
                    os.makedirs(plot_path, exist_ok=True)
                    sc.settings.figdir = plot_path
                    print(f"--- Saving to {plot_path}")

                    for condition in tmp.obs[
                        splitter].value_counts().index.tolist():
                        cond_sliced = tmp[
                                      tmp.obs[splitter] == condition, :]
                        interest = [e for e in list(cond_sliced.obs)
                                    if "sc_score" in e and
                                    scenario in e and str(on_raw) in e] + \
                                   ["condition",
                                    "controlstatus",
                                    "phase",
                                    cluster_identifier]
                        for e in [0.2, 0.4, 1.0, 2]:
                            sc.pl.draw_graph(cond_sliced,
                                             color=interest,
                                             save=condition + f"_{e}",
                                             cmap=cmap,
                                             palette=cond_palette,
                                             vmin=0,
                                             vmax=e)
                        print(f"--- Plotted overview for {condition}.")

            ##################################################################
            # 1.5 Plot genes dotplot & violin plots
            ##################################################################
            sc.settings.figdir = plot_path
            for i, geneset in enumerate(geneset2path):
                print(f"--- Plotting gene dotplot {geneset}.")
                # 1.1 Read the custom genelist from config
                try:
                    df = pd.read_csv(geneset2path[geneset])
                except:
                    df = pd.read_table(geneset2path[geneset])
                # 1.2 Generate a boolean variable in the df
                geneset_list = df.iloc[:, 0].tolist()
                # 1.3  Generate the genelist object
                geneset_list = remove_duplicates(
                    adjust_list(geneset_list, tmp.var['gene_ids-1']))

                # Make dict
                gene_dict = {i: geneset_list}
                genes = gene_dict[i]
                if len(genes) < 1:
                    continue

                for celltype in ["all_cells"] + tmp.obs[
                    cluster_identifier].value_counts().index.tolist():
                    plot_path = f"{root}/{status}/{scenario}/{celltype}/"
                    os.makedirs(plot_path, exist_ok=True)
                    sc.settings.figdir = plot_path

                    # Slice to cell type of
                    if celltype != "all_cells":
                        tmp_sliced = tmp[tmp.obs[cluster_identifier]
                                         == celltype, :]
                    else:
                        tmp_sliced = tmp

                    for bl in [True, False]:
                        # Plot gene expression
                        sc.pl.matrixplot(tmp_sliced,
                                         genes,
                                         groupby="controlstatus",
                                         save=f"CONT_{geneset}_scaled_{bl}",
                                         cmap="Reds",
                                         vmin=0,
                                         vmax=0.3,
                                         layer="scaled",
                                         dendrogram=True,
                                         use_raw=bl,
                                         )
                        sc.pl.dotplot(tmp_sliced,
                                      genes,
                                      groupby="controlstatus",
                                      save=f"CONT_{geneset}_scaled_{bl}",
                                      cmap="Reds",
                                      vmin=0,
                                      vmax=0.3,
                                      layer="scaled",
                                      dendrogram=True,
                                      use_raw=bl,
                                      )
                        sc.pl.dotplot(tmp_sliced,
                                      genes,
                                      groupby="controlstatus",
                                      save=f"overview_{geneset}_raw_{bl}",
                                      cmap=cmap,
                                      standard_scale="var",
                                      dendrogram=True,
                                      use_raw=bl,
                                      )
                        sc.pl.dotplot(tmp_sliced,
                                      genes,
                                      groupby="controlstatus",
                                      save=f"overview_{geneset}_scaled_{bl}",
                                      cmap=cmap,
                                      layer="scaled",
                                      dendrogram=True,
                                      use_raw=bl,
                                      )
                        sc.pl.dotplot(tmp_sliced,
                                      genes,
                                      groupby="controlstatus",
                                      save=f"CONT_{geneset}_scaled_{bl}",
                                      cmap="Reds",
                                      vmin=0,
                                      layer="scaled",
                                      dendrogram=True,
                                      use_raw=bl,
                                      )

                        sc.pl.matrixplot(tmp_sliced,
                                         genes,
                                         groupby="controlstatus",
                                         save=f"CONT_{geneset}_raw_{bl}",
                                         cmap="Reds",
                                         vmin=0,
                                         standard_scale="var",
                                         dendrogram=True,
                                         use_raw=bl,
                                         )
                        sc.pl.dotplot(tmp_sliced,
                                      genes,
                                      groupby="controlstatus",
                                      save=f"CONT_{geneset}_raw_{bl}",
                                      cmap="Reds",
                                      vmin=0,
                                      dendrogram=True,
                                      use_raw=bl,
                                      )
                    # Plot the score as a violin
                    sc.pl.violin(tmp_sliced,
                                 keys=interest_sub,
                                 groupby="controlstatus",
                                 rotation=90,
                                 save=f"_{geneset}_{celltype}",
                                 palette=[GROUP_COLS[e] for e in order],
                                 stripplot=False,
                                 order=order,
                                 inner="box")
    # END OF FUNCTION


def integration_plot(adata, dataset_name="", keys_to_plot=False,
                     outputpath="", umap=False, scores=False,
                     rerank=False):
    """

    General function for plotting of integrated scRNA-seq data.

    :param adata: AnnData object
    :param dataset_name: str
    :param keys_to_plot: list, should be in adata.obs
    :param outputpath: str
    :param UMAP: bool, whether to repeat time consuming plotting
    :param scores: bool, whether to score (time consuming)
    :param rerank: bool, whether export DGE (time consuming)

    :return: outputfolder with plots

    """

    ########################################################################
    # 0. Path, settings
    ########################################################################
    root_dir = outputpath + f"/INTEGR/{dataset_name}/"
    sc.set_figure_params(scanpy=True, fontsize=5)
    os.makedirs(f"{outputpath}/QC", exist_ok=True)
    infofile = f"{outputpath}/QC/metrics.txt"
    info = f"{adata.shape}"
    for e in [e for e in list(adata.obs) if "condition" in e
                                            or "sample" in e or
                                            "control" in e]:
        info += f"{adata.obs[e].value_counts()}\n"

    for e in [e for e in list(adata.obs) if len(adata.obs[e].value_counts())
                                            < 50]:
        info += f"{adata.obs[e].value_counts()}\n"
    fi = open(infofile, "w")
    fi.write(info)
    fi.close()

    ########################################################################
    # 1.1 Plot scores (on key of interest)
    ########################################################################
    if "composite" in list(adata.obs) and "D21" not in outputpath and \
            scores:
        score_dir = outputpath + "/SCORES/"
        os.makedirs(score_dir, exist_ok=True)
        for e in [False]:
            print(f"--- Scoring on raw {e}")
            outdir = f"{score_dir}/raw_{e}/"
            sc.settings.figdir = outdir
            os.makedirs(outdir, exist_ok=True)
            plot_scores(adata, cluster_identifier="composite",
                        on_raw=e,
                        interest= ["sc_score_26-gene_signature",
                                   "sc_score_98-gene_signature"],)

    ########################################################################
    # 1.2 QC Plot, metrics
    ########################################################################
    print("--- QC plot")
    qc_plot(adata, sample_id=dataset_name, group_by="sample",
            outputpath=outputpath, foldername="QC")

    ########################################################################
    # 1.3 Overview plot
    ########################################################################
    sc.settings.figdir = root_dir
    key_genes = [e for e in GENERAL_MARKERS if e in adata.var.index]
    geneset_adjusted = reduce_to_genes_expressed(adata, genesets_lineage)
    genesets_plots_order = reduce_to_genes_expressed(adata, genesets_plots)
    genesets_plots_big_order = reduce_to_genes_expressed(adata, genesets_plots_big)
    for gset in [geneset_adjusted, genesets_plots_order, genesets_plots_big_order]:
        for cellype in gset:
            for gene in gset[cellype]:
                if gene not in key_genes:
                    key_genes.append(gene)

    keys_of_interest = []
    for ci in [e for e in list(adata.obs) if "leiden" in e]:
        umap_vars.append(ci)
        keys_of_interest.append(ci)

    us = [e for e in umap_vars if e in list(adata.obs)]

    ########################################################################
    # Define general parameters & plot on UMAP
    ########################################################################
    interest = key_genes + us
    if umap:
        if not os.path.isfile(f"{root_dir}draw_graph_fa_{dataset_name}"
                              f"_scaled-cut.pdf"):
            # Plot them on UMAP
            sc.pl.umap(adata,
                       color=interest,
                       vmin=-4, vmax=4,
                       layer="scaled",
                       cmap=cmap,
                       save=f"_{dataset_name}"
                            f"_scaled-cut")

            sc.pl.umap(adata,
                       color=interest,
                       layer="scaled",
                       cmap=cmap,
                       save=f"_{dataset_name}"
                            f"_scaled")

            # Plot them on graph
            sc.pl.draw_graph(adata,
                             color=interest,
                             save=f"_{dataset_name}",
                             layer="scaled")

            sc.pl.draw_graph(adata,
                             color=interest,
                             save=f"_{dataset_name}"
                                  f"_scaled-cut",
                             layer="scaled",
                             cmap=cmap,
                             vmin=-4, vmax=4)
            print("--- Saved overview plots.")

    if not keys_to_plot:
        keys_to_plot = [e for e in list(adata.obs) if
                        ("leiden" in e or
                         "composite" in e) and
                        "counts" not in e]
        keys_to_plot.sort()

    ########################################################################
    # 2. Cluster-key specific plots
    ########################################################################
    print(f"--- Plots for: {keys_to_plot}")
    for cluster_identifier in keys_to_plot:
        sc.settings.figdir = f"{root_dir}{cluster_identifier}"
        os.makedirs(f"{root_dir}{cluster_identifier}", exist_ok=True)

        if "Zonation_score" not in cluster_identifier:
            if rerank:
                sc.tl.rank_genes_groups(adata, cluster_identifier,
                                        method="wilcoxon",
                                        key_added="test")
                genes = pd.DataFrame(adata.uns["test"]['names'])
                csv_path = f"{root_dir}{cluster_identifier}/" \
                           f"gene_ranks_wilcoxon_allvsall.csv"
                genes.to_csv(csv_path)

        # RAW GENE EXPRESSION PLOT
        for i, gset in enumerate([genesets_lineage, genesets_plots]):
            for e in gset:
                sc.pl.draw_graph(adata,
                                 color=gset[e],
                                 cmap=cmap,
                                 save=f"_{i}_{e}_scaled")
        if rerank:
            csv_path = f"{root_dir}{cluster_identifier}/" \
                       f"gene_ranks_wilcoxon.csv"
            rank_key = f"wilcoxon_{cluster_identifier}"
            # Rank genes
            sc.tl.rank_genes_groups(adata, cluster_identifier,
                                    method="wilcoxon",
                                    key_added=rank_key,
                                    use_raw=True)

            # Export the differentially expressed genes
            genes = pd.DataFrame(adata.uns[rank_key]['names'])
            genes.to_csv(csv_path)

            print(f"--- UMAP for {cluster_identifier}")
            sc.pl.umap(adata,
                       color=cluster_identifier,
                       save=cluster_identifier,
                       cmap="rocket"
                       )
            sc.pl.draw_graph(adata,
                             color=cluster_identifier,
                             save=cluster_identifier,
                             cmap="rocket")

        if "Zonation_score" not in cluster_identifier:
            umap_cells_custom_color(adata, cluster_identifier)
            sc.pl.umap(adata,
                       color=cluster_identifier,
                       save=cluster_identifier,
                       )
            sc.pl.draw_graph(adata,
                             color=cluster_identifier,
                             save=cluster_identifier)
        else:
            return "Done"

        ########################################################################
        # 3. Cell fractions
        ########################################################################
        plot_path = outputpath + f"/INTEGR/{dataset_name}/{cluster_identifier}" \
                                 f"/PROPROTIONS/"
        os.makedirs(plot_path, exist_ok=True)
        sc.settings.figdir = plot_path
        plot_diff_props(adata, cluster_identifier=cluster_identifier,
                        plot_path=plot_path)

        plot_cellproportions_per_cluster(adata, cluster_identifier=cluster_identifier,
                                         plot_path=plot_path, dataset_name=dataset_name)

        ########################################################################
        # 4. Marker genes
        ########################################################################
        print(f"--- Marker genes {cluster_identifier}")
        for i, gset in enumerate([geneset_adjusted, genesets_plots_order,
                                  genesets_plots_big_order]):
            sc.pl.dotplot(adata, gset, cluster_identifier,
                          dendrogram=True,
                          cmap=cmap,
                          layer='scaled',
                          save=f"{cluster_identifier}"
                               f"_KEYMARKERS_{i}")
            sc.pl.matrixplot(adata, gset,
                             cluster_identifier,
                             dendrogram=True,
                             cmap=cmap,
                             layer='scaled',
                             save=f"{cluster_identifier}"
                                  f"_KEYMARKERS_{i}")
        # END OF FUNCTION


def qc_plot(adata, sample_id="", group_by=False, outputpath="", foldername=""):
    """

    Plotting the pre/post filtered QC params.

    :param adata: AnnData object
    :param sample_id: str
    :param group_by: bool
    :param outputpath: str
    :param foldername: str

    :return: folder containing QC plots

    """

    ########################################################################
    # 0. Path, settings
    ########################################################################
    DR_path = f"{outputpath}{foldername}"
    os.makedirs(DR_path, exist_ok=True)
    sc.settings.figdir = DR_path
    sc.pl.highest_expr_genes(adata, n_top=10, palette=palettediff2)
    sc.settings.set_figure_params(fontsize=11)

    ########################################################################
    # Assing sample name as observation
    ########################################################################
    adata.obs["sample"] = sample_id

    vals = [e for e in ['n_genes_by_counts',
                        'total_counts',
                        'pct_mito', 'pct_ribo',
                        'pct_lnc', 'pct_mir',
                        "predicted_doublet",
                        "doublet_score"] if e in
            list(adata.obs)]

    ########################################################################
    # Generate QC plots.
    ########################################################################
    if group_by:
        sc.pl.violin(adata, vals,
                     jitter=0.4,
                     groupby=group_by,
                     save=f"_{sample_id}",
                     palette=palettediff2,
                     show=False,
                     stripplot=False)
        sc.pl.scatter(adata,
                      x='total_counts',
                      y='pct_mito',
                      save=f"_pct_counts_mt_{sample_id}",
                      color=group_by,
                      palette=palettediff2,
                      frameon=False)
        sc.pl.scatter(adata,
                      x='total_counts',
                      y='pct_ribo',
                      save=f"_pct_counts_rb_{sample_id}",
                      color=group_by,
                      palette=palettediff2)
        sc.pl.scatter(adata,
                      x='total_counts',
                      y='n_genes_by_counts',
                      save=f"_n_genes_by_counts_{sample_id}",
                      color=group_by,
                      palette=palettediff2)
        sc.pl.scatter(adata,
                      x='pct_mito',
                      y='pct_ribo',
                      color=group_by,
                      palette=palettediff2)
        sc.pl.scatter(adata,
                      x='total_counts',
                      y='doublet_score',
                      color=group_by,
                      palette=palettediff2)
    else:
        for elem in vals:
            sc.pl.violin(adata,
                         elem,
                         jitter=0.4,
                         save=f"_{elem}_{sample_id}",
                         palette=palettediff2,
                         stripplot=False)
        sc.pl.scatter(adata,
                      x='total_counts',
                      y='pct_mito',
                      save=f"_pct_counts_mt_{sample_id}",
                      palette=palettediff2)
        sc.pl.scatter(adata,
                      x='total_counts',
                      y='pct_ribo',
                      save=f"_pct_counts_mt_{sample_id}",
                      palette=palettediff2)
        sc.pl.scatter(adata,
                      x='total_counts',
                      y='n_genes_by_counts',
                      save=f"_n_genes_by_counts_{sample_id}",
                      palette=palettediff2)
        sc.pl.scatter(adata,
                      x='pct_mito',
                      y='pct_ribo',
                      palette=palettediff2)
    # END OF FUNCTION
# END OF FILE
