"""

Processing and plotting functions for CellPhoneDB
outputs (cell-cell interactions in scrnaseq analysis)
@author: Anja Hess
@date: 2022-OCT-01

"""

import sys
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pyvenn.venn import venn, generate_petal_labels
import holoviews as hv
from holoviews import opts, dim
from bokeh.plotting import output_file

sns.set_context("paper")
sns.set_style("whitegrid",
              {'axes.grid': False})
script_path = str(os.path.dirname(os.path.abspath(__file__)))
maindir = script_path.split("/tools/plotting")[0]
src_path = f"{maindir}/src"
bin_path = f"{maindir}/bin"
sys.path.insert(0, script_path)
sys.path.insert(0, maindir)
sys.path.insert(0, src_path)
sys.path.insert(0, bin_path)
sys.path.insert(0, os.path.abspath('../plotting'))
sys.path.insert(0, os.path.abspath(f'../'))
sys.path.insert(0, os.path.abspath(f'../src'))
cwd = os.getcwd()
absFilePath = os.path.abspath(__file__)
script_path = str(os.path.dirname(os.path.abspath(__file__)))

############################################################################
# 1. Helper dictionaries
############################################################################
diff2controlcondition = {
    "PA": "CTRL-PA",
    "OA500": "CTRL-OA",
    "TGFB1": "CTRL-TGFB1"
}

CELL_HYPERCATS = {'HB1': '11_Hepatocytes', 'HB2': '2_Hepatocytes',
                  'FH1': "3_Hepatocytes", 'FH2': '4_Hepatocytes',
                  'cAH': '4_Hepatocytes', 'AH': '5_Hepatocytes',
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
                  'Smooth muscle cells': '9999_Fibroblasts',
                  'SMCs': '9999_Fibroblasts',
                  }

cell2abbr = {"HB1": "HB1", "Hepatic stellate cells": "HSCs",
             "HB2": "HB2", "Fibroblasts": "FIBs",
             "Cholangiocytes": "CHOLs", "FH1": "FH1",
             "Smooth muscle cells": "SMCs",
             "Activated hepatic stellate cells": "HSCs-ACT",
             "FH2": "cAH", "Ductal cells": "DCs",
             "AH": "AH"}

CELL_COLS = {"epato": ["#2d435b", "#436488", "#537ca9", "#bcd4e6",
                       "#bccae6", "#c3bce6", "#bccae6", "#e4e1f4",
                       "#9184d1", "#ffb6c1", "#d6dfea", "#ebb5b3"],
             "AH": ["#2d435b", "#436488", "#537ca9", "#bcd4e6",
                    "#bccae6", "#c3bce6", "#bccae6", "#e4e1f4",
                    "#9184d1", "#ffb6c1", "#d6dfea", "#ebb5b3"],
             "HB1": ["#477b80", "#bcd4e6", "#bccae6", "#c3bce6",
                     "#bccae6", "#e4e1f4", "#9184d1", "#ffb6c1",
                     "#d6dfea", "#ebb5b3"],
             "HB2": ["#bcd4e6"],
             "FH": ["#6b9ac8", "#96a2ae", "#3e74a8"],
             "yoepithelial": ["#477b80", "darkgreen", "#8b4513"],
             "epatoblast": ["#ffb6c1", "#d6dfea",
                            "#ebb5b3"],
             "holangio": ["#d56763", "#a02f2b", "#c73a35",
                          "#d25f5b", "#dd8481", "#e7aaa8",
                          "#f2d0ce"], "CHOLs": ["#a02f2b"],
             "ibroblast": ["#fcd2a1", "#965304", "#c56d06",
                           "#f58707", "#f89e35", "#fab564",
                           "#fbcb94", "#fde2c3"],
             "stellate": ["#fcd2a1", "#965304", "#c56d06",
                          "#f58707", "#f89e35", "#fab564",
                          "#fbcb94", "#fde2c3"],
             "HSCs": ["#fcd2a1", "#965304", "#c56d06",
                      "#f58707", "#f89e35", "#fab564",
                      "#fbcb94", "#fde2c3"],
             "FIBs": ["black"],
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

GROUP_COLS_VENN = {
    "CTRL-PA": "lightgrey", "TGFB1": "#d56763",
    "PA": "#477b80",
    "OA500": "#fcd2a1",
}

cmap = sns.diverging_palette(220, 20, as_cmap=True)


def assign_me_a_color(cell, colors_used=[]):
    """
    To return color for a cell type

    :param cell: str
    :param colors_used: list
    :return: str

    """
    counter = 0

    for celltype in CELL_COLS:
        if celltype in cell:
            color = CELL_COLS[celltype][0]
            while color in colors_used:
                counter += 1
                color = CELL_COLS[celltype][counter]

            return color

    # END OF FUNCTION


def heatmap(path_to_heatmap, vmaxes={"count": 20, "fraction_all_counts": 0.018},
            percent_id="fraction_all_counts",
            fixed_order=False, title=""):
    """

    Plotting heatmap from cellphonedb output

    :param path_to_heatmap: str
    :param vmaxes: dict
    :param percent_id: str
    :param fixed_order: bool
    :param: title: str
    :return: plots
    """

    # Load data
    df = pd.read_csv(path_to_heatmap,
                     sep="\t")
    df = df.sort_values("SOURCE")
    df[percent_id] = df["count"] / df["count"].sum()
    header = df.SOURCE.unique()
    columns = len(header)

    # Plot relative and absolut interaction numbers
    for e in [percent_id, "count"]:
        # Path settings
        save_file = path_to_heatmap.rsplit(".txt", 1)[0] \
                    + f"_{e}_{fixed_order}_{title}.pdf"
        save_file = save_file.replace(title + "/", "overview/")

        if fixed_order:
            hmp_df = df.pivot("SOURCE", "TARGET", e)
            sns.heatmap(hmp_df,
                        vmax=vmaxes[e],
                        cmap=cmap,
                        xticklabels=True,
                        yticklabels=True)
            plt.xticks(
                fontsize=4,
                fontweight='normal')
            plt.yticks(fontsize=4,
                       fontweight='normal')
        else:
            counts = df[e].tolist()
            list_of_counts = [counts[x:x + columns] for x in
                              range(0, len(counts), columns)]
            count_array = np.array([np.array(xi)
                                    for xi in list_of_counts])
            hmp_df = pd.DataFrame(data=count_array[0:, 0:],  # values
                                  index=header,  # 1st column as index
                                  columns=header)
            sns.clustermap(hmp_df,
                           vmax=vmaxes[e],
                           cmap=cmap)
        plt.title(f"Interactions {title} ({e})")
        plt.savefig(save_file)
        plt.close()
    # END OF FUNCTION


def chord(df, root, id="", select=False, adjust=True, multi=False):
    """

    Holoviews chord diagram for cell-cell interactions from
    CellPhoneDB output.

    :param df: pandas dataframe
    :param root: str
    :param id: str
    :param select: bool
    :param adjust: bool
    :param multi: bool

    :return: Chord diagrams are generated, chord object is returned.

    """

    ############################################################################
    # 0. Define, dfs, names, and dirs
    ############################################################################
    df.columns = ["source", "target", "value"]
    if multi:
        df_dir = root + f"/normalized_interactions_{id}.csv"
    else:
        df_dir = root + f"/normalized_interactions.csv"
        basic_dir = root + "/interactions.csv"

    ############################################################################
    # 1. Simplify labels
    ############################################################################
    for e in cell2abbr:
        if e == select:
            select = cell2abbr[e]

    ############################################################################
    # 2. Normalize by total interactions
    ############################################################################
    if not os.path.isfile(df_dir):
        new_df = df.copy()
        if adjust:
            new_df["value"] = new_df["value"] / \
                              (new_df["value"].sum() / 100)
            new_df["source"] = new_df["source"
            ].str.split("_", expand=True)[0]
            new_df["target"] = new_df["target"
            ].str.split("_", expand=True)[0]
            new_df["source"] = new_df["source"].map(cell2abbr)
            new_df["target"] = new_df["target"].map(cell2abbr)
        new_df.to_csv(df_dir)
        df.to_csv(basic_dir)
        return ""

    if os.path.isfile(df_dir):
        new_df = pd.read_csv(df_dir,
                             index_col=0)

    ############################################################################
    # 3. Select cell of interest for plotting
    ############################################################################
    new_df["to_top"] = new_df["source"].isin([select])
    new_df = new_df.sort_values("to_top")
    # Sort for interactions affecting the cell of interest
    df1 = new_df.loc[(new_df["target"].isin([select]))]
    df2 = new_df.loc[(new_df["source"].isin([select]))]
    # CRITICAL! Drop duplicates
    new_df = pd.concat([df1,
                        df2]).drop_duplicates().reset_index(
        drop=True)
    # Exclude cells without interactions
    new_df = new_df[new_df["value"] > 0]

    ############################################################################
    # 4. Assign group numbers to the cell types (for compatibility)
    ############################################################################
    node_info = pd.DataFrame(new_df["target"].value_counts())
    node_info["name"] = node_info.index
    node_info["hypercat"] = node_info["name"].map(CELL_HYPERCATS
                                                  ).fillna("Unknown")

    ############################################################################
    # 5.  Sort df & convert to holoviews object
    ############################################################################
    node_info.sort_values(by="hypercat",
                          ascending=False,
                          inplace=True)
    # 5.1 Add the numerical index
    index = pd.RangeIndex(start=0,
                          stop=len(node_info),
                          step=1,
                          name="index")
    node_info.index = index
    # 5.2 Assign each group an int id
    node_info["group"] = index
    # Optional:
    # self_id = node_info[node_info["name"] == select]["target"].values[0]
    # prev_id = node_info[node_info["name"] == select]["group"].values[0]
    # node_info['group'] = node_info['group'].replace(prev_id, self_id)
    # 5.3 Create a dict with the info
    sample2condition = node_info.set_index("name"
                                           )["group"].to_dict()
    # 5.4 Convert to an hv object
    nodes = hv.Dataset(node_info, 'index')
    new_df["source"] = new_df["source"].map(sample2condition).astype(
        'category')
    new_df["target"] = new_df["target"].map(sample2condition).astype(
        'category')

    ############################################################################
    # 5. Add each cell type's color
    ############################################################################
    colors = []
    for cell in sample2condition:
        col = assign_me_a_color(cell, colors)
        if col == None:
            print(cell, " without specified color.")
            col = "lightgrey"
        colors.append(col)

    ############################################################################
    # 6. In a second loop grey out to avoid color overlap
    ############################################################################
    for i, cell in enumerate(sample2condition):
        if select:
            if cell != select:
                if (i % 2) == 0:
                    colors[i] = "lightgrey"
                else:
                    colors[i] = "#dfdfdf"
    hv.extension("matplotlib")
    hv.output(fig='pdf', size=250)

    ############################################################################
    # 7. Plot the chord
    ############################################################################
    c = hv.Chord((new_df, nodes))
    c.opts(opts.Chord(cmap=colors,
                      node_size=1,
                      edge_linewidth=10,
                      edge_cmap=colors,
                      edge_color=dim('source').str(),
                      labels='name',
                      node_color=dim('index').str())
           ).opts(title=id)

    # Save as pdf
    hv.save(c, root + f"/chord_overview_{id}.pdf")
    return c

    # END OF FUNCTION


def receptor_dotplot(path_to_heatmap, id="", filter_pval=True,
                     heatmap=True):
    """

    Dot-plot from CellPhoneDB outputs, other plotting functions.

    :param path_to_heatmap: str
    :param id: str
    :param filter_pval: Only interactions with p < 0.05
    :param heatmap: bool

    :return: plots are generated.

    """

    ############################################################################
    # 0. Define paths
    ############################################################################
    result_dir = f"{path_to_heatmap.rsplit('/', 1)[0]}" \
                 f"/results"

    chord_dir = f"{path_to_heatmap.rsplit('/', 1)[0].rsplit('/', 1)[0]}" \
                f"overview/{id}.pdf"

    print(f"--- Saving to {result_dir}")
    os.makedirs(result_dir, exist_ok=True)

    ############################################################################
    # 1. Load the cellphoneDB output
    ############################################################################
    print(f"--- Reading p values for {id}")
    df = pd.read_csv(path_to_heatmap, sep="\t", index_col=0)
    pval_df = pd.read_csv(f"{path_to_heatmap.rsplit('/', 1)[0]}/pvalues.txt",
                          sep="\t", index_col=0)

    ############################################################################
    # 2. Generate chord diagrams
    ############################################################################
    chord_df = pd.read_csv(f"{path_to_heatmap.rsplit('/', 1)[0]}"
                           f"/count_network.txt", sep="\t")
    chord_df_bu = chord_df
    common = []
    for celltype in chord_df_bu["SOURCE"
    ].value_counts().index.tolist():
        fig = chord(chord_df, result_dir,
                    select=celltype,
                    id=f"{celltype}_{id}")
        if type(common) == list:
            common = fig
        else:
            common = common + fig
    hv.save(common, chord_dir)

    ############################################################################
    # 2. Generate more human readable cell pair names
    ############################################################################
    cell_pairs = [e.split('|')[0].split('_')[0] + '|' +
                  e.split('|')[1].split('_')[0]
                  for e in df.columns.tolist() if "|" in e]

    # remove unneeded columns
    for i, frame in enumerate([df, pval_df]):
        for column in frame.columns:
            if "|" not in column and "interacting_pair" not in column:
                del frame[column]

    df.columns = ['interacting_pair'] + cell_pairs
    pval_df.columns = ['interacting_pair'] + cell_pairs

    ############################################################################
    # 3. P-value filtering (remove pairs where no individual p-val < 0.05)
    ############################################################################
    if filter_pval:
        print("--- Filtering for pairs with min.1 p-val < 0.05.")
        pairs = [e for e in pval_df.columns.tolist() if "|" in e]
        len_pairs = len(pairs)
        pval_df['row_min'] = pval_df.min(axis=1)
        alpha = 0.05
        # if the smallest value is above alpha, there is no sign. interaction
        pval_df = pval_df[pval_df['row_min'] < alpha]
        del pval_df['row_min']
    print(f"--- Exporting p-value <{alpha} dataframe.")
    # Save this df
    pval_df.to_csv(f"{result_dir}/p_vals_signi.csv")

    ############################################################################
    # 4. Transform to sns long format & merge
    ############################################################################
    pval_df = pval_df.melt(id_vars=["interacting_pair"],
                           var_name="cell-cell_pair",
                           value_name=f"pval")

    df = df.melt(id_vars=["interacting_pair"],
                 var_name="cell-cell_pair",
                 value_name=f"mean")
    merged = pd.merge(pval_df, df, how='left')
    merged["sample"] = id

    ############################################################################
    # 5. Group
    ############################################################################
    print("--- Defining groups of receptor-ligand pairs.")
    merged["group"] = merged["interacting_pair"].str.split("_", 1,
                                                           expand=True)[0]
    merged.to_csv(f"{result_dir}/merged_min1signi.csv")
    # Even if selection was disabled here filter significant interactions
    merged_signi = merged[merged["pval"] < alpha]
    merged_signi.to_csv(f"{result_dir}/merged_signi.csv")
    counts = pd.DataFrame(merged_signi["group"].value_counts(),
                          index=None)
    counts["sample"] = id
    counts.to_csv(f"{result_dir}/sign_interact_per_group.csv")
    counts["name"] = counts.index

    ############################################################################
    # 6. Plot interactions
    ############################################################################
    sns.scatterplot(data=counts, x="name", y="group", palette="Blues")
    sns.despine(offset=5, trim=True)
    plt.grid(False)
    plt.xticks(rotation=90, fontweight='normal')
    plt.ylim(0)
    plt.title(f"""Partners with significant interactions {id}""")
    plt.tight_layout()
    plot_path = f"{result_dir}/sign_ia_by_count_{id}.pdf"
    plt.savefig(plot_path)
    plt.close()

    # END OF FUNCTION


def venn_from_dict(dct, outpath=""):
    """

    Multi group Venn diagram

    :param dct: python dictionary
    :param outpath: str

    :return: Venn plots are generated

    """

    summary_dct = {}

    ############################################################################
    # 0. Merge controls
    ############################################################################
    for val in dct:
        if "CTRL" in val:
            if "CTRL" not in summary_dct:
                summary_dct.update({"CTRL": dct[val]})
            else:
                # APPEND TO OTHER CTRL VALS
                summary_dct["CTRL"].update({
                    val for val in dct[val]
                    if val not in summary_dct["CTRL"]})
        else:
            summary_dct[val] = dct[val]
    cols = list(GROUP_COLS_VENN.values())

    ############################################################################
    # 1. Plot Venn
    ############################################################################
    venn(summary_dct, cmap=cols)
    petal_labels = generate_petal_labels(summary_dct.values(),
                                         fmt="{size}",
                                         len_petal_show=30)
    plt.title("Overview")
    plt.savefig(outpath)
    plt.close()
    venn(summary_dct, cmap=cols)
    petal_labels = generate_petal_labels(summary_dct.values(),
                                         fmt="{size}",
                                         len_petal_show=0)
    plt.title("Overview_label")
    plt.savefig(outpath.replace(".pdf", "_labels.pdf"))
    plt.close()

    ############################################################################
    # 2. Plot Venn for each condition separately
    ############################################################################
    for pair in diff2controlcondition:
        smalldict = {k: dct[k] for k in (pair, diff2controlcondition[pair])}
        venn(smalldict,
             cmap=list(GROUP_COLS.values()))
        petal_labels = generate_petal_labels(smalldict.values(),
                                             len_petal_show=0)
        plt.title(pair)
        plt.savefig(outpath.replace(".pdf",
                                    f"_{pair}.pdf"))
        plt.close()

    # END OF FUNCTION


def venn_summary(path_to_results, vars=[]):
    """

    Wrapper function for generating a venn diagram

    :param path_to_results: str
    :param vars: list

    :return: Start diagram

    """

    ############################################################################
    # 0. Define paths
    ############################################################################
    result_dir = f"{path_to_results.rsplit('/', 1)[0]}/overview/"
    os.makedirs(result_dir, exist_ok=True)
    result_dir = f"{path_to_results.rsplit('/', 1)[0]}/overview/VENN/"
    print(f"--- Saving to {result_dir}")
    os.makedirs(result_dir, exist_ok=True)

    ############################################################################
    # 1. Merge dataframes of interest
    ############################################################################
    for venn_df in ["p_vals_signi"]:
        initial_df = []
        venn_dict = {}
        for i, e in enumerate(vars):
            print(f"-- Adding {e}")
            df = pd.read_csv(f"{path_to_results}/{e}/results/{venn_df}.csv",
                             header=0)
            ####################################################################
            # 1.1 Get interactions specific for each condition/group
            ####################################################################
            if type(initial_df) == list:
                initial_df = df
                try:
                    initial_df.columns = ["group", "count", "sample"]
                    vals = initial_df["group"].values.tolist()
                except:
                    initial_df["comp"] = initial_df["interacting_pair"]
                    vals = initial_df["interacting_pair"].values.tolist()
                    initial_frac = pd.DataFrame(initial_df["interacting_pair"])
                    initial_frac["count"] = 1
                venn_dict.update({e: set(vals)})
            else:
                # Merg depending on which scheme you work with (group/pairs)
                try:
                    df.columns = ["group", "count", "sample"]
                    vals = df["group"].values.tolist()
                except:
                    vals = df["interacting_pair"].values.tolist()
                    frac = pd.DataFrame(df["interacting_pair"])
                    frac["count"] = 1
                venn_dict.update({e: set(vals)})
                try:
                    merge = pd.merge(initial_df,
                                     df,
                                     on="group",
                                     indicator=e,
                                     how="outer",
                                     suffixes=["_" + vars[i - 1], "_" + e])
                    initial_df = merge
                except:
                    merge = pd.merge(initial_frac, frac,
                                     on="interacting_pair",
                                     indicator=e,
                                     how="outer",
                                     suffixes=["_" + vars[i - 1], "_" + e])
                    initial_frac = merge
                    continue

        ########################################################################
        # 2. SUMMARY VENN
        ########################################################################
        venn_from_dict(venn_dict, outpath=f"{result_dir}/venn_{venn_df}.pdf")

        ########################################################################
        # 3. Export results
        ########################################################################
        # Add 0 where no count detected, then calc percentage of all counts
        for col in merge.columns:
            if "count" in col:
                merge[f"pct_total_{col}"] = merge[col] / \
                                            (merge[col].sum(axis=0) / 100)
                merge[col].fillna(0, inplace=True)
                merge[f"pct_total_{col}"].fillna(0, inplace=True)
            if "sample" in col:
                merge[col] = col.split("sample_")[1]
        save_merged = f"{result_dir}/counts_by_group_" \
                      f"{e}_{venn_df}_merged.csv"
        merge.to_csv(save_merged)
    # END OF FUNCTION


def interaction_summary(path_to_results, vars=[], p_val_filter=True):
    """

    Comparing number of interactions (counts) per group.

    :param path_to_results: str
    :param vars: list
    :param p_val_filter: bool

    :return: plots are generated.

    """

    ############################################################################
    # 0. Define paths
    ############################################################################

    initial_df = []
    savedir = path_to_results + "/overview/"
    merged_df = f"{savedir}merged_interactions.csv"
    os.makedirs(savedir, exist_ok=True)

    ############################################################################
    # 1. Merge dataframes of all conditions
    ############################################################################
    if not os.path.isfile(merged_df):
        for i, e in enumerate(vars):
            df = pd.read_csv(f"{path_to_results}{e}/results/merged_min1signi.csv",
                             header=0, index_col=0)
            sample = df["sample"].value_counts().index.tolist()[0]
            if "CTRL" in sample or "50" in sample or "51" in  sample or "62" \
                    in sample or "63" in sample:
                df["condition"] = "Control"
            else:
                df["condition"] = "Treated"
            df["sample"] = e
            if "-" in e:
                df["sample_group"] = e.split("-")[1]
                if e.split("-")[1] == "OA":
                    df["sample_group"] = "OA500"
            else:
                df["sample_group"] = e
            if type(initial_df) == list:
                initial_df = df
            else:
                merged = initial_df.append(df)
                initial_df = merged
        merged.to_csv(merged_df)

        if p_val_filter:
            merged = merged[merged["pval"] < 0.05]
            merged = merged.sort_values(by="interacting_pair")
        merged.sort_values(by="cell-cell_pair", inplace=True)
        merged["class"] = merged["cell-cell_pair"].str.rsplit('|', 1,
                                                              expand=True)[0]
        merged.to_csv(merged_df.replace('.csv', '_significant.csv'))

    ############################################################################
    # 2. Load df
    ############################################################################
    sns.set(font_scale=1)
    sns.set_style("whitegrid", {'axes.grid': False})
    if os.path.isfile(merged_df):
        print("--- Merged df exists. Loading.")
        all_df = pd.read_csv(merged_df, index_col=0)
        if "signi" not in all_df.columns:
            all_df['signi'] = np.where(all_df['pval'] < 0.05, 'yes', 'no')
            all_df["unique_interaction"] = all_df["interacting_pair"] + "_" \
                                           + all_df["cell-cell_pair"]
            all_df.to_csv(merged_df)

        #########################################################################
        # 3. Prepare receptor-ligand plots groups
        #########################################################################
        savedir = path_to_results + "/overview/GROUPS/"
        os.makedirs(savedir, exist_ok=True)

        # If looking at specific sets you can select with "|" operation
        interest = "TNFSF10|SPP1|PDGFB|NOTCH|CXCL2|_VEGF|TGFB1_AR|IL|COL"

        if interest:
            savedir = path_to_results + f"/overview/GROUPS/interest/"
            os.makedirs(savedir, exist_ok=True)

        #########################################################################
        # 4. Generate a dataframe for matrix-like plots (cells x pairs)
        #########################################################################
        if not os.path.isfile(savedir + "all_incl_missing_ia.csv"):
            # Add non-interaction value
            conds = all_df["sample"].value_counts().index.tolist()
            for ui in all_df["unique_interaction"
            ].value_counts().index.tolist():
                ui_df = all_df[all_df["unique_interaction"] == ui]
                for cond in conds:
                    if cond not in ui_df["sample"
                    ].value_counts().index.tolist():
                        rlp = ui.rsplit("_", 1)[0]
                        group = rlp.rsplit("_", 1)[0]
                        ccp = ui.rsplit("_", 1)[1]
                        if "CTRL" in cond:
                            trt = "Control"
                        else:
                            trt = "Treated"
                        ls = [rlp, ccp, 1, 0, cond, group, trt,
                              "other", "no", ui]
                        row = pd.Series(ls, index=all_df.columns)
                        all_df = all_df.append(row, ignore_index=True)
            all_df.sort_values(by=["signi", "mean"],
                               inplace=True)
            all_df.to_csv(savedir + "all_incl_missing_ia.csv")

        if os.path.isfile(savedir + "all_incl_missing_ia.csv"):
            all_df = pd.read_csv(savedir + "all_incl_missing_ia.csv")
            all_df.sort_values(by=["mean"], inplace=True)

        #########################################################################
        # 5. Now iterate through groups and mark induced pairs
        #########################################################################
        for cl in all_df["group"].value_counts().index.tolist():
            if interest:
                group_df = all_df[all_df["unique_interaction"
                ].str.contains(interest)]
                group_df.to_csv(savedir + "interest.csv")

                # 5.1 Remove non-significant (p < 0.05) interactions
                for ui in group_df["cell-cell_pair"
                ].value_counts().index.tolist():
                    if len(group_df[group_df["cell-cell_pair"] == ui][
                               "signi"].value_counts().index.tolist()) == 1 and \
                            group_df[group_df["cell-cell_pair"] == ui]["signi"
                            ].value_counts().index.tolist()[0] == "no":
                        group_df = group_df[group_df["cell-cell_pair"] != ui]
                    else:
                        # 5.2 Select induced pairs
                        info_dct = {}
                        for cond in group_df[group_df["cell-cell_pair"]
                                             == ui]["condition"
                        ].value_counts().index.tolist():
                            # 5.2.1 Generate an info df to help selct
                            info_df = group_df[(group_df["cell-cell_pair"]
                                                == ui)
                                               & (group_df["condition"]
                                                  == cond)]["signi"
                            ].value_counts()
                            # 5.2.2 Fill it with the significant interactions
                            if "yes" in info_df:
                                info_dct.update({cond: info_df["yes"]})
                            else:
                                info_dct.update({cond: 0})

                        # 5.3 Now take out reduced/similar interactions
                        if info_dct["Control"] >= info_dct["Treated"]:
                            group_df = group_df[group_df[
                                    "cell-cell_pair"] != ui]
                        delta = info_dct["Treated"] - info_dct["Control"]

                        # 5.4 Threshold for representation in plot (size of change)
                        if delta < 2:
                            group_df = group_df[
                                group_df["cell-cell_pair"] != ui]
            else:
                group_df = all_df[all_df["group"] == cl]

            #######################################################################
            # 6. Plot
            #######################################################################
            custom_order = ['TNF_NOTCH1', 'TNFRSF11B_TNFSF10', 'CXCL2_DPP4',
                            'ICAM1_AREG', 'RAreceptor_RXRA_atRetinoicAcid_byALDH1A1',
                            'SPP1_PTGER4', 'SPP1_CD44', 'NOTCH1_JAG1', 'NOTCH1_DLK1',
                            'NRP2_VEGFA', 'NRP1_VEGFA', 'PDGFR_complex_PDGFD',
                            'PDGFRB_PDGFD', 'PDGFB_ADGRV1', 'PDGFB_PDGFR_complex',
                            'PDGFB_PDGFRB', 'PDGFB_PDGFRA','TGFB1_AR']

            df_mapping = pd.DataFrame({'order': custom_order,})
            sort_mapping = df_mapping.reset_index().set_index('order')
            group_df['size_num'] = group_df["interacting_pair"
            ].map(sort_mapping['index'])
            group_df.sort_values(['size_num', "cell-cell_pair"], inplace=True)
            for e in ["sample"]:
                if e == "sample_group":
                    order_plot = ["OA500", "PA", "TGFB1"]
                else:
                    order_plot = ["CTRL-PA", "CTRL-OA", "CTRL-TGFB1", "PA",
                                  "OA500", "TGFB1"]

                p = sns.relplot(data=group_df,
                                y="interacting_pair",
                                x="cell-cell_pair",
                                hue="signi",
                                size="mean",
                                col=e,
                                col_wrap=3,
                                edgecolor="black",
                                rasterized=True,
                                sizes=(20, 130),
                                linewidth=0,
                                height=3.2,
                                aspect=1.8,
                                col_order=order_plot,
                                palette=["lightblue",
                                         "coral"],
                                facet_kws={'sharey': True,
                                           'sharex': True}
                                )
                plt.title(f"{cl} : all interactions")
                plt.xticks(rotation=90, fontsize=3)
                plt.yticks(fontsize=3)
                sns.despine(offset=5, trim=True)
                for ax in p.axes.flatten():
                    ax.set_xticklabels(
                        ax.get_xticklabels(), rotation=90,
                        ha='right', rotation_mode='anchor')
                    ax.tick_params(axis="y")
                plt.savefig(savedir + f"/GROUP_{e}_{cl}.pdf",
                            dpi=200, bbox_inches="tight")
                plt.close()
                del group_df

            if interest:
                break  # needed only once

    ############################################################################
    # 7. Another venn diagram
    ############################################################################
    merged = pd.read_csv(merged_df.replace('.csv', '_significant.csv'),
                         index_col=0)
    if "signi" not in merged.columns:
        merged['signi'] = np.where(merged['pval'] < 0.05, 'yes', 'no')
        merged["unique_interaction"] = merged["interacting_pair"] + "_" \
                                       + merged["cell-cell_pair"]
        merged.to_csv(merged_df.replace('.csv', '_significant.csv'))
    venn_dict = {}
    for group in merged["sample"]:
        new_df = merged[merged["sample"] == group]
        venn_dict.update({group: set(
            new_df["unique_interaction"].values.tolist())})

    ############################################################################
    # 8. Generate information on condition-exclusive interactions
    ############################################################################
    os.makedirs(path_to_results + f"/overview/VENN_UNIQUE/", exist_ok=True)
    df = []
    names = []
    all_ctrl_vals = [venn_dict[e] for e in venn_dict if "CTRL" in e]
    all_ctrl_set = set().union(*all_ctrl_vals)

    for cond_set in venn_dict:
        #########################################################################
        # 8.3 Interactions unique for treatments
        #########################################################################
        names.append(cond_set)
        all_other_vals = [venn_dict[e] for e in venn_dict if e != cond_set]
        all_other_set = set().union(*all_other_vals)
        intersection = venn_dict[cond_set].intersection(all_other_set)
        exclusive = venn_dict[cond_set] - intersection
        exclusive = list(exclusive)
        exclusive.sort()
        if type(df) == list:
            df = pd.DataFrame.from_dict({cond_set: exclusive})
        else:
            print(cond_set, len(exclusive))
            orig = pd.DataFrame.from_dict({cond_set: exclusive})
            df = pd.concat([df, orig], ignore_index=True, axis=1)

        #########################################################################
        # 8.2 Interactions shared between treatments
        #########################################################################
        if "CTRL" not in cond_set:
            for other_treatment in [e for e in venn_dict if
                                    "CTRL" not in e
                                    and e != cond_set]:
                ia_term = f"{cond_set}_and_{other_treatment}"
                names.append(ia_term)
                # 3.2.1 Shared between the two treatments:
                other_treatment_vals = [venn_dict[e] for e in
                                        venn_dict if e == other_treatment]
                other_set = set().union(*other_treatment_vals)
                intersection = venn_dict[cond_set].intersection(other_set)
                # 3.2.2 Of those, everything that is not in CTRLs
                exclusive = intersection - all_ctrl_set
                exclusive = list(exclusive)
                exclusive.sort()
                if type(df) == list:
                    df = pd.DataFrame.from_dict({f"{cond_set}_{other_treatment}"
                                                 : exclusive})
                else:
                    orig = pd.DataFrame.from_dict({f"{cond_set}_{other_treatment}"
                                                   : exclusive})
                    df = pd.concat([df, orig], ignore_index=True, axis=1)

    ############################################################################
    # 9. Export information on condition-exclusive interactions
    ############################################################################
    df.columns = names
    df.to_csv(f"{path_to_results}/overview/VENN_UNIQUE/unique_interactions.csv")
    venn_from_dict(venn_dict, f"{path_to_results}/overview/VENN_UNIQUE/"
                   f"interaction_venn.pdf")
    # END OF FUNCTION


def stacked_barplot(path_to_df, title=""):
    """

    To generate a stacked barplot displaying the fractions of interactions
    for varying conditions

    :param path_to_df: str
    :param title: str
    :return: the barplots are generated

    """

    ############################################################################
    # 0. Define paths, load data, prepare dataframe
    ############################################################################
    df = pd.read_csv(path_to_df, sep="\t", index_col=0)
    save_file = path_to_df.rsplit(".txt", 1)[0] + f"props_{title}.pdf"
    save_file = save_file.replace(title + "/", "overview/")
    df["name"] = df.index
    df = df.sort_values("name")
    colors = []
    for e in df["name"].values:
        col = assign_me_a_color(e, colors_used=colors)
        colors.append(col)
    del df["name"]
    df["fraction"] = df["all_sum"] / df["all_sum"].sum()
    del df["all_sum"]
    df = df.transpose()

    ############################################################################
    # 1. Plot the bar plot
    ############################################################################
    df.plot.bar(stacked=True,
                figsize=(2, 6),
                width=1,
                linewidth=0.01,
                edgecolor="black",
                rot=45,
                color=colors
                )
    plt.tight_layout()
    plt.title(title)
    plt.savefig(save_file)
    plt.close()

    # END OF FUNCTION


def find_changing_interactions(path_to_results, vars=[]):
    """

    Comparing number of relative interactions (counts) per group.

    :param path_to_results: str
    :param vars: list
    :return: csvs are created

    """

    ############################################################################
    # 0. Define paths, load data, prepare dataframe
    ############################################################################
    savedir = path_to_results + "/overview/"
    merged_df_dir = f"{savedir}normalized_interactions.csv"
    raw_df_dir = f"{savedir}interactions.csv"
    delta_df_dir = f"{savedir}delta_interactions.csv"
    savedir_delta = path_to_results + "/overview/changes/"
    os.makedirs(savedir, exist_ok=True)
    os.makedirs(savedir_delta, exist_ok=True)
    merged_df = []

    if not os.path.isfile(raw_df_dir):
        for i, e in enumerate(vars):
            df = pd.read_table(f"{path_to_results}{e}/count_network.txt",
                               header=0)
            df[f"value_{e}"] = df["count"]
            df["source_target"] = df["SOURCE"] + "_" + df["TARGET"]
            del df["SOURCE"]
            del df["count"]
            del df["TARGET"]
            if type(merged_df) == list:
                merged_df = df
            else:
                merged_df = pd.merge(merged_df, df, on="source_target",
                                     how='outer')
        merged_df.set_index("source_target", inplace=True)
        merged_df.to_csv(raw_df_dir)

    ############################################################################
    # 1. Merge dataframes of all conditions
    ############################################################################
    if not os.path.isfile(merged_df_dir):
        for i, e in enumerate(vars):
            df = pd.read_csv(f"{path_to_results}{e}/results/"
                             f"normalized_interactions.csv",
                             header=0, index_col=0)
            df[f"value_{e}"] = df["value"]
            df["source_target"] = df["source"] + "_" + df["target"]
            del df["source"]
            del df["value"]
            del df["target"]
            if type(merged_df) == list:
                merged_df = df
            else:
                # merge
                merged_df = pd.merge(merged_df,
                                     df,
                                     on="source_target",
                                     how='outer')
        merged_df.set_index("source_target", inplace=True)
        merged_df.to_csv(merged_df_dir)
        delta_df = merged_df
        for e in [e for e in vars if "CTRL" not in e]:
            ctrl = [i for i in vars if "CTRL" in i and e[:2] in i][0]
            print(f"-- Adding {e} vs {ctrl}")
            delta_df[f"{e}_vs_CTRL"] = delta_df[f"value_{e}"] - \
                                       delta_df[f"value_{ctrl}"]
        for e in delta_df.columns:
            if "vs" not in e:
                del delta_df[e]
        delta_df.to_csv(delta_df_dir)

    ############################################################################
    # 2. Plot clustermap
    ############################################################################
    if os.path.isfile(delta_df_dir):
        delta_df = pd.read_csv(delta_df_dir, header=0)
        for i, df in enumerate([delta_df_dir, merged_df_dir]):
            if os.path.isfile(df):
                sns.set(font_scale=0.5)
                n_df = pd.read_csv(df, index_col=0)
                sns.clustermap(n_df, cmap="magma", rasterized=True)
                plt.savefig(f"{savedir}changing_interactions_{i}.pdf", dpi=200)

        ########################################################################
        # 3. Plot chord (enforced / reduced)
        ########################################################################
        for comparison in delta_df.columns:
            if "vs" not in comparison:
                continue
            # Select comparison
            delta_df_sub = delta_df[["source_target", comparison]]
            # Prepare for chord
            delta_df_sub["SOURCE"] = delta_df_sub["source_target"].str.split(
                    "_", 1, expand=True)[0]
            delta_df_sub["TARGET"] = \
                delta_df_sub["source_target"].str.split( "_", 1, expand=True)[1]
            del delta_df_sub["source_target"]
            delta_df_sub = delta_df_sub.reindex(columns=["SOURCE", "TARGET",
                                                         comparison])
            savedir_delta_comp = path_to_results + f"/overview/changes/{comparison}/"
            os.makedirs(savedir_delta_comp, exist_ok=True)

            # 3.1 Select enforced or reduced interactions
            for variable in ["reduced", "enforced"]:
                print(f"--- {variable}")
                chord_dir = f"{savedir_delta_comp}" \
                            f"_{comparison}" \
                            f"_{variable}.pdf"
                if variable == "enforced":
                    delta_df_sub_var = \
                        delta_df_sub[delta_df_sub[comparison] > 0]
                else:
                    delta_df_sub_var = \
                        delta_df_sub[delta_df_sub[comparison] < 0]
                    delta_df_sub_var[comparison] = \
                        delta_df_sub_var[comparison].abs()
                common = []
                celltypes = delta_df_sub_var["SOURCE"].value_counts().index.tolist()
                celltypes.sort()

                # 3.2 Plot for each cell type
                for celltype in celltypes:
                    fig = chord(delta_df_sub_var, savedir_delta_comp,
                                select=celltype, adjust=False,
                                id=f"{celltype}_{comparison}_{variable}",
                                multi=True)
                    if type(common) == list:
                        common = fig
                    else:
                        common = common + fig
                hv.save(common, chord_dir)
    # END OF FUNCTION

############################################################################
# 3. START ANALYSIS
############################################################################

coi = ["CTRL-OA", "OA500", "PA", "CTRL-PA", "CTRL-TGFB1", "TGFB1", "CTRL"]
root = os.getcwd() + "/"
os.makedirs(root + "overview", exist_ok=True)
for elem in coi:
    receptor_dotplot(f"{root}/{elem}/means.txt", id=elem)
    stacked_barplot(f"{root}/{elem}/interaction_count.txt", title=elem)
    heatmap(f"{root}/{elem}/count_network.txt", fixed_order=True,
            title=elem)
venn_summary(root, vars=coi)
find_changing_interactions(root, vars=coi)
interaction_summary(root, vars=coi)
# END OF SCRIPT.
