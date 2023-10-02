"""

Python script for score heatmaps
@author: Anja Hess
@date: 2023-JUN-10

"""
import os.path
import matplotlib
import matplotlib.pyplot as plt
import scanpy as sc
import seaborn as sns
import pandas as pd
import glob
matplotlib.rcParams['figure.figsize'] = [6, 6]
cmap = "Spectral_r"
scores = ['sc_score_98-gene_signature_all_False',
          'sc_score_26-gene_signature_all_False']



def plotlabel(ax, xvar, yvar, label, pval=False):
    """

    """

    if label < 1:
        label = "↓"
    else:
        if label > 1:
            label = "↑"
        else:
            if label == 1:
                label = "="

    if pval:
        if label >= 0.05:
            label = "     $\it{ns}$"
        else:
            if label > 0.01:
                label = "     *".format(label)
            else:
                if label > 0.001:
                    label = "     **".format(label)
                else:
                    label = "     ***".format(label)
            #label = "     $\it{p}$ = "+"{:.2E}".format(label)
    ax.text(xvar, yvar, label,
            fontsize=12)

def score_mtx_plot(path_to_mtx="", cluster_identifier='composite_simple',
                   condition_identifier="controlstatus", cell_obs="composite",
                   pval_col="mann-whitney_pval", bool_col='p<0.05'):
    """

    :param path_to_mtx: str
    :param file_id: str
    :param cluster_identifier: obs slot
    :return: Matrix plot for score visualization across cell types /
    conditions

    Normalization see: https://github.com/holtzy/The-Python-Graph-Gallery/
    blob/master/src/notebooks/94-use-normalization-on-seaborn-heatmap.ipynb

    """
    GROUP_COLS = {"CTRL": "lightgrey",
                  "CTRL-OA": "lightgrey",
                  "CTRLOA": "lightgrey",
                  "CTRL-PA": "lightgrey",
                  "CTRL-TGFB1": "lightgrey",
                  "PA": "#477b80",
                  "TGFB1": "#d56763",
                  "OA500": "#fcd2a1",
                  }

    df = pd.read_csv(path_to_mtx)
    order_to_cond = {
        "detailed_condition": ["CTRL-OA", "OA500", "CTRL-PA", "PA", "CTRL-TGFB1", "TGFB1"] ,
        "controlstatus": ["CTRL", "OA500", "PA", "TGFB1"]
    }
    diff2controlcondition = {
        "OA500": "CTRL-OA",
        "PA": "CTRL-PA",
        "TGFB1": "CTRL-TGFB1"
    }

    # Loop through comparison modes

    for condition_obs in order_to_cond:
        order = order_to_cond[condition_obs]
        for score in df.columns:
            if "score" not in score:
                continue
            score_id = score.rsplit("sc_score_")[1].rsplit("_")[0]
            means_file = path_to_mtx.rsplit('/',1)[0]+ f"/means_{score}.csv"

            ################################################################
            # Generate p-values (pairwise comparison)
            ################################################################
            pvals = []
            for condition in diff2controlcondition:
                ctrl = diff2controlcondition[condition]
                comp = f"{condition}_by{ctrl}"
                cond_vals = df.loc[(df[condition_obs] == condition)
                                   ][score].values.tolist()
                ctrl_vals = df.loc[(df[condition_obs] == ctrl)

                                   ][score].values.tolist()
                try:
                    stats, p_value = scipy_stats.mannwhitneyu(cond_vals,
                                                              ctrl_vals)
                    pvals.append(["Bulk", comp, p_value])
                    print("bulk", comp, p_value)
                except:
                    continue

                for celltype in df[cell_obs].unique():
                    cond_vals = df.loc[(df[condition_obs] == condition)
                                       & (df[cell_obs] == celltype)
                                       ][score].values.tolist()
                    ctrl_vals = df.loc[(df[condition_obs] == ctrl)
                                       & (df[cell_obs] == celltype)
                                       ][score].values.tolist()
                    try:
                        stats, p_value = scipy_stats.mannwhitneyu(cond_vals,
                                                                  ctrl_vals)
                    except:
                        continue
                    pvals.append([celltype, comp, p_value])

            ################################################################
            # 1. Heatmap (calc mean scores per cell type and BULK)
            ################################################################
            df_mean = df.groupby([cell_obs, condition_obs]
                                 )[score].mean()

            col_list = [GROUP_COLS[e] for e in order]
            sns.violinplot(data=df, y=score, x=condition_obs,
                           order=order,
                           palette=col_list)
            plt.title(score_id)
            plt.savefig(path_to_mtx.rsplit('/',1)[0] +
                        f"/vio_{score_id}.pdf")
            plt.show()
            plt.close()

            df_mean_bulk = df.groupby(condition_obs)[score].mean()
            df_mean_bulk = pd.concat([df_mean_bulk], keys=['Bulk'], names=[cell_obs])
            df_mean = pd.concat([df_mean, df_mean_bulk], axis=0)
            df_mean.to_csv(means_file)
            df_mean = pd.read_csv(means_file)
            print(df_mean)
            df_h = df_mean.pivot(cell_obs, condition_obs,
                                 values=score)
            ################################################################
            # 1.1 Normalized by control
            ################################################################
            df_norm = df_h
            del df_norm["CTRL"]
            df_norm.dropna(inplace=True)

            order_norm = []
            for condition in diff2controlcondition:
                ctrl = diff2controlcondition[condition]
                df_norm[f"{condition}_by{ctrl}"] = (df_norm[condition]
                                                    / df_norm[ctrl])
                order_norm.append(f"{condition}_by{ctrl}")

            df_norm.to_csv((path_to_mtx.rsplit('/',1)[0] +
                            f"/heatmap_source_by-ctrl_{score_id}.csv"))
            for e in df_norm.columns:
                if e not in order_norm:
                    del df_norm[e]
            ################################################################
            # Bubble scatter (long format)
            ################################################################
            df_norm = pd.DataFrame(df_norm)
            norm_var = "Score treatment / Score control"
            df_norm["composite"] = df_norm.index
            df_long = df_norm.melt(value_vars=df_norm.columns.tolist(),
                                   id_vars=cell_obs,
                                   value_name=norm_var)

            # Map in the pvalues.
            df_long[pval_col] = ""
            for seq in pvals:
                idx = df_long.index[
                    (df_long[cell_obs] == seq[0])
                    & (df_long[condition_obs] == seq[1])
                    ].tolist()
                df_long.loc[idx, pval_col] = seq[2]
                # END OF LOOP
            df_long[pval_col] = df_long[pval_col].astype(float)
            df_long['p<0.05'] = np.where(
                df_long[pval_col] < 0.05, "*", "ns")
            print(df_long[pval_col])
            df_long["label"] = df_long[pval_col]
            df_long.to_csv(path_to_mtx.rsplit('/',1)[0] +
                        f"/bubble-source_{score_id}.csv",)

            ax = sns.scatterplot(data=df_long, y=cell_obs, x=condition_obs,
                                 hue=norm_var, palette="Spectral_r",
                                 s=400, linewidth=1, edgecolor="black",
                                 center=1)
            plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5), ncol=1)

            df_long.apply(lambda x: plotlabel(ax, x[condition_obs],
                                              x[cell_obs], x[norm_var]),
                          axis=1)

            plt.xticks(rotation=90)
            plt.grid(axis='y')
            plt.title(score_id)

            plt.savefig(path_to_mtx.rsplit('/',1)[0] +
                        f"/bubble-by-ctrl_{score_id}.pdf",
                        bbox_inches="tight")
            plt.show()
            plt.close()
            continue
            #exit()

            ################################################################
            # 1.1 Normalized by row ( = Cell type)
            ################################################################
            df_norm_row = df_h.apply(lambda x: (x - x.mean()) / x.std(),
                                     axis=1)
            sns.heatmap(df_norm_row, cmap=cmap,
                        cbar_kws={'label': f"normalized \n {score_id}"})
            plt.title(score)
            plt.tight_layout()
            plt.savefig(path_to_mtx.rsplit('/',1)[0]+f"/{score_id}_heatmap.pdf")
            plt.close()

            sns.barplot(data=df,  x=condition_obs, y=score,
                        palette="Paired", order=order)
            plt.savefig(path_to_mtx.replace(".csv",
                                            f"{score_id}_bulk.pdf"))
            plt.close()
            sns.catplot(
                data=df, x=condition_obs, y=score,
                col_wrap=2,
                kind="bar", col="composite", aspect=.7,
                palette="Paired", order=order
            )
            plt.savefig(path_to_mtx.replace(".csv",
                                            f"{score_id}.pdf"))
            plt.close()
            kruskal_groups = []
            val_dict = {}
            for e in df[condition_obs].unique():
                if e == "CTRL":
                    continue
                sliced = df.loc[df[condition_obs] == e]
                print(len(sliced), e)
                kruskal_groups.append(
                    sliced[score].values.tolist())
                val_dict.update({e: sliced[score
                ].values.tolist()})

            if len(kruskal_groups) > 1:
                stats, p_value = scipy_stats.kruskal(*kruskal_groups)
                print(stats, p_value)

                if p_value < 0.05:
                    p_adjust_test = 'bonferroni'
                    kruskal_groups_for_posthoc = np.array(kruskal_groups)
                    if len(kruskal_groups) < 3:
                            results = sp.posthoc_conover([kruskal_groups[0],
                                                          kruskal_groups[1]],
                                                         p_adjust=p_adjust_test)
                    else:
                        results = sp.posthoc_conover(kruskal_groups_for_posthoc,
                                                     p_adjust=p_adjust_test)
                    results.columns = list(val_dict.keys())
                    results.index = list(val_dict.keys())
                    print(results)
                    results.to_csv(path_to_mtx.replace(".csv", f"{score}.csv"))
        # ENF OF FUNCTION


############################################################################
# START
############################################################################
dir = "./"
for file in glob.glob(dir+"*sourcefile.csv"):
    score_mtx_plot(path_to_mtx=file)

def score_mtx_plot(path_to_mtx="", cluster_identifier='composite_simple',
                   condition_identifier="controlstatus",):
    """

    :param path_to_mtx: str
    :param file_id: str
    :param cluster_identifier: obs slot
    :return: Matrix plot for score visualization across cell types /
    conditions

    Normalization see: https://github.com/holtzy/The-Python-Graph-Gallery/
    blob/master/src/notebooks/94-use-normalization-on-seaborn-heatmap.ipynb

    """

    ###################################################################
    # 1. Load the query data
    ###################################################################
    for score in scores:
        save_df = f"{score}_means.csv"

        ################################################################
        # 2. Make the df
        ################################################################

        if not os.path.isfile(save_df):
            adata = sc.read_h5ad(path_to_mtx)
            df = adata.obs.groupby([cluster_identifier,
                                    condition_identifier]
                                   )[score].mean()
            df.to_csv(f"{score}_means.csv")

        if os.path.isfile(save_df):
            df = pd.read_csv(save_df)
            ############################################################
            # 3. Heatmap
            ############################################################
            df = df.pivot(cluster_identifier, condition_identifier,
                          values=score)
            sns.heatmap(df, cmap=cmap,
                        cbar_kws={'label': score})
            plt.title(score)
            plt.tight_layout()
            plt.savefig(f"{score}.pdf")
            plt.savefig(f"{score}.jpg")
            plt.close()

            ############################################################
            # 4. Heatmap -  Normalize by row (cell type)
            ############################################################
            df_norm_row = df.apply(lambda x: (x - x.mean()) / x.std(),
                                   axis=1)
            sns.heatmap(df_norm_row, cmap=cmap,
                        cbar_kws={'label': f"normalized \n {score}"})
            plt.title(score)
            plt.tight_layout()
            plt.savefig(f"{score}_norm.pdf")
            plt.savefig(f"{score}_norm.jpg")
            plt.close()
    # ENF OF FUNCTION


############################################################################
# START
############################################################################
dir = "./"
file_id = "tmp_file_all_genes.h5ad"
score_mtx_plot(path_to_mtx=f"{dir}/{file_id}")

# END OF SCRIPT
