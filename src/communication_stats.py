"""
Statistical evaluation of cell-cell interactions
across conditions

@author: Anja Hess
@date: 2023-MAY-01

"""

import scipy.stats as scipy_stats
import pandas as pd
import numpy as np

root = "./"

df = pd.read_csv(f"{root}/normalized_interactions.csv")

info = f"interacting_pair\tcondition" \
       f"\tcondition-control\tstats" \
       f"\tpval\tsignificant\n"

for e in df["source_target"].unique():
    # Select value for statistics
    val_dict = {}
    kruskal_groups = []

    sliced = df.loc[df['source_target'] == e]
    ctrl_vals = [sliced["value_50"].values[0],
                 sliced["value_51"].values[0],
                 sliced["value_62"].values[0],
                 sliced["value_63"].values[0],
                 sliced["value_SM-L3XWE"].values[0],
                 sliced["value_SM-L3XWF"].values[0],
                sliced["value_SM-L3XWG"].values[0],
                 sliced["value_SM-L3XWH"].values[0]]
    tgfb1_vals = [sliced["value_52"].values[0],
                  sliced["value_53"].values[0]]
    pa_vals = [sliced["value_64"].values[0],
               sliced["value_65"].values[0]]
    oa_vals = [sliced["value_SM-L3XWI"].values[0],
               sliced["value_SM-L3XWJ"].values[0]]

    # 1. TGFB1
    stats, p_value = scipy_stats.mannwhitneyu(ctrl_vals,
                                              tgfb1_vals,
                                              method="exact")
    if p_value < 0.05:
        signi = True
    else:
        signi = False
    info += f"{e}\tTGFB1\t" \
            f"{np.average(tgfb1_vals) - np.average(ctrl_vals)}" \
            f"\t{stats}\t{p_value}\t{signi}\n"

    # 2. PA
    stats, p_value = scipy_stats.mannwhitneyu(ctrl_vals,
                                              pa_vals,
                                              method="exact")
    if p_value < 0.05:
        signi = True
    else:
        signi = False

    info += f"{e}\tPA\t" \
            f"{np.average(pa_vals) - np.average(ctrl_vals)}" \
            f"\t{stats}\t{p_value}\t{signi}\n"

    # 3. OA
    stats, p_value = scipy_stats.mannwhitneyu(ctrl_vals, oa_vals)
    if p_value < 0.05:
        signi = True
    else:
        signi = False
    info += f"{e}\tOA\t" \
            f"{np.average(oa_vals) - np.average(ctrl_vals)}" \
            f"\t{stats}\t{p_value}\t{signi}\n"

fi = open(f"{root}mann_whitney_u_stats.tsv", "w")
fi.write(info)
fi.close()

# END OF SCRIPT
