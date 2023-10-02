"""

Python script for upset plot from cpdb output
@author: Anja Hess
@date: 2023-JUN-09

"""

import pandas as pd
import matplotlib.pyplot as plt
from upsetplot import UpSet, from_indicators, from_contents
from collections import defaultdict

def upset(path_to_csv=""):
    """

    :param path_to_csv: str
    :return: Upset plot

    """

    df = pd.read_csv(path_to_csv)
    df = df[["sample", "unique_interaction"]]

    # Combine all controls
    rename = {}
    for ctr in df["sample"].unique():
        if "CTRL" in ctr and ctr != "CTRL":
            rename.update({ctr: "CTRL"})
    df.replace(rename, inplace=True)

    ####################################################################
    # Sum data for Fig. 4f
    ####################################################################
    df_dict = df.groupby(['sample'])[
        "unique_interaction"].apply(list).to_dict()
    # remove dups:
    for cond in df_dict:
        df_dict[cond] = list(dict.fromkeys(df_dict[cond]))

    input = from_contents(df_dict)
    upset = UpSet(input, show_counts='{:d}', orientation='vertical')
    upset.plot()
    plt.savefig(
        f"{path_to_csv.replace('.csv', '_summary.pdf')}")
    plt.show()
    plt.close()
    exit()

    ####################################################################
    # Sum data for Fig. 4e
    ####################################################################
    df_final = pd.crosstab(df["unique_interaction"],
                           df["sample"]
                           ).ne(0).rename_axis(index='thing',
                                               columns=None)

    # get the ones NOT in ctrl
    df_final = df_final[df_final["CTRL"] != True]
    df_final = df_final[df_final["CTRL-OA"] != True]
    df_final = df_final[df_final["CTRL-PA"] != True]
    df_final = df_final[df_final["CTRL-TGFB1"] != True]
    df_final = df_final.transpose()
    print(df_final)

    # Order is CTRL, PA, OA, TGFB1
    df_final.to_csv(f"{path_to_csv.replace('.csv', '_source.csv')}")
    UpSet(from_indicators(
        lambda df: df.select_dtypes(bool),
        data=df_final),
        min_subset_size=1,
        show_counts=True).plot()
    plt.savefig(f"{path_to_csv.replace('.csv', '.pdf')}")
    plt.show()
    plt.close()
    # ENF OF FUNCTION


############################################################################
# START
############################################################################
dir = "./"
file_id = "merged_interactions_significant.csv"
upset(path_to_csv=f"{dir}{file_id}")

# END OF SCRIPT
