"""

Python script for intersection of pathway enrichment
@author: Anja Hess
@date: 2023-JUN-09

"""

import seaborn as sns
import pandas as pd
import os
cmap = sns.diverging_palette(220, 20,
                             as_cmap=True)


def intersect_pws(path_to_csvs="", file_id=""):
    """

    :param path_to_mtx: str
    :param file_id: str
    :param cluster_identifier: obs slot
    :return: Statistics for abundance

    """

    file = open(f"{path_to_csvs}intersect.csv", 'w')

    # Iterate through all cell types (=folders)
    for subdir, dirs, files in os.walk(path_to_csvs):
        for cell_type in dirs:
            print(cell_type)

            ###################################################################
            # 1. Load the query data
            ###################################################################
            df = pd.read_csv(f"{path_to_csvs}{cell_type}/{file_id}", index_col=0)

            set_oa = df[df["Condition"] == "OA500"]["Term"].values.tolist()
            set_pa = df[df["Condition"] == "PA"]["Term"].values.tolist()
            set_tgfb = df[df["Condition"] == "TGFB1"]["Term"].values.tolist()

            set_dict = {"tgf": set_tgfb,
                        "oa": set_oa,
                        "pa": set_pa}

            ###################################################################
            # 2. Make result file
            ###################################################################

            for i, e in enumerate(set_dict):
                cond_set = set(set_dict[e])
                all_other_vals = [set_dict[i] for i in set_dict if i != e]
                all_other_set = set().union(*all_other_vals)
                intersection = cond_set.intersection(all_other_set)
                exclusive = cond_set - intersection
                exclusive = list(exclusive)
                intersection = list(intersection)

                # Only needed once
                if i == 0:
                    file.write(f'shared_{cell_type}\t')
                    for inters in intersection:
                        file.write(f'{inters}\t')
                    file.write(f'\n')

                file.write(f'exclusive_{cell_type}_{e}\t')
                for excl in exclusive:
                    file.write(f'{excl}\t')
                file.write(f'\n')
        file.close()
        break
    # ENF OF FUNCTION


############################################################################
# START
############################################################################

# 1. OS-HLOs
dir = "./"
file_id = "sourcefile_WikiPathway_2021_Human_Adjusted P-value_induced.csv"
intersect_pws(path_to_csvs=dir, file_id=file_id)

# END OF SCRIPT
