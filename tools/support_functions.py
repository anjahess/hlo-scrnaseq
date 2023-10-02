"""

Support functions

@author: Anja Hess
@date: 2023-MAY-01

"""
import os
from shutil import copyfile
import pandas as pd
import scanpy as sc
from src.sc_constants import sample2condition, \
    sample2culture, sample2controlstatus, hep2loc, \
    sample2detailedcondition, CELL_COLS
from src.data_parsing import export
from tools.annotate import anno


def clean(adata, clust_path=""):
    """

    In case there are more var slots than needed.

    :param adata: AnnData object
    :param clust_path:
    :return: AnnData object

    """

    if "gene_ids" not in list(adata.var):
        for i in [e for e in list(adata.var)
                  if "gene_ids" in e]:
            print(adata.var[i])
            if i != "gene_ids-1-0" and i != "gene":
                del adata.var[i]

        adata.var["gene_ids"] = adata.var[
            "gene_ids-1-0"]
        del adata.var["gene_ids-1-0"]
    adata.write(clust_path)
    print("--- Wrote cleaned")
    # END OF FUNCTION


def simplify_labels(adata):
    """

    To unite subgroups of similar annotations.

    :param adata: AnnData object

    :return: AnnData object

    """
    for cluster_identifier in [
        e for e in list(adata.obs) if
        "leiden" in e and "simple" not in e
        and "2" in e]:
        try:
            adata.obs[cluster_identifier + "_simple"] = \
                adata.obs[
                    cluster_identifier
                ].str.split("_", expand=True)[0]
        except:
            continue

    # Manual annotation of stromal -> ductal
    adata.obs[cluster_identifier + "_simple"] \
        = adata.obs[cluster_identifier + "_simple"
                      ].replace(regex='Stromal cells',
                  value='Ductal cells')
    return adata
    # END OF FUNCTION


def meta_from_dict(adata):
    """

    Parse meta infos.

    :param adata: AnnData object
    :return: modified  AnnData object

    """
    adata.obs["condition"] = adata.obs["sample"] \
        .map(sample2condition).astype(
        'category')
    adata.obs["culture"] = adata.obs["sample"] \
        .map(sample2culture).astype(
        'category')
    adata.obs["controlstatus"] = adata.obs["sample"] \
        .map(sample2controlstatus).astype(
        'category')
    adata.obs["detailed_condition"] = adata.obs["sample"] \
        .map(sample2detailedcondition).astype(
        'category')
    # END OF FUNCTION



def zone_score(adata, cluster_identifier=""):
    """

    Classify liver zone score based on DGE
    described in MacParland, S. A. et al.
    Single cell RNA sequencing of human liver reveals
    distinct intrahepatic macrophage populations.
    Nat Commun 9, 4383 (2018).

    :param adata:
    :param cluster_identifier:
    :return: Score 1-6 (0-5)
    """

    print("--- Adding obs Zonation_score")
    adata.obs["zone"] = adata.obs[
        cluster_identifier
    ].str.rsplit("_", 1, expand=True)[0]
    adata.obs["zone"] = adata.obs[
        "zone"].str.rsplit("_ns", 1, expand=True)[0]
    adata.obs["Zonation_score"] = adata.obs["zone"] \
        .map(hep2loc)
    return adata
    # END OF FUNCTION


def split_adata_by_condition(adata_path, variable="",
                             variable_dict=False,
                             clean_up=False, cells=False,
                             re_anno=False):
    """

    To subset into individual h5ads
    Best on qc-filtered, raw data.

    :param adata_path:
    :param variable: the adata.obs slot to split on

    :return: each value in the adata.obs slot gets
    its own .h5ad.

    """

    adata = sc.read_h5ad(adata_path)
    root = adata_path.rsplit('/', 1)[0]

    # Map
    if variable_dict:
        adata.obs[variable] = adata.obs[
            "sample"].map(variable_dict).astype(
            'category')

    print(adata.obs[variable].value_counts())
    if "detailed_condition" not in list(
            adata.obs) and variable == \
            "detailed_condition":
        adata.obs["detailed_condition"] = \
            adata.obs["sample"].map(
                sample2detailedcondition
            ).astype('category')
        adata.write_h5ad(adata_path)
        print(adata.obs["detailed_condition"])
        print("--- Wrote")
        exit()

    splitted_files = []

    # Split & save each condition
    outdir = f"{root}/SUBSETS/"
    os.makedirs(outdir, exist_ok=True)
    outdir = f"{root}/SUBSETS/{variable}/"
    os.makedirs(outdir, exist_ok=True)

    for cond in adata.obs[
        variable].value_counts().index.tolist():
        split_name = adata_path.rsplit(
            '.h5ad', 1)[0] + f"_{cond}.h5ad"
        split_name = outdir + split_name.rsplit("/",
                                                1)[1]
        if cells:
            cont = False
            for cell in cells:
                if cell in cond:
                    print(f"{cond} in target: {cells}")
                    cont = True
                    break
            if not cont:
                continue

        if not os.path.isfile(split_name):
            print(f"--- Splitting {cond}")
            keep = (adata.obs[variable].isin([cond]))
            adata_split = adata[keep, :]
            print("--- Split ctrl:", adata_split.obs[
                variable].value_counts().index.tolist())

            if clean_up:
                print("--- CLEAN UP")
                adata_split.obs.drop(columns=[
                    e for e in list(adata_split.obs)
                    if "leiden" in e], inplace=True)
                for e in [e for e in
                          list(adata_split.uns)
                          if "leiden" in e]:
                    del adata_split.uns[e]

            adata_split.write_h5ad(split_name)

        if os.path.isfile(split_name):
            print(f"--- Split version exists: {cond}")
            adata_split = sc.read_h5ad(split_name)
            if re_anno:
                # Annotate
                print("Annotating")
                adata_split, keys_to_plot = anno(adata_split,
                                                 outputpath=
                                                 f"{root}/"
                                                 f"SUBSETS/"
                                                 f"{cond}/",
                                                 dataset_name=
                                                 "ANNOTATION",
                                                 )
                print("--- Wrote annotated.")
                adata_split.write_h5ad(split_name)

            # Save to condition-specific METADATA folder
            metadata_dir = f"{root}/SUBSETS/" \
                           f"{variable}/METADATA/"
            os.makedirs(metadata_dir,
                        exist_ok=True)
            metadata_dir = f"{root}/SUBSETS/" \
                           f"{variable}/METADATA/{cond}/"
            os.makedirs(metadata_dir,
                        exist_ok=True)
            export(adata_split,
                   outputpath=f"{metadata_dir}")

        splitted_files.append(split_name)
    print("All splitted.")
    return splitted_files
    # END OF FUNCTION


def create_composite(adata_path, reintegrate=[],
                     re_key=[], in_key=""):
    """

    To merge in sub-clustered ids

    :param adata: AnnData object
    :return: AnnData object

    """

    adata = sc.read_h5ad(adata_path)
    for e in ["new", "composite","composite_simple"]:
        if e in list(adata.obs):
            del adata.obs[e]

    # 1. Iterate through cell subsets
    for i, e in enumerate(reintegrate):
        print(f"--- Reintegrating {i}, {e}")

        prefixed = [filename for
                    filename in
                    os.listdir(adata_path.rsplit(
                        '/', 1)[0] + f"/SUBSETS/{in_key}/")
                    if filename.startswith(
                adata_path.rsplit(
                    '/', 1)[1].rsplit(".h5ad")[0] + f"_{e}")]

        if len(prefixed) > 1:
            print("--- Warning multiple files match!")
            exit()

        reference_data = adata_path.rsplit(
            '/', 1)[0] + f"/SUBSETS/{in_key}/" + prefixed[0]

        ref = sc.read_h5ad(reference_data)
        # 1.1 Generate the new key in both sets
        ref.obs["new"] = ref.obs[re_key[i]].astype(str)
        # So no category dilemma
        for e in list(ref.obs):
            if "new" not in e:
                del ref.obs[e]

        # 1.2 Merge both sets (yields new_x,y)
        adata.obs = pd.merge(adata.obs,
                             ref.obs,
                             on="barcode",
                             how="left")

        print(adata.obs["new"])
        print(adata.obs[in_key])
        print(list(adata.obs))
        print(list(ref.obs))

        # 1.3 Combine first to fill nans
        if i == 0:
            adata.obs["composite"] = adata.obs[
                "new"].combine_first(
                adata.obs[in_key]).astype("category")
            del adata.obs["new"]
            print(adata.obs["composite"].value_counts())
        else:
            adata.obs["composite"] = adata.obs[
                "new"].combine_first(
                adata.obs["composite"]
            ).astype("category")

    adata.obs["composite_simple"] = adata.obs[
        "composite"].str.split(
        "_", expand=True)[0]

    print(adata.obs["composite"
                    "_simple"].value_counts())
    print(adata.obs["composite"].value_counts())

    if "new" in list(adata.obs):
        del adata.obs["new"]
    adata.write_h5ad(adata_path)
    print("Saved composite.")
    return adata
    # END OF FUNCTION


def split_celltype_by_condition(adata, sample2condition,
                                old_celltype="",
                                condition_column=""):
    """

    :param sample2condition dict
    :param adata: h5ad object
    :param old_celltype: str, column in adata.obs
    with cell type labels to be splitted
    :param condition_column: str, column in adata.obs
    which will split the cells

    :return:

    """

    adata_str = adata
    adata = sc.read_h5ad(adata)
    adata.obs["condition"] = adata.obs[
        "sample"].map(sample2condition).astype(
        'category')
    print(adata.obs[old_celltype])
    print(adata.obs[condition_column])
    adata.obs[f"{old_celltype}_by_cond"] = \
        adata.obs[old_celltype
        ].astype(str).replace("_PA_FULL_", "") + \
        adata.obs[
            condition_column].astype(str)
    print(adata.obs[f"{old_celltype}_by_cond"])
    adata.write_h5ad(adata_str)

    return f"{old_celltype}_by_cond"
    # END OF FUNCTION


def save_script(script_title, output_path):
    """

    :param output_path: path to save the py script to
    :return:

    """

    SCRIPT_PATH = str(os.path.dirname(
        os.path.abspath(__file__)))
    os.makedirs(output_path,
                exist_ok=True)
    scripts = [os.path.abspath(SCRIPT_PATH +
                               "/" +
                               script_title)]
    for scriptfile in scripts:
        script_title = scriptfile.split("/")[-1]
        copyfile(scriptfile, f"{output_path}/{script_title}")
    # END OF FUNCTION
# END OF SCRIPT
