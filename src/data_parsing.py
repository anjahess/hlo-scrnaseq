"""

Module containing basic functions for loading and parsing scRNAseq data.

@author: Anja Hess
@date: 2022-OCT-01

"""

import scanpy as sc
import pandas as pd
import os
from sys import exit
from cell_qc import easy_scrublet

def load_10X(path_to_matrix, sample_name=False, scrublet=False):
    """

    Routine for 10X data loading into AnnData.

    :param sample_name: str
    :param path_to_matrix: str, path/to/10X/outs/filtered_feature_bc_matrix
    :return: AnnData object

    """

    save_file = f"" \
                f"{path_to_matrix.rsplit('/filtered_feature_bc_matrix', 1)[0]}" \
                f"/{sample_name}.h5ad"

    print(f"--- Saving to {save_file}")

    if not os.path.isfile(save_file):
        adata = sc.read_10x_mtx(path_to_matrix +
                                "/filtered_feature_bc_matrix",
                                var_names='gene_symbols',
                                cache=True,
                                gex_only=True)
        if sample_name:
            adata.obs['sample'] = sample_name

        if scrublet:
            scrub_file = f"{path_to_matrix}/" \
                         f"{sample_name}_doublet-scores.csv"

            if os.path.isfile(scrub_file):
                print("--- Parsing scrublet info.")
                meta = pd.read_csv(scrub_file,
                                   header=0)
                del meta["Unnamed: 0"]
                adata.obs["barcode"] = adata.obs.index.tolist()
                merged = pd.merge(adata.obs,
                                  meta,
                                  how='left',
                                  on=["barcode", "barcode"])
                merged.set_index("barcode", inplace=True)
                adata.obs = merged
                print(adata)
            else:
                print("--- Warning: no scrublet info.")
        adata.write_h5ad(save_file)
    return adata
    # END OF FUNCTION


def export(adata, method="", outputpath=""):
    """

    Export metadata, eg. for CellPhoneDB
    :param adata: AnnData object
    :param method:sre
    :param outputpath:str
    :return:

    """

    if type(adata) == str:
        adata = sc.read_h5ad(adata)

    observations = list(adata.obs)

    outputpath = f"{outputpath}/METADATA"
    os.makedirs(outputpath, exist_ok=True)
    for elem in observations:
        if "leiden" in elem or \
                "composite" in elem:
            meta_file = f"{outputpath}/" \
                        f"metadata_{method}_{elem}.txt"
            if not os.path.isfile(meta_file):
                metadata = adata.obs[elem]
                metadata.to_csv(meta_file,
                                sep="\t",
                                header=["cell_type"],
                                index_label="Cell")

    print(f"--- Exported metadata to {outputpath}")
    # END OF FUNCTION


def find_cell(cells, target_cells, all=False):
    """

    To find all cell clusters belonging to a desired group

    :param cells: list of all cluster names
    :param root_cell: str or list
    :param all: bool
    :return:

    """
    hits = []
    if type(target_cells) == str:
        target_cells = [target_cells]

    for candidate in cells:
        for target in target_cells:
            if target in candidate:
                if all:
                    hits.append(candidate)
                else:
                    hit_cell = candidate
                    return hit_cell
    if not all:
        if not hits:
            return False
    return hits
    # END OF FUNCTION


def merge_sc(sample2path_dict, outputpath="", data_path=False,
             overwrite=False, mode="inner", scrublet_analysis=True):
    """

    :param file_name: name of the output file
    :param data_path:
    :param outputpath:
    :param run_id:
    :param sample2path_dict: sample_id:"PATH/TO/SAMPLE"
    :param scrublet_analysis: bool
    :return: adata object, merged if more than one input sample matrix

    """

    common_matrix = []
    batch_categories = []

    #########################################################################
    # 1. Set the outputpath
    #########################################################################
    save_file = outputpath + "-merged.h5ad"

    if os.path.isfile(save_file) and not overwrite:
        print(f"--- Merged file already created: {save_file}")
        return save_file
    else:
        print(f"--- Not found. Creating it.")

    #########################################################################
    # 2. Perform doublet detection if not done yet
    #########################################################################
    if scrublet_analysis:
        print("--- Performing doublet detection")
        for e in sample2path_dict:
            file_path = sample2path_dict[e]
            easy_scrublet(file_path,
                          sample_name=e,
                          outputpath=file_path)

    #########################################################################
    # 3. Load the h5ad
    #########################################################################
    for study in sample2path_dict:
        print(f"--- Adding {study}")
        batch_categories.append(study)
        if data_path:
            if not scrublet_analysis:
                if ".h5ad" in sample2path_dict[study]:
                    # Simply use the h5ad in dict if present
                    h5ad_file = sample2path_dict[study]
                    h5ad_folder = sample2path_dict[study].rsplit("/", 1)[0]
                else:
                    h5ad_folder = sample2path_dict[study]
                    h5ad_file = f'{h5ad_folder}/{study}.h5ad'
            else:
                h5ad_folder = sample2path_dict[study]
        else:
            h5ad_file = sample2path_dict[study]

        scrub_file = f'{h5ad_folder}/{study}_doublet-scores.csv'
        meta_file = f'{h5ad_folder}/{study}_meta.csv'

        ######################################################################
        # 4. Now parse doublet info from the prev. generated csv
        ######################################################################
        if os.path.isfile(h5ad_file):
            single_matrix = sc.read_h5ad(h5ad_file)
            try:
                test = single_matrix.obs["predicted_doublet"]
            except:
                print("--- No scrublet info found. Looking for csv.")
                if os.path.isfile(scrub_file):
                    meta = pd.read_csv(scrub_file,
                                       header=0)
                    del meta["Unnamed: 0"]
                    single_matrix.obs["barcode"] = \
                        single_matrix.obs.index.tolist()
                    merged = pd.merge(single_matrix.obs,
                                      meta,
                                      how='left',
                                      on=["barcode", "barcode"])
                    merged.set_index("barcode",
                                     inplace=True)
                    merged.to_csv(meta_file)
                    single_matrix.obs = merged
                else:
                    print("--- Continuing without scrublet info.")
                    exit()

            if type(common_matrix) == list:  # First matrix will be main matrix
                common_matrix = single_matrix
                del single_matrix
            else:
                single_matrix.obs["sample"] = study
                common_matrix.obs_names_make_unique()
                common_matrix.var_names_make_unique()
                common_matrix = common_matrix.concatenate(single_matrix,
                                                          join=mode)
                del single_matrix
        else:
            print(f"--- Did not find {h5ad_file}."
                  f"--- Creating and saving to {h5ad_file}")
            print(f"--- Loading as mtx from {h5ad_folder}")

            single_matrix = load_10X(f"{h5ad_folder}",
                                     sample_name=study)
            single_matrix.write_h5ad(h5ad_file)
            print(f"--- Wrote to {h5ad_file}")

    sample_no = common_matrix.obs[
        "sample"].value_counts().index.tolist()
    print(common_matrix.obs["sample"])

    if len(sample_no) < 2:
        print("--- Error: less than two samples. "
              "... Merging did not work")
        exit()

    common_matrix.write_h5ad(save_file)
    print(f"Merging complete. Saved to {save_file}")

    del common_matrix
    return save_file
    # END OF FUNCTION

# END OF SCRIPT
