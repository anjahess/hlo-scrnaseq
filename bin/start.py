"""

MAIN SCRIPT TO START SCRNA-SEQ ANALYSIS
@author: Anja Hess
@date: 2022-OCT-03

# 1. For splitting a subset (of QC, PP file):
start --split COM all_scrub_clean_0.5

# 2. For subsets of subset of QC/PP data:
# 2.1 OS-HLOs:
start.py --harmony -nopp \
all_scrub_clean_0.5_qc_dr_OS.h5ad OS_harmony

# 2.2 FOR ULA-HLOs:
start.py --harmony -nopp \
all_scrub_clean_0.5_qc_dr_Day21.h5ad D21

"""

import os
import sys
import argparse
import scanpy as sc

path_delim = "/"
script_path = str(os.path.dirname
                  (os.path.abspath(__file__)))
maindir = script_path.split(path_delim + "bin")[0]
src_path = f"{maindir}{path_delim}src"
bin_path = f"{maindir}{path_delim}bin"
bash_path = f"{maindir}{path_delim}bash"
sys.path.insert(0, script_path)
sys.path.insert(0, maindir)
sys.path.insert(0, src_path)
sys.path.insert(0, bin_path)
sys.path.insert(0, os.path.abspath('.'))
sys.path.insert(0, os.path.abspath('..' + path_delim))
sys.path.insert(0, os.path.abspath('..' + path_delim + 'src'))
cwd = os.getcwd()
absFilePath = os.path.abspath(__file__)
script_path = str(os.path.dirname
                  (os.path.abspath(__file__)))

from src.sc_constants import MIN_CELLS, \
    MITO_CUTOFF, COUNT_LOWER, \
    DOUBLET_SCORE, NMARKERS, TARGET_CHOL, \
    TARGET_HEPATO, RES_CHOL, RES_HEPATO, \
    MAIN_ANNO, ANNO_DICT_HEPATO, ANNO_DICT_CHOL, \
    GENES_BY_COUNT_FRAME, sample2condition
from src.core_functions import easy_dr, easy_pp, \
    harmony, paga, easy_clust, GO_terms, hts_easy
from src.cell_qc import easy_qc, easy_scrublet
from src.imputation_control import plot_imp_norm
from src.data_parsing import merge_sc, load_10X
from tools.plotting.plots import integration_plot
from tools.annotate import anno, pathways,\
    hepato_zonescoring
from tools.support_functions import \
    split_adata_by_condition, meta_from_dict, \
    simplify_labels, create_composite, clean

#########################################################################
# Parse the folder with .h5ads / cellranger outs.
# This folder should contain all samples you wish
# to integrate into ONE .h5ad object
#########################################################################
parser = argparse.ArgumentParser(description=
                                 'Analyse 10X')
parser.add_argument('path',
                    metavar='P',
                    type=str,
                    nargs='+',
                    help='Path to folder')
parser.add_argument('name',
                    metavar='N',
                    type=str,
                    nargs='+',
                    help='Name of composite')
parser.add_argument('-r', '--raw',
                    default=False,
                    required=False,
                    dest='raw',
                    action='store_true')
parser.add_argument('-harm', '--harmony',
                    default=False,
                    required=False,
                    dest='harm',
                    action='store_true')
parser.add_argument('-sp', '--split',
                    default=False,
                    required=False,
                    dest='split',
                    action='store_true')
parser.add_argument('-nopp', '--nopp',
                    default=False,
                    required=False,
                    dest='nopp',
                    action='store_true')
parser.add_argument('-noclust', '--noclust',
                    default=False,
                    required=False,
                    dest='noclust',
                    action='store_true')
parser.add_argument('-ss', '--singlesample',
                    default=False,
                    required=False,
                    dest='single_sample',
                    action='store_true')

args = parser.parse_args()
screen_folder = args.path[0]
mergename = args.name[0]

# If working with single 10X outputs
if args.raw and args.single_sample:
    # Here no need to merge so get the h5ad right away
    easy_scrublet(screen_folder,
                  sample_name=mergename,
                  outputpath=screen_folder)
    load_10X(screen_folder,
             sample_name=mergename,
             scrublet=True)
    screen_folder = screen_folder + f"/{mergename}.h5ad"


#########################################################################
# IMPORTANT SETTINGS (DEFAULT INPUT: .h5ad files)
#########################################################################
outpath = f"{screen_folder.rsplit('/', 1)[0]}" \
          f"/processing/"
os.makedirs(outpath, exist_ok=True)
matrixpath = f"{outpath}/data/raw/datasets/"
all_data_path = f"{outpath}/data/"
overwrite = False
cpdb = True
re_anno = False
re_clust = False
GO = True
clean_up = False

#########################################################################
# 1. SAMPLE MERGING (OPTIONAL)
#########################################################################
print(f""" 
        --- STARTING ANALSYSIS ON {mergename} ---
        """)

merge_dict = {}

if ".h5ad" in screen_folder:
    # Option 1: single .h5ad file (merging done)
    print("--- ! Input is an .h5ad. Assuming "
          "all samples to be merged.")
    id = screen_folder.rsplit(".h5ad", 1)[0]
    if "/" in id:
        id = id.rsplit("/", 1)[1]
    merge_dict.update({id: screen_folder})
else:
    print("No h5ad, loading from mtx")
    merge_str = ''
    merged_sample_path = ''
    # Iterate through the folder
    for dirpath, dirnames, file_names \
            in os.walk(screen_folder):
        if args.raw and not args.single_sample:
            # Find cellranger mtx folder
            for elem in dirnames:
                if "matrix" not in elem:
                    id = elem
                    merge_dict.update({
                        id: f"{dirpath}/{elem}" \
                            f"/filtered_feature"
                            f"_bc_matrix/"
                    })
                    merge_str += f"{id}-"
        else:
            # Find Anndata files
            for elem in file_names:
                if ".h5ad" in elem and "merged" \
                        not in elem:
                    id = elem.rsplit(".h5ad",
                                     1)[0].split(
                        "_", 1)[0]
                    merge_dict.update(
                        {id: f"{dirpath}/{elem}"})
                    merge_str += f"{id}-"
                if 'merged.h5ad' in elem \
                        and mergename in elem:
                    merged_sample_path = elem

        break

    merge_str = merge_str[:-1]
    if mergename:
        merge_str = mergename

    if args.noclust:
        # No need to merge/clust
        merge_dict = {mergename: screen_folder}
    else:
        #################################################################
        # Merge to new file
        #################################################################
        if merged_sample_path != '':
            print(f"--- Found merged file: "
                  f"{merged_sample_path}")
            merge_str = merged_sample_path
            merged_sample_path = f"{outpath}" \
                                 f"/" \
                                 f"{merged_sample_path}"
        else:
            if args.raw and not args.single_sample:
                print(f"""--- Mode: RAW (from raw 
                    cellranger out)""")
                de_novo_folder = os.getcwd() + \
                                 f"/merged_data/" \
                                 f"{merge_str}.h5ad"

                merged_sample_path = merge_sc(merge_dict,
                                              outputpath=
                                              de_novo_folder,
                                              mode='inner',
                                              scrublet_analysis=
                                              True,
                                              data_path=
                                              os.path.abspath
                                              (all_data_path))
            else:
                print("""--- Mode: PROCESSED (merging .h5ads)
                     """)
                merged_sample_path = merge_sc(merge_dict,
                                              outputpath=
                                              outpath + merge_str,
                                              mode='inner',
                                              scrublet_analysis=
                                              False,
                                              data_path=
                                              os.path.abspath
                                              (all_data_path))

        merge_dict = {merge_str: merged_sample_path}
        # END OF LOOP


#########################################################################
# 2. START ACTUAL ANALYSIS
#########################################################################

def analyse_per_sample(sample2matrix,
                       group_by=None):
    """
    Analyse an AnnData scRNA-seq object.
    :param sample2matrix: name:path-to-h5ad
    :param group_by: condition variable
    :return: tons of plots
    """

    for sample in sample2matrix:
        #################################################################
        # 2.1 PATH SETTINGS
        #################################################################
        # 2.1.1 File, name, check if existing
        sample_h5ad_path = merge_dict[sample]
        if ".h5ad" in sample:
            sample = sample.split(".h5ad", 1)[0]
            if "-merged" in sample:
                sample = sample.split("-merged")[0]
        if not os.path.isfile(sample_h5ad_path):
            print(f"--- File {sample} not found!")
            continue

        # 2.1.2 Directory for results
        if args.harm:
            outputfolder = f"results/{sample}_Harmony"
        else:
            outputfolder = f"results/{sample}"

        outputpath = f"{outpath}{outputfolder}/" \
                     f"{COUNT_LOWER}_" \
                     f"{GENES_BY_COUNT_FRAME[0]}_" \
                     f"{GENES_BY_COUNT_FRAME[1]}" \
                     f"_maxcount_{MITO_CUTOFF}_doubl" \
                     f"{DOUBLET_SCORE}_{MIN_CELLS}/"
        print(f"--- results @ {outputpath}")

        qc_path = f'{outputpath}{sample}_qc.h5ad'
        dr_path = f"{qc_path.rsplit('.h5ad')[0]}" \
                  f"_dr.h5ad"
        clust_path = f"{qc_path.rsplit('.h5ad')[0]}" \
                     f"_clust.h5ad"
        if "clust" in sample:
            clust_path = sample_h5ad_path

        if args.nopp:
            # Will make the qc/dr exist
            qc_path = sample_h5ad_path
            dr_path = sample_h5ad_path

        if args.noclust:
            clust_path = sample_h5ad_path

        harmony_path = f"{qc_path.rsplit('.h5ad')[0]}" \
                       f"_harmony.h5ad"

        ####################################################################
        # 2.2 QC
        ####################################################################
        if "clust" not in sample:
            if not os.path.isfile(qc_path):
                print(f"--- QC for {mergename}.")
                qc_path = easy_qc(sample_h5ad_path,
                                  run_id=sample,
                                  outputpath=
                                  outputpath,
                                  group_by=group_by)
            else:
                print(f"--- Using existing QC file.")

        ####################################################################
        # 2.3 PP, DR, CLUST
        ####################################################################
        if not os.path.isfile(clust_path):
            ################################################################
            # 2.3.1 PP, DR
            ################################################################
            if not os.path.isfile(dr_path):
                print(f"--- PP, DR for {mergename}.")
                # 1. PP
                adata = sc.read_h5ad(qc_path)
                adata = easy_pp(adata)

                # 2. DR
                easy_dr(adata)
                adata.write_h5ad(dr_path)

            if os.path.isfile(dr_path):
                ############################################################
                # Optional: Split commonly PPed data
                ############################################################
                if args.split:
                    print("--- Splitting PP file.")
                    split_adata_by_condition(
                        dr_path,
                        variable=
                        "condition", # "culture",
                        variable_dict=
                        sample2condition)
                    # sample2culture)
                    exit()
                ############################################################
                # Optional: Batch correct on PP data
                ############################################################
                if args.harm:
                    if not os.path.isfile(
                            harmony_path):
                        print(f""" 
                                --- HARMONY for: 
                                    {dr_path}
                                --- SAVING to: 
                                    {harmony_path}
                                """)
                        adata = sc.read_h5ad(dr_path)
                        harmony(adata)
                        adata.write_h5ad(harmony_path)
                        print(f"--- Wrote harmony to"
                              f" {harmony_path}")

                    if os.path.isfile(harmony_path):
                        print(""
                              "--- Loading harmony.")
                        dr_path = harmony_path

                ############################################################
                # 2.3.2 CLUSTERING
                ############################################################
                adata = sc.read_h5ad(dr_path)
                adata = easy_clust(
                    adata,
                    outputpath=outputpath,
                    require_all=True)
                print("--- Clustering complete.")
                paga(adata)
                print("--- PAGA complete.")
                adata.write_h5ad(clust_path)

        if os.path.isfile(clust_path):
            print(f"--- Found clustered file.")
            adata = sc.read_h5ad(clust_path)
            print(adata)

            ################################################################
            # 2.4 Optional: Sub-cluster cell types
            ################################################################
            subcluster = False
            cells_to_reanno = ["Hepatocytes",
                               "Cholangiocytes"]
            targets = [TARGET_HEPATO,
                       TARGET_CHOL]

            if subcluster:
                print(f"--- Subanalysis of "
                      f"{cells_to_reanno}")

                # 2.4.1 Split cell types
                subcells = split_adata_by_condition(
                    clust_path,
                    variable=MAIN_ANNO,
                    clean_up=True,
                    cells=cells_to_reanno,
                    re_anno=False)

                subcells.sort()

                for subcell in subcells:
                    print(f"--- Subcell {subcell}")
                    adata = sc.read_h5ad(subcell)
                    cell_id = subcell.rsplit(
                        'clust_', 1)[1
                    ].rsplit('.h5ad', 1)[0]
                    subout = f"{outputpath}" \
                             f"SUBSETS/{cell_id}/"
                    os.makedirs(subout,
                                exist_ok=True)
                    sc.settings.figdir = \
                        subout

                    if "Hepatocyte" in subcell:
                        target = TARGET_HEPATO
                        res = RES_HEPATO
                        anno_dict = ANNO_DICT_HEPATO

                    if "Cholangiocyte" in subcell:
                        target = TARGET_CHOL
                        res = RES_CHOL
                        anno_dict = ANNO_DICT_CHOL

                    ##########################################################
                    # 2.4.2 Sub-cluster + anno
                    ##########################################################
                    redo = False
                    if target not in list(adata.obs) \
                            or redo:
                        print(f"""
                        --- Sub-clustering 
                        {cell_id}
                        """)
                        re_clust_sub = False

                        if re_clust_sub:
                            adata = easy_clust(
                                adata,
                                outputpath=subout,
                                keys_to_plot=[res],
                                require_all=True)
                            adata.write_h5ad(subcell)

                        adata, keys_to_plot = anno(
                            adata,
                            keys_to_plot=[res],
                            outputpath=subout,
                            nmarkers=NMARKERS,
                            anno_dict=anno_dict)
                        simplify_labels(adata) # Re-name stromal
                        adata.write_h5ad(subcell)
                        print(f"--- Wrote to {subcell}.")

                    #########################################################
                    # 2.4.3 Hepatocyte Zone Scores
                    #########################################################
                    hepato_zonescoring(
                       subcell=subcell,
                       outputpath=outputpath,
                       target=target)

                    #########################################################
                    # 2.4.4 Pathways
                    #########################################################
                    pw_analysis = False
                    if pw_analysis:
                        adata, keys_to_plot = \
                            pathways(
                                adata,
                                cluster_identifier=target,
                                keys_to_plot=[res],
                                outputpath=subout,
                                nmarkers=NMARKERS)
                    ########################################################
                    # 2.4.5 Plot
                    ########################################################
                    replot = True
                    keys = [res, target,
                           target.rsplit("_simple")[0]]
                    if replot:
                        adata = sc.read_h5ad(subcell)
                        integration_plot(adata,
                                         keys_to_plot=keys,
                                         dataset_name="ALL",
                                         outputpath=subout)
                        keep = (adata.obs["controlstatus"]
                                == "CTRL")
                        adata_split = adata[keep, :]
                        integration_plot(adata_split,
                                         keys_to_plot=keys,
                                         dataset_name="CTRL",
                                         outputpath=subout)
                    print("Plotted.")

            ################################################################
            # 2.5 Append metadata
            ################################################################
            if "detailed_condition" not in \
                    list(adata.obs):
                meta_from_dict(adata)
                adata.write_h5ad(clust_path)
                print("--- Wrote meta.")

            ################################################################
            # 2.6 Optional: Clean up
            ################################################################
            clean_up = False
            if clean_up:
                clean(adata, clust_path)

            ################################################################
            # 2.7 Optional: Create composite from subclustered annotations
            ################################################################
            composite = False
            redo = False
            if composite:
                if "composite" not in list(adata.obs) \
                        or redo:
                    adata = create_composite(
                        clust_path,
                        reintegrate=cells_to_reanno,
                        # whom to sub-analyse
                        re_key=targets, # of sub-cell
                        in_key=MAIN_ANNO, # of adata
                        )
                    adata.write_h5ad(clust_path)
            final_cluster_identifier = "composite"

            ###############################################################
            # 2.8 Optional: CellPhoneDB analysis
            ###############################################################
            cpdb = False
            cpdb_cluster_identifier = \
                final_cluster_identifier + "_simple"

            if cpdb:
                splitted_files = \
                    split_adata_by_condition(
                    clust_path,
                    variable="sample")
                #"detailed_condition"

                for subfile in splitted_files:
                    condition_id = subfile.rsplit(
                        ".h5ad", 1)[0].rsplit("_", 1)[1]
                    cpdb_out_dir = \
                        f"{os.getcwd()}/" \
                        f"{subfile.rsplit('/', 1)[0]}" \
                        f"/{condition_id}"
                    print(f"--- CPDB on {condition_id}")

                    if not os.path.isfile(
                            f"{cpdb_out_dir}"
                            f"/count_network.txt"):
                        print(f"Missing "
                              f"{cpdb_out_dir}"
                              f"/count_network.txt")

                        # Initiate job submission
                        script_loc = f"{bash_path}" \
                                     f"/cpdb.sh"
                        key = "metadata__"
                        meta_dir = \
                            f"{subfile.rsplit('/', 1)[0]}" \
                            f"/METADATA/{condition_id}/" \
                            f"METADATA/" \
                            f"{key}" \
                            f"{cpdb_cluster_identifier}" \
                            f".txt"
                        command = os.popen(
                            f"sbatch {script_loc} {subfile} "
                            f"{meta_dir} {cpdb_out_dir}")
                        print(f"sbatch {script_loc} {subfile} "
                              f"{meta_dir} {cpdb_out_dir}")
                        command.close()

            ################################################################
            # 2.9 Annotate cell identity
            ################################################################
            re_anno = False
            if re_anno:
                print(f""" 
                        --- ANNO for {mergename} ---
                        """)

                # Annotate based on control samples
                adata, keys_to_plot = anno(adata,
                                           outputpath=
                                           outputpath,
                                           dataset_name=
                                           "ANNOTATION",
                                           rerank=rerank
                                           )
                adata.write_h5ad(clust_path)
                print("--- Wrote annotated.")
                ############################################################
                # 2.9.1 Smplify labels
                ############################################################
                simplify_labels(adata)
                print(adata.obs["sample"].value_counts())
                adata.write(clust_path)
                print("--- Wrote simplified labels.")

            ################################################################
            # 2.10 Plotting
            ################################################################
            print(f""" 
                    --- Plots/DGE/PATHWAYS/TI for {mergename} ---
                    """)
            # Plot for only controls and full each
            for subset in [False, "CTRL"]:
                if subset == "CTRL":
                    print(f"--- Analyzing {subset} only.")
                    keep = (adata.obs["controlstatus"]
                            == subset)
                    adata_split = adata[keep, :]
                    integration_plot(adata_split,
                                     dataset_name=
                                     subset,
                                     outputpath=
                                     outputpath,
                                     scores=False,
                                     umap=False,
                                     keys_to_plot=[
                                         "leiden_0.12sctype_200",
                                         final_cluster_identifier,
                                         f"{final_cluster_identifier}"
                                         f"_simple"
                                     ])

                else:
                    print(f"--- Ananlyzing ALL conditions.")
                    adata_split = adata
                    plot = True
                    if plot:
                        integration_plot(adata_split,
                                         dataset_name="ALL",
                                         outputpath=outputpath,
                                         keys_to_plot=[
                                             "leiden_0.12sctype_200",
                                             f"{final_cluster_identifier}"
                                             f"_simple",
                                             final_cluster_identifier
                                         ],
                                         umap=False,
                                         scores=True)
                    ########################################################
                    # 2.11 DGE/GO (on ALL samples)
                    ########################################################
                    GO = False
                    if GO:
                        if len(adata.obs[
                                   "condition"
                               ].value_counts().index.tolist()) > 1:
                            for e in [e for e in list(
                                    adata_split.obs)
                                      if final_cluster_identifier in e]:
                                GO_terms(clust_path,
                                         id=e,
                                         outputpath=outputpath,
                                         cluster_identifier=e)

                    ########################################################
                    # 2.12 Trajectory inference
                    ########################################################
                    for lineage in [
                        "Hepatic", "Stellate"
                    ]:
                        plot_imp_norm(adata_path=clust_path,
                                 cluster_identifier=
                                 cpdb_cluster_identifier,
                                 nn=30, ncomp=30,
                                 num_wp=9000,
                                 lineage=lineage)

                        hts_easy(adata_path=clust_path,
                                 cluster_identifier=
                                 cpdb_cluster_identifier,
                                 nn=30, ncomp=30,
                                 num_wp=9000,
                                 lineage=lineage)
            # END OF ANALYSIS
    print(f"___ FINISHED SAMPLE {sample} and "
          f"saved to {outputpath}{sample}")
    # END OF SAMPLE LOOP


############################################################################
# 3. Execute!
############################################################################
analyse_per_sample(merge_dict, group_by="sample")

# END OF SCRIPT
