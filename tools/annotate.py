"""

Module to annotate cell clusters.

@author: Anja Hess
@date: 2023-MAY-01

"""
import scanpy as sc
import pandas as pd
import gseapy
import os
import sys

path_delim = "/"
script_path = str(os.path.dirname(os.path.abspath(__file__)))
maindir = script_path.split(path_delim + "bin")[0]
src_path = maindir + path_delim + "src"
bin_path = maindir + path_delim + "bin"
sys.path.insert(0, script_path)
sys.path.insert(0, maindir)
sys.path.insert(0, src_path)
sys.path.insert(0, bin_path)
sys.path.insert(0, os.path.abspath('.'))
sys.path.insert(0, os.path.abspath('..' + path_delim))
sys.path.insert(0, os.path.abspath('..' + path_delim + 'src'))
from src.sc_constants import annoset2path, RECLUSTER, STATSTEST, \
    keys_to_plot, SEARCHTERMS, NMARKERS, CELLMARKER_PATH
from tools.plotting.plots import enrichment_barplot, integration_plot


def hepato_zonescoring(subcell="", outputpath="",
                       target="",
                       zone_identifier=
                       "leiden_1.02macparlandtop100_200"):
    """

    Marker gene-based (MacParland et al., 2018)
    information on hepatocyte zone identity.

    :param subcell: str
    :param outputpath: str
    :param target: str
    :param zone_identifier: str

    :return:

    """

    ########################################################################
    # 0. Imports
    ########################################################################
    from src.core_functions import easy_clust
    from tools.support_functions import split_adata_by_condition, zone_score

    ########################################################################
    # Zone score @ MacParland if AH-annotated Heps
    ########################################################################
    if "Hepato" in subcell:
        hepatos_to_reanno = ["AH"]
        print(f"--- Subanalysis of hepatos {hepatos_to_reanno}")

        # Split
        subcells = split_adata_by_condition(subcell,
                                            variable=target,
                                            clean_up=True,
                                            cells=hepatos_to_reanno,
                                            re_anno=False)
        for subcell in subcells:
            if "AH" in subcell:
                adata = sc.read_h5ad(subcell)
                cell_id = subcell.rsplit(
                    'clust_', 1)[1].rsplit('.h5ad', 1)[0]

                subout = f"{outputpath}SUBSETS/{cell_id}/"
                os.makedirs(subout, exist_ok=True)
                sc.settings.figdir = subout

                if zone_identifier not in list(adata.obs):
                    adata = easy_clust(adata,
                                       outputpath=subout,
                                       keys_to_plot=["leiden_1.0"],
                                       require_all=False)
                    adata.write_h5ad(subcell)

                    adata, keys_to_plot = anno(adata,
                                               keys_to_plot=[
                                                   "leiden_1.0"],
                                               outputpath=
                                               subout,
                                               nmarkers=NMARKERS,
                                               anno_dict=
                                               {"macparlandtop100":
                                                    f"{CELLMARKER_PATH}"
                                                    f"macparland_top.csv",
                                                })
                    adata.write_h5ad(subcell)
                ############################################################
                # Zone score
                ############################################################
                adata = zone_score(adata,
                                   cluster_identifier=
                                   zone_identifier)
                adata.write_h5ad(subcell)

                print("--- Added zonation score and wrote file.")
                ############################################################
                # Plot individually
                ############################################################
                replot = True
                if replot:
                    integration_plot(adata,
                                     dataset_name="ALL",
                                     keys_to_plot=["Zonation_score"],
                                     outputpath=subout)

                    keep = (adata.obs["controlstatus"]
                            == "CTRL")
                    adata_split = adata[keep, :]
                    integration_plot(adata_split,
                                     dataset_name="CTRL",
                                     keys_to_plot=["Zonation_score"],
                                     outputpath=subout)
    # END OF FUNCTION


def adjust_dict(dictionary, list_valid):
    myDict = {key: [e for e in dictionary[key]
                    if e in list_valid]
              for key in dictionary
              }
    return myDict


def reduce_to_genes_expressed(adata, dictionary=""):
    """

    :param adata: h5ad
    :param dictionary: dict

    :return: geneset_adjusted: dict

    """

    try:
        geneset_adjusted = adjust_dict(
            dictionary, adata.var['gene_ids'])
    except:
        try:
            geneset_adjusted = adjust_dict(
                dictionary, adata.var['gene_ids-1'])
        except:
            try:
                geneset_adjusted = adjust_dict(
                    dictionary, adata.var['gene_ids-0'])
            except:
                try:
                    geneset_adjusted = adjust_dict(
                        dictionary, adata.var_names)
                except:
                    print(adata.var_names)
                    print(f"--- None of the "
                          f"genes expressed: {dictionary}.")
                    return False
    return geneset_adjusted
    # END OF FUNCTION


def pathways(adata, outputpath="", rerank=True,
             nmarkers=NMARKERS,
             keys_to_plot=keys_to_plot,
             dataset_name="ANNOTATION_PW",
             cluster_identifier="",
             save=False):
    """

    Pathway enrichment

    :param adata: h5ad
    :param outputpath: str
    :param rerank: bool
    :param nmarkers: int
    :param keys_to_plot: list
    :param dataset_name: str
    :param cluster_identifier: str

    :return: enrichment plots

    """

    ########################################################################
    # 0. Directories
    ########################################################################
    print(f"... Pathway enrichment for {cluster_identifier}")
    os.makedirs(f"{outputpath}/INTEGR/{dataset_name}/", exist_ok=True)

    ########################################################################
    # Rank genes & export
    ########################################################################
    stats_key = f"{STATSTEST}_{cluster_identifier}"
    if stats_key + "_CTRLS" not in list(adata.uns):
        keep = (adata.obs["controlstatus"] == "CTRL")
        adata_split = adata[keep, :]
        # Use controls for enrichment only
        sc.tl.rank_genes_groups(adata_split,
                                cluster_identifier,
                                method=STATSTEST,
                                key_added=stats_key + "_CTRLS"
                                )
        adata.uns[stats_key + "_CTRLS"] = adata_split.uns[stats_key + "_CTRLS"]
        if save:
            adata.write_h5ad(save)

    genes = pd.DataFrame(adata.uns[stats_key + "_CTRLS"]['names'])
    csv_path = f"{outputpath}/INTEGR/{dataset_name}/" \
               f"gene_ranks_{cluster_identifier}_CTRL_ONYL.csv"
    genes.to_csv(csv_path)

    ########################################################################
    # Enrich
    ########################################################################
    for cl in adata.obs[cluster_identifier] \
            .cat.categories.tolist():
        print(f"--- Pathways for {cl}")
        for term in SEARCHTERMS:
            os.makedirs(outputpath + f"/{term}", exist_ok=True)
            glist = sc.get.rank_genes_groups_df(adata,
                                                pval_cutoff=0.05,
                                                log2fc_min=2,
                                                group=cl,
                                                key=stats_key + "_CTRLS")[
                'names'].tolist()
            glist = [e for e in glist if "MT-" not in e and "RPS"
                     not in e and "RPL" not in e]
            enr_res = gseapy.enrichr(gene_list=glist,
                                     organism='Human',
                                     gene_sets=term,
                                     description='pathway',
                                     cutoff=0.05,
                                     outdir=outputpath + f"/{term}")
            enr_res.results.to_csv(f"{outputpath}/{term}/"
                                   f"{term}_{cl}"
                                   f"_{cl}_PVAL.csv")
            enrichment_barplot(f"{outputpath}/{term}/{term}_{cl}_{cl}_PVAL.csv",
                               outputpath=f"{outputpath}/{term}/",
                               title=f"{cl}_{cl}_PVAL",
                               cutoff_col="P-value")
    # END OF FUNCTION


def anno(adata, outputpath="", rerank=True,
         nmarkers=NMARKERS,
         keys_to_plot=keys_to_plot,
         dataset_name="ANNOTATION",
         anno_dict=annoset2path,
         save=False):
    """

    Function to assign cell type identity.

    :param adata_path:
    :param dataset_name:
    :param keys_to_plot:
    :param outputpath:

    :return: adata.obs slot with ID

    """

    ########################################################################
    # 0. Directories
    ########################################################################
    os.makedirs(f"{outputpath}/INTEGR/", exist_ok=True)
    os.makedirs(f"{outputpath}/INTEGR/{dataset_name}", exist_ok=True)

    ########################################################################
    # 1. Define annotations and start anno
    ########################################################################
    for cluster_identifier in keys_to_plot:
        print(f"--- Anno {cluster_identifier}")

        # Define cluster-specific gene rank key
        stats_key = f"{STATSTEST}_{cluster_identifier}"

        if cluster_identifier not in list(adata.obs):
            print(f"Adding {cluster_identifier}")
            sc.tl.leiden(adata, resolution=
                         float(cluster_identifier.rsplit("_")[1]),
                         key_added=cluster_identifier)

        if stats_key not in adata.uns:
            print("--- Reranking")
            sc.tl.rank_genes_groups(adata, cluster_identifier,
                                    method=STATSTEST,
                                    key_added=stats_key)
        if len(adata.obs[cluster_identifier].value_counts()
                       .index.tolist()) < 2:
            print("--- Too little clusters")
            continue

        ####################################################################
        # Base annotation on control samples
        ####################################################################
        if rerank:
            keep = (adata.obs["controlstatus"] == "CTRL")
            adata_split = adata[keep, :]
            if stats_key + "_CTRLS" not in list(adata.uns):
                len_all = len(adata.obs[
                                  cluster_identifier
                              ].value_counts().index.tolist())
                len_ctrl = len(adata_split.obs[
                                   cluster_identifier
                               ].value_counts().index.tolist())
                print(f"--- {len_ctrl}/{len_all} clusters in CTRL.")
                sc.tl.rank_genes_groups(adata_split,
                                        cluster_identifier,
                                        method=STATSTEST,
                                        key_added=stats_key + "_CTRLS")
                genes = pd.DataFrame(adata_split.uns[
                                         stats_key + "_CTRLS"]['names'])
                csv_path = f"{outputpath}/INTEGR/{dataset_name}/" \
                           f"gene_ranks_{cluster_identifier}_CTRL_ONYL.csv"
                genes.to_csv(csv_path)
                adata.uns[stats_key + "_CTRLS"] = adata_split.uns[stats_key
                                                                  + "_CTRLS"]
                if save:
                    adata.write_h5ad(save)
        else:
            adata_split = adata

        ####################################################################
        # 1. Iterate through marker gene sets
        ####################################################################
        for marker_gene_set in anno_dict:
            n = nmarkers[0]
            anno_prefix = f"{cluster_identifier}" \
                          f"2{marker_gene_set}_{n}"
            plot_path = outputpath + f"/INTEGR/" \
                                     f"{dataset_name}/" \
                                     f"{anno_prefix}/"
            observations = list(adata.obs)  # not to top!
            if f"{cluster_identifier}2{marker_gene_set}" \
                    not in observations or RECLUSTER:
                print(f"""
                        --- Annotation: {anno_prefix}
                        """)
                anno_path = anno_dict[marker_gene_set]
                os.makedirs(plot_path, exist_ok=True)
                sc.settings.figdir = plot_path
                if rerank:
                    adata = celltype_annotation(adata, anno_path,
                                                cluster_identifier=
                                                cluster_identifier,
                                                outputpath=plot_path,
                                                prefix=anno_prefix,
                                                nmarkers=n,
                                                control_set=adata_split,
                                                stats_key=stats_key)
                else:
                    adata = celltype_annotation(adata, anno_path,
                                                cluster_identifier=
                                                cluster_identifier,
                                                outputpath=plot_path,
                                                prefix=anno_prefix,
                                                nmarkers=n,
                                                stats_key=stats_key)

    ########################################################################
    # 2. Add the annotations to the keys to plot
    ########################################################################
    novel_keys = []
    for cluster_identifier in keys_to_plot:
        if "SCN" not in cluster_identifier:
            for marker_gene_set in anno_dict:
                novel_keys.append(f"{cluster_identifier}"
                                  f"2{marker_gene_set}")
    keys_to_plot = sorted(keys_to_plot + novel_keys)
    return adata, keys_to_plot
    # END OF FUNCTION


def parse_annos(path_to_csv):
    """

    Function to read marker gene lists from csvs.

    :param path_to_csv: str
    :return: python dictionary

    """

    gene_dict = {}
    df = pd.read_csv(path_to_csv)

    if "sctype" not in path_to_csv:
        # Exclude mus musculus
        df = df[df['species'] != "Mm"]
        # Parse
        for idx, row in df.iterrows():
            if row['cell type'] not in gene_dict:
                gene_dict.update({
                    row['cell type']:
                        [row['official gene symbol']]})
            else:
                gene_dict[row['cell type']
                ].append(row['official gene symbol'])
    else:
        for idx, row in df.iterrows():
            if row['cell type'].replace("/", "-") not in gene_dict:
                gene_dict.update({
                    row['cell type'].replace("/", "-"):
                        row['official gene symbol'].split(",")})
            else:
                gene_dict[row['cell type'].replace("/", "-")] \
                    += row['official gene symbol'].split(",")

    return gene_dict
    # END OF FUNCTION

def celltype_annotation(adata,
                        path_to_csv,
                        prefix="",
                        cluster_identifier="",
                        outputpath="",
                        nmarkers=0,
                        cutoff=0.05,
                        organism="Human",
                        control_set=False,
                        stats_key="",
                        pathways=False):
    """

    Based on tutorial:
    https://nbisweden.github.io/workshop-scRNAseq/
    labs/compiled/scanpy/scanpy_06_celltype.html

    :param adata: AnnData object
    :param path_to_csv: str
    :param cluster_identifier: str
    :param stats_key: str
    :param outputpath: str

    :return: AnnData object (modified)

    """

    ########################################################################
    # 1. Read marker genes
    ########################################################################
    gene_dict = parse_annos(path_to_csv)

    ########################################################################
    # 2. Check if data present
    ########################################################################
    if control_set:
        if stats_key+"_CTRLS" not in list(control_set.uns) \
                or cluster_identifier not in list(adata.obs):
            print("--- Missing values")
            print(stats_key, list(control_set.uns), )
            print(cluster_identifier, list(adata.obs))
            exit()

    ########################################################################
    # 3. Create empty prediction dictionary
    ########################################################################
    pred = {}

    ########################################################################
    # 4. Iterate trough clusters and annotate
    ########################################################################
    for cl in adata.obs[cluster_identifier] \
            .cat.categories.tolist():
        try:
            glist = sc.get.rank_genes_groups_df(control_set, group=cl,
                                                key=stats_key + "_CTRLS")[
                'names'].tolist()
        except:
            # If a cluster is only present in treatments,
            # get gene ranks from the full dataset
            print(f"--- {cl} not in control.")
            glist = sc.get.rank_genes_groups_df(adata, group=cl,
                                                key=stats_key)[
                'names'].tolist()

        ####################################################################
        # Cell type analysis
        ####################################################################
        enr_res = gseapy.enrichr(gene_list=glist[:nmarkers],
                                 organism=organism,
                                 gene_sets=gene_dict,
                                 background=adata.X.shape[1],
                                 description='pathway',
                                 cutoff=1,outdir=outputpath)

        if enr_res.results.shape[0] == 0:
            # Case of no result
            pred.update({cl: "Unknown"})

        else:
            ##################################################################
            # Results: Sort by adjusted p-value
            ##################################################################
            enr_res.results.sort_values(by="Adjusted P-value",
                                        axis=0, ascending=True,
                                        inplace=True)

            # CAVE ERROR HERE DUE TO INDEX NUMERICAL !
            # CHANGE INDEX TO ANY STR !
            enr_res.results.index = enr_res.results["Term"]
            enr_res.results["Len_overlap"] = enr_res.results[
                "Overlap"].apply(lambda x: int(x.split("/")[0]))

            ##################################################################
            # Select significant enrichments
            ##################################################################
            significant_subset = enr_res.results.loc[
                enr_res.results['Adjusted P-value'] < cutoff]

            if len(significant_subset) == 0:
                # 2.1 The case of no significant results
                enr_res.results.sort_values(by=["P-value", "Len_overlap"],
                                            axis=0,
                                            ascending=[True, False],
                                            inplace=True)

                # Case of similar n genes decide by lowest p val.
                if enr_res.results['Len_overlap'][0] == \
                            enr_res.results['Len_overlap'][1]:
                    enr_res.results = enr_res.results[:2]
                    enr_res.results.sort_values(by=["P-value"],
                                                axis=0,
                                                ascending=[True],
                                                inplace=True)
                pred_clust = f"{enr_res.results['Term'][0]}_ns"
                genes = enr_res.results['Genes'][0]

                # For ESCs require at least one pluripotency gene
                if enr_res.results['Term'][0] == "Embryonic stem cells":
                    if "SOX2" not in genes and "NANOG" not in genes \
                            and "KLF4" not in genes and "POU5F1" not in genes:
                        try:
                            pred_clust = f"Premature " \
                                         f"{enr_res.results['Term'][1]}" \
                                         f"-{enr_res.results['Term'][2]}" \
                                         f"_ns"
                        except:
                            pred_clust = "Unknown"

            else:
                # 2.2 Rank significant results and rank by
                # number of markers
                # If all with adjutsed p value below 0.05,
                # sort by n overlapping genes
                significant_subset = significant_subset.sort_values(
                    by=["Len_overlap"], axis=0, ascending=[False])
                pred_clust = significant_subset['Term'][0]

                genes = significant_subset['Genes'][0]
                # For embryonic stem cells we require at
                # least one of the following: NANOG,
                if pred_clust == "Embryonic stem cells":
                    if "SOX2" not in genes  and "NANOG" not in genes \
                            and "KLF4" not in genes and "POU5F1" not in genes:
                        try:
                            pred_clust = f"Premature " \
                                         f"{enr_res.results['Term'][1]}"
                        except:
                            pred_clust = "Unknown"

            pred_cluster = f"{pred_clust}_{cl}"
            print(f"--- Prediction: {pred_cluster}")
            ##################################################################
            # Export and update
            ##################################################################
            pred.update({cl: pred_cluster})
            enr_res.results.to_csv(f"{outputpath}{pred_clust}_{cl}.csv")

        ######################################################################
        # 2. Pathway analysis (functional)
        ######################################################################
        if pathways:
            for term in SEARCHTERMS:
                os.makedirs(outputpath + f"/{term}", exist_ok=True)
                glist = sc.get.rank_genes_groups_df(control_set,
                                                    pval_cutoff=0.05,
                                                    log2fc_min=1,
                                                    group=cl,
                                                    key=stats_key
                                                        + "_CTRLS")[
                    'names'].tolist()
                glist = [e for e in glist if "MT-" not in e and
                         "RPS" not in e and "RPL" not in e]
                enr_res = gseapy.enrichr(gene_list=glist,
                                         organism=organism,
                                         gene_sets=term,
                                         description='pathway',
                                         cutoff=cutoff,
                                         outdir=outputpath
                                                + f"/{term}")
                print(enr_res.results["Term"][:10])
                enr_res.results.to_csv(f"{outputpath}/{term}/{term}_{pred_clust}"
                                       f"_{cl}_PVAL.csv")
                if len(enr_res.results["Term"]) == 0:
                    print("--- To little hits.")
                else:
                    enrichment_barplot(f"{outputpath}/{term}/{term}_{pred_clust}"
                                       f"_{cl}_PVAL.csv",
                                       outputpath=f"{outputpath}/{term}/",
                                       title=f"{pred_clust}_{cl}_PVAL",
                                       cutoff_col="P-value")

    ########################################################################
    # 3. FINAL: Parse annotation
    ########################################################################
    print("--- Annotating prediction")
    prediction = [pred[x] for x in adata.obs[cluster_identifier]]
    adata.obs[prefix] = prediction
    adata.obs[prefix].astype('category')
    print(adata.obs[prefix].value_counts())
    df_numbers = adata.obs[prefix].value_counts()

    try:
        sc.pl.umap(adata, color=["condition", "sample", prefix],
                   save=f"_{prefix}")
        sc.pl.draw_graph(adata, color=["condition", "sample", prefix],
                         save=f"_{prefix}")
    except:
        print("--- Projections not available")

    ########################################################################
    # Export the to differential genes to csv
    ########################################################################
    genes = pd.DataFrame(adata.uns[stats_key]['names'])
    genes.rename(columns=pred, inplace=True)
    csv_path_new = f"{outputpath}gene_ranks{prefix}.csv"
    genes.to_csv(csv_path_new)
    df_numbers.to_csv(f"{outputpath}cell_numbers_{prefix}.csv")
    return adata
    # END OF FUNCTION

# END OF SCRIPT
