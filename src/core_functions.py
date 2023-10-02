"""

Core functions for scRNA-seq analysis.

@author: Anja Hess
@date:  2022-11-19

"""
import os
import palantir
import pickle
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import gseapy
import scanpy as sc
import scanpy.external as sce
import pandas as pd
import seaborn as sns
from sklearn.preprocessing import StandardScaler

from sc_constants import FIBROSIS_MARKERS, INFLAMMATION_MARKERS, \
    GROUP_COLS, STATSTEST, SEARCHTERMS, palettediff2, geneset2path, \
    cmap, keys_to_plot, sample2detailedcondition, diff2controlcondition, \
    cell_cycle_file, geneset2path_fibrosis, geneset2path_inflammation
from tools.lazy_functions import remove_duplicates, adjust_list
from data_parsing import find_cell

sns.set_style("whitegrid", {'axes.grid': False})


def split_celltype_by_condition(adata, old_celltype="",
                                condition_column=""):
    """

    :param adata: h5ad object
    :param old_celltype: str, column in adata.obs
    with cell type labels to be splitted
    :param condition_column: str, column in
    adata.obs which will split the cells
    :return:

    """

    adata.obs[f"{old_celltype}_by_cond"] = adata.obs
    [old_celltype].replace(np.NaN, 'Is Null value')
    # END OF FUNCTION


def score_geneset(adata_path, save_file="", on_raw=True, per_cluster=True,
                  cluster_identifier="", ext_csv="", pairwise=False,
                  interest=False):
    """

    Score Gene lists.

    :param adata_path:
    :param save_file:
    :param on_raw:
    :return:

    """

    #####################################################################
    # 0. Load data
    #####################################################################
    if type(adata_path) == str:
        adata = sc.read_h5ad(adata_path)
    else:
        adata = adata_path

    #####################################################################
    # 1. Iterate through gene sets
    #####################################################################
    for i, geneset in enumerate(geneset2path):
        if interest:
            if geneset not in interest:
                continue
        # Read the custom genelist from config
        try:
            df = pd.read_csv(geneset2path[geneset],
                             header=None)
        except:
            df = pd.read_table(geneset2path[geneset],
                               header=None)

        # Generate a boolean variable in the df
        geneset_list = df.iloc[:, 0].tolist()
        if "gene_ids" not in list(adata.var):
            try:
                adata.var["gene_ids"] = adata.var.index
            except:
                print("failed to append gene_ids")
        geneset_list = remove_duplicates(geneset_list)

        # Generate the gene list object
        print(f"--- Loading {geneset} with"
              f"{len(geneset_list)} genes.")

        # Reduce to genes in adata.var
        geneset_list = adjust_list(geneset_list, adata.var[
            'gene_ids'])
        print(f"--- Reduced to {len(geneset_list)}"
              f"genes expressed.")

        #################################################################
        # 1.1 Control for bias: calculate & exclude top 5% genes
        #################################################################
        if "total_counts" not in adata.var:
            sc.pp.calculate_qc_metrics(adata, inplace=True)
        targets_total_counts = adata[
                               :, geneset_list
                               ].var["total_counts"].sort_values(
            ascending=False)
        minus_top5 = [e for e in geneset_list if e not in list(
            targets_total_counts.index[: int(5 *
                                             (len(geneset_list)
                                              / 100))])]
        scenarios = {"all": geneset_list,
                     "minus_top5": minus_top5}

        #################################################################
        # 1.2 Score
        #################################################################
        for scenario in scenarios:
            print(f"--- Scoring {scenario}, pairwise: {pairwise}")
            gene_list = scenarios[scenario]

            #############################################################
            # 1.2.1 Separate scoring for control vs. treated
            #############################################################
            if pairwise:
                for dc in [e for e in adata.obs[
                    "detailed_condition"].value_counts().index.tolist()
                           if "CTRL" not in e]:

                    save_plots = save_file.rsplit("/", 1)[0] \
                                 + "/scores/"
                    os.makedirs(save_plots, exist_ok=True)
                    print(f"--- saving to {save_plots}")
                    sc.settings.figdir = save_plots

                    ref = diff2controlcondition[dc]
                    print(f"--- Scoring {dc} & {ref}.")
                    target = [dc, ref]
                    keep = (adata.obs[
                                "detailed_condition"].isin(target))
                    adata_split = adata[keep, :]
                    print(adata_split.obs["condition"].value_counts())

                    sc.tl.score_genes(adata_split,
                                      gene_list,
                                      ctrl_size=len(gene_list),
                                      score_name=
                                      f'sc_score_'
                                      f'{geneset}_'
                                      f'{scenario}_'
                                      f'{on_raw}_pairwise',
                                      copy=False,
                                      use_raw=on_raw)

                    # Get max score
                    max_val = max(adata_split.obs[f'sc_score_'
                                                  f'{geneset}'
                                                  f'_{scenario}'])
                    for cond in target:
                        keep = (adata_split.obs[
                                    "detailed_condition"] == cond)
                        adata_subsplit = adata_split[keep, :]
                        sc.pl.draw_graph(adata_subsplit,
                                         color=f'sc_score_'
                                               f'{geneset}_{scenario}',
                                         save=f'{cond}_sc_score_'
                                              f'{geneset}_'
                                              f'{scenario}'
                                              f'_{cond}',
                                         cmap="Reds",
                                         vmax=max_val,
                                         vmin=0,
                                         palette=palettediff2,
                                         )
                    print(f"--- Scored {geneset}.")

            #################################################################
            # 1.2.2 Default (all)
            #################################################################
            else:
                sc.tl.score_genes(adata,
                                  gene_list,
                                  ctrl_size=len(gene_list),
                                  score_name=
                                  f'sc_score_'
                                  f'{geneset}_{scenario}_{on_raw}',
                                  copy=False,
                                  use_raw=on_raw)
                print(f"--- Scored {geneset}")

    print(f"--- Scoring complete.")
    return adata
    # END OF FUNCTION


def cell_cycle_scoring(adata, outputpath=""):
    """

    Cell cycle scoring following a tutorial:
    https://nbisweden.github.io/workshop-scRNAseq/
    labs/compiled/scanpy/scanpy_01_qc.html

    :param adata:
    :param outputpath:
    :return:
    """
    cell_cycle_genes = [x.strip() for x in
                        open(cell_cycle_file)]
    cell_cycle_genes = [x for x in cell_cycle_genes
                        if x in adata.var_names]

    # Split into 2 lists
    s_genes = cell_cycle_genes[:43]
    g2m_genes = cell_cycle_genes[43:]
    sc.tl.score_genes_cell_cycle(adata,
                                 s_genes=s_genes,
                                 g2m_genes=g2m_genes)
    return adata
    # END OF FUNCTION


def rank_genes(adata, cluster_identifier="", csv_path="",
               rank_key=""):
    """

    DGE per cluster etc.

    :param adata: AnnData object
    :param cluster_identifier: str
    :param csv_path: str
    :param rank_key: str


    """

    # 1. Define cells with > 1 count
    cell_counts = adata.obs[
        cluster_identifier].value_counts()
    cells = adata.obs[
        cluster_identifier].value_counts().index.tolist()
    valid_groups = [
        cell for cell in cells if cell_counts[
                                      cell] > 1]

    # 2. Rank genes
    sc.tl.rank_genes_groups(adata,
                            cluster_identifier,
                            method=STATSTEST,
                            key_added=rank_key,
                            use_raw=True,
                            groups=valid_groups)

    # 3. Export the differentially expressed genes
    genes = pd.DataFrame(adata.uns[
                             rank_key]['names'])
    print("--- Exported to ", csv_path)
    genes.to_csv(csv_path)

    return adata
    # END OF FUNCTION


def easy_pp(adata, n_genes=5000):
    """

    Preprocessing.

    :param adata: AnnData object
    :param n_genes: int
    :return: AnnData with cell cycle scores etc.

    """

    # PP
    print("--- PP")
    print("   ", adata.X.shape)
    adata.raw = adata
    sc.pp.normalize_total(adata,
                          exclude_highly_expressed=True,
                          max_fraction=0.1)
    sc.pp.log1p(adata)
    print("--- Cell cycle scoring")
    adata = cell_cycle_scoring(adata)

    print("--- HVG")
    sc.pp.highly_variable_genes(adata,
                                n_top_genes=n_genes)
    adata = adata[:, adata.var.highly_variable]
    print("   ", adata.X.shape)

    print("--- Scaling")
    adata.layers['scaled'] = sc.pp.scale(adata,
                                         copy=True,
                                         max_value=5).X
    return adata
    # END OF FUNCTION


def harmony(adata):
    """

    Batch correction

    :param adata:
    :return: harmony projection

    """

    print("--- PCA")
    sc.tl.pca(adata)
    sce.pp.harmony_integrate(adata,
                             "sample")

    print("--- Calculating UMAP with harmony rep")
    sc.pp.neighbors(adata,
                    n_pcs=50,
                    use_rep="X_pca_harmony")
    sc.tl.umap(adata)
    # END OF FUNCTION


def easy_dr(adata):
    """

    :param adata: AnnData
    :return: AnnData with PCA, neighbors, UMAP.

    """
    print("--- PCA")
    sc.tl.pca(adata)

    print("--- Neighbors")
    sc.pp.neighbors(adata,
                    n_pcs=50)
    print("--- UMAP")
    sc.tl.umap(adata)
    # END OF FUNCTION


def paga(adata):
    """

    Initially termed paga, for graph generation
    :param adata: AnnData object
    :return:

    """

    sc.tl.draw_graph(adata)
    # END OF FUNCTION


def easy_clust(adata, outputpath="", require_all=True,
               keys_to_plot=keys_to_plot):
    """

    :param adata: AnnData object
    :param outputpath: str
    :param require_all: bool (if cluster needs to contain cells
    from each replicate)

    :return: AnnData with leiden clusters

    """

    print("--- Clustering.")
    observations = list(adata.obs)

    for cluster_identifier in keys_to_plot:
        # Set path
        plot_path = f"{outputpath}/INTEGR/" \
                    f"{cluster_identifier}"
        os.makedirs(plot_path,
                    exist_ok=True)
        csv_path = f"{plot_path}/gene_ranks_" \
                   f"{cluster_identifier}.csv"

        # Double check for ranking, add if missing
        if cluster_identifier not in observations:
            print(f"--- Missing {cluster_identifier}")

            # Get leiden resolution to calculate
            alg = cluster_identifier.rsplit("leiden_", 1)[1]
            alg_num = float(alg)
            rank_key = f"{STATSTEST}_{cluster_identifier}"

            # Calculate Leiden clusters
            if alg not in observations:
                print(f"--- Clustering {alg_num}.")
                if "neighbors" not in list(adata.uns):
                    easy_dr(adata)
                sc.tl.leiden(adata,
                             resolution=alg_num,
                             key_added=cluster_identifier)

                # Check if present in all replicates:
                if require_all:
                    n_samples = len(adata.obs[
                                        "sample"
                                    ].value_counts(
                    ).index.tolist())

                    for cl in adata.obs[
                        cluster_identifier
                    ].value_counts().index.tolist():
                        slice = adata[
                                (adata.obs[cluster_identifier]
                                 == cl), :]
                        n_samples_cl = len(
                            slice.obs["sample"
                            ].value_counts().index.tolist())

                        if n_samples_cl != n_samples:
                            print(f"--- Warning: {cl} "
                                  f"not in all samples.")
                            print(
                                adata.obs[
                                    cluster_identifier
                                ].value_counts())
                            keep = (adata.obs[
                                        cluster_identifier]
                                    != cl)
                            adata = adata[keep, :]

                # Rank genes (Cluster-specific marker genes)
                if rank_key not in list(adata.uns):
                    print(f"--- Ranking {alg}, adding {rank_key}.")
                    rank_genes(adata,
                               cluster_identifier=cluster_identifier,
                               csv_path=csv_path,
                               rank_key=rank_key)

                    # Export DGE table
                    genes = pd.DataFrame(
                        adata.uns[rank_key]['names'])
                    genes.to_csv(csv_path)
                    print("--- Exported gene ranks.")
    return adata
    # END OF FUNCTION


def subset_and_mitoclear(adata, celltype="",
                         cluster_identifier=""):
    """

    Subset for cluster-wise DGE.

    :param adata: AnnData object
    :param celltype: str, which cell to keep
    :param cluster_identifier: key to find cell above
    :return: subset of AnnData object

    """

    # 1. Subset the cluster
    cl1 = adata[adata.obs[cluster_identifier]
                == celltype, :]
    # 2. Convert to raw
    cl1 = cl1.raw.to_adata()

    # 3. Remove mito and ribo genes
    mito_genes = cl1.var_names.str.startswith(
        "MT-")
    ribo_genes = cl1.var_names.str.startswith(
        ("RPS", "RPL"))
    remove = np.add(mito_genes,
                    ribo_genes)
    keep = np.invert(remove)
    cl1 = cl1[:, keep]

    return cl1
    # END OF FUNCTION


def export_dge(adata, key="", outfile="", append=False):
    """

    Export ranked gene tables.

    :param adata: AnnData object
    :param key: str
    :param outfile: str
    :param append: bool (whether to add to existing file)
    :return: creates csv files.

    """

    common_results = []

    for rel in ["names",
                "pvals_adj",
                "logfoldchanges",
                "scores"]:

        df = pd.DataFrame(
            adata.uns[key][rel])
        cols = list(df.columns)
        new_cols = [e + f"_{rel}" for e in cols]
        df.columns = new_cols

        if type(common_results) == list:
            common_results = df
        else:
            common_results = pd.concat(
                [common_results, df],
                axis=1)

    if os.path.isfile(outfile) and append:
        data_new = pd.read_csv(outfile,
                               index_col=0)
        big_df = pd.concat([
            data_new,
            common_results.reset_index(drop=True)],
            axis=1)
        big_df.to_csv(outfile)
    else:
        common_results.to_csv(outfile)

    print(f"--- Exported DGE to {outfile}")
    # END OF FUNCTION


def dge(adata, variable="", rank_key="", outdir=""):
    """

    Differential gene expression.

    :param adata: AnnData object, 1 celltype
    :param variable: str, in adata.obs defining
    the conditions
    :param rank_key: how to name the wilcoxon
    result section
    :return: DGE directory with genes.

    """

    ############################################################################
    # 1. Find conditions in cell type
    ############################################################################

    conditions_here = adata.obs[
        "condition"
    ].value_counts().index.tolist()
    if len(conditions_here) < 2:
        diffs_to_check = [conditions_here[0]]
    else:
        diffs_to_check = conditions_here

    ############################################################################
    # 2. Iterate through them to compare pairwise to their controls
    ############################################################################
    for dc in [e for e in diffs_to_check
               if "CTRL" not in e]:
        try:
            dc_key = f"{dc}_vs_{diff2controlcondition[dc]}"
        except:
            print(f"{dc} not in dictionary {conditions_here}.")
        ref = diff2controlcondition[dc]

        print(f"--- Ranking {dc} vs. {diff2controlcondition[dc]}.")
        sc.tl.rank_genes_groups(adata, variable,
                                groups=[dc, ref],
                                method=STATSTEST,
                                reference=ref,
                                key_added=dc_key)

        values_df = sc.pl._tools._get_values_to_plot(
            adata,
            "scores",
            INFLAMMATION_MARKERS + FIBROSIS_MARKERS,
            key=dc_key,
        )

        var_names = {"Inflammation": INFLAMMATION_MARKERS,
                     "Fibrosis": FIBROSIS_MARKERS}
        e = "score"

        sc.pl.matrixplot(
            adata,
            var_names,
            groupby=variable,
            values_df=values_df,
            colorbar_title=e,
            vmin=-10,
            vmax=10,
            save=f"DGE_{dc_key}_NORM.pdf",
            cmap=cmap,
        )

        export_dge(adata, key=dc_key,
                   outfile=f"{outdir}pairwise.csv",
                   append=True)

        sc.pl.matrixplot(
            adata,
            var_names,
            groupby=variable,
            values_df=values_df,
            colorbar_title=e,
            vmin=-20,
            vmax=20,
            save=f"DGE_{dc_key}_NORM20.pdf",
            cmap=cmap,
        )

        sc.pl.matrixplot(
            adata,
            var_names,
            groupby=variable,
            values_df=values_df,
            colorbar_title=e,
            save=f"DGE_{dc_key}.pdf",
            cmap=cmap,
        )
        sc.pl.matrixplot(
            adata,
            var_names,
            groupby=variable,
            values_df=values_df,
            save=f"DGE_{dc_key}_scaled.pdf",
            vmin=-(np.amax(values_df.values)),
            vmax=np.amax(values_df.values),
            cmap=cmap,
        )

    print("--- Ranked all conditions pairwise.")

    # 3. Also do a general ranking test
    print("--- Comparing vs rest.")
    try:
        sc.tl.rank_genes_groups(adata,
                                variable,
                                method=STATSTEST,
                                key_added=rank_key)
        export_dge(adata,
                   key=rank_key,
                   outfile=f"{outdir}vs_rest.csv")
    except:
        print("Too little cells")

    return adata
    # END OF FUNCTION


def GO_terms(adata_path, cluster_identifier="", search_terms=SEARCHTERMS,
             id="", pval_enr="TOP200", cutoff=0.05, outputpath="",
             pairwise=True, variable="sample2detailedcondition",
             overwrite=False):
    """

    Pathway enrichment for DGE

    :param adata: AnnData object
    :param cluster_identifier: str
    :param search_terms: list
    :param id: str
    :param cutoff: float
    :param outputpath: str
    :return:

    """

    # Import to avoid circular
    from tools.plotting.plots import circular_barplot

    ############################################################################
    # 1. Load data, map conditions, define cell types
    ############################################################################
    adata = sc.read_h5ad(adata_path)
    if pairwise:
        id = f"{id}_pairwise"

    ############################################################################
    # 2. Find the cluster identifier to be used
    ############################################################################
    pot_ids = list(adata.obs)
    for e in pot_ids:
        if cluster_identifier in e and "counts" not in e:
            cluster_identifier = e
            break

    rank_key = f"{STATSTEST}_{cluster_identifier}"
    print("--- Found cluster identifier:",
          cluster_identifier)

    # Make output directories
    os.makedirs(f"{outputpath}/ENRICHMENT_{pval_enr}/",
                exist_ok=True)
    os.makedirs(f"{outputpath}/ENRICHMENT_FILES/",
                exist_ok=True)
    root = f"{outputpath}/ENRICHMENT_{pval_enr}/{id}/"
    os.makedirs(root,
                exist_ok=True)
    sc.settings.figdir = root

    ############################################################################
    # 3. Map condition
    ############################################################################
    if variable not in list(adata.obs):
        adata.obs[variable] = adata.obs["sample"].map(
            sample2detailedcondition).astype('category')
    cells = adata.obs[
        cluster_identifier] \
        .value_counts().index.tolist()

    ############################################################################
    # 4. Test each cell cluster separately
    ############################################################################
    for celltype in cells:
        # Define paths
        celltype_id = celltype.replace(" ", "-")
        outpath = f"{root}{celltype_id}/"
        celltype_file = f"{outputpath}/ENRICHMENT_FILES/" \
                        f"{celltype_id}.h5ad"
        os.makedirs(outpath,
                    exist_ok=True)
        os.makedirs(f"{outputpath}/ENRICHMENT_FILES/",
                    exist_ok=True)
        sc.settings.figdir = outpath
        print(f"--- DGE for cluster {celltype}")

        if not os.path.isfile(celltype_file) or overwrite:
            # 4.1 Slice the cell type
            print(f"--- Preparing cell type {celltype}")
            cl1 = subset_and_mitoclear(adata,
                                       cluster_identifier=
                                       cluster_identifier,
                                       celltype=celltype)
            ####################################################################
            # 4.2 Rank genes (general)
            ####################################################################
            cl1 = dge(cl1,
                      variable=variable,
                      rank_key=rank_key,
                      outdir=outpath)
            cl1.write_h5ad(celltype_file)
            print("--- Wrote cell type file.")

        if os.path.isfile(celltype_file):
            print("--- Found file.")
            cl1 = sc.read_h5ad(celltype_file)
            ####################################################################
            # 4.3 Rank genes (pairwise comparison)
            ####################################################################
            for comparison in [e for e in list(cl1.uns)
                               if "vs" in e or e == rank_key]:
                if rank_key not in comparison:
                    comparison_path = f"{outpath}{comparison}"
                    os.makedirs(comparison_path,
                                exist_ok=True)
                    sc.settings.figdir = comparison_path
                    # Export DGE
                    if not os.path.isfile(f"{comparison_path}"
                                          f"pairwise-full.csv"):
                        export_dge(cl1,
                                   key=comparison,
                                   outfile=f"{comparison_path}"
                                           f"pairwise-full.csv",
                                   append=True)

            ####################################################################
            # 4.4 Enrich the found DEG (pathway enrichment)
            ####################################################################
            for term in search_terms:
                for bl in ["induced", "reduced"]:
                    csv_res_file = f"{outpath}/enr_res_" \
                                   f"{celltype_id}_{term}_" \
                                   f"{comparison}_{bl}.csv"
                    redo = False

                    # Split by pos/neg l2fc
                    if bl == "induced":
                        # only genes with pos l2fc
                        minlog2fc = 0
                        maxlog2fc = 10000
                    else:
                        # only genes with neg l2fc
                        minlog2fc = -10000
                        maxlog2fc = 0

                    # Select top genes of the up/down list
                    cut = False
                    if type(pval_enr) == str and "TOP" in pval_enr:
                        len_glist = int(pval_enr.rsplit('TOP')[1])
                        cut = True

                    if not os.path.isfile(csv_res_file) or redo:
                        print(f"--- Querying {term} ({comparison}) {bl}")
                        glist = sc.get.rank_genes_groups_df(
                            cl1,
                            group=dc,
                            log2fc_max=maxlog2fc,
                            log2fc_min=minlog2fc,
                            key=comparison)[
                            'names']
                        try:
                            # If there are byte formatting issues
                            glist = [e.decode('UTF-8') for e in glist]
                        except:
                            glist = sc.get.rank_genes_groups_df(
                                cl1,
                                group=dc,
                                log2fc_max=maxlog2fc,
                                log2fc_min=minlog2fc,
                                key=comparison)[
                                'names'].squeeze().str.strip().tolist()

                        if cut:
                            glist = glist[:len_glist]

                        if len(glist) == 0:
                            print(f"--- No genes for {bl}")
                            continue
                        else:
                            print(f"---- Enriching {len(glist)} genes")

                        # Enrich
                        enr_res = gseapy.enrichr(gene_list=glist,
                                                 organism='Human',
                                                 gene_sets=term,
                                                 description='pathway',
                                                 cutoff=cutoff,
                                                 outdir=comparison_path)

                        ######################################################
                        # Append to diff dataframe
                        ######################################################
                        if bl == "induced":
                            enr_res.res2d['UP_DW'] = "UP"
                        else:
                            enr_res.res2d['UP_DW'] = "DOWN"

                        # 3.3 Export the information to csv
                        enr_res.results.to_csv(csv_res_file)

                    if os.path.isfile(csv_res_file):
                        enr_res = pd.read_csv(csv_res_file,
                                              index_col=0)
                        if bl == "induced":
                            enr_res['UP_DW'] = "UP"
                        else:
                            enr_res['UP_DW'] = "DOWN"

            ##################################################################
            # 4.5 Meta-plot for each pathway (all conditions)
            ##################################################################
            for term in search_terms:
                for bl in ["induced", "reduced"]:
                    for sorter in ["Combined Score"]:
                        # Initiate df for each sorter
                        dc_go_frame = []
                        ######################################################
                        # For each comparison plot the sign. pathways (p<0.05)
                        ######################################################
                        for comparison in [e for e in list(cl1.uns)
                                           if "vs" in e]:
                            dc = comparison.rsplit("_vs_", 1)[0]
                            res_file = f"{outpath}/enr_res_" \
                                       f"{celltype_id}_{term}" \
                                       f"_{comparison}_{bl}.csv"
                            print(f"--- Found {comparison}")

                            # Filter significant results!
                            try:
                                df = pd.read_csv(res_file,
                                                 index_col=0)
                            except:
                                print("--- No res.")
                                continue
                            df = df.loc[df["Adjusted P-value"]
                                        < cutoff]
                            new_df = pd.DataFrame()
                            new_df["Term"] = df["Term"]
                            new_df["Adjusted P-value"] = df["Adjusted P-value"]
                            new_df["Condition"] = dc
                            new_df["Combined Score"] = df["Combined Score"]
                            ascending = False

                            # Sort the sign. results by score!
                            new_df.sort_values(by=sorter,
                                               ascending=ascending,
                                               inplace=True)

                            if type(dc_go_frame) == list:
                                dc_go_frame = new_df
                            else:
                                dc_go_frame = dc_go_frame.append(
                                    new_df,
                                    ignore_index=True)

                        if type(dc_go_frame) != list:
                            dc_go_frame.to_csv(f"{outpath}sourcefile_"
                                               f"{term}_{sorter}_{bl}.csv")

                            # Circular barplot
                            circular_barplot(dc_go_frame,
                                             outfolder=outpath,
                                             id=f"{term}_{sorter}_{bl}",
                                             title=f"{term} {celltype_id} {bl}",
                                             value="Adjusted P-value",
                                             name="Term",
                                             group="Condition",
                                             top=15)
                            print("--- Circular barpolot saved.")
        del cl1
        # END OF CELL CLUSTER LOOP
    # END OF FUNCTION


def plot_gene_trends(gene_trends, genes=None, length=6):
    """

    For palantir trajectory inference.

    :param: gene_trends: Results
    :param genes: list
    :param length: plotting n genes

    """

    ############################################################################
    # 1. Get branches and corresponding genes
    ############################################################################
    branches = list(gene_trends.keys())
    colors = pd.Series(palettediff2[:length], index=branches)
    if genes is None:
        genes = gene_trends[branches[0]]["trends"].index

    ############################################################################
    # 2. Set up figure
    ############################################################################
    fig = plt.figure(figsize=[7, 3 * len(genes)])
    for i, gene in enumerate(genes):
        ax = fig.add_subplot(len(genes), 1, i + 1)
        for branch in branches:
            trends = gene_trends[branch]["trends"]
            stds = gene_trends[branch]["std"]
            ax.plot(
                trends.columns, trends.loc[gene, :],
                color=colors[branch],
                label=branch
            )
            ax.set_xticks([0, 1])
            ax.fill_between(
                trends.columns,
                trends.loc[gene, :] - stds.loc[gene, :],
                trends.loc[gene, :] + stds.loc[gene, :],
                alpha=0.1,
                color=colors[branch],
            )
            ax.set_title(gene)
        # Add legend
        if i == 0:
            ax.legend()
    sns.despine()
    # END OF FUNCTION


def plot_pal_res(pr_res, tsne, s=2):
    """

    *** Taken from original palantir module to reduce pdf
    size *** (Otherwise to heavy for local PC to open)

    """
    ############################################################################
    # 1. Set up figure
    ############################################################################
    n_branches = \
        pr_res.branch_probs.shape[1]
    n_cols = 6
    n_rows = int(np.ceil(n_branches / n_cols))
    fig = plt.figure(figsize=
                     [2 * n_cols, 2 * (n_rows + 2)])
    gs = plt.GridSpec(
        n_rows + 2,
        n_cols,
        height_ratios=
        np.append([0.75, 0.75],
                  np.repeat(1, n_rows))
    )
    cmap = matplotlib.cm.plasma

    ############################################################################
    # 2. Pseudotime plot
    ############################################################################
    ax = plt.subplot(gs[0:2, 1:3])
    c = pr_res.pseudotime[tsne.index]
    ax.scatter(tsne.iloc[:, 0],
               tsne.iloc[:, 1],
               s=s,
               cmap=matplotlib.cm.plasma,
               c=c,
               rasterized=True)
    normalize = matplotlib.colors.Normalize(
        vmin=np.min(c),
        vmax=np.max(c))
    cax, _ = matplotlib.colorbar.make_axes(ax)
    cbar = matplotlib.colorbar.ColorbarBase(cax,
                                            norm=normalize,
                                            cmap=cmap)
    ax.set_axis_off()
    ax.set_title("Pseudotime")

    ############################################################################
    # 3.  Entropy plot
    ############################################################################
    ax = plt.subplot(gs[0:2, 3:5])
    c = pr_res.entropy[tsne.index]
    ax.scatter(tsne.iloc[:, 0],
               tsne.iloc[:, 1],
               s=s,
               cmap=matplotlib.cm.plasma,
               c=c,
               rasterized=True)
    normalize = matplotlib.colors.Normalize(
        vmin=np.min(c),
        vmax=np.max(c))
    cax, _ = matplotlib.colorbar.make_axes(ax)
    cbar = matplotlib.colorbar.ColorbarBase(cax,
                                            norm=normalize,
                                            cmap=cmap)
    ax.set_axis_off()
    ax.set_title("Differentiation potential")

    ############################################################################
    # 4. Terminal branches
    ############################################################################
    for i, branch in enumerate(pr_res.branch_probs.columns):
        row = int(np.floor(i / n_cols))
        ax = plt.subplot(gs[row + 2, np.remainder(i, n_cols)])
        c = pr_res.branch_probs.loc[tsne.index, branch]
        ax.scatter(
            tsne.iloc[:, 0],
            tsne.iloc[:, 1],
            s=s,
            cmap=matplotlib.cm.plasma,
            c=c,
            rasterized=True
        )
        normalize = matplotlib.colors.Normalize(
            vmin=np.min(c),
            vmax=np.max(c))
        cax, _ = matplotlib.colorbar.make_axes(ax)
        cbar = matplotlib.colorbar.ColorbarBase(cax,
                                                norm=normalize,
                                                cmap=cmap)
        ax.set_axis_off()
        ax.set_title(branch, fontsize=10)
    # END OF FUNCTION EDITED FROM PALANTIR MODULE


def get_fig(fig=None, ax=None, figsize=[4, 4]):
    """

    fills in any missing axis or figure with the currently active one

    :param ax: matplotlib Axis object
    :param fig: matplotlib Figure object

    """

    if not fig:
        fig = plt.figure(figsize=figsize)
    if not ax:
        ax = plt.gca()
    return fig, ax
    # END OF FUNCTION

def highlight_cells(plot_tsne, cells, fig=None, ax=None):
    """

    Taken from original palantir to reduce pdf size
    Function to highlight specific cells on the tSNE map

    """
    fig, ax = get_fig(fig=fig,
                      ax=ax)
    tsne = plot_tsne.copy()
    tsne.columns = ['x', 'y']
    ax.scatter(tsne["x"], tsne["y"],
               s=5,
               color="lightgrey",
               rasterized=True)
    ax.scatter(tsne.loc[cells, "x"],
               tsne.loc[cells, "y"],
               rasterized=True,
               s=20)
    ax.set_axis_off()
    return fig, ax
    # END OF FUNCTION FROM PALANTIR MODULE


def plot_gene_trend_heatmaps(gene_trends, outformat="jpeg",
                             save_file="", matrix=False):
    """
    Plot the gene trends on heatmap: a
    heatmap is generated or each branch
    ** Taken and adjusted from palantir to fix the multi-gene
    issue **

    :param gene_trends: Results of the compute_marker_trends function
    :param save_file: str
    :param matrix: bool

    """

    ############################################################################
    # 1. Define cell names & get genes
    ############################################################################
    branches = list(gene_trends.keys())

    # Get the genes from any index
    for e in gene_trends:
        if matrix:
            genes = gene_trends[e][
                "mat"].index.tolist()
            break
        else:
            genes = gene_trends[e][
                "trends"].index.tolist()
            break

    ############################################################################
    # 2. Save the ordered genes copy-friendly
    ############################################################################
    infofile = f"{save_file.replace('.pdf', '_LEGEND.txt')}"
    fi = open(infofile, "w")
    fi.writelines(f"{gene}\n" for gene in genes)
    fi.close()

    # Subdivide to avoid numpy error raised if too many genes
    n = 18
    final = [genes[i * n:(i + 1) * n]
             for i in range((len(genes) + n - 1) // n)]

    # Now plot for each chunk of gene list
    for no, genes in enumerate(final):
        # Plot height
        height = 0.7 * len(genes) * len(branches)

        #  Set up plot
        fig = plt.figure(figsize=[7, height])
        for i, branch in enumerate(branches):
            ax = fig.add_subplot(len(branches), 1, i + 1)
            # Standardize the matrix
            # Reduce df to genes
            if matrix:
                trend_df = gene_trends[branch]["mat"]
                trend_df = trend_df[trend_df.index.isin(genes)]
                mat = trend_df
            else:
                trend_df = gene_trends[branch]["trends"]
                trend_df = trend_df[trend_df.index.isin(genes)]

                mat = trend_df  # gene_trends[branch]["trends"]
                mat = pd.DataFrame(
                    StandardScaler().fit_transform(mat.T).T,
                    index=mat.index,
                    columns=mat.columns,
                )
            sns.heatmap(mat,
                        xticklabels=False,
                        ax=ax,
                        cmap=matplotlib.cm.Spectral_r,
                        rasterized=True,
                        vmin=-3,
                        vmax=3)
            ax.set_title(branch,
                         fontsize=12)
        save_file_finl = save_file.replace(".pdf",
                                           f"_{no}.{outformat}")
        plt.savefig(save_file_finl, dpi=100)
        plt.close()
    # END OF FUNCTION


def translate_cells(cells=[], adata="",
                    cluster_identifier="",
                    obs_name="controlstatus"):
    """

    Get condition and cluster
    by cell barcode index

    """

    transl = []

    for i, cell in enumerate(cells):
        transl.append(
            f"{adata.obs.loc[cell][obs_name]}"
            f"_{adata.obs.loc[cell][cluster_identifier]}"
            f"_{i}")
    return transl
    # END OF FUNCTION


def export_pr(pr_res, root_dir="",
              nn="", ncomp="", num_wp="",
              use_start=""):
    """

    To export palantir results.

    :param palantir_results: palantir obj
    :param root_dir: str
    :param nn: int
    :param ncomp: int
    :param num_wp: int
    :param use_start: bool

    """
    pr_res.branch_probs.to_csv(f"{root_dir}pal_"
                               f"{nn}_{ncomp}_{num_wp}"
                               f"_use{use_start}"
                               f"_BRANCHPROBS.csv")
    pr_res.pseudotime.to_csv(f"{root_dir}pal_"
                             f"{nn}_{ncomp}_{num_wp}"
                             f"_use{use_start}"
                             f"_PSEUDOTIME.csv")
    pr_res.entropy.to_csv(f"{root_dir}pal_"
                          f"{nn}_{ncomp}_{num_wp}"
                          f"_use{use_start}"
                          f"_ENTROPY.csv")
    pr_res.waypoints.to_frame(index=False,
                              name='waypoint'
                              ).to_csv(f"{root_dir}pal_"
                                       f"{nn}_{ncomp}_{num_wp}"
                                       f"_use{use_start}"
                                       f"_WAYPOINTS.csv")
    print("--- Exported results to csv.")
    # END OF FUNCTION


def stacked_barplot_pal(df, root_dir="", title="",
                        interest="condition"):
    """

    To plot cell proportions representing the
    pseudotime-ordered terminal cells from
    the palantir outputs.

    :param path_to_df:
    :param title:
    :return:

    """

    save_file = f"{root_dir}props_{title}.pdf"
    df = df.sort_values(interest)
    df = pd.DataFrame(df[interest].value_counts())
    df.columns = ["all_sum"]

    colors = []
    for e in df.index.values:
        col = GROUP_COLS[e]
        colors.append(col)

    df["fraction"] = df["all_sum"] / df["all_sum"].sum()
    del df["all_sum"]

    df = df.transpose()

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


def get_high_branchprob_cells(pr_res, terminal_cells=[],
                              root_dir="",
                              transl=[],
                              top_n=100, adata="",
                              cluster_identifier=""):
    """

    To find cells with top n high branching probs

    :param pr_res: palantir obj
    :param terminal_cells: list
    :param: transl: list
    :param top_n: int
    :param adata: AnnData obj
    :param cluster_identifier: str

    """

    # Initiate the df for cells with hightest branch probs
    high_probs_df = pd.DataFrame(columns=["branch",
                                          "branch_id",
                                          "cell_id",
                                          "branchprob",
                                          "condition_clust"
                                          "condition",
                                          cluster_identifier])
    # Embedding for highlighting cells
    embedding = pd.DataFrame(adata.obsm[f"X_draw_graph_fa"],
                             index=adata.obs_names,
                             columns=["x", "y"])

    # Iterate through terminal cells
    for i, tc in enumerate(terminal_cells):
        # Translate branch name
        tc_name = transl[i]
        # Sort by highest pranch probs for this branch
        sorted_frame = pr_res.branch_probs. \
            sort_values(by=tc,
                        ascending=False)
        # Select top_n cells
        hit_cells = sorted_frame[
                    :top_n].index.tolist()
        # Get their branch probs
        hit_probs = sorted_frame[
                    :top_n][tc].tolist()
        print(f"Top {top_n} cells :", hit_cells[:1], " ...")

        # Translate to condition
        hit_cells_hr = translate_cells(hit_cells,
                                       adata=adata,
                                       cluster_identifier=
                                       cluster_identifier)
        # Fill the dataframe
        high_probs_df["cell_id"] = hit_cells
        high_probs_df["branchprob"] = hit_probs
        high_probs_df["branch"] = tc
        high_probs_df["branch_id"] = tc_name
        high_probs_df["condition_clust"] = hit_cells_hr

        high_probs_df["condition"] = \
            high_probs_df["condition_clust"].str.split(
                "_", 1, expand=True)[0]

        high_probs_df[cluster_identifier] = \
            high_probs_df["condition_clust"].str.split(
                "_", 1, expand=True)[1]

        # Export the dataframe
        high_probs_df.to_csv(f"{root_dir}"
                             f"top_{top_n}"
                             f"_probs_{tc_name}")

        # Plot the proportions
        stacked_barplot_pal(high_probs_df,
                            root_dir=root_dir,
                            title=f"{tc_name}"
                                  f"_top_{top_n}")

        # Highlight cells on embedding
        print("--- Plotting cells on embedding. ")
        highlight_cells(embedding, hit_cells)
        plt.savefig(f"{root_dir}HIGHBRANCHPROB_"
                    f"top_{top_n}"
                    f"_{tc_name}.pdf",
                    dpi=400)
        plt.close()
    # END OF FUNCTION HIGH BRANCHPROBS

def run_magc_imputation(data, dm_res, n_steps=3):
    """
    **TAKEN FROM PALANTIR PACKAGE TO REDUCE STEPS

    Run MAGIC imputation

    :param dm_res: Diffusion map results from run_diffusion_maps
    :param n_steps: Number of steps in the diffusion operator
    :return: Imputed data matrix
    """
    print(data.obsm["X_draw_graph_fa"])
    if type(data) is sc.AnnData:
        data = pd.DataFrame(data.X.todense(),
                            index=data.obs_names,
                            columns=data.var_names)

    T_steps = dm_res["T"]
    print(T_steps)
    exit()
    n_steps = 1
    T_steps = dm_res["T"] ** n_steps
    print(T_steps)

    n_steps = 3
    T_steps = dm_res["T"] ** n_steps
    print(T_steps)


    imputed_data = pd.DataFrame(
        np.dot(T_steps.todense(), data), index=data.index, columns=data.columns
    )
    return imputed_data
    # END OF FUNCTION

def hts_easy(adata_path, cluster_identifier="", nn=20,
             ncomp=20, num_wp=1200, use_start=False,
             pca_proj="X_draw_graph_fa", hvg=False,
             lineage=False, id=False,
             lineage_info={"Hepatic": ["FH", "HB", "AH",
                                       "Hepatocyte","DC",
                                       "Cholangiocytes",
                                       "Ductal"],
                           "Stellate": ["ellate", "ibro", "uscle"]},
             impute_data=False
             ):
    """

    Palantir implementation, see for https://nbviewer.org/github/
    dpeerlab/Palantir/blob/master/notebooks/
    Palantir_sample_notebook.ipynb

    :param adata_path: str
    :param cluster_identifier: str
    :param nn: int
    :param ncomp: int
    :param num_wp: int
    :param use_start: bool
    :param pca_proj: str
    :param hvg: bool
    :param lineage: bool
    :param id: bool / str
    :param lineage_info: python dictionary
    ;param impute_data: bool
    :return: trajectory inference analysis

    """

    from tools.annotate import reduce_to_genes_expressed
    sns.set_style('ticks')
    matplotlib.rcParams['figure.figsize'] = [4, 4]
    matplotlib.rcParams['figure.dpi'] = 500
    matplotlib.rcParams['image.cmap'] = 'Spectral_r'

    ############################################################################
    # 1. Define cell names & paths
    ############################################################################
    target_lineage = lineage
    suffix = adata_path.rsplit(".h5ad", 1)[0]
    if lineage:
        save_file = f"{suffix}_palantir_{hvg}_{pca_proj}_{ncomp}_" \
                    f"{target_lineage}_lineage.h5ad"
    else:
        save_file = f"{suffix}_palantir_{hvg}" \
                    f"_{pca_proj}_{ncomp}.h5ad"

    # Directories
    save_folder = f"{suffix.rsplit('/', 1)[0]}/PALANTIR/"
    os.makedirs(save_folder,
                exist_ok=True)
    base_folder = f"{suffix.rsplit('/', 1)[0]}/PALANTIR/FILES/"
    os.makedirs(base_folder,
                exist_ok=True)
    if id:
        root_dir = f"{suffix.rsplit('/', 1)[0]}/PALANTIR/" \
                   f"{pca_proj}_{hvg}_{ncomp}_{nn}" \
                   f"_{num_wp}_{target_lineage}_{id}/"
    else:
        root_dir = f"{suffix.rsplit('/', 1)[0]}/PALANTIR/" \
                   f"{pca_proj}_{hvg}_{ncomp}_{nn}" \
                   f"_{num_wp}_{target_lineage}/"
    os.makedirs(root_dir, exist_ok=True)
    sc.settings.figdir = root_dir
    pal_res_path = f'{root_dir}palres_{pca_proj}_' \
                   f'{ncomp}_{nn}_{num_wp}' \
                   f'_{target_lineage}.pickle'
    magic_key = f"MAGIC_imputed_data_X_draw_graph_fa_{nn}_{lineage}"

    ############################################################################
    # 2.  Load normalized, log-transformed, hvg-red data
    ############################################################################
    if os.path.isfile(save_file):
        adata = sc.read_h5ad(save_file)
    else:
        adata = sc.read_h5ad(adata_path)

        if pca_proj not in list(adata.obsm):
            del adata.obsm['X_diffmap']
            sc.pp.neighbors(adata,
                            n_pcs=50)
            sc.tl.draw_graph(adata)
            adata.write_h5ad(adata_path)

        if lineage:
            cells = adata.obs[
                cluster_identifier
            ].value_counts().index.tolist()
            targets = find_cell(cells,
                                target_cells=
                                lineage_info[
                                    target_lineage],
                                all=True)
            keep = (adata.obs[cluster_identifier
                    ].isin(targets))
            adata = adata[keep, :]
            adata.write_h5ad(save_file)
            for ob in [cluster_identifier, "sample", "phase",
                       "condition", "controlstatus"]:
                try:
                    sc.pl.embedding(adata,
                                    basis=pca_proj,
                                    color=ob,
                                    save=ob)
                except:
                    print(f"--- Couldn't draw {ob}")
        if hvg:
            print("--- HVG")
            sc.pp.highly_variable_genes(adata,
                                        n_top_genes=hvg,
                                        flavor=
                                        "cell_ranger")
            adata = adata[:, adata.var.highly_variable]
            adata.write_h5ad(save_file)

    ############################################################################
    # 3.  Define projections
    ############################################################################

    pca_projections = pd.DataFrame(adata.obsm[pca_proj],
                                   index=
                                   adata.obs_names)

    embedding = pd.DataFrame(adata.obsm[
                                 f"X_draw_graph_fa"],
                             index=adata.obs_names,
                             columns=["x", "y"])

    ############################################################################
    # 4.  Define early cell based on marker gene expression
    ############################################################################
    if target_lineage == "Stellate":
        afp_pos = adata[adata[:, 'IGFBP2'].X > 0, :]
        afp_pos_cebpa_pos = afp_pos[afp_pos[:,
                                    'DDIT4'].X > 0, :]

        # Define highest afp
        max_val = afp_pos_cebpa_pos[:, 'IGFBP2'].X.max()
        early_cell = afp_pos_cebpa_pos[
                     afp_pos_cebpa_pos[
                     :, 'IGFBP2'].X == max_val, :
                     ].obs.index.tolist()[0]
    else: # Hepatocyte lineage
        afp_pos = adata[adata[:, 'AFP'].X > 0, :]
        afp_pos_cebpa_pos = afp_pos[afp_pos[:,
                                    'CEBPA'].X > 0, :]

        # Define highest afp
        max_val = afp_pos_cebpa_pos[:, 'AFP'].X.max()
        early_cell = afp_pos_cebpa_pos[
                     afp_pos_cebpa_pos[
                     :, 'AFP'].X == max_val, :
                     ].obs.index.tolist()[0]

    print(f"--- Early cell {early_cell}")

    ############################################################################
    # 5.  Save metrics
    ############################################################################

    infofile = f"{root_dir}metrics.txt"
    info = f"{adata.shape}, " \
           f"{adata.obs['condition'].value_counts()}" \
           f"{adata.obs['composite'].value_counts()}" \
           f"{adata.obs['composite_simple'].value_counts()}"
    fi = open(infofile, "w")
    fi.write(info)
    fi.close()

    if magic_key not in list(adata.layers):
        print("--- Imputation missing")
        dm_res_path = f"dms_{lineage}"
        if not os.path.isfile(dm_res_path):
            dm_res = palantir.utils.run_diffusion_maps(
                pca_projections,
                n_components=ncomp)

            with open(dm_res_path, 'wb') as file:
                print("--- Preparing DM for MAGIC")
                pickle.dump(dm_res, file)

        if os.path.isfile(dm_res_path):
            with open(dm_res_path, "rb") \
                    as file2:
                dm_res = pickle.load(file2)
        ########################################################################
        # 6.   MAGIC imputation
        ########################################################################
        if impute_data:
            adata.layers[magic_key] = \
                palantir.utils.run_magic_imputation(adata,
                                                    dm_res)
            adata.write_h5ad(save_file)
            print("--- Wrote imputed data.")

    if not os.path.isfile(pal_res_path):
        print(f"--- Not found: {pal_res_path}")
        save_folder = f"{root_dir}pal_" \
                      f"{nn}_{ncomp}_{num_wp}" \
                      f"_use{use_start}_{impute_data}/DMs/"
        os.makedirs(save_folder,
                    exist_ok=True)
        sc.settings.figdir = save_folder

        # DM
        dm_res = palantir.utils.run_diffusion_maps(
            pca_projections,
            n_components=ncomp)
        palantir.plot.plot_diffusion_components(embedding,
                                                dm_res)
        plt.savefig(f"{save_folder}{ncomp}_DM.png",
                    dpi=600)
        plt.close()

        # MS
        ms_data = palantir.utils.determine_multiscale_space(
            dm_res)

        ########################################################################
        # 5. Palantir results
        ########################################################################
        pr_res = palantir.core.run_palantir(ms_data,
                                            early_cell,
                                            knn=nn,
                                            num_waypoints=num_wp,
                                            terminal_states=None,
                                            use_early_cell_as_start=use_start,
                                            n_jobs=30)
        # Save
        try:
            with open(pal_res_path, 'wb') as file:
                pickle.dump(pr_res, file)
        except:
            print("--- Pickle failed.")

    if os.path.isfile(pal_res_path):
        # Load
        with open(pal_res_path, "rb") \
                as file2:
            pr_res = pickle.load(file2)

        # export
        export_pr(pr_res, root_dir=root_dir, nn=nn, ncomp=ncomp,
                  num_wp=num_wp, use_start=use_start)

        ########################################################################
        # 6. FIND CELLS OF INTEREST
        ########################################################################
        terminal_cells = pr_res.branch_probs.columns
        transl = translate_cells(terminal_cells,
                                 adata=adata,
                                 cluster_identifier=
                                 cluster_identifier)

        # Analyse cells with highest branching probs
        get_high_branchprob_cells(pr_res,
                                  top_n=100,
                                  root_dir=root_dir,
                                  terminal_cells=terminal_cells,
                                  transl=transl,
                                  adata=adata,
                                  cluster_identifier=
                                  cluster_identifier)
        terminal_states = pd.Series(transl, index=terminal_cells)
        terminal_states.to_csv(f"{root_dir}pal_"
                               f"_{nn}_{ncomp}_{num_wp}"
                               f"_use{use_start}"
                               f"_terminals.csv")

        cells = terminal_states.index.tolist()
        cells.append(early_cell)
        pr_res.branch_probs.columns = terminal_states[
            pr_res.branch_probs.columns]

        ########################################################################
        # 7. General visualization
        ########################################################################
        visualize = False
        if visualize:
            for e in ["draw_graph_fa"]:
                save_folder = f"{root_dir}OBSERVATIONS/"
                os.makedirs(save_folder, exist_ok=True)
                sc.settings.figdir = save_folder

                # OBSERVATIONS
                for ob in [cluster_identifier, "sample", "phase",
                           "condition", "controlstatus"]:
                    if ob in list(adata.obs):
                        sc.pl.embedding(adata,
                                        basis=e,
                                        color=ob,
                                        save=ob)

                # TERMINAL AND START CELLS
                print("--- Plotting cells on embedding. ")
                highlight_cells(embedding, cells)
                plt.savefig(f"{save_folder}pl_res_cells_{e}.pdf",
                            dpi=400)
                plt.close()

                # PALANTIR RESULTS
                plot_pal_res(pr_res, embedding)
                plt.savefig(f"{save_folder}pal-res_{e}_{nn}_"
                            f"{ncomp}_{pca_proj}_{num_wp}.pdf",
                            dpi=400)
                plt.close()

        ########################################################################
        # 8. Gene trends for fibrosis/inflammation-specif genes
        ########################################################################
        save_folder_top = f"{root_dir}GENE_TRENDS_IMP-{impute_data}/"
        os.makedirs(save_folder, exist_ok=True)
        genes = []

        geneset = True
        if geneset:
            genesets_sort = [geneset2path_inflammation, geneset2path_fibrosis]
            for i, g in enumerate(genesets_sort):
                save_folder = f"{save_folder_top}{i}/"
                genetrends_path = f'{root_dir}palres_{pca_proj}_' \
                                  f'{ncomp}_{nn}_{num_wp}' \
                                  f'_{target_lineage}_genetrends_{i}_' \
                                  f'{impute_data}.pickle'
                os.makedirs(save_folder, exist_ok=True)
                sc.settings.figdir = save_folder
                trend_id = "genesets"
                for i, geneset in enumerate(g):
                    print(f"--- Trends for {geneset}")
                    # 1.1 Read the custom genelist from config
                    try:
                        df = pd.read_csv(geneset2path[geneset])
                    except:
                        df = pd.read_table(geneset2path[geneset])
                    # 1.2 Generate a boolean variable in the df
                    geneset_list = df.iloc[:, 0].tolist()
                    # 1.3  Generate the genelist object
                    geneset_list = remove_duplicates \
                        (adjust_list(geneset_list, adata.var[
                            'gene_ids']))

                    # Make dict
                    gene_dict = {i: geneset_list}
                    gene_dict = reduce_to_genes_expressed(
                        adata, gene_dict)
                    gene_lst = gene_dict[i]

                    genes += [e for e in gene_lst if e not in genes]

                    # Save metrics individually
                    infofile = f"{save_folder}{geneset}.txt"
                    fi = open(infofile, "w")
                    fi.writelines(f"{gene}\n" for gene in gene_lst)
                    fi.close()
                    # END OF LOOP

                # Save metrics of ALL genes
                infofile = f"{save_folder}ALL_GENES.txt"
                fi = open(infofile, "w")
                fi.writelines(f"{gene}\n" for gene in genes)
                fi.close()

                ################################################################
                # 8.1 Now calc trends, save, and sort based on final expr.
                ################################################################
                # Calc and save gene trends
                if not os.path.isfile(genetrends_path):
                    print("Imputation retrieval")
                    imp_df = pd.DataFrame(adata[:, genes].layers[magic_key],
                                          index=adata.obs_names,
                                          columns=genes)
                    print("Calculating trends")
                    gene_trends = palantir.presults.compute_gene_trends(
                        pr_res, imp_df.loc[:, genes])
                    with open(genetrends_path, 'wb') as file:
                        pickle.dump(gene_trends, file)

                # Load genetrends
                if os.path.isfile(genetrends_path):
                    with open(genetrends_path, "rb") as file2:
                        gene_trends = pickle.load(file2)

                # Plot
                print("Plotting imputed / non-imp.")
                sc.pl.draw_graph(adata,
                                 color=genes,
                                 layer="scaled",
                                 save=f"scaled_n={len(genes)}")
                sc.pl.draw_graph(adata,
                                 color=genes,
                                 layer=f"MAGIC_imputed_data_X_draw_graph_fa_30_{lineage}",
                                 save=f"imputed_n={len(genes)}")

                # Iterate through each terminal state
                cells_of_interest = transl

                heatmaps = False
                if heatmaps:
                    for coi in cells_of_interest:
                        for celltype in transl:
                            trend_df = gene_trends[celltype]["trends"]
                            print(trend_df)
                            # Add the matrix
                            mat = pd.DataFrame(
                                StandardScaler().fit_transform(trend_df.T).T,
                                index=trend_df.index,
                                columns=trend_df.columns,
                            )
                            # add to sorted dict
                            gene_trends[celltype]["mat"] = mat


                            # For the case the cell is the cell of interest
                            if coi in celltype:
                                print(f"--- Found {coi}, sorting after its trends.")
                                #Sort descencing
                                mat = mat.sort_values(by=[mat.columns[-1],
                                                          mat.columns[-50],
                                                          mat.columns[-100],
                                                          mat.columns[0]],
                                                      ascending=[False, False,
                                                                 False, True])
                                gene_trends[celltype]["mat"] = mat
                                sorter = gene_trends[celltype]["mat"].index.tolist()

                        # Now sort the others in same order
                        for celltype in transl:
                            if coi not in celltype:
                                if sorter:
                                    print("Sorting")
                                    gene_trends[celltype]["mat"] = gene_trends[
                                        celltype]["mat"].reindex(sorter)
                        ############################################################
                        # 8.2 Visualize
                        ############################################################
                        print(f"{save_folder}{trend_id}_overview_{coi}.pdf",)
                        print("GENETRENDS HEATMAP")
                        plot_gene_trend_heatmaps(gene_trends,
                                                 save_file=
                                                 f"{save_folder}"
                                                 f"{trend_id}_overview_{coi}.pdf",
                                                 matrix=True)
                        # END OF LOOP
    # END OF FUNCTION
# END OF SCRIPT