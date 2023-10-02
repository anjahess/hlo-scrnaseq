"""

Python script for the detection of rare cell types
10X scRNAseq data.
(Cluster Independent Algorithm for the identification
of markers of RAre cell types)
https://github.com/ScialdoneLab/CIARA_python

@author: Anja Hess
@date: 2023-MAY-01

"""

import scanpy as sc
import os
import seaborn as sns
from ciara_python import get_background_full, ciara

sc.settings.verbosity = 3
sc.logging.print_header()
sc.settings.file_format_figs = 'svg'
sc.settings.dpi_save = 10
sc.settings.verbosity = 0
sc.settings.autosave = True
sc.settings.set_figure_params(dpi=300,
                              fontsize=14,
                              figsize=[3, 3])
cmap = sns.diverging_palette(220, 20, as_cmap=True)
from matplotlib.pyplot import rc_context


def run_ciara(path_to_mtx="", file_id="", n_top=30):
    """

    :param path_to_mtx: str
    :param file_id: str
    :param n_top: int, how many genes to show
    :return:

    """
    ########################################################################
    # Define names and dirs
    ########################################################################
    file_name = file_id.rsplit(".h5ad")[0]
    outdir = f"{file_id}/"
    outpath = f"{file_id}/"
    save_file = f"{outpath}/{file_id}_ciara.h5ad"
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(outpath, exist_ok=True)
    sc.settings.figdir = outpath

    print(f"--- Searching rare celltypes for {file_name}.")
    if not os.path.isfile(save_file):
        ####################################################################
        # 1. Load the query data
        ####################################################################
        if ".h5ad" in file_id:
            adata = sc.read_h5ad(path_to_mtx)
            if "OS" in file_id:
                # Analyse undisturbed (controls) only
                keep = (adata.obs["controlstatus"]
                        == "CTRL")
                qDat = adata[keep, :]
                print(qDat.obs["controlstatus"].value_counts())
        else:
            adata = sc.read_10x_mtx(path_to_mtx +
                                   "/filtered_feature_bc_matrix/",
                                   var_names='gene_symbols',
                                   cache=True,
                                   gex_only=True)

        ####################################################################
        # 2. Get background
        ####################################################################
        get_background_full(adata,
                            threshold=1,
                            n_cells_low=2,
                            n_cells_high=5)
        ####################################################################
        # 3. Run ciara
        ####################################################################
        ciara(adata, n_cores=4,
              p_value=0.001,
              approximation=True,
              local_region=1)

        top_markers = adata.var.nsmallest(n_top,
                                          ["CIARA_p_value"])

        with rc_context({'figure.figsize': (3, 3)}):
            sc.pl.umap(adata,
                       color=top_markers.index.tolist(),
                       save=f"_{file_name}_{n_top}")
        import numpy as np
        np.savetxt(f"{outpath}ciara_topmarkers.csv",
        top_markers.index.tolist(),
                   delimiter=",",
                   fmt='%s')
        adata.write_h5ad(save_file)
        print("--- Done and wrote.")
    else:
        print("Ciara annotated version exists.")
        adata = sc.read_h5ad(save_file)
        top_markers = adata.var.nsmallest(n_top,
                                          ["CIARA_p_value"])

        import numpy as np
        np.savetxt(f"{outpath}{file_id}_topmarkers.csv",
        top_markers.index.tolist(),
                   delimiter=",",
                   fmt='%s')

        with rc_context({'figure.figsize': (3, 3)}):
            sc.pl.umap(adata,
                       color=top_markers.index.tolist(),
                       save=f"_{file_name}_{n_top}")
    # ENF OF FUNCTION


############################################################################
# START THE ANALYSIS
############################################################################

# 1. ULA-HLOs
dir = "$HOME/"
file_id = "all_scrub_clean_0.5_qc_dr_Day21_qc_clust.h5ad"
run_ciara(path_to_mtx=f"{dir}/{file_id}",
          file_id=file_id)

# 2. OS-HLOs
file_id = "all_scrub_clean_0.5_qc_dr_OS_qc_clust.h5ad"
run_ciara(path_to_mtx=f"{dir}/{file_id}",
          file_id=file_id)

# END OF SCRIPT