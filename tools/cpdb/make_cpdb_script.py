root_file = "all_scrub_clean_0.5_qc_dr_OS_qc_clust_"
script_loc = "$HOME/work/applications/scrnaseq/cpdb/cpdb_plot.sh"
key = "metadata__composite_simple"

for cond in ["CTRL", "CTRL-OA", "CTRL-TGFB1", "CTRL-PA", "OA500", "PA", "TGFB1"]:
    print(f"sbatch {script_loc} {root_file}{cond}.h5ad METADATA/{cond}/METADATA/{key}.txt {cond}")