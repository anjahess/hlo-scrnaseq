#!/bin/bash

## Usage
# 1. Make sure there is the METADATA folder in your h5ad dir
# 2. Run: cpdb cpdb.sh META.txt OUTDIR
# 3. For multiple files: use make_cpdb_script.py
# 4. Then use custom scripts for downstream analysis: python3.6 cpdb_heatmap.py

conda_dir=$HOME/work/applications/miniconda3/
source $conda_dir"bin/activate" cpdb

counts=$1
meta=$2
out=$3

echo "Running cellphonedb on " $counts "meta:" $meta "SAVE:" $out
cellphonedb method statistical_analysis $meta $counts --counts-data gene_name --output-path $out
cellphonedb plot heatmap_plot --pvalues-path $out/pvalues.txt --output-path $out $meta && cellphonedb plot dot_plot --means-path $out/means.txt --pvalues-path $out/pvalues.txt --output-path $out

