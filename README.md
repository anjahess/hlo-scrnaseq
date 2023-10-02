Stratification of Human Liver Organoid Models for MASLD
===========

Author: Anja Hess

Date: 2023-MAY-27

## 1. Description
Command line interface and scripts to reproduce results from
the scRNA-seq study of human liver organoid models of MASLD. For publication details see below.
The main start script located in bin points to a folder containing ALL 10X samples of interest
or an .h5ad file. 

## 2. Before you start

**Check out our [Jupyter notebook](hlo_masld.ipynb) guiding through the analysis!**

These scripts are written for **python => 3.6**. 

Jobs were run in a **Unix** environment (required memory **~ 40-60 GB**).
Please make sure you installed dependencies for the individual packages (e.g. scanpy, cpdb).

## 3. Usage
    cd ./bin
    python3 start.py path/to/FOLDER/ run_name
    python3 start.py file.h5ad run_name

## 4. Required folder architecture:

### 4.1 For 10X cellranger outs:
    FOLDER
    -- SAMPLE1
        -- .mtx.gz
        -- .tsv.gz
        -- .tsv.gz
    -- SAMPLE2
        -- ...

### 4.2 For individual AnnData objects to be merged
    FOLDER
    -- sample1.h5ad
    -- sample2.h5ad


### 4.3 For an already concatenated AnnData object

    python3 start.py file.h5ad run_name

## 5. Data access

Raw sequence data for scRNA-seq performed in this study are publicly available on GEO under the accession GSE207889. 

Previously published scRNA-seq data that were used for analyses:
MacParland et al.: GSE115469, Wesley et al.: E-MTAB-8210 (accessed 2023-MAY-28 through https://collections.cellatlas.io/liver-development)

## 6. Citation

**Hess, A. et al. Single cell transcriptomic landscapes of human liver organoids stratify models of non-alcoholic fatty liver disease. 2022.07.19.500693 Preprint at https://doi.org/10.1101/2022.07.19.500693 (2023).**