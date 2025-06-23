# Single-Cell-Gene-Networks-in-P4HA1-Knockout-Mouse-Mammary-Glands

# Code/Network_Analysis
Pipeline for reconstructing gene-regulatory networks (GRNs) in two mouse genotypes (5Ht and 6Ho).

## Pre-processing

### Notebook: 
- PreprocessHVGs5_no_ptprc_adgre1.ipynb

### Purpose:

- Load raw counts & metadata

- Filter low-quality cells/genes, select HVGs

- Normalise & log-transform

#### Input: data/counts.tsv, data/metadata.tsv

#### Output: data/processed/processed_{strain}.h5ad

## Step 1 – GRN inference (20× replicates)

### Notebooks:

- pySCENIC_step1_GRN_5ht_basal_no_ptprc_adgre1.ipynb

- pySCENIC_step1_GRN_6ho_basal_no_ptprc_adgre1.ipynb

### Purpose: 

- Run GRNBoost2 twenty times to infer gene network and reduce stochasticity

#### Input: processed .h5ad

#### Output: results/step1/run*/adjacencies.tsv

## Step 2 – Motif enrichment

### Notebooks:

- pySCENIC_step2_5ht_basal_no_ptprc_adgre1.ipynb

- pySCENIC_step2_6ho_basal_no_ptprc_adgre1.ipynb

### Purpose:

- Convert adjacencies → putative regulons using motif databases

#### Input: step-1 outputs + data/motifs/* + data/tfs.txt

#### Output: regulons_{rep}.pkl (per replicate)

## Step 2b – QC / Statistics

### Notebooks: 

- Get_stats_5Ht_basal_step2.ipynb

- Get_stats_6ho_basal_step2.ipynb

### Purpose: 

- Summarise TFs, targets, and replicate frequency for each preliminary regulon

## Step 2c – Aggregation

### Notebook

- Agg_regulons_step2_basal_no_ptprc_adgre1.ipynb

### Purpose: 

- Merge 20 replicates into consensus regulons

#### Output: Two Pickle files containing aggregated regulons from all runs in the two mouse types 

## Step 3 – AUCell

### Notebooks:

- pySCENIC_step3_Agg_5ht_basal_no_ptprc_adgre1.ipynb

- pySCENIC_step3_Agg_6ho_basal_no_ptprc_adgre1.ipynb

### Purpose: 

- Compute regulon AUC per cell → activity matrices

#### Output: Two CSV files containing AUC matrices for the two mouse types 

