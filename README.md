# Single-Cell-Gene-Networks-in-P4HA1-Knockout-Mouse-Mammary-Glands

# Code/Network_Analysis
Pipeline for reconstructing gene-regulatory networks (GRNs) in two mouse genotypes (5Ht and 6Ho).

---

## Pre-processing

**Notebook**
- `PreprocessHVGs5_no_ptprc_adgre1.ipynb`

**Purpose**

- Load raw counts & metadata

- Filter low-quality cells/genes, select HVGs

- Normalise & log-transform

**Input**

- `Data/seurat.integrated.5Ht_6Ho.counts.tar.xz` (must be de-compressed to get csv file first)
- `data/seurat.integrated.5Ht_6Ho.metadata.csv`

**Output**   

- data/processed/processed_{strain}.h5ad

---

## Step 1 – GRN inference (20× replicates)

**Notebooks**:

- `pySCENIC_step1_GRN_5ht_basal_no_ptprc_adgre1.ipynb`

- `pySCENIC_step1_GRN_6ho_basal_no_ptprc_adgre1.ipynb`

**Purpose**

- Run GRNBoost2 twenty times to infer gene network and reduce stochasticity

**Input**   

- processed .h5ad

**Output**   

- results/step1/run*/adjacencies.tsv

---

## Step 2 – Motif enrichment

**Notebooks**

- `pySCENIC_step2_5ht_basal_no_ptprc_adgre1.ipynb`

- `pySCENIC_step2_6ho_basal_no_ptprc_adgre1.ipynb`

**Purpose**

- Convert adjacencies → putative regulons using motif databases

**Input**

-
-

**Output**  

- 

## Step 2b – QC / Statistics

**Notebooks**

- `Get_stats_5Ht_basal_step2.ipynb`

- `Get_stats_6ho_basal_step2.ipynb`

**Purpose**

- Summarise TFs, targets, and replicate frequency for each preliminary regulon

## Step 2c – Aggregation

**Notebooks**

- Agg_regulons_step2_basal_no_ptprc_adgre1.ipynb

**Purpose**

- Merge 20 replicates into consensus regulons

**Output** 

- Two Pickle files containing aggregated regulons from all runs in the two mouse types 

---

## Step 3 – AUCell

**Notebooks**

- `pySCENIC_step3_Agg_5ht_basal_no_ptprc_adgre1.ipynb`

- `pySCENIC_step3_Agg_6ho_basal_no_ptprc_adgre1.ipynb`

**Purpose**

- Compute regulon AUC per cell → activity matrices

**Output**

- Two CSV files containing AUC matrices for the two mouse types 

# Code/Cluster_Analysis

Pipeline for **cluster-level** comparison of regulon activity and gene expression between the two mouse genotypes (**5Ht** vs **6Ho**).

---

## Prepare AUC matrices for clustering  
**Notebook**

- `Prepare_AUCMtx_for_cluster_analysis_no_ptprc_adgre1.ipynb`  

**Purpose**

- Scale the AUC matrices (generated in *Network Analysis Step 3*) and subset for common regulons between 5Ht and 6Ho for **Cluster 3.0** and **TreeView**.
- Reorders the 6Ho AUC matrix to match the regulon order of the clustered 5Ht AUC matrix  

**Input**  

- AUC matrices generated in *Network Analysis Step 3*
- .CDT file of the clustered 5Ht AUC matrix

**Output**  

- Scaled AUC matrices (.csv files) for 5Ht and 6Ho

---

## Step 1 – Differential gene expression per sub-cluster  

**Scripts**
- `DGE_for_subcluster/find_markerGenes_subclusters_5ht.R`  
- `DGE_for_subcluster/find_markerGenes_subclusters_6ho.R`  

**Purpose**  
- Perform DGE for each subcluster in 5Ht and 6Ho.
- For similar subclusters, DGE is peformed between subcluster in 5Ht vs 6Ho
- For unique subclusters, DGE is performed between subcluster vs all other cells in that mouse type.

**Input**  
- Raw or normalised counts matrix (`data/counts.tsv`).  
- `metadata` with sub-cluster labels.  

**Output**  
- `dge_results/5ht_subcluster*_dge.tsv`  
- `dge_results/6ho_subcluster*_dge.tsv`

---
## Step 1 – UMAP & centroid analysis 

**Script**  
- `UMAP_centroid_analysis_reg_act_plots.R`  

**Purpose**  
- Run UMAP on the AUC matrices.  
- Plot UMAPs coloured by sub-cluster.  
- Compute Euclidean distances between the centroids of all sub-clusters (5Ht ↔ 6Ho).  

**Input**  
- Scaled AUC matrices (from *Prepare AUC matrices for clustering*)  
- Cell lists for subclusters

**Output**  
- `figures/UMAP_regulon_activity.pdf`  
- `results/centroid_distance.tsv`

---

## Step 3 – Identify significant regulons between genotypes  
**Notebook**  
- `Basal_find_sig_regulon_Wilcox.ipynb`  

**Purpose**  
- Append sub-cluster labels to the per-cell AUC matrix.  
- For each sub-cluster, run Wilcoxon rank-sum tests (5Ht vs 6Ho) on every regulon.  
- Keep regulons with **BH-adjusted p < 0.05** *and* **|ΔAUC| ≥ 0.15**.  

**Input**  
- Merged AUC matrix from Step 0.  
- Sub-cluster labels.  

**Output**  
- `signif_regulons_per_subcluster.tsv`

---

## Step 4 – Retrieve differentially-expressed target genes (DETGs) for significantly different regulons  
**Notebook**  
- `Get_differential_regulons_Basal_fc15_p005.ipynb`  

**Purpose**  
- Intersect significant regulons (Step 3) with DGE lists (Step 2) to obtain **DETGs** per regulon, per sub-cluster.  

**Input**  
- `signif_regulons_per_subcluster.tsv`  
- `dge_results/*_dge.tsv`  

**Output**  
- `detg_lists/<subcluster>_detg.tsv`

---

## Step 5 – Summarise DETG statistics  
**Notebook**  
- `Get_Stats_basal_fc15_p005.ipynb`  

**Purpose**  
- Process DETG counts per regulon and per sub-cluster.  
- Produce summary tables and publication-ready plots.  

**Input**  
- `detg_lists/` (from Step 4)  

**Output**  
- `detg_stats_summary.tsv`  
- `figures/DETGs_summary_plots.pdf`

