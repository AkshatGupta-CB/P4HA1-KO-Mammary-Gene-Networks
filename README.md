# Single-Cell-Gene-Networks-in-P4HA1-Knockout-Mouse-Mammary-Glands

# Code/Network_Analysis
Pipeline for reconstructing gene-regulatory networks (GRNs) in two mouse groups (5Ht and 6Ho).
Please ensure all input and output paths are changed to match your needs.

Versions of key packages have been listed in `packages.txt`

---

## Cell annotation by reference dataset

**R Notebook**
- `Cell_annotation_by_reference_dataset.Rmd`

**Purpose**

- Use an external reference dataset ([GSE216542](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE216542)) for cell type annotation.

- Cross validate with the cell type label of our dataset to ensure more reliable identification of basal epithelial.


**Input**
- Reference: `GSE216542_RNA_Metadata_Final.csv.gz`, `GSE216542_RNA_Counts_Final.rds.gz` ([download link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE216542))
- Target: `Data/seurat.integrated.5Ht_6Ho.counts.tar.xz` (must be de-compressed), `Data/seurat.integrated.5Ht_6Ho.metadata.csv`


**Output**   
- `5Ht_with_annotation_by_ref_data.csv`
- `6Ho_with_annotation_by_ref_data.csv`

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
- `Data/seurat.integrated.5Ht_6Ho.metadata.csv`

**Output**   
- `NetworkData_HVGs_basal_5ht6ho_without_PTPRC_Adgre1.h5ad`

---

## Step 1 – GRN inference (20× replicates)

**Notebooks**:

- `pySCENIC_step1_GRN_5ht_basal_no_ptprc_adgre1.ipynb`
- `pySCENIC_step1_GRN_6ho_basal_no_ptprc_adgre1.ipynb`

**Purpose**

- Run GRNBoost2 twenty times to infer gene network and reduce stochasticity

**Input**   
- `NetworkData_HVGs_basal_5ht6ho_without_PTPRC_Adgre1.h5ad`
- List of Transcription Factors in Mus Musculus. Also provided here - `Data/allTFs_mm.txt`

**Output**   
- Twenty CSV files (for each mouse group) containing adjancencies inferred by GRNBoost2

---

## Step 2 – Motif enrichment

**Notebooks**

- `pySCENIC_step2_5ht_basal_no_ptprc_adgre1.ipynb`
- `pySCENIC_step2_6ho_basal_no_ptprc_adgre1.ipynb`

**Purpose**

- Convert adjacencies → putative regulons using motif databases

**Input**
- Feather files for motif enrichment - can be download here: https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc_v10_clust/gene_based/
- Twenty CSV files containing adjancencies inferred by GRNBoost2
- List of Transcription Factors in Mus Musculus. Also provided here - `Data/allTFs_mm.txt`

**Output**  
- Twenty Pickle Files (for each mouse group) containing regulons

## Step 2b – QC / Statistics

**Notebooks**

- `Get_stats_5Ht_basal_step2.ipynb`
- `Get_stats_6ho_basal_step2.ipynb`

**Purpose**

- Summarise TFs, targets, and replicate frequency for each preliminary regulon

## Step 2c – Aggregation of regulons across all runs

**Notebooks**

- Agg_regulons_step2_basal_no_ptprc_adgre1.ipynb

**Purpose**

- Aggregare 20 replicates into consensus regulons

**Input** 
- Twenty Pickle Files (for each mouse group) containing regulons

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

## Step 1 – Identify significant regulons between the two mouse groups 

**Notebook**  
- `Basal_find_sig_regulon_Wilcox.ipynb`  

**Purpose**  
- Append sub-cluster labels to the per-cell AUC matrix.  
- For each sub-cluster, run Wilcoxon rank-sum tests (5Ht vs 6Ho) on every regulon.  
- Keep regulons with **BH-adjusted p < 0.05** *and* **|ΔAUC| ≥ 0.15**.  

**Input**  
- Scaled AUC matrices (from *Prepare AUC matrices for clustering*)  
- List of cells in each subcluster/Sub-cluster labels (.txt)  

**Output**  
- AUC matrices (.csv) with subcluster information
- Significantly different regulons (.txt) for each subcluster
  
---

## Step 2 – Differential gene expression per sub-cluster  

**Scripts**
- `DGE_for_subcluster/find_markerGenes_subclusters_5ht.R`  
- `DGE_for_subcluster/find_markerGenes_subclusters_6ho.R`  

**Purpose**  
- Perform DGE for each subcluster in 5Ht and 6Ho
- For similar subclusters, DGE is peformed between subcluster in 5Ht vs 6Ho
- For unique subclusters, DGE is performed between subcluster vs all other cells in that mouse type

**Input**  
- `Data/seurat.integrated.5Ht_6Ho.counts.tar.xz` (must be de-compressed to get csv file first)
- `data/seurat.integrated.5Ht_6Ho.metadata.csv`
- `metadata` with cell information of cells used in network analysis. Can be obtained from the `obs` layer in `NetworkData_HVGs_basal_5ht6ho_without_PTPRC_Adgre1.h5ad`
-  AUC matrices (.csv) with subcluster information (from **Identify significant regulons between the two mouse groups**)

**Output**  
- Differentially expressed genes (.txt) for each subcluster in 5Ht and 6Ho

---

## Step 3 – Identify differentially-expressed target genes (DETGs) for significantly different regulons

**Notebook**  
- `Get_differential_regulons_Basal_fc15_p005.ipynb`  

**Purpose**  
- Identify DETGs of significant regulons (Step 1) with DGE lists (Step 2) to obtain **DETGs** per regulon, per sub-cluster

**Input**  
- Significantly different regulons (.txt) for each subcluster (Step 1)
- Differentially expressed genes (.txt) for each subcluster in 5Ht and 6Ho
  
**Output**  
- DETGs for significant regulons (.csv) in each subcluster of the 5Ht and 6Ho mice

---

## Step 5 – Summarise DETG statistics  
**Notebook**  
- `Get_Stats_basal_fc15_p005.ipynb`  

**Purpose**  
- Process and summarize DETG results from step 4.  

**Input**  
- DETGs for significant regulons (.csv) in each subcluster of the 5Ht and 6Ho mice (from Step 4)  

**Output**  
- Statistics and Summary of DETGs identfified in each regulon in each subcluster of the 5Ht and 6Ho mice

---

## UMAP & centroid analysis 

**Script**  
- `UMAP_centroid_analysis_reg_act_plots.R`  

**Purpose**  
- Run UMAP on the AUC matrices.  
- Plot UMAPs coloured by sub-cluster.  
- Compute Euclidean distances between the centroids of all sub-clusters (5Ht ↔ 6Ho)

**Input**  
- AUC matrices (.csv) with subcluster information (From **Identify significant regulons between the two mouse groups**) 

**Output**  
- UMAP plots of AUC matrices with subcluster annotation
- Distances between centroids of pairs of subclusters (5Ht ↔ 6Ho)

---

# Code/Figures

Code to generate all figures is present in this folder. The file name includes the figures it generates. 

---

##### For any questions, please contact akshatg2@alumni.cmu.edu.
