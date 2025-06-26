rm(list=ls())

### Find Marker genes for subclusters
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(dplyr)
################################################################################################################################################################
# Data preparation
################################################################################################################################################################
### read original data ###
data.dir <- '../../../OriginalData/'

results.dir <- '../../../Results/Results_no_ptprc_adgre1/Cluster_Analysis_avg_link_gene_cell_corr/'

output.dir <- paste0(results.dir,"DEG_results_fc15/")

output.dir.pval.cutoff <- paste0(results.dir,"DEG_p005/DEG_results_fc15_p005/")

#output.dir.pval.adj.cutoff <- paste0(results.dir,"p_val_adj_cutoff/DEG_results_fc15/")

AUC.mtx.data.dir<- paste0(results.dir,'AUC_mtx_with_clusters/')

## Load the meta data into a data frame
meta_data <- read.csv(paste0(data.dir,'seurat.integrated.5Ht_6Ho.metadata.csv'),row.names=1)
## Load count data
count <- read.csv(paste0(data.dir,'seurat.integrated.5Ht_6Ho.counts.csv'), header = TRUE,row.names = 1,check.names = FALSE)
# colnames(count) <- gsub(".", "-", colnames(count), fixed = TRUE)

### read subcluster info ###
files <- list.files(AUC.mtx.data.dir, pattern = "_AUC_with_clusts.csv", full.names = TRUE)
#Skip basal for now since there is an error
#files <- c(files[2:3],files[5:length(files)]) 

all_data <- lapply(files, function(file) {
  data <- read.csv(file, row.names = NULL)[, 1:2] 
  type <- gsub("_AUC_with_clusts.csv", "", basename(file)) # Extract type 
  data[, 2] <- paste(type, data[, 2], sep = "_") 
  colnames(data) <- c("cell_id", "subcluster") # Set appropriate column names
  return(data)
})
# Combine all data frames into one
combined_data <- bind_rows(all_data)

### add 2 cols to seurat object: sub_info, cell_type ###
meta_data$cell_id <- rownames(meta_data)
# Merge
meta_data <- merge(meta_data, combined_data, by = "cell_id", all.x = TRUE)
# Replace NA in subcluster with 'undecided'
meta_data$subcluster[is.na(meta_data$subcluster)] <- "undecided"
# Rename 'subcluster' to 'sub_info'
names(meta_data)[names(meta_data) == "subcluster"] <- "sub_info"

meta_data[["cell_type"]] <- ifelse(meta_data[["sub_info"]] == "undecided", 
                                   "undecided", 
                                   sub("^[^_]*_([^_]*)_.*$", "\\1", meta_data[["sub_info"]]))


meta_data <- meta_data[match(colnames(count), meta_data$cell_id), ]
row.names(meta_data) <- meta_data$cell_id
# create seurat object
seuratObj <- CreateSeuratObject(counts = count, project = 'project', meta.data = meta_data,min.cells = 5,min.features = 300)
seuratObj <- subset(seuratObj, subset = nFeature_RNA > 300 & nFeature_RNA < 5000)

sub_info_table <- table(seuratObj@meta.data$sub_info)

# Output the table to a text file
write.table(sub_info_table, file = file.path(output.dir, "sub_info_distribution.txt"), sep = "\t", quote = FALSE, col.names = NA)

seuratObj <- NormalizeData(seuratObj)

##################Subset basal cells with No_ptprc_no_adgre1############################
################################################################
metadata_No_ptprc_no_adgre1 = read.csv('//wsl.localhost/Ubuntu/home/akshatgupta/wwLab/wwlab-grn/MMG_Analysis/Processed_Data/metaData_basal_5ht6ho_without_PTPRC_Adgre1_filtered.csv')

cells_to_keep_No_ptprc_no_adgre1 = metadata_No_ptprc_no_adgre1$cell
obj_No_ptprc_no_adgre1 = subset(seuratObj,cells = cells_to_keep_No_ptprc_no_adgre1)
################################################################################################################################################################
# Find Markers
################################################################################################################################################################
# Function to perform analysis and output results
perform_analysis <- function(unique_type, seuratObj, output_path,pval_cutoff = FALSE,pval_adj_cutoff = FALSE) {
  # Create output folder named by type
  type_dir <- file.path(output_path, unique_type)
  dir.create(type_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Subset cells of the specific type from Seurat object
  subset_cells <- subset(seuratObj, subset = sub_info == unique_type)
  
  unique_type_cellType <- sub("^[^_]*_([^_]*)_.*$", "\\1", unique_type)
  
  # # Condition 1: Compare subclusters in 5Ht vs the corresponding subclusters in 6Ho
  corresponding_5Ht_type <- gsub("6ho", "5ht", unique_type)
  if (corresponding_5Ht_type %in% seuratObj$sub_info) {
     subset_6Ho_cells <- subset(seuratObj, sub_info == corresponding_5Ht_type)
     Idents(seuratObj) <- "sub_info"
     markers_cross_ho_ht <- FindMarkers(seuratObj, ident.1 = unique_type, ident.2 = corresponding_5Ht_type, only.pos = FALSE, min.pct = 0.15, min.diff.pct = 0.15, test.use = "MAST",logfc.threshold = 0.585)
     markers_cross_ho_ht <- data.frame(gene = rownames(markers_cross_ho_ht), markers_cross_ho_ht)
     markers_cross_ho_ht <- markers_cross_ho_ht[, c("gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")]
     if(pval_cutoff == TRUE){
       markers_cross_ho_ht <- markers_cross_ho_ht[markers_cross_ho_ht$p_val < 0.05, ] 
     }else if(pval_adj_cutoff == TRUE){
       markers_cross_ho_ht<- markers_cross_ho_ht[markers_cross_ho_ht$p_val_adj < 0.05, ]
     }
     write.table(markers_cross_ho_ht, file = file.path(type_dir, paste(unique_type, "5Ht_6Ho.txt", sep = "_")), sep = "\t", row.names = FALSE, col.names = TRUE)
   }
  
}

unique_types <- unique(combined_data$subcluster[grepl("6Ho", combined_data$subcluster, ignore.case = TRUE)])

#Apply p_value cutoff
for (type in unique_types) {
  perform_analysis(type, obj_No_ptprc_no_adgre1, output.dir.pval.cutoff,pval_cutoff = TRUE)
}
