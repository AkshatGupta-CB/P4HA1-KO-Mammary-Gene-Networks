rm(list=ls())

### Find Marker genes for subclusters
#library(nichenetr) # Please update to v2.0.4
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(dplyr)
################################################################################################################################################################
# Data preparation
################################################################################################################################################################
### read original data ###
data.path <- 'data/'
output_path <- "Results_061424/DEG_results_fc15/"
## Load the meta data into a data frame
meta_data <- read.csv(paste0(data.path,'seuratMD.csv'),row.names=1)
## Load count data
count <- read.csv(paste0(data.path,'seuratCounts.csv'), row.names=1)
colnames(count) <- gsub(".", "-", colnames(count), fixed = TRUE)

### read subcluster info ###
files <- list.files(data.path, pattern = "_auc_clusts.csv", full.names = TRUE)
#Skip basal for now since there is an error
#files <- c(files[2:3],files[5:length(files)]) 
all_data <- lapply(files, function(file) {
  data <- read.csv(file, row.names = NULL)[, 1:2] 
  type <- gsub("_auc_clusts.csv", "", basename(file)) # Extract type 
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

# because meta_data from Akshat's result has different cell order, that why we need to do -
# reorder mata_data to match count
meta_data <- meta_data[match(colnames(count), meta_data$cell_id), ]

# create seurat object
seuratObj <- CreateSeuratObject(counts = count, project = 'project', meta.data = meta_data,min.cells = 3,min.features = 200)
seuratObj <- subset(seuratObj, subset = nFeature_RNA > 200 & nFeature_RNA < 5000)

# if need to filter the cells undecided, uncommment below
#seuratObj <- subset(x = seuratObj, subset = (sub_info != "undecided"))

sub_info_table <- table(seuratObj[["sub_info"]])

# Output the table to a text file
write.table(sub_info_table, file = file.path(output_path, "sub_info_distribution.txt"), sep = "\t", quote = FALSE, col.names = NA)

seuratObj <- NormalizeData(seuratObj)
################################################################################################################################################################
# Find Markers
################################################################################################################################################################
# Function to perform analysis and output results
perform_analysis <- function(unique_type, seuratObj, output_path) {
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
     write.table(markers_cross_ho_ht, file = file.path(type_dir, paste(unique_type, "5Ht_6Ho.txt", sep = "_")), sep = "\t", row.names = FALSE, col.names = TRUE)
   }

  # Condition 2: Marker genes for this subcluster vs subclusters in same cell type (basal, etc.)
  #For unique subclusters
  if ("cell_type" %in% colnames(seuratObj@meta.data)) {

    same_type_cells <- subset(seuratObj, cell_type == unique_type_cellType & orig.ident == "6Ho")
    Idents(same_type_cells) <- "sub_info"
    
    markers_same_type <- FindMarkers(same_type_cells, ident.1 = unique_type,only.pos = FALSE, min.pct = 0.15, min.diff.pct = 0.15, test.use = "MAST",logfc.threshold = 0.585)
      #}
    markers_same_type <- data.frame(gene = rownames(markers_same_type), markers_same_type)
    markers_same_type <- markers_same_type[, c("gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")]
    write.table(markers_same_type, file = file.path(type_dir, paste(unique_type, "same_celltype.txt", sep = "_")), sep = "\t", row.names = FALSE, col.names = TRUE)
    }
}

unique_types <- unique(combined_data$subcluster[grepl("6Ho", combined_data$subcluster, ignore.case = TRUE)])


for (type in unique_types) {
  perform_analysis(type, seuratObj, output_path)
}
