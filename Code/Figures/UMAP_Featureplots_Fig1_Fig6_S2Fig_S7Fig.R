library(Seurat)
library(patchwork)
library(ggplot2)
library(readxl)

###################
#Fig 1A
#Load the metadata
rm(list=ls())
metaData <- read_xlsx(path = "5Ht_6Ho_files_whole-population_12212023_complete_clustering/5Ht_6Ho_files_whole-population_12212023_complete_clustering/meta_data_unified_v5_AG.xlsx")
counts <- read.csv(file = "seurat.integrated.5Ht_6Ho.counts.csv",header = TRUE)
# Set the first row as column names
row.names(counts) <- counts[,1]

# Remove the first row
counts <- counts[,-1 ]

obj <- CreateSeuratObject(counts = counts, project = 'project', meta.data = metaData,min.cells = 3,min.features = 200)
obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 5000)
metaData <- obj@meta.data

obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
obj <- RunUMAP(obj, dims = 1:15, reduction = "pca",reduction.name = "UMAP")

#Fig1A
Idents(obj) <- as.character(metaData$sample)
p <- DimPlot(obj, reduction = "UMAP",group.by = "ident") + theme(legend.position = 'none' ,legend.text = element_text(size = 8))+ theme(panel.border = element_rect(colour = "black", fill=NA, size=1),axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20),axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))
ggsave('../MAnuscript/Figures/Integrated_UMAP.png', plot = p, width = 10, height = 8)

#Fig1B
# Define cluster-to-color mapping
cluster_colors <- c(
  "Epithelial-Basal" = "tomato",
  "Epithelial-Luminal" = "orange",
  "ILC2" = 'seagreen',#updated
  "T" = "purple",
  "NK" = "steelblue",#updated
  "B" = "violetred4",
  "Macrophages" = "gold3",#updated
  "Monocytes" = "springgreen2",#updated
  "Neutrophils" = "cyan",
  "DC" = "aquamarine",
  "Endothelial" = "violet",
  "Collagen" = "dodgerblue1",#updated
  "MAST" = "magenta",
  "Other" = "lightgrey"  # Set "Other" to grey
)

# Generate the UMAP plot with the custom colors
Idents(obj) <- as.character(metaData$Reclassified_cluster)
p <- DimPlot(obj, reduction = "UMAP", group.by = "ident") +
  scale_color_manual(values = cluster_colors) +  # Apply custom color mapping
  theme(
    legend.position = 'none',
    legend.text = element_text(size = 8),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20)
  )

# Generate the labeled plot
labelled_plot <- LabelClusters(plot = p, id = "ident", repel = TRUE,fontface = "bold") 
  
# Display the plot
print(labelled_plot)
ggsave('Figures/Integrated_UMAP_unified_cluster_annot.png', plot = labelled_plot, width = 10, height = 8)

#S2 Fig
p <- DimPlot(obj, reduction = "UMAP",group.by = "ident",split.by = "ident") + theme(legend.position = 'none' ,legend.text = element_text(size = 8))+ theme(panel.border = element_rect(colour = "black", fill=NA, size=1),axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20),axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))
ggsave('../MAnuscript/Figures/Split_Integrated_UMAP.png', plot = p, width = 16, height = 8)
##############################################################################################
#S7 Fig
#Load the metadata
rm(list=ls())
data.dir = 'MMG_Analysis'
metadata.path = file.path(data.dir,'OriginalData/seurat.integrated.5Ht_6Ho.metadata.csv')
counts.path = file.path(data.dir,'OriginalData/seurat.integrated.5Ht_6Ho.counts.csv')
metadata_No_ptprc_no_adgre1.path = file.path(data.dir,'Processed_Data/metaData_basal_5ht6ho_without_PTPRC_Adgre1_filtered.csv')
figure.output.dir= 'Manuscript/Figures/'

metaData <- read.csv(file = metadata.path,row.names = 1)

#Load the counts
counts <- read.csv(file = counts.path,header = TRUE,row.names = 1,check.names = FALSE)


#Create the Seurat Object 
obj <- CreateSeuratObject(counts = counts, project = 'project', meta.data = metaData,min.cells = 3,min.features = 200)
obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 5000)
metaData <- obj@meta.data

#Normalize the data
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
obj <- RunUMAP(obj, dims = 1:15, reduction = "pca",reduction.name = "UMAP")


#Subset the data to only get basal cells that dont have Ptprc or Adgre1
metadata_No_ptprc_no_adgre1 = read.csv(metadata_No_ptprc_no_adgre1.path)
cells_to_keep_No_ptprc_no_adgre1 = metadata_No_ptprc_no_adgre1$cell
obj_No_ptprc_no_adgre1 = subset(obj,cells = cells_to_keep_No_ptprc_no_adgre1)

#####################################
#Read cell lists for clusters
#####################################
#5Ht
#5Ht cell list dir
cluster.cell.5ht.dir = file.path(data.dir,"Results/Results_no_ptprc_adgre1/Cluster3/gene_cell_corr/Average_linkage/5Ht/Reorder_cells_reg_paper/cluster_cell_list")
cluster.cell.6ho.dir = file.path(data.dir,"Results/Results_no_ptprc_adgre1/Cluster3/gene_cell_corr/Average_linkage/6Ho/Reorder_cells_reg_paper/cluster_cell_list")


#5Ht basal
ht5_basal_S1 <- readLines(paste0(cluster.cell.5ht.dir,'/S1_489cells_5Ht_list.txt'))
ht5_basal_S2_all <- readLines(paste0(cluster.cell.5ht.dir,'/S2-all_483cells_5Ht_list.txt'))
ht5_basal_S3 <- readLines(paste0(cluster.cell.5ht.dir,'/Hif_45cells_5Ht_list.txt'))
ht5_basal_S2_3 <- readLines(paste0(cluster.cell.5ht.dir,'/S2-3_42cells_5Ht_list.txt'))
ht5_basal_U1 <- readLines(paste0(cluster.cell.5ht.dir,'/U1_wt_5ht_29cells.txt'))
ht5_S2_1_2_3 <- setdiff(ht5_basal_S2_all, c(ht5_basal_S3, ht5_basal_U1))

#6Ho basal
ho6_basal_S1 <- readLines(paste0(cluster.cell.6ho.dir,'/S1_301cells_6Ho_list.txt'))
ho6_basal_S2_all <- readLines(paste0(cluster.cell.6ho.dir,'/S2-all_215cells_6Ho_list.txt'))
ho6_basal_S3 <- readLines(paste0(cluster.cell.6ho.dir,'/Hif_59cells_6Ho_list.txt'))
ho6_S2 <- setdiff(ho6_basal_S2_all, ho6_basal_S3)

##################################################################################
#CLUSTER 1
values_array_basal_5ht_S1 <- vector(mode = "character", length = nrow(metadata_No_ptprc_no_adgre1))
values_array_basal_6ho_S1 <- vector(mode = "character", length = nrow(metadata_No_ptprc_no_adgre1))
values_array_basal_both_S1 <- vector(mode = "character", length = nrow(metadata_No_ptprc_no_adgre1))
# Iterate over each row in the dataframe meta
for (i in 1:nrow(metadata_No_ptprc_no_adgre1)) {
  # Get the sample_id from the current row
  sample_id <- metadata_No_ptprc_no_adgre1$cell[i]
  
  # Initialize a variable to store the value
  value_1 <- "not_of_interest"  # Default value if the sample_id is not found in any list
  value_2 <- "not_of_interest" 
  value_3 <- "not_of_interest"
  if (sample_id %in% ht5_basal_S1) {
    
    value_1 <- paste("5Ht_basal_clust_S1-",length(ht5_basal_S1))
    value_3 <- paste("5Ht_basal_clust_S1-",length(ht5_basal_S1))
  }else if (sample_id %in% ho6_basal_S1) {
    
    value_2 <- paste("6ho_basal_clust_S1-",length(ho6_basal_S1))  
    value_3 <- paste("6ho_basal_clust_S1-",length(ho6_basal_S1))  
    
  }
  
  # Append the value to the array
  values_array_basal_5ht_S1[i] <- value_1
  values_array_basal_6ho_S1[i] <- value_2
  values_array_basal_both_S1[i] <- value_3
}

#CLUSTER 2
values_array_basal_5ht_S2 <- vector(mode = "character", length = nrow(metadata_No_ptprc_no_adgre1))
values_array_basal_6ho_S2 <- vector(mode = "character", length = nrow(metadata_No_ptprc_no_adgre1))
values_array_basal_both_S2 <- vector(mode = "character", length = nrow(metadata_No_ptprc_no_adgre1))
# Iterate over each row in the dataframe meta
for (i in 1:nrow(metadata_No_ptprc_no_adgre1)) {
  # Get the sample_id from the current row
  sample_id <- metadata_No_ptprc_no_adgre1$cell[i]
  
  # Initialize a variable to store the value
  value_1 <- "not_of_interest"  # Default value if the sample_id is not found in any list
  value_2 <- "not_of_interest" 
  value_3 <- "not_of_interest"
  if (sample_id %in% ht5_S2_1_2_3) {
    
    value_1 <- paste("5Ht_basal_clust_S2_1_2_3-",length(ht5_S2_1_2_3))
    value_3 <- paste("5Ht_basal_clust_S2_1_2_3-",length(ht5_S2_1_2_3))
  }else if (sample_id %in% ho6_S2) {
    
    value_2 <- paste("6ho_basal_clust_S2-",length(ho6_S2))  
    value_3 <- paste("6ho_basal_clust_S2-",length(ho6_S2))  
    
  }
  
  # Append the value to the array
  values_array_basal_5ht_S2[i] <- value_1
  values_array_basal_6ho_S2[i] <- value_2
  values_array_basal_both_S2[i] <- value_3
}

#CLUSTER Hif - S3
values_array_basal_5ht_S3 <- vector(mode = "character", length = nrow(metadata_No_ptprc_no_adgre1))
values_array_basal_6ho_S3 <- vector(mode = "character", length = nrow(metadata_No_ptprc_no_adgre1))
values_array_basal_both_S3 <- vector(mode = "character", length = nrow(metadata_No_ptprc_no_adgre1))
# Iterate over each row in the dataframe meta
for (i in 1:nrow(metadata_No_ptprc_no_adgre1)) {
  # Get the sample_id from the current row
  sample_id <- metadata_No_ptprc_no_adgre1$cell[i]
  
  # Initialize a variable to store the value
  value_1 <- "not_of_interest"  # Default value if the sample_id is not found in any list
  value_2 <- "not_of_interest" 
  value_3 <- "not_of_interest"
  if (sample_id %in% ht5_basal_S3) {
    
    value_1 <- paste("5Ht_basal_clust_S3-",length(ht5_basal_S3))
    value_3 <- paste("5Ht_basal_clust_S3-",length(ht5_basal_S3))
  }else if (sample_id %in% ho6_basal_S3) {
    
    value_2 <- paste("6ho_basal_clust_S3-",length(ho6_basal_S3))  
    value_3 <- paste("6ho_basal_clust_S3-",length(ho6_basal_S3))  
    
  }
  
  # Append the value to the array
  values_array_basal_5ht_S3[i] <- value_1
  values_array_basal_6ho_S3[i] <- value_2
  values_array_basal_both_S3[i] <- value_3
}

#CLUSTER Unique
values_array_basal_5ht_U1 <- vector(mode = "character", length = nrow(metadata_No_ptprc_no_adgre1))

# Iterate over each row in the dataframe meta
for (i in 1:nrow(metadata_No_ptprc_no_adgre1)) {
  # Get the sample_id from the current row
  sample_id <- metadata_No_ptprc_no_adgre1$cell[i]
  cell <- metadata_No_ptprc_no_adgre1$cluster1[i]
  
  # Initialize a variable to store the value
  value_1 <- "not_of_interest"  # Default value if the sample_id is not found in any list
  
  
  if (sample_id %in% ht5_basal_U1) {
    
    value_1 <- paste("5Ht_basal_clust_U1_wt-",length(ht5_basal_U1))
  }
  # Append the value to the array
  values_array_basal_5ht_U1[i] <- value_1
  
}

ht5_color <- 'red'
ho6_color <- 'blue'
others_color <- 'grey'





#S1
#5ht
Idents(obj_No_ptprc_no_adgre1) <- values_array_basal_5ht_S1
p<-DimPlot(obj_No_ptprc_no_adgre1, reduction = "UMAP", group.by = "ident",order = c('5Ht_basal_clust_S1- 489','6ho_basal_clust_S2- 156','not_of_interest'), cols = c(others_color,ht5_color)) +  theme(legend.position = 'none' ,legend.text = element_text(size=8))+ theme(panel.border = element_rect(colour = "black", fill=NA, size=1),axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20),axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))+xlim(-5, 5) + ylim(-12, -7)
#print(p)
ggsave(file.path(figure.output.dir,'Fig3_UMAP/S1_left_1_basal.png'), plot = p, width = 8, height = 6)

#6ho
Idents(obj_No_ptprc_no_adgre1) <- values_array_basal_6ho_S1
p<-DimPlot(obj_No_ptprc_no_adgre1, reduction = "UMAP",group.by = "ident",order = c('6ho_basal_clust_S1- 301','not_of_interest') ,cols = c(others_color,ho6_color)) + theme(legend.text = element_text(size=8))+theme(legend.position = 'none' ,legend.text = element_text(size=8))+ theme(panel.border = element_rect(colour = "black", fill=NA, size=1),axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20),axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))+xlim(-5, 5) + ylim(-12, -7)
#print(p)
ggsave(file.path(figure.output.dir,'Fig3_UMAP/S1_cent_1_basal.png'), plot = p, width = 8, height = 6)

#Both
Idents(obj_No_ptprc_no_adgre1) <- values_array_basal_both_S1
p<-DimPlot(obj_No_ptprc_no_adgre1, reduction = "UMAP",group.by = "ident",order = c('5Ht_basal_clust_S1- 489','6ho_basal_clust_S1- 301','not_of_interest'),cols = c(others_color,ho6_color,ht5_color))  + theme(legend.text = element_text(size=8))+theme(legend.position = 'none' ,legend.text = element_text(size=8))+ theme(panel.border = element_rect(colour = "black", fill=NA, size=1),axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20),axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))+xlim(-5, 5) + ylim(-12, -7)
#print(p)
ggsave(file.path(figure.output.dir,'Fig3_UMAP/S1_right_1_basal.png'), plot = p, width = 8, height = 6)

#S2 
#5Ht
Idents(obj_No_ptprc_no_adgre1) <- values_array_basal_5ht_S2
p<-DimPlot(obj_No_ptprc_no_adgre1, reduction = "UMAP", group.by = "ident",order = c('5Ht_basal_clust_S2_1_2_3- 409','not_of_interest'),cols = c(others_color,ht5_color))  + theme(legend.text = element_text(size=8))+theme(legend.position = 'none' ,legend.text = element_text(size=8))+ theme(panel.border = element_rect(colour = "black", fill=NA, size=1),axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20),axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))+xlim(-5, 5) + ylim(-12, -7)
#print(p)
ggsave(file.path(figure.output.dir,'Fig3_UMAP/S2_left_1_basal.png'), plot = p, width = 8, height = 6)

#6Ho
Idents(obj_No_ptprc_no_adgre1) <- values_array_basal_6ho_S2
p<-DimPlot(obj_No_ptprc_no_adgre1, reduction = "UMAP",group.by = "ident",order = c('6ho_basal_clust_S2- 156','not_of_interest'),cols = c(others_color,ho6_color))  + theme(legend.text = element_text(size=8))+theme(legend.position = 'none' ,legend.text = element_text(size=8))+ theme(panel.border = element_rect(colour = "black", fill=NA, size=1),axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20),axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))+xlim(-5, 5) + ylim(-12, -7)
#print(p)
ggsave(file.path(figure.output.dir,'Fig3_UMAP/S2_cent_1_basal.png'), plot = p, width = 8, height = 6)

#Both
Idents(obj_No_ptprc_no_adgre1) <- values_array_basal_both_S2
p<-DimPlot(obj_No_ptprc_no_adgre1, reduction = "UMAP",group.by = "ident",order = c('5Ht_basal_clust_S2_1_2_3- 409','6ho_basal_clust_S2- 156','not_of_interest'),cols = c(others_color,ho6_color,ht5_color))  + theme(legend.text = element_text(size=8))+theme(legend.position = 'none' ,legend.text = element_text(size=8))+ theme(panel.border = element_rect(colour = "black", fill=NA, size=1),axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20),axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))+xlim(-5, 5) + ylim(-12, -7)
#print(p)
ggsave(file.path(figure.output.dir,'Fig3_UMAP/S2_right_1_basal.png'), plot = p, width = 8, height = 6)

#S3
#5Ht
Idents(obj_No_ptprc_no_adgre1) <- values_array_basal_5ht_S3
p<-DimPlot(obj_No_ptprc_no_adgre1, reduction = "UMAP", group.by = "ident",order = c('5Ht_basal_clust_S3- 45','not_of_interest'), cols = c(others_color,ht5_color))  + theme(legend.text = element_text(size=8))+theme(legend.position = 'none' ,legend.text = element_text(size=8))+ theme(panel.border = element_rect(colour = "black", fill=NA, size=1),axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20),axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))+xlim(-5, 5) + ylim(-12, -7)
#print(p)
ggsave(file.path(figure.output.dir,'Fig3_UMAP/S3_left_1_basal.png'), plot = p, width = 8, height = 6)

#6Ho
Idents(obj_No_ptprc_no_adgre1) <- values_array_basal_6ho_S3
p<-DimPlot(obj_No_ptprc_no_adgre1, reduction = "UMAP",group.by = "ident",order = c('6ho_basal_clust_S3- 59','not_of_interest'),cols = c(others_color,ho6_color))  + theme(legend.text = element_text(size=8))+theme(legend.position = 'none' ,legend.text = element_text(size=8))+ theme(panel.border = element_rect(colour = "black", fill=NA, size=1),axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20),axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))+xlim(-5, 5) + ylim(-12, -7)
#print(p)
ggsave(file.path(figure.output.dir,'Fig3_UMAP/S3_cent_1_basal.png'), plot = p, width = 8, height = 6)

#Both
Idents(obj_No_ptprc_no_adgre1) <- values_array_basal_both_S3
p<-DimPlot(obj_No_ptprc_no_adgre1, reduction = "UMAP",group.by = "ident",order = c('5Ht_basal_clust_S3- 45','6ho_basal_clust_S3- 59','not_of_interest'),cols = c(others_color,ho6_color,ht5_color)) + theme(legend.text = element_text(size=8))+theme(legend.position = 'none' ,legend.text = element_text(size=8))+ theme(panel.border = element_rect(colour = "black", fill=NA, size=1),axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20),axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))+xlim(-5, 5) + ylim(-12, -7)
#print(p)
ggsave(file.path(figure.output.dir,'Fig3_UMAP/S3_right_1_basal.png'), plot = p, width = 8, height = 6)

#Unique - U1_wt
Idents(obj_No_ptprc_no_adgre1) <- values_array_basal_5ht_U1
p<-DimPlot(obj_No_ptprc_no_adgre1, reduction = "UMAP", group.by = "ident",order = c('5Ht_basal_clust_U1_wt- 29','not_of_interest'), cols = c(others_color,ht5_color))+ theme(legend.text = element_text(size=8))+theme(legend.position = 'none' ,legend.text = element_text(size=8))+ theme(panel.border = element_rect(colour = "black", fill=NA, size=1),axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20),axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))+xlim(-5, 5) + ylim(-12, -7)
print(p)
ggsave(file.path(figure.output.dir,'Fig3_UMAP/U1_wt_basal.png'), plot = p, width = 8, height = 6)
#dev.off()

##############################################################################
##############################################################################
#FEATURE PLOTS - Fig 6
##############################################################################
##############################################################################
# Create a named list of regulon target‐gene vectors
Inflammatory <- list(
  Fli1 = c(
    "Sgms1", "Fcer1g", "C5ar1", "Pld4", "Slc11a1", "Ccl24", "Stab1", "Rbpj",
    "Ctsc", "C1qa", "Adam8", "Ednrb", "Nrros", "Ccl7", "Alox5ap", "Ccl9",
    "Nfkbid", "Itgam", "Cxcl2", "Prcp", "Clec10a", "Pf4", "Syk", "Cd163",
    "Ly86", "Ccl6", "Cybb", "Csf1r", "Ptafr", "Sirpa", "Adam17", "Gpx1"
  )
)

for (gene in names(Inflammatory)[1]) {
  dir_name <- paste(figure.output.dir,'FeaturePlots/Inflammatory/',gene,'/',sep = '')
  dir.create(dir_name, showWarnings = FALSE,recursive = TRUE)
  genes_to_plot = c(Inflammatory[[gene]],gene)
  genes_to_plot <- unique(genes_to_plot)
  for (g in  genes_to_plot) {
    p<- FeaturePlot(obj_No_ptprc_no_adgre1, features = g,order = TRUE,pt.size = 1.5)+ theme(panel.border = element_rect(colour = "black", fill=NA, size=1.5),axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20),axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), plot.title = element_text(size = 35, face = "plain"))+
      scale_colour_gradientn(colours = c("grey", "darkblue"),limits = c(0, 2.5)) +xlim(-5, 5) + ylim(-12, -7)  
    save_path <- paste(dir_name,g,'.png',sep = '')
    ggsave(save_path, plot = p, width = 8, height = 6)
  }

}

