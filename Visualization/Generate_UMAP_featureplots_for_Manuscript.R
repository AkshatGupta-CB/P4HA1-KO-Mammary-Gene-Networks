#Code for generating UMAP and feature plots for the manuscript.

library(Seurat)
library(patchwork)
library(ggplot2)
library(readxl)

#Load the metadata
rm(list=ls())
metaData <- read.csv(file = "seuratMD.csv",header = TRUE)
counts <- read.csv(file = "seuratCounts.csv",header = TRUE)


#5Ht basal
ht5_basal_S1 <- readLines("New_cell_lists_clusters/Basal/5ht/5ht_basal_S1.txt")
ht5_basal_S2 <- readLines("New_cell_lists_clusters/Basal/5ht/5ht_basal_S2.txt")
ht5_basal_S3 <- readLines("New_cell_lists_clusters/Basal/5ht/5ht_basal_S3.txt")
ht5_basal_U1 <- readLines("New_cell_lists_clusters/Basal/5ht/5ht_basal_U1_wt.txt")
ht5_basal_S4 <- readLines("New_cell_lists_clusters/Basal/5ht/5ht_basal_S4.txt")

#6Ho luminal
ho6_basal_S1 <- readLines("New_cell_lists_clusters/Basal/6ho/6ho_basal_S1.txt")
ho6_basal_S2 <- readLines("New_cell_lists_clusters/Basal/6ho/6ho_basal_S2.txt")
ho6_basal_S3 <- readLines("New_cell_lists_clusters/Basal/6ho/6ho_basal_S3.txt")
ho6_basal_U1 <- readLines("New_cell_lists_clusters/Basal/6ho/6ho_basal_U1_ko.txt")
ho6_basal_S4 <- readLines("New_cell_lists_clusters/Basal/6ho/6ho_basal_S4.txt")

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

#######################################
#UPDATE: SUBSET FOR BASAL CELLS AND GENERATE PLOTS W.R.T. BASAL CELLS ONLY
######################################
obj <- subset(obj,subset = cluster1 == 'Mammary epithelial cells-Basal')
metaData <- obj@meta.data



values_array <- vector(mode = "character", length = nrow(metaData))
cells_array <- vector(mode = "character", length = nrow(metaData))
cell_cycle_array <- vector(mode = "character", length = nrow(metaData))

# Iterate over each row in the dataframe meta
for (i in 1:nrow(metaData)) {
  # Get the sample_id from the current row
  sample_id <- metaData$cell[i]
  cell <- metaData$cluster1[i]
  
  # Initialize a variable to store the value
  value <- "not_of_interest"  # Default value if the sample_id is not found in any list
  if (cell == "Mammary epithelial cells-Basal") {
    cells_array[i] <- cell
  }else if (cell == "Mammary epithelial cells-Luminal") {
    cells_array[i] <- cell
  }else if (cell == "Monocytes-macrophages") {
    cells_array[i] <- cell
  }else{
    cells_array[i] <- "not_of_interest"
  }
  
}#


Idents(obj) <- as.character(metaData$sample)
p <- DimPlot(obj, reduction = "UMAP",group.by = "ident") + theme(legend.position = 'none' ,legend.text = element_text(size = 8))+ theme(panel.border = element_rect(colour = "black", fill=NA, size=1),axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20),axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))
ggsave('../MAnuscript/Figures/Integrated_UMAP.png', plot = p, width = 10, height = 8)

p <- DimPlot(obj, reduction = "UMAP",group.by = "ident",split.by = "ident") + theme(legend.position = 'none' ,legend.text = element_text(size = 8))+ theme(panel.border = element_rect(colour = "black", fill=NA, size=1),axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20),axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))
ggsave('../MAnuscript/Figures/Split_Integrated_UMAP.png', plot = p, width = 16, height = 8)

Idents(obj) <- cells_array#as.character(metaData$cluster1)
p<- DimPlot(obj, reduction = "UMAP",group.by = "ident")+ theme(legend.position = 'none' ,legend.text = element_text(size = 8))+ theme(panel.border = element_rect(colour = "black", fill=NA, size=1),axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20),axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))
ggsave('../MAnuscript/Figures/Integrated_Annotated_CELL_UMAP.png', plot = p, width = 10, height = 8)

#Idents(obj) <- values_array
#DimPlot(obj, reduction = "UMAP",group.by = "ident")




#CLUSTER 1
values_array_basal_5ht_S1 <- vector(mode = "character", length = nrow(metaData))
values_array_basal_6ho_S1 <- vector(mode = "character", length = nrow(metaData))
values_array_basal_both_S1 <- vector(mode = "character", length = nrow(metaData))
# Iterate over each row in the dataframe meta
for (i in 1:nrow(metaData)) {
  # Get the sample_id from the current row
  sample_id <- metaData$cell[i]
  
  # Initialize a variable to store the value
  value_1 <- "not_of_interest"  # Default value if the sample_id is not found in any list
  value_2 <- "not_of_interest" 
  value_3 <- "not_of_interest"
  if (sample_id %in% ht5_basal_S1) {
    
    value_1 <- paste("5Ht_basal_clust_S1-",length(ht5_basal_S1))#"5Ht_luminal_clust_5"  
    value_3 <- paste("5Ht_basal_clust_S1-",length(ht5_basal_S1))#"5Ht_luminal_clust_5"
  }else if (sample_id %in% ho6_basal_S1) {
    
    value_2 <- paste("6ho_basal_clust_S1-",length(ho6_basal_S1))#"6ho_luminal_clust_5"  
    value_3 <- paste("6ho_basal_clust_S1-",length(ho6_basal_S1))  
    
  }
  
  # Append the value to the array
  values_array_basal_5ht_S1[i] <- value_1
  values_array_basal_6ho_S1[i] <- value_2
  values_array_basal_both_S1[i] <- value_3
}

#CLUSTER 2
values_array_basal_5ht_S2 <- vector(mode = "character", length = nrow(metaData))
values_array_basal_6ho_S2 <- vector(mode = "character", length = nrow(metaData))
values_array_basal_both_S2 <- vector(mode = "character", length = nrow(metaData))
# Iterate over each row in the dataframe meta
for (i in 1:nrow(metaData)) {
  # Get the sample_id from the current row
  sample_id <- metaData$cell[i]
  
  # Initialize a variable to store the value
  value_1 <- "not_of_interest"  # Default value if the sample_id is not found in any list
  value_2 <- "not_of_interest" 
  value_3 <- "not_of_interest"
  if (sample_id %in% ht5_basal_S2) {
    
    value_1 <- paste("5Ht_basal_clust_S2-",length(ht5_basal_S2))#"5Ht_luminal_clust_5"  
    value_3 <- paste("5Ht_basal_clust_S2-",length(ht5_basal_S2))#"5Ht_luminal_clust_5"
  }else if (sample_id %in% ho6_basal_S2) {
    
    value_2 <- paste("6ho_basal_clust_S2-",length(ho6_basal_S2))#"6ho_luminal_clust_5"  
    value_3 <- paste("6ho_basal_clust_S2-",length(ho6_basal_S2))  
    
  }
  
  # Append the value to the array
  values_array_basal_5ht_S2[i] <- value_1
  values_array_basal_6ho_S2[i] <- value_2
  values_array_basal_both_S2[i] <- value_3
}

#CLUSTER 3
values_array_basal_5ht_S3 <- vector(mode = "character", length = nrow(metaData))
values_array_basal_6ho_S3 <- vector(mode = "character", length = nrow(metaData))
values_array_basal_both_S3 <- vector(mode = "character", length = nrow(metaData))
# Iterate over each row in the dataframe meta
for (i in 1:nrow(metaData)) {
  # Get the sample_id from the current row
  sample_id <- metaData$cell[i]
  
  # Initialize a variable to store the value
  value_1 <- "not_of_interest"  # Default value if the sample_id is not found in any list
  value_2 <- "not_of_interest" 
  value_3 <- "not_of_interest"
  if (sample_id %in% ht5_basal_S3) {
    
    value_1 <- paste("5Ht_basal_clust_S3-",length(ht5_basal_S3))#"5Ht_luminal_clust_5"  
    value_3 <- paste("5Ht_basal_clust_S3-",length(ht5_basal_S3))#"5Ht_luminal_clust_5"
  }else if (sample_id %in% ho6_basal_S3) {
    
    value_2 <- paste("6ho_basal_clust_S3-",length(ho6_basal_S3))#"6ho_luminal_clust_5"  
    value_3 <- paste("6ho_basal_clust_S3-",length(ho6_basal_S3))  
    
  }
  
  # Append the value to the array
  values_array_basal_5ht_S3[i] <- value_1
  values_array_basal_6ho_S3[i] <- value_2
  values_array_basal_both_S3[i] <- value_3
}

#CLUSTER 4
values_array_basal_5ht_S4 <- vector(mode = "character", length = nrow(metaData))
values_array_basal_6ho_S4 <- vector(mode = "character", length = nrow(metaData))
values_array_basal_both_S4 <- vector(mode = "character", length = nrow(metaData))
# Iterate over each row in the dataframe meta
for (i in 1:nrow(metaData)) {
  # Get the sample_id from the current row
  sample_id <- metaData$cell[i]
  
  # Initialize a variable to store the value
  value_1 <- "not_of_interest"  # Default value if the sample_id is not found in any list
  value_2 <- "not_of_interest" 
  value_3 <- "not_of_interest"
  if (sample_id %in% ht5_basal_S4) {
    
    value_1 <- paste("5Ht_basal_clust_S4-",length(ht5_basal_S4))#"5Ht_luminal_clust_5"  
    value_3 <- paste("5Ht_basal_clust_S4-",length(ht5_basal_S4))#"5Ht_luminal_clust_5"
  }else if (sample_id %in% ho6_basal_S4) {
    
    value_2 <- paste("6ho_basal_clust_S4-",length(ho6_basal_S4))#"6ho_luminal_clust_5"  
    value_3 <- paste("6ho_basal_clust_S4-",length(ho6_basal_S4))  
    
  }
  
  # Append the value to the array
  values_array_basal_5ht_S4[i] <- value_1
  values_array_basal_6ho_S4[i] <- value_2
  values_array_basal_both_S4[i] <- value_3
}

#CLUSTER Unique
values_array_basal_5ht_U1 <- vector(mode = "character", length = nrow(metaData))
values_array_basal_6ho_U1 <- vector(mode = "character", length = nrow(metaData))
# Iterate over each row in the dataframe meta
for (i in 1:nrow(metaData)) {
  # Get the sample_id from the current row
  sample_id <- metaData$cell[i]
  cell <- metaData$cluster1[i]
  
  # Initialize a variable to store the value
  value_1 <- "not_of_interest"  # Default value if the sample_id is not found in any list
  value_2 <- "not_of_interest" 
  value_3 <- "not_of_interest"
  if (sample_id %in% ht5_basal_U1) {
    
    value_1 <- paste("5Ht_basal_clust_U1_wt-",length(ht5_basal_U1))#"5Ht_luminal_clust_5"  
  }else if (sample_id %in% ho6_basal_U1) {
    
    value_2 <- paste("6ho_basal_clust_U1_ko-",length(ho6_basal_U1))#"6ho_luminal_clust_5"  
  }
  
  # Append the value to the array
  values_array_basal_5ht_U1[i] <- value_1
  values_array_basal_6ho_U1[i] <- value_2
  
}


ht5_color <- 'red'
ho6_color <- 'blue'
others_color <- 'grey'

# order subset to be plotted last (on top)
# levels(Idents(obj))
# "5Ht_basal_clust_S1- 490" "not_of_interest"     
# We swap order, also need to swap cols
pdf("UMAP_plots_v6.pdf",paper = "USr")
Idents(obj) <- values_array_basal_5ht_S1
p<-DimPlot(obj, reduction = "UMAP", group.by = "ident", order = c("5Ht_basal_clust_S1- 490","not_of_interest"), cols = c(others_color,ht5_color)) + theme(legend.position = 'none' ,legend.text = element_text(size = 8))+ theme(panel.border = element_rect(colour = "black", fill=NA, size=1),axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20),axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))+xlim(-5, 5) + ylim(-12, -7) 
#print(p)
ggsave('../MAnuscript/Figures/S2_left_1_basal.png', plot = p, width = 8, height = 6)


Idents(obj) <- values_array_basal_6ho_S1
p<-DimPlot(obj, reduction = "UMAP",group.by = "ident",order = c("6ho_basal_clust_S1- 208","not_of_interest"),cols = c(others_color,ho6_color)) + theme(legend.position = 'none',legend.text = element_text(size=8))+ theme(panel.border = element_rect(colour = "black", fill=NA, size=1),axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20),axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))+xlim(-5, 5) + ylim(-12, -7)
#print(p)
ggsave('../MAnuscript/Figures/S2_cent_1_basal.png', plot = p, width = 8, height = 6)


Idents(obj) <- values_array_basal_both_S1
p<-DimPlot(obj, reduction = "UMAP",group.by = "ident",order = c("5Ht_basal_clust_S1- 490", "6ho_basal_clust_S1- 208","not_of_interest") ,cols = c(others_color,ho6_color, ht5_color)) + theme(legend.position = 'none' ,legend.text = element_text(size=8))+ theme(panel.border = element_rect(colour = "black", fill=NA, size=1),axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20),axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))+xlim(-5, 5) + ylim(-12, -7)
#print(p)
ggsave('../MAnuscript/Figures/S2_right_1_basal.png', plot = p, width = 8, height = 6)


Idents(obj) <- values_array_basal_5ht_S2
p<-DimPlot(obj, reduction = "UMAP",group.by = "ident",order = c("5Ht_basal_clust_S2- 30", "not_of_interest") ,cols = c(others_color,ht5_color)) + theme(legend.position = 'none' ,legend.text = element_text(size=8))+ theme(panel.border = element_rect(colour = "black", fill=NA, size=1),axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20),axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))+xlim(-5, 5) + ylim(-12, -7)
#print(p)
ggsave('../MAnuscript/Figures/S2_left_2_basal.png', plot = p, width = 8, height = 6)

#DimPlot(obj, reduction = "UMAP",group.by = "ident",order = ht5_basal_S2 ,cols = c(ht5_color,others_color),split.by = "ident") + theme(legend.text = element_text(size=5))

Idents(obj) <- values_array_basal_6ho_S2
p<-DimPlot(obj, reduction = "UMAP",group.by = "ident",order = c("6ho_basal_clust_S2- 23","not_of_interest"),cols = c(others_color,ho6_color)) + theme(legend.position = 'none' ,legend.text = element_text(size=8))+ theme(panel.border = element_rect(colour = "black", fill=NA, size=1),axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20),axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))+xlim(-5, 5) + ylim(-12, -7)
#print(p)
ggsave('../MAnuscript/Figures/S2_cent_2_basal.png', plot = p, width = 8, height = 6)

#DimPlot(obj, reduction = "UMAP",group.by = "ident",order = ho6_basal_S2 ,cols = c(ho6_color,others_color),split.by = "ident") + theme(legend.text = element_text(size=5))

Idents(obj) <- values_array_basal_both_S2
p<-DimPlot(obj, reduction = "UMAP",group.by = "ident",order = c("5Ht_basal_clust_S2- 30","6ho_basal_clust_S2- 23","not_of_interest") ,cols = c(others_color,ho6_color,ht5_color)) + theme(legend.position = 'none' ,legend.text = element_text(size=8))+ theme(panel.border = element_rect(colour = "black", fill=NA, size=1),axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20),axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))+xlim(-5, 5) + ylim(-12, -7)
#print(p)
ggsave('../MAnuscript/Figures/S2_right_2_basal.png', plot = p, width = 8, height = 6)


Idents(obj) <- values_array_basal_5ht_S3
p<-DimPlot(obj, reduction = "UMAP",group.by = "ident",order = c("5Ht_basal_clust_S3- 163", "not_of_interest") ,cols = c(others_color,ht5_color)) + theme(legend.position = 'none' ,legend.text = element_text(size=8))+ theme(panel.border = element_rect(colour = "black", fill=NA, size=1),axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20),axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))+xlim(-5, 5) + ylim(-12, -7)
#print(p)
ggsave('../MAnuscript/Figures/S2_left_3_basal.png', plot = p, width = 8, height = 6)

#DimPlot(obj, reduction = "UMAP",group.by = "ident",order = ht5_basal_S3 ,cols = c(ht5_color,others_color),split.by = "ident") + theme(legend.text = element_text(size=5))

Idents(obj) <- values_array_basal_6ho_S3
p<-DimPlot(obj, reduction = "UMAP",group.by = "ident",order = c("6ho_basal_clust_S3- 60","not_of_interest"),cols = c(others_color,ho6_color)) + theme(legend.position = 'none' ,legend.text = element_text(size=8))+ theme(panel.border = element_rect(colour = "black", fill=NA, size=1),axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20),axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))+xlim(-5, 5) + ylim(-12, -7)
#print(p)
ggsave('../MAnuscript/Figures/S2_cent_3_basal.png', plot = p, width = 8, height = 6)


# DimPlot(obj, reduction = "UMAP",group.by = "ident",order = ho6_basal_S3 ,cols = c(ho6_color,others_color),split.by = "ident") + theme(legend.text = element_text(size=5))

Idents(obj) <- values_array_basal_both_S3
p<-DimPlot(obj, reduction = "UMAP",group.by = "ident",order = c("5Ht_basal_clust_S3- 163","6ho_basal_clust_S3- 60","not_of_interest") ,cols = c(others_color,ho6_color,ht5_color)) + theme(legend.position = 'none' ,legend.text = element_text(size=8))+ theme(panel.border = element_rect(colour = "black", fill=NA, size=1),axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20),axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))+xlim(-5, 5) + ylim(-12, -7)
#print(p)
ggsave('../MAnuscript/Figures/S2_right_3_basal.png', plot = p, width = 8, height = 6)


Idents(obj) <- values_array_basal_5ht_S4
p<-DimPlot(obj, reduction = "UMAP",group.by = "ident",order = c("5Ht_basal_clust_S4- 357","not_of_interest") ,cols = c(others_color,ht5_color)) + theme(legend.position = 'none' ,legend.text = element_text(size=5))+ theme(panel.border = element_rect(colour = "black", fill=NA, size=1),axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20),axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))+xlim(-5, 5) + ylim(-12, -7)
#print(p)
ggsave('../MAnuscript/Figures/S2_left_4_basal.png', plot = p, width = 8, height = 6)

#DimPlot(obj, reduction = "UMAP",group.by = "ident",order = ht5_basal_S4 ,cols = c(ht5_color,others_color),split.by = "ident") + theme(legend.text = element_text(size=5))

Idents(obj) <- values_array_basal_6ho_S4
p<-DimPlot(obj, reduction = "UMAP",group.by = "ident",order = c("6ho_basal_clust_S4- 173","not_of_interest"),cols = c(others_color,ho6_color)) + theme(legend.position = 'none' ,legend.text = element_text(size=5))+ theme(panel.border = element_rect(colour = "black", fill=NA, size=1),axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20),axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))+xlim(-5, 5) + ylim(-12, -7)
#print(p)
ggsave('../MAnuscript/Figures/S2_cent_4_basal.png', plot = p, width = 8, height = 6)

#DimPlot(obj, reduction = "UMAP",group.by = "ident",order = ho6_basal_S4 ,cols = c(ho6_color,others_color),split.by = "ident") + theme(legend.text = element_text(size=5))

Idents(obj) <- values_array_basal_both_S4
p<-DimPlot(obj, reduction = "UMAP",group.by = "ident",order = c("5Ht_basal_clust_S4- 357","6ho_basal_clust_S4- 173","not_of_interest") ,cols = c(others_color,ho6_color,ht5_color)) + theme(legend.position = 'none' ,legend.text = element_text(size=5))+ theme(panel.border = element_rect(colour = "black", fill=NA, size=1),axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20),axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))+xlim(-5, 5) + ylim(-12, -7)
#print(p)
ggsave('../MAnuscript/Figures/S2_right_4_basal.png', plot = p, width = 8, height = 6)



Idents(obj) <- values_array_basal_5ht_U1
p<-DimPlot(obj, reduction = "UMAP",group.by = "ident",order = c("5Ht_basal_clust_U1_wt- 129","not_of_interest") ,cols = c(others_color,ht5_color)) + theme(legend.position = 'none' ,legend.text = element_text(size=5))+ theme(panel.border = element_rect(colour = "black", fill=NA, size=1),axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20),axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))+xlim(-5, 5) + ylim(-12, -7)
#print(p)
ggsave('../MAnuscript/Figures/S2_left_5_basal.png', plot = p, width = 8, height = 6)

#DimPlot(obj, reduction = "UMAP",group.by = "ident",order = ht5_basal_U1 ,cols = c(ht5_color,others_color),split.by = "ident") + theme(legend.text = element_text(size=5))

Idents(obj) <- values_array_basal_6ho_U1
p<-DimPlot(obj, reduction = "UMAP",group.by = "ident",order = c("6ho_basal_clust_U1_ko- 113","not_of_interest"),cols = c(others_color,ho6_color)) + theme(legend.position = 'none' ,legend.text = element_text(size=5))+ theme(panel.border = element_rect(colour = "black", fill=NA, size=1),axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20),axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))+xlim(-5, 5) + ylim(-12, -7)
#print(p)
ggsave('../MAnuscript/Figures/S2_right_5_basal.png', plot = p, width = 8, height = 6)
dev.off()
#DimPlot(obj, reduction = "UMAP",group.by = "ident",order = ho6_basal_U1 ,cols = c(ho6_color,others_color),split.by = "ident") + theme(legend.text = element_text(size=5))




###############################FEATURE PLOTS FINAL#######################################
path<- "C:/Users/gupta/OneDrive/Desktop/SSAnalysis/MAnuscript/Figures/Feature_plots_final/"
genes_5ht_u1_reg_hif <- c('Arid3a','Egln2','Hmox1','Nfkb1')

reg_stem_cell = c('Ikzf1','Rel','Fli1','Mafb','Cited2','Zeb2','Nrp2','Fn1','Rbpj'
                  ,'Mef2c','Ncoa3','Nrp1','Ptprc','Sema4a','Zfp36')

fli_DETG <- c("Fli1",'Cd14','Nfkbid','Fcgr2b','Syk','Ccr2','Cybb','Csf1r','Itgb2')

gene = 'Fli1'
for (genes in fli_DETG) {
  #dir_name <- paste(path,sep = '')
  #dir.create(dir_name, showWarnings = FALSE)
  for (gene in genes) {
    p<- FeaturePlot(obj, features = gene,order = TRUE,pt.size = 1.5)+ theme(panel.border = element_rect(colour = "black", fill=NA, size=1.5),axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20),axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), plot.title = element_text(size = 35, face = "plain"))+
      scale_colour_gradientn(colours = c("grey", "darkblue"),limits = c(0, 2.5)) +xlim(-5, 5) + ylim(-12, -7) 
    save_path <- paste(path,gene,'.png',sep = '')
    ggsave(save_path, plot = p, width = 8, height = 6)
  }
}
