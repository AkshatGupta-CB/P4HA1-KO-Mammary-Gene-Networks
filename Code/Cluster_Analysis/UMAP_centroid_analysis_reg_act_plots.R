library(Seurat)
library(patchwork)
library(ggplot2)
library(readxl)
library(dplyr)
library(umap)
#Load the AUC Matrix
AUC_mtx_5ht = read.csv('Results/Results_no_ptprc_adgre1/Cluster_Analysis_avg_link_gene_cell_corr/AUC_mtx_with_clusters/5ht_basal_AUC_with_clusts.csv')

AUC_mtx_6ho = read.csv('Results/Results_no_ptprc_adgre1/Cluster_Analysis_avg_link_gene_cell_corr/AUC_mtx_with_clusters/6ho_basal_AUC_with_clusts.csv')

#Add cluster information
clust_5ht = paste0('5Ht_',AUC_mtx_5ht$Cluster)
clust_6ho = paste0('6Ho_',AUC_mtx_6ho$Cluster)

AUC_mtx_5ht$Cluster = clust_5ht
AUC_mtx_6ho$Cluster = clust_6ho


#Combine the matrices
combined_AUC_activity = rbind(AUC_mtx_5ht,AUC_mtx_6ho) 
metadata = combined_AUC_activity[c('Cell','Cluster')]
row.names(combined_AUC_activity)= combined_AUC_activity$Cell
combined_AUC_activity$Cluster=NULL
combined_AUC_activity$Cell=NULL
group <- ifelse(grepl("5Ht", metadata$Cluster), "5Ht", "6Ho")

metadata$group <- group

#Scale the matrix
scaled_matrix = scale(combined_AUC_activity)

# Perform PCA
pca_result <- prcomp(scaled_matrix, center = TRUE, scale. = TRUE)

# Extract the top 15 PCs
pca_coords <- as.data.frame(pca_result$x[, 1:15])  # Take first 15 PCs
pca_coords$Group <- metadata$Cluster  # Add the cluster information 

# Run UMAP on the PCA-reduced data
set.seed(2)                # R-level RNG
umap_config$random_state <- 2  # C++ side RNG
umap_config <- umap.defaults
umap_config$n_neighbors <- 30  
umap_config$min_dist <- 0.3    
umap_config$n_components <- 2  # 2D UMAP
umap_config$metric <- "euclidean"

umap_result <- umap(pca_coords[,1:15], config = umap_config)

# Convert UMAP results to dataframe
umap_df <- as.data.frame(umap_result$layout)
colnames(umap_df) <- c("UMAP_1", "UMAP_2")
umap_df$Group <- metadata$group   # Add group information back


# Plot UMAP for each group
ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = Group)) +
  geom_point(alpha = 0.6) +
  #facet_wrap(~Group) +  # Create separate plots for each group
  theme_minimal() +
  ggtitle("UMAP of Scaled Regulon Activities by Group")


umap_result <- umap(pca_coords[, 1:15], config = umap_config)  # Exclude cluster column

# Convert UMAP results to dataframe
umap_df <- as.data.frame(umap_result$layout)
colnames(umap_df) <- c("UMAP_1", "UMAP_2")
umap_df$Cluster <- pca_coords$Group  # Add cluster information back

# Define unique subclusters
similar_clusters <- c('S1','S2','S3')

# Loop through subclusters and create separate plots
for (subcluster in similar_clusters) {
  # Create a new column marking the selected subcluster pair
  group_vec <- c()
  for (clust in umap_df$Cluster) {
    if (grepl(subcluster, clust)) {
      group_vec <- c(group_vec, clust)
    } else {
      group_vec <- c(group_vec, "Others")
    }
  }
  
  ht5_clust <- paste0("5Ht_", subcluster)  # Generate full cluster names
  ho6_clust <- paste0("6Ho_", subcluster)
  
  umap_df$Group <- factor(group_vec, levels = c("Others", ht5_clust, ho6_clust))  # Order groups
  
  # Dynamically map colors
  color_mapping <- setNames(c("grey", "red", "blue"), c("Others", ht5_clust, ho6_clust))
  
  p <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = Group)) +
    # plot grey “Others” first
    geom_point(data = subset(umap_df, Group == "Others"), alpha = 1) +
    # then highlight the subcluster points
    geom_point(data = subset(umap_df, Group != "Others"), alpha = 1) +
    # use your custom colour mapping
    scale_color_manual(values = color_mapping) +
    # minimal base theme + customized elements
    theme_minimal() +
    theme(
      panel.grid       = element_blank(),                        # remove all grid lines
      panel.border   = element_rect(colour = "black", fill = NA, size = 1.5),
      axis.ticks        = element_line(color = "black", size = 1),
      axis.ticks.length = unit(0.1, "cm"),
      axis.text.x    = element_text(size = 20),
      axis.text.y    = element_text(size = 20),
      axis.title.x   = element_text(size = 20),
      axis.title.y   = element_text(size = 20),
      plot.title     = element_text(size = 35, face = "plain"),
      legend.position  = "none"        # remove the legend
    ) +xlim(-6, 4) + ylim(-5, 5) 
  
  # Display the plot
  print(p)
  ggsave(paste0('../MAnuscript/PLOS_Genetics_sub_2/Figures/Regulon_activityv3_plot_',subcluster,'.png'), plot = p, width = 8, height = 6)
}

# Define unique subclusters
unique_clusters <- c('U1')

# Loop through subclusters and create separate plots
for (subcluster in unique_clusters) {
  # Create a new column marking the selected subcluster pair
  group_vec <- c()
  for (clust in umap_df$Cluster) {
    if (grepl(subcluster, clust)) {
      group_vec <- c(group_vec, clust)
    } else {
      group_vec <- c(group_vec, "Others")
    }
  }
  
  # Dynamically map colors
  if (subcluster == 'U1') {
    c = '5Ht_U1'
    umap_df$Group <- factor(group_vec, levels = c("Others", c))
    color_mapping <- setNames(c("grey", "red"), c("Others",c))
  }
  
  p <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = Group)) +
    # plot grey “Others” first
    geom_point(data = subset(umap_df, Group == "Others"), alpha = 1) +
    # then highlight the subcluster points
    geom_point(data = subset(umap_df, Group != "Others"), alpha = 1) +
    # use your custom colour mapping
    scale_color_manual(values = color_mapping) +
    # minimal base theme + customized elements
    theme_minimal() +
    theme(
      panel.grid       = element_blank(),                        # remove all grid lines
      panel.border   = element_rect(colour = "black", fill = NA, size = 1.5),
      axis.ticks        = element_line(color = "black", size = 1),
      axis.ticks.length = unit(0.1, "cm"),
      axis.text.x    = element_text(size = 20),
      axis.text.y    = element_text(size = 20),
      axis.title.x   = element_text(size = 20),
      axis.title.y   = element_text(size = 20),
      plot.title     = element_text(size = 35, face = "plain"),
      legend.position  = "none"        # remove the legend
    ) +xlim(-6, 4) + ylim(-5, 5) 
  
  # Display the plot
  print(p)
  ggsave(paste0('../MAnuscript/PLOS_Genetics_sub_2/Figures/Regulon_activityv3_plot_',subcluster,'.png'), plot = p, width = 8, height = 6)
}

################################################################################
#Centroid analysis
################################################################################
# Step 1: Calculate centroids for each subcluster
centroids <- umap_df %>%
  group_by(Cluster) %>%
  summarize(UMAP_1_mean = mean(UMAP_1), UMAP_2_mean = mean(UMAP_2))

# Step 2: Separate centroids for 5Ht and 6Ho
centroids_5Ht <- centroids %>% filter(grepl("5Ht", Cluster))
centroids_6Ho <- centroids %>% filter(grepl("6Ho", Cluster))

# Step 3: Compute pairwise Euclidean distances
distance_matrix <- outer(
  1:nrow(centroids_5Ht),
  1:nrow(centroids_6Ho),
  Vectorize(function(i, j) {
    sqrt(
      (centroids_5Ht$UMAP_1_mean[i] - centroids_6Ho$UMAP_1_mean[j])^2 +
        (centroids_5Ht$UMAP_2_mean[i] - centroids_6Ho$UMAP_2_mean[j])^2
    )
  })
)

# Step 4: Create a dataframe for the distance matrix
rownames(distance_matrix) <- centroids_5Ht$Cluster
colnames(distance_matrix) <- centroids_6Ho$Cluster
distance_df <- as.data.frame(distance_matrix)

# Step 5: Display the result
library(openxlsx)
write.xlsx(distance_df,'../MAnuscript/PLOS_Genetics_sub_2/Centroid_distances_basal_clust_reg_act.xlsx',rowNames=TRUE)



















