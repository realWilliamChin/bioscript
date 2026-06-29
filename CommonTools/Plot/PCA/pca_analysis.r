#!/usr/bin/env Rscript
# Enhanced PCA Analysis with Automatic Grouping
# Usage: Rscript pca_analysis.R <input_file> <working_dir> [options]
#Rscript pca_analysis.R test.txt D:\\shen3projectOthersList\\gongshi20250530spatial\\DZOE2024040507_刘明-罗书航老师_XMSH_202505_23937_qV1eNfNU\\1.downsample\\split_sample\\PCAexp
#source("pca_analysis.R test.txt D:\\shen3projectOthersList\\gongshi20250530spatial\\DZOE2024040507_刘明-罗书航老师_XMSH_202505_23937_qV1eNfNU\\1.downsample\\split_sample\\PCAexp")

#cd C:\Program Files\R\R-4.3.1\bin
#Rscript pca_analysis.R data.txt D:\\shen3projectOthersList\\gongshiA20250610PCA TRUE testtest.txt 3 TRUE TRUE
#Rscript --vanilla D:\shen3projectOthersList\gongshiA20250610PCA\pca_analysis.R data.txt D:\shen3projectOthersList\gongshiA20250610PCA
#Rscript D:\shen3projectOthersList\gongshiA20250610PCA\pca_analysis.R data.txt D:\shen3projectOthersList\gongshiA20250610PCA


#R4.3.1
#install.packages("optparse")#if missing this package
#install.packages("factoextra")#if missing this package
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(optparse)
  library(ggforce)
  library(factoextra)  # For clustering
  library(cluster)     # For pam clustering
})

# ------------------------- Command Line Arguments -------------------------
option_list <- list(
  make_option(c("--scale"), action="store_true", default=TRUE,
              help="Scale variables before PCA"),
  make_option(c("--output"), type="character", default="PCA_plot",
              help="Base name for output files"),
  #make_option(c("--k-groups"), type="integer", default=NULL,
  make_option(c("--k-groups"), type="integer", default=3,
              help="Number of groups to identify (NULL for automatic determination)"),
  make_option(c("--show-ellipses"), action="store_true", default=TRUE,
              help="Show confidence ellipses for groups"),
  make_option(c("--show-cluster-centers"), action="store_true", default=TRUE,
              help="Show cluster centers on plot")
)

parser <- OptionParser(usage = "%prog <input_file> <working_dir> [options]", 
                       option_list = option_list)
args <- parse_args(parser, positional_arguments = 2)

input_file <- args$args[1]
working_dir <- args$args[2]
# input_file <-"data.txt"
#input_file <-"merged_exp.txt"
#input_file <-"test.txt"
#working_dir <- "D:\\shen3projectOthersList\\gongshiA20250610PCA"

#working_dir <- "D:\\shen3projectOthersList\\gongshi20250530spatial\\DZOE2024040507_刘明-罗书航老师_XMSH_202505_23937_qV1eNfNU\\1.downsample\\split_sample"
opts <- args$options

# ------------------------- Data Loading and Preprocessing -------------------------
#setwd("D:\\shen3projectOthersList\\gongshi20250530spatial\\DZOE2024040507_刘明-罗书航老师_XMSH_202505_23937_qV1eNfNU\\1.downsample\\split_sample")
setwd(working_dir)
data <- read.table(input_file, row.names = 1, check.names = FALSE,head=TRUE)
#data <- read.delim(input_file, row.names = 1, check.names = FALSE)

convert_to_numeric <- function(df) {
  for (col in names(df)) {
    converted_col <- as.numeric(as.character(df[[col]]))
    if (any(is.na(converted_col))) {
      warning(paste("Column", col, "contains non-numeric values and will not be converted."))
    } else {
      df[[col]] <- converted_col
    }
  }
  return(df)
}
data <- convert_to_numeric(data)
data_filtered <- data[rowSums(data) > 0, ]

# ------------------------- PCA Calculation -------------------------
pca_result <- prcomp(t(data_filtered), scale. = opts$scale)
#pca_result <- prcomp(t(data_filtered), scale. = TRUE)
pca_df <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  Sample = colnames(data_filtered)
)

# ------------------------- Automatic Grouping -------------------------
# Determine optimal number of clusters if not specified
if (is.null(opts$`k-groups`)) {
#if (is.null(3)) {
  gap_stat <- clusGap(pca_result$x[,1:2], FUN = pam, K.max = 5, B = 50)
  opts$`k-groups` <- maxSE(gap_stat$Tab[, "gap"], gap_stat$Tab[, "SE.sim"])
  message("Automatically determined optimal number of groups: ", opts$`k-groups`)
}

# Perform clustering
clusters <- pam(pca_result$x[,1:2], k = opts$`k-groups`)
#clusters <- pam(pca_result$x[,1:2], k = 3)
pca_df$Cluster <- as.factor(clusters$clustering)

# Add cluster centers to the data
centers_df <- data.frame(
  PC1 = clusters$medoids[,1],
  PC2 = clusters$medoids[,2],
  #Cluster = as.factor(1:3),
  #Sample = paste("Center", 1:3)
  Cluster = as.factor(1:opts$`k-groups`),
  Sample = paste("Center", 1:opts$`k-groups`)
)

# ------------------------- Generate Plot -------------------------
p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Cluster, label = Sample)) +
  geom_point(size = 3, alpha = 0.7) +
  xlim(range(pca_df$PC1)*1.2) +  # 扩展20%范围
  ylim(range(pca_df$PC2)*1.2) +
  coord_cartesian(expand=TRUE) +   # 确保边缘不被裁剪
  geom_text_repel(size = 3, max.overlaps = 20, show.legend = FALSE) +
  labs(
    x = paste0("PC1 (", round(100 * pca_result$sdev[1]^2/sum(pca_result$sdev^2), 1), "% Variance)"),
    y = paste0("PC2 (", round(100 * pca_result$sdev[2]^2/sum(pca_result$sdev^2), 1), "% Variance)"),
    title = "PCA with Automatic Grouping",
    color = "Cluster",
  ) +
  theme_minimal()

# Add cluster centers if requested
if (opts$`show-cluster-centers`) {
  p <- p + 
    geom_point(data = centers_df, size = 5, shape = 8, show.legend = FALSE) +
    geom_text_repel(data = centers_df, aes(label = Cluster), 
                   color = "black", size = 5, fontface = "bold")
}

# Add ellipses if requested
if (opts$`show-ellipses`) {
  p <- p + geom_mark_ellipse(aes(group = Cluster, label = Cluster), 
                            alpha = 0.1, show.legend = FALSE)
}

# Save plot
plot_file <- paste0(opts$output, ".png")
ggsave(plot_file, p, width = 20, height = 20, dpi = 300)
#plot_file <- paste0("output", ".png")
#ggsave("plot_file.png", p, width = 10, height = 8, dpi = 300)
# ------------------------- Save Cluster Information -------------------------
# Combine sample coordinates with cluster assignments
cluster_output <- pca_df %>%
  mutate(
    Cluster = as.integer(Cluster),
    Distance_to_Center = sqrt((PC1 - centers_df$PC1[Cluster])^2 + 
                             (PC2 - centers_df$PC2[Cluster])^2)
  )

cluster_output <- cluster_output[order(cluster_output$Cluster, cluster_output$Sample), ]

# Save cluster information
#opts$output<-"outputtest"
cluster_file <- paste0(opts$output, "_cluster_info.csv")
#cluster_file <- paste0("klustoutput", "_cluster_info.csv")
write.csv(cluster_output, cluster_file, row.names = FALSE)

# Save cluster centers
centers_file <- paste0(opts$output, "_cluster_centers.csv")
#centers_file <- paste0("kcluster", "_cluster_centers.csv")

write.csv(centers_df, centers_file, row.names = FALSE)

# ------------------------- Output Summary -------------------------
message("\nPCA with automatic grouping completed!")
message("Number of clusters identified: ", opts$`k-groups`)
message("\nOutput files:")
message("- PCA plot: ", plot_file)
message("- Cluster information: ", cluster_file)
message("- Cluster centers: ", centers_file)

# Print cluster summary
cat("\nCluster summary:\n")
print(cluster_output %>% 
  group_by(Cluster) %>% 
  summarise(
    Samples = n(),
    PC1_center = mean(PC1),
    PC2_center = mean(PC2),
    Avg_Distance = mean(Distance_to_Center)
  ))