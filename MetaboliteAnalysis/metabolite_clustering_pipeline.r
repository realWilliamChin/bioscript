#!/usr/bin/env Rscript

library(ggplot2)
library(ggdendro)
library(readxl)
library(showtext)
library(factoextra)
library(cluster)
library(dplyr)
library(FactoMineR)
library(ggplotify)
library(grid)
library(NbClust)
library(writexl)
library(openxlsx)
library(tibble)
library(optparse)

set.seed(123)

# ========== 命令行参数 ==========
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = "metabolite_expression.xlsx",
              help = "代谢物丰度表，第一列为Sample_ID，其余列为代谢物丰度，不包含分组列 [默认: %default]"),
  make_option(c("-s", "--samples_described"), type = "character", default = NULL,
              help = "样本分组描述文件(可选)，需含Sample_ID和Group_ID两列;不提供则视为无分组"),
  make_option(c("-c", "--compound_def"), type = "character", default = NULL,
              help = "化合物分类注释表(可选)，需含化合物名(Compound_ID或Compound_name)、Class、SubClass列;不提供则热图不加分类注释"),
  make_option(c("-o", "--outdir"), type = "character", default = "clustering_result",
              help = "结果输出目录 [默认: %default]"),
  make_option(c("-k", "--k"), type = "integer", default = NULL,
              help = "手动指定聚类数;不指定时用NbClust全指标投票自动决定"),
  make_option(c("--min_k"), type = "integer", default = 2,
              help = "NbClust最小聚类数 [默认: %default]"),
  make_option(c("--max_k"), type = "integer", default = 10,
              help = "NbClust最大聚类数 [默认: %default]"),
  make_option(c("--pca_cutoff"), type = "double", default = 70,
              help = "PCA累计解释度阈值，用于NbClust降维 [默认: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

# 导入通用热图脚本
source("/home/colddata/qinqiang/script/CommonTools/Plot/Heatmap/heatmap_1.r")

testDir <- opt$outdir
dir.create(testDir, showWarnings = FALSE, recursive = TRUE)

# ========== 读取数据（保持原变量名：df） ==========
df_raw <- read_excel(opt$input)
df_raw <- as.data.frame(df_raw)
if (!("Sample_ID" %in% colnames(df_raw))) {
  colnames(df_raw)[1] <- "Sample_ID"
}

rownames(df_raw) <- df_raw$Sample_ID
df_expr <- df_raw[, -which(colnames(df_raw) == "Sample_ID"), drop = FALSE]

# 读取分组信息（可选）：提供samples_described则用其分组，否则视为无分组
has_group <- !is.null(opt$samples_described)
if (has_group) {
  sample_info <- read.table(opt$samples_described, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  if (!all(c("Sample_ID", "Group_ID") %in% colnames(sample_info))) {
    stop("samples_described文件必须包含Sample_ID和Group_ID两列")
  }
  sample_info <- sample_info[match(rownames(df_expr), sample_info$Sample_ID), ]
  if (any(is.na(sample_info$Sample_ID))) {
    stop("samples_described中的Sample_ID与input中的Sample_ID不完全匹配")
  }
  df <- data.frame(Group_ID = sample_info$Group_ID, df_expr, check.names = FALSE)
  rownames(df) <- rownames(df_expr)
} else {
  df <- df_expr
}

if (has_group) {
  clustering_df <- df %>% select(-Group_ID)
} else {
  clustering_df <- df
}

# 转换为数值型
for (e in colnames(clustering_df)) {
  clustering_df[[e]] <- as.numeric(clustering_df[[e]])
}

# 数据标准化
clustering_df_scale <- scale(clustering_df)

max_ncp <- min(nrow(clustering_df_scale), ncol(clustering_df_scale))
pca_for_clust <- PCA(clustering_df_scale, graph = FALSE, ncp = max_ncp)
var_contrib <- pca_for_clust$eig[, 2]
cum_contrib <- cumsum(var_contrib)
n_pc <- which(cum_contrib > opt$pca_cutoff)[1]
if (is.na(n_pc)) n_pc <- length(cum_contrib)
n_pc <- max(n_pc, 2)
pca_scores <- pca_for_clust$ind$coord[, 1:n_pc]

# 确定聚类数：指定--k时直接使用，否则用NbClust全指标投票
if (!is.null(opt$k)) {
  bestK <- as.integer(opt$k)
  if (bestK < 2) stop("--k 必须 >= 2")
  used_nbclust <- FALSE
  cat("使用手动指定的聚类数 K =", bestK, "\n")
} else {
  used_nbclust <- TRUE
  # 用PCA分数做NbClust投票
  nb_ward <- NbClust(pca_scores, min.nc = opt$min_k, max.nc = opt$max_k, method = "ward.D2", index = "all")

  index_names <- colnames(nb_ward$Best.nc)
  opt_k_vals <- as.integer(nb_ward$Best.nc[1, ])

  index_k_table <- data.frame(
    Index = index_names,
    Opt_K = opt_k_vals,
    row.names = NULL
  )

  # 定义核心指标清单（保持原变量名）
  core_index <- c("Hartigan", "CH", "Silhouette", "DB", "Tracew", "Ball")
  core_k_table <- subset(index_k_table, Index %in% core_index)

  # 统计每个K获得的票数
  vote_count <- table(index_k_table$Opt_K)

  bestK <- as.integer(names(which.max(vote_count)))
  if (length(bestK) == 0 || bestK == 0) bestK <- 2
}

hc <- hclust(dist(clustering_df_scale), method = "ward.D2")
ward_group <- cutree(hc, k = bestK)
ward_group <- as.factor(ward_group)

clustering_df.pca <- PCA(clustering_df_scale, graph = F)

add_square_coord <- function(grp_pca) {
  gb <- ggplot_build(grp_pca)

  all_x <- c()
  all_y <- c()

  for (i in seq_along(gb$data)) {
    if ("x" %in% names(gb$data[[i]])) {
      all_x <- c(all_x, gb$data[[i]]$x)
    }
    if ("y" %in% names(gb$data[[i]])) {
      all_y <- c(all_y, gb$data[[i]]$y)
    }
  }

  lim <- max(abs(c(all_x, all_y)), na.rm = TRUE) * 1.1

  grp_pca + coord_fixed(
    ratio = 1,
    xlim = c(-lim, lim),
    ylim = c(-lim, lim)
  )
}

# 按Cluster分组
grp_pca <- fviz_pca_ind(
  clustering_df.pca,
  mean.point = F,
  label = "none",
  habillage = ward_group,
  palette = "lancet",
  addEllipses = TRUE,
  ellipse.type = "t",
  ellipse.level = 0.95,
  ggtheme = theme_bw(base_size = 14)
)
grp_pca <- add_square_coord(grp_pca)
ggsave(file.path(testDir, "01_Cluster_PCA_ward_by_cluster.png"), grp_pca, width = 8, height = 8, dpi = 300)

# 按Group分组（仅在提供分组信息时输出）
if (has_group) {
  grp_pca <- fviz_pca_ind(
    clustering_df.pca,
    mean.point = F,
    label = "none",
    habillage = df$Group_ID,
    palette = "lancet",
    addEllipses = TRUE,
    ellipse.type = "t",
    ellipse.level = 0.95,
    ggtheme = theme_bw(base_size = 14)
  )
  grp_pca <- add_square_coord(grp_pca)
  ggsave(file.path(testDir, "02_Cluster_PCA_ward_by_group.png"), grp_pca, width = 8, height = 8, dpi = 300)
}

# ========== 热图 - 增加Compound_def注释（可选） ==========
# 提供compound_def时读取化合物分类作为行注释，否则热图不加行注释
if (!is.null(opt$compound_def)) {
  compound_def <- read_excel(opt$compound_def)
  compound_def <- as.data.frame(compound_def)

  # 兼容不同Compound_def表头：Compound_ID/Compound_name，Subclass/SubClass
  compound_id_col <- if ("Compound_ID" %in% colnames(compound_def)) "Compound_ID" else "Compound_name"
  subclass_col <- if ("Subclass" %in% colnames(compound_def)) "Subclass" else "SubClass"
  if (!(compound_id_col %in% colnames(compound_def))) {
    stop("Compound_def中必须包含Compound_ID或Compound_name列")
  }
  if (!("Class" %in% colnames(compound_def))) {
    stop("Compound_def中必须包含Class列")
  }
  if (!(subclass_col %in% colnames(compound_def))) {
    stop("Compound_def中必须包含Subclass或SubClass列")
  }

  # 确保化合物顺序一致
  compound_order <- colnames(clustering_df_scale)
  compound_anno <- compound_def[match(compound_order, compound_def[[compound_id_col]]), ]
  rownames(compound_anno) <- compound_order

  # 未匹配到注释的化合物标记为Unknown，避免因空列报错
  compound_anno$Class[is.na(compound_anno$Class)] <- "Unknown"
  compound_anno[[subclass_col]][is.na(compound_anno[[subclass_col]])] <- "Unknown"

  # 行注释：代谢物分类
  row_anno <- data.frame(
    Class = compound_anno$Class,
    Subclass = compound_anno[[subclass_col]]
  )
  rownames(row_anno) <- compound_order
} else {
  row_anno <- NA  # 无分类注释
}

# 列注释：样本信息（注意：annotation_col的行名 = 数据的列名）
# 有分组时增加Group色条，无分组时只标Cluster
if (has_group) {
  sample_anno <- data.frame(
    Cluster = factor(paste0("C", ward_group)),
    Group = df$Group_ID
  )
} else {
  sample_anno <- data.frame(
    Cluster = factor(paste0("C", ward_group))
  )
}
rownames(sample_anno) <- rownames(clustering_df_scale)  # 样本名

# 转置矩阵（行=代谢物，列=样本）
metabolite_matrix <- t(clustering_df_scale)

# 版本1：行按原始顺序（不聚类）
smart_heatmap(
  matrix_data = metabolite_matrix,
  annotation_row = row_anno,
  annotation_col = sample_anno,
  filename = file.path(testDir, "03_Heatmap_original_order.png"),
  main = paste0("Cluster number (ward.D2) = ", bestK),
  scale = "row",
  cluster_rows = FALSE,              # 不聚类行，保持原始顺序
  cluster_cols = TRUE,
  clustering_method = "ward.D2",
  clustering_distance_cols = "euclidean",
  treeheight_col = 40,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 7,
  fontsize_col = 7,
  color = colorRampPalette(c("#4575B4","white","#D73027"))(100),
  autoset_image_specification = TRUE
)

# 版本2：行聚类（按相似度排序）
smart_heatmap(
  matrix_data = metabolite_matrix,
  annotation_row = row_anno,
  annotation_col = sample_anno,
  filename = file.path(testDir, "03_Heatmap_row_clustered.png"),
  main = paste0("Cluster number (ward.D2) = ", bestK),
  scale = "row",
  cluster_rows = TRUE,               # 聚类行
  cluster_cols = TRUE,
  clustering_method = "ward.D2",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  treeheight_row = 30,
  treeheight_col = 40,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 7,
  fontsize_col = 7,
  color = colorRampPalette(c("#4575B4","white","#D73027"))(100),
  autoset_image_specification = TRUE
)

# ========== 输出Excel - 保持原文件名 ==========
# NbClust投票结果（仅在自动模式下输出）
if (used_nbclust) {
  vote_count_df <- as.data.frame(vote_count)
  colnames(vote_count_df) <- c('Cluster_k', 'Vote_Num')
  write_xlsx(vote_count_df, file.path(testDir, "04_Clustering_full_index.xlsx"))

  colnames(core_k_table) <- c('Index_Name', 'Recommend_K')
  write_xlsx(core_k_table, file.path(testDir, "04_Clustering_core_index.xlsx"))
}

df$Cluster <- paste0("C", ward_group)
df_result <- rownames_to_column(df, var = "Sample_ID")
df_result$.orig_order <- seq_len(nrow(df_result))  # 记录输入原始顺序
df_result <- df_result[order(df_result$Cluster, df_result$.orig_order), ]
df_result$.orig_order <- NULL
if (has_group) {
  df_result <- df_result %>%
    select(Sample_ID, Group_ID, Cluster, everything())
} else {
  df_result <- df_result %>%
    select(Sample_ID, Cluster, everything())
}
write_xlsx(df_result, file.path(testDir, "04_SampleID_to_clusterID.xlsx"))

cat("===== 分析完成 =====\n")
print(paste("最佳聚类数K =", bestK))
