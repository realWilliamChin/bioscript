#!/usr/bin/env Rscript
# Created Time  : 2026/06/29
# Author        : William GoGo
# ClusteringClinical 第 3 步: 在 Ward 聚类的基础上画标准化热图
# 行=指标, 列=样本, 列顶部标 cluster 颜色条。

suppressPackageStartupMessages({
    library(optparse)
    library(readxl)
    library(pheatmap)
    library(grid)
    library(dplyr)
})

option_list <- list(
    make_option(c("-i", "--input"), type = "character",
                help = "输入数据 xlsx/csv/tsv"),
    make_option(c("-o", "--output"), type = "character",
                help = "输出 png 路径"),
    make_option(c("--sheet"), type = "character", default = NULL),
    make_option(c("--id_col"), type = "character", default = NULL),
    make_option(c("--features"), type = "character",
                help = "参与聚类的列名, 逗号分隔 (必须)"),
    make_option(c("--feature_aliases"), type = "character", default = NULL,
                help = "用于显示的特征英文名, 逗号分隔; 顺序对应 --features"),
    make_option(c("--best_k"), type = "integer", default = 3),
    make_option(c("--width"), type = "integer", default = 3200),
    make_option(c("--height"), type = "integer", default = 2000),
    make_option(c("--dpi"), type = "integer", default = 300)
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$input) || is.null(opt$output) || is.null(opt$features)) {
    print_help(OptionParser(option_list = option_list))
    stop("必须提供 --input / --output / --features")
}

dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)
set.seed(123)

# 读
ext <- tolower(tools::file_ext(opt$input))
if (ext %in% c("xlsx", "xls")) {
    df <- as.data.frame(read_excel(opt$input, sheet = opt$sheet))
} else {
    df <- read.csv(opt$input, sep = if (ext == "tsv") "\t" else ",",
                   header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
}
if (!is.null(opt$id_col)) {
    rownames(df) <- df[[opt$id_col]]
}

feats <- trimws(strsplit(opt$features, ",")[[1]])
clustering_df <- df[, feats, drop = FALSE]
for (col in colnames(clustering_df)) clustering_df[[col]] <- as.numeric(clustering_df[[col]])
clustering_df_scale <- scale(clustering_df)

if (!is.null(opt$feature_aliases)) {
    aliases <- trimws(strsplit(opt$feature_aliases, ",")[[1]])
    colnames(clustering_df_scale) <- aliases
}

dist_mat <- dist(clustering_df_scale, method = "euclidean")
hc <- hclust(dist_mat, method = "ward.D2")
cluster_result <- cutree(hc, k = opt$best_k)

annotation_col <- data.frame(Cluster = factor(paste0("C", cluster_result)))
rownames(annotation_col) <- rownames(clustering_df_scale)

p <- pheatmap(
    t(clustering_df_scale),
    scale = "none",
    clustering_method = "ward.D2",
    clustering_distance_cols = "euclidean",
    clustering_distance_rows = "euclidean",
    annotation_col = annotation_col,
    treeheight_col = 40, treeheight_row = 30,
    show_rownames = TRUE, show_colnames = TRUE,
    silent = TRUE, fontsize = 12,
    color = colorRampPalette(c("#4575B4", "white", "#D73027"))(100),
    main = paste0("heatmap | best_k = ", opt$best_k))

png(opt$output, width = opt$width, height = opt$height, res = opt$dpi)
grid.rect(gp = gpar(fill = "gray95", col = NA))
print(p, newpage = FALSE)
dev.off()
message("保存 → ", opt$output)
