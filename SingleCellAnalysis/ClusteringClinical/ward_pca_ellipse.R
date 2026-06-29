#!/usr/bin/env Rscript
# Created Time  : 2026/06/29
# Author        : William GoGo
# ClusteringClinical 第 2 步: Ward 层次聚类 + PCA 椭圆图
# 在 NbClust 选定的 K 上, 先做 Ward 聚类得到 cluster 标签, 再做 PCA 可视化。
# 支持按多个分组列(如 cluster / 年龄 / 疾病 / 性别)分别上色画图。

suppressPackageStartupMessages({
    library(optparse)
    library(readxl)
    library(writexl)
    library(dplyr)
    library(FactoMineR)
    library(factoextra)
    library(tibble)
})

option_list <- list(
    make_option(c("-i", "--input"), type = "character",
                help = "输入数据 xlsx/csv/tsv"),
    make_option(c("-o", "--output_dir"), type = "character",
                help = "输出目录"),
    make_option(c("--sheet"), type = "character", default = NULL),
    make_option(c("--id_col"), type = "character", default = NULL),
    make_option(c("--features"), type = "character",
                help = "参与聚类的列名, 逗号分隔 (必须)"),
    make_option(c("--feature_aliases"), type = "character", default = NULL,
                help = "用于显示的特征英文名, 逗号分隔; 顺序对应 --features"),
    make_option(c("--best_k"), type = "integer", default = 3),
    make_option(c("--color_by"), type = "character", default = "cluster",
                help = "椭圆图分组方式, 逗号分隔: cluster + 数据中其他列名"),
    make_option(c("--output_table"), type = "character", default = NULL,
                help = "保存样本 → cluster 标签的 xlsx 路径")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$input) || is.null(opt$output_dir) || is.null(opt$features)) {
    print_help(OptionParser(option_list = option_list))
    stop("必须提供 --input / --output_dir / --features")
}

dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(123)

# 读数据
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

# Ward 聚类
hc <- hclust(dist(clustering_df_scale), method = "ward.D2")
ward_group <- as.factor(cutree(hc, k = opt$best_k))
df$Cluster <- ward_group

# PCA
pca_obj <- PCA(clustering_df_scale, graph = FALSE)

# 用别名替换 PCA 显示的变量名
if (!is.null(opt$feature_aliases)) {
    aliases <- trimws(strsplit(opt$feature_aliases, ",")[[1]])
    colnames(clustering_df_scale) <- aliases
}

# 多种着色方式画图
color_keys <- trimws(strsplit(opt$color_by, ",")[[1]])
for (key in color_keys) {
    habillage <- if (key == "cluster") ward_group else df[[key]]
    if (is.null(habillage)) {
        warning("跳过 ", key, ": 数据中无此列")
        next
    }
    p <- fviz_pca_ind(
        pca_obj,
        mean.point = FALSE, label = "none",
        habillage = habillage,
        palette = "lancet",
        addEllipses = TRUE,
        ellipse.type = "t",
        ellipse.level = 0.95,
        ggtheme = theme_bw(base_size = 14))
    out_png <- file.path(opt$output_dir,
                         paste0("02_PCA_ward_by_", key, ".png"))
    ggplot2::ggsave(out_png, p, width = 9, height = 7, dpi = 300)
    message("保存 ", out_png)
}

# 保存样本 → cluster
if (!is.null(opt$output_table)) {
    out_df <- df %>% rownames_to_column(var = "SampleID")
    write_xlsx(out_df, opt$output_table)
    message("保存 sample→cluster 表 → ", opt$output_table)
}
