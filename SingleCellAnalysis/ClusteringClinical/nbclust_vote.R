#!/usr/bin/env Rscript
# Created Time  : 2026/06/29
# Author        : William GoGo
# ClusteringClinical 第 1 步: NbClust 多指标投票, 推荐最佳 K
# 输入数值数据 (行=样本, 列=指标), 输出:
#   - <out_dir>/01_NbClust_full_index.xlsx     全部指标 → 推荐 K 投票表
#   - <out_dir>/01_NbClust_core_index.xlsx     核心指标 (CH/Silhouette/DB ... ) 推荐 K

suppressPackageStartupMessages({
    library(optparse)
    library(readxl)
    library(writexl)
    library(dplyr)
    library(NbClust)
})

option_list <- list(
    make_option(c("-i", "--input"), type = "character",
                help = "输入数据 xlsx/csv/tsv, 行=样本, 列=指标"),
    make_option(c("-o", "--output_dir"), type = "character",
                help = "输出目录"),
    make_option(c("--sheet"), type = "character", default = NULL,
                help = "xlsx 的 sheet 名"),
    make_option(c("--id_col"), type = "character", default = NULL,
                help = "用作样本 ID 的列名 (会设为 rownames 然后去掉)"),
    make_option(c("--features"), type = "character", default = NULL,
                help = "参与聚类的列名, 逗号分隔; 默认全部数值列"),
    make_option(c("--min_nc"), type = "integer", default = 2),
    make_option(c("--max_nc"), type = "integer", default = 10),
    make_option(c("--method"), type = "character", default = "ward.D2")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$input) || is.null(opt$output_dir)) {
    print_help(OptionParser(option_list = option_list))
    stop("必须提供 --input 与 --output_dir")
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
    df <- df[, setdiff(colnames(df), opt$id_col), drop = FALSE]
}

if (!is.null(opt$features)) {
    feats <- trimws(strsplit(opt$features, ",")[[1]])
    df <- df[, feats, drop = FALSE]
}

# 强制数值
for (col in colnames(df)) df[[col]] <- as.numeric(df[[col]])
df_scale <- scale(df)

message("NbClust 跑 method=", opt$method, " min=", opt$min_nc, " max=", opt$max_nc)
nb <- NbClust(df_scale, min.nc = opt$min_nc, max.nc = opt$max_nc,
              method = opt$method, index = "all")

index_names <- colnames(nb$Best.nc)
opt_k_vals <- as.integer(nb$Best.nc[1, ])
index_k_table <- data.frame(Index = index_names,
                            Opt_K = opt_k_vals,
                            row.names = NULL)

core_index <- c("Hartigan", "CH", "Silhouette", "DB", "TraceW", "Ball")
core_k_table <- subset(index_k_table, Index %in% core_index)
colnames(core_k_table) <- c("Index_Name", "Recommend_K")

vote_count <- as.data.frame(table(index_k_table$Opt_K))
colnames(vote_count) <- c("Cluster_k", "Vote_Num")

write_xlsx(vote_count, file.path(opt$output_dir, "01_NbClust_full_index.xlsx"))
write_xlsx(core_k_table, file.path(opt$output_dir, "01_NbClust_core_index.xlsx"))

message("最常被推荐的 K = ",
        as.character(vote_count$Cluster_k[which.max(vote_count$Vote_Num)]))
message("结果 → ", opt$output_dir)
