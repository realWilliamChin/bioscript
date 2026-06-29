#!/usr/bin/env Rscript
# Created Time  : 2026/06/29
# Author        : William GoGo
# Pipeline 第 16 步: SingleR 注释 (默认 MonacoImmuneData)
# 在 11 步导出的 Seurat 对象上做注释, 新增 SingleR_Label/Score 列, 保存为新 rds 并出统计表。

suppressPackageStartupMessages({
    library(optparse)
    library(Seurat)
    library(SingleR)
    library(celldex)
    library(SummarizedExperiment)
    library(dplyr)
    library(tibble)
})

option_list <- list(
    make_option(c("-i", "--input_rds"), type = "character",
                help = "11 步输出的 Seurat rds"),
    make_option(c("-o", "--output_rds"), type = "character",
                help = "带 SingleR 注释的新 rds 路径"),
    make_option(c("--output_dir"), type = "character",
                help = "csv 统计表输出目录"),
    make_option(c("--ref"), type = "character", default = "MonacoImmuneData",
                help = "celldex 参考集 [default %default]; 也支持 HumanPrimaryCellAtlasData 等"),
    make_option(c("--label_col"), type = "character", default = "label.fine",
                help = "SingleR 使用的 label 列 [default %default]"),
    make_option(c("--cluster_col"), type = "character", default = "leiden_r025",
                help = "做 cluster × annot 交叉表用的列"),
    make_option(c("--sample_col"), type = "character", default = "sampleID")
)
opt <- parse_args(OptionParser(option_list = option_list))
for (k in c("input_rds", "output_rds", "output_dir")) {
    if (is.null(opt[[k]])) {
        print_help(OptionParser(option_list = option_list))
        stop("必须提供 --", k)
    }
}

dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(opt$output_rds), recursive = TRUE, showWarnings = FALSE)

# 加载参考
message("加载 celldex::", opt$ref)
ref_fn <- get(opt$ref, envir = asNamespace("celldex"))
ref <- ref_fn(legacy = TRUE)
labels <- ref[[opt$label_col]]

seurat_obj <- readRDS(opt$input_rds)
seurat_obj <- NormalizeData(seurat_obj,
                            normalization.method = "LogNormalize",
                            scale.factor = 1e4, verbose = FALSE)

sce <- as.SingleCellExperiment(seurat_obj, assay = "RNA")
pred <- SingleR(test = sce, ref = ref, labels = labels,
                genes = "de", de.method = "wilcox", de.n = 50,
                aggr.ref = FALSE)
seurat_obj$SingleR_Label <- pred$labels
seurat_obj$SingleR_Score <- pred$scores

saveRDS(seurat_obj, file = opt$output_rds)
message("保存 ", opt$output_rds)

# 统计表
cluster_x_label <- as.data.frame.matrix(
    table(seurat_obj[[opt$cluster_col]][, 1], seurat_obj$SingleR_Label))
cluster_x_label <- rownames_to_column(cluster_x_label, "clusterID")
write.csv(cluster_x_label,
          file.path(opt$output_dir,
                    paste0(opt$cluster_col, "_clusters_singleRAnnot.csv")),
          row.names = FALSE)

sample_x_label <- as.data.frame.matrix(
    table(seurat_obj[[opt$sample_col]][, 1], seurat_obj$SingleR_Label))
sample_x_label <- rownames_to_column(sample_x_label, opt$sample_col)
write.csv(sample_x_label,
          file.path(opt$output_dir, "sample_singleRAnnot.csv"),
          row.names = FALSE)

prop_by_sample <- prop.table(
    t(table(seurat_obj[[opt$sample_col]][, 1], seurat_obj$SingleR_Label)),
    margin = 2) * 100
prop_df <- as.data.frame.matrix(round(prop_by_sample, 2))
prop_df <- rownames_to_column(prop_df, "cellType")
write.csv(prop_df,
          file.path(opt$output_dir, "cellProp_by_sample.csv"),
          row.names = FALSE)

message("SingleR 注释完成 → ", opt$output_rds)
