#!/usr/bin/env Rscript
# Created Time  : 2026/06/29
# Author        : William GoGo
# Pipeline 第 11 步: 把 03 步导出的 h5ad 转成 Seurat rds, 给下游 Seurat/CellChat/monocle3 用。

suppressPackageStartupMessages({
    library(optparse)
    library(Seurat)
    library(reticulate)
})

option_list <- list(
    make_option(c("-i", "--input_h5ad"), type = "character",
                help = "scanpy 流程输出的 h5ad (例: adataFiltrd_celltypist.h5ad)"),
    make_option(c("-o", "--output_rds"), type = "character",
                help = "输出的 Seurat rds 路径"),
    make_option(c("--python"), type = "character",
                default = "~/miniconda3/envs/test_scanpy_python311/bin/python",
                help = "reticulate 使用的 python 解释器 [default %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$input_h5ad) || is.null(opt$output_rds)) {
    print_help(OptionParser(option_list = option_list))
    stop("必须提供 --input_h5ad 与 --output_rds")
}

use_python(opt$python, required = TRUE)
ad <- import("anndata")

message("读取 ", opt$input_h5ad)
adata <- ad$read_h5ad(opt$input_h5ad)

# Seurat 需要 基因×细胞 的矩阵
counts <- t(adata$raw$X)
cell_meta <- as.data.frame(adata$obs)
gene_meta <- as.data.frame(adata$var)
rownames(counts) <- rownames(gene_meta)
colnames(counts) <- rownames(cell_meta)

seurat_obj <- CreateSeuratObject(counts = counts, meta.data = cell_meta)
message("Seurat 对象: ", paste(dim(seurat_obj@meta.data), collapse = " x "))

dir.create(dirname(opt$output_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(seurat_obj, file = opt$output_rds)
message("保存 → ", opt$output_rds)
