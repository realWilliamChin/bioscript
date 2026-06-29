#!/usr/bin/env Rscript
# Created Time  : 2026/06/29
# Author        : William GoGo
# Pipeline 第 13 步: 给每个 group 构建 CellChat 对象, 保存为 <group>_cellchat.rds。
# 后续 14 步用这些 rds 做组对比较。

suppressPackageStartupMessages({
    library(optparse)
    library(Seurat)
    library(CellChat)
    library(dplyr)
    library(stringr)
})

option_list <- list(
    make_option(c("-i", "--input_rds"), type = "character",
                help = "11 步输出的 Seurat rds"),
    make_option(c("-o", "--output_dir"), type = "character",
                help = "存放 <group>_cellchat.rds 的目录"),
    make_option(c("--annot_tsv"), type = "character",
                help = "两列 (cluster,label) 的注释表, tab/csv 分隔"),
    make_option(c("--cluster_col"), type = "character", default = "leiden_r025",
                help = "用作 cluster 的列名 [default %default]"),
    make_option(c("--group_col"), type = "character", default = "group"),
    make_option(c("--target_clusters"), type = "character", default = NULL,
                help = "要保留的 cluster ID, 逗号分隔; 默认全部"),
    make_option(c("--target_groups"), type = "character", default = NULL,
                help = "要构建 CellChat 的 group 列表, 逗号分隔; 默认 group 列所有值"),
    make_option(c("--db"), type = "character", default = "human",
                help = "CellChatDB 选择 human/mouse [default %default]"),
    make_option(c("--min_cells"), type = "integer", default = 10)
)
opt <- parse_args(OptionParser(option_list = option_list))
for (k in c("input_rds", "output_dir", "annot_tsv")) {
    if (is.null(opt[[k]])) {
        print_help(OptionParser(option_list = option_list))
        stop("必须提供 --", k)
    }
}

# utils
script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- script_args[grep("^--file=", script_args)]
script_dir <- if (length(file_arg)) {
    normalizePath(dirname(sub("^--file=", "", file_arg[1])))
} else {
    getwd()
}
source(file.path(script_dir, "..", "utils", "seurat_helpers.R"))

dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

seurat_obj <- readRDS(opt$input_rds)

annot_df <- read.table(opt$annot_tsv, header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE,
                       comment.char = "", check.names = FALSE)
if (ncol(annot_df) < 2) {
    annot_df <- read.csv(opt$annot_tsv, header = TRUE,
                         stringsAsFactors = FALSE, check.names = FALSE)
}
cluster2annot_map <- cluster2annot_mapping(annot_df)

target_clusters <- if (!is.null(opt$target_clusters))
    trimws(strsplit(opt$target_clusters, ",")[[1]]) else NULL
target_groups <- if (!is.null(opt$target_groups)) {
    trimws(strsplit(opt$target_groups, ",")[[1]])
} else {
    unique(as.character(seurat_obj[[opt$group_col]][, 1]))
}

CellChatDB <- if (opt$db == "human") CellChatDB.human else CellChatDB.mouse

for (grp in target_groups) {
    message("==== group: ", grp, " ====")
    cells_keep <- rownames(seurat_obj@meta.data)[
        as.character(seurat_obj[[opt$group_col]][, 1]) == grp]
    if (!is.null(target_clusters)) {
        cluster_vals <- as.character(seurat_obj[[opt$cluster_col]][, 1])
        cells_keep <- intersect(cells_keep,
                                rownames(seurat_obj@meta.data)[
                                    cluster_vals %in% target_clusters])
    }
    grpObj <- subset(seurat_obj, cells = cells_keep)
    grpObj <- add_cluster_annot(grpObj, opt$cluster_col, cluster2annot_map,
                                new_col = "r025_annot")
    grpObj$r025_annot <- factor(grpObj$r025_annot)

    data.input <- GetAssayData(grpObj, assay = "RNA", layer = "counts")
    meta <- data.frame(labels = grpObj$r025_annot,
                       row.names = colnames(grpObj))

    cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
    cellchat@DB <- CellChatDB
    cellchat <- subsetData(cellchat)
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    cellchat <- computeCommunProb(cellchat)
    cellchat <- filterCommunication(cellchat, min.cells = opt$min_cells)
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)

    out_file <- file.path(opt$output_dir, paste0(grp, "_cellchat.rds"))
    saveRDS(cellchat, file = out_file)
    message("保存 ", out_file)
}
