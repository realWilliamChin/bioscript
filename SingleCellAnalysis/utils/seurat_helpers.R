#!/usr/bin/env Rscript
# Created Time  : 2026/06/29
# Author        : William GoGo
# Seurat 流程通用工具函数。通过 source('utils/seurat_helpers.R') 在其他 R 脚本中调用。

suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
})


#' 从 cluster→annotation 的两列数据框生成命名向量, 用于 AddMetaData。
#' @param df 两列数据框: 第 1 列 cluster ID, 第 2 列 annotation 名
cluster2annot_mapping <- function(df) {
    df <- as.data.frame(df)
    annot_map <- df[[2]]
    names(annot_map) <- df[[1]]
    return(annot_map)
}


#' 把 cluster 注释附加到 Seurat 对象的 meta.data 上
#' @param seuratObj Seurat 对象
#' @param cluster_col 用于查表的 cluster 列名 (例: "leiden_r025")
#' @param annot_map 命名向量 (cluster2annot_mapping 输出)
#' @param new_col 新 meta 列名 (默认 "annot")
add_cluster_annot <- function(seuratObj, cluster_col, annot_map, new_col = "annot") {
    clusters <- as.character(seuratObj[[cluster_col]][, 1])
    annot_result <- annot_map[clusters]
    names(annot_result) <- colnames(seuratObj)
    AddMetaData(seuratObj, metadata = annot_result, col.name = new_col)
}


#' 按 group 列子集 Seurat 对象, 必要时做 NormalizeData
#' @param seuratObj Seurat 对象
#' @param target_groups 字符向量, 要保留的 group 值
#' @param group_col group 所在 meta 列 (默认 "group")
subset_groups <- function(seuratObj, target_groups, group_col = "group") {
    cells_keep <- rownames(seuratObj@meta.data)[
        as.character(seuratObj[[group_col]][, 1]) %in% target_groups]
    grpsObj <- subset(seuratObj, cells = cells_keep)
    Idents(grpsObj) <- group_col

    if (is.null(grpsObj[["RNA"]]$data) || nrow(grpsObj[["RNA"]]$data) == 0) {
        message("data 层为空, 进行 LogNormalize 归一化 ...")
        grpsObj <- NormalizeData(grpsObj,
                                 normalization.method = "LogNormalize",
                                 scale.factor = 1e4, verbose = FALSE)
    }
    return(grpsObj)
}


#' 在两个 group 之间跑 FindMarkers, 同时返回合并后的 Seurat 子集 (供后续 heatmap 用)
#' @param seuratObj 已经做过 NormalizeData 的 Seurat 对象, meta 中有 group 列
#' @param cond1 参考组
#' @param cond2 比较组
#' @return list(degs = data.frame, seuratObj = 合并子集)
find_markers_pair <- function(seuratObj, cond1, cond2,
                              min_pct = 0.1, logfc_threshold = 0.25,
                              test_use = "wilcox") {
    degs <- FindMarkers(object = seuratObj,
                        ident.1 = cond1, ident.2 = cond2,
                        min.pct = min_pct,
                        logfc.threshold = logfc_threshold,
                        test.use = test_use,
                        only.pos = FALSE, verbose = FALSE)

    cells1 <- rownames(seuratObj@meta.data)[seuratObj$group == cond1]
    cells2 <- rownames(seuratObj@meta.data)[seuratObj$group == cond2]
    comb_seuratObj <- subset(seuratObj, cells = c(cells1, cells2))
    list(degs = degs, seuratObj = comb_seuratObj)
}


#' 在 DEG 表上添加 direction (Up/Down/Not significant) 和 significance 标签
annotate_deg_direction <- function(degs,
                                   padj_col = "p_val_adj",
                                   logfc_col = "avg_log2FC",
                                   padj_cutoff = 0.05,
                                   logfc_cutoff = 0.5) {
    degs %>%
        mutate(
            direction = case_when(
                !!sym(padj_col) < padj_cutoff & !!sym(logfc_col) >  logfc_cutoff ~ "Up",
                !!sym(padj_col) < padj_cutoff & !!sym(logfc_col) < -logfc_cutoff ~ "Down",
                TRUE ~ "Not significant"
            ),
            significance = ifelse(!!sym(padj_col) < padj_cutoff,
                                  "Significant", "Not significant")
        )
}
