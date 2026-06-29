#!/usr/bin/env Rscript
# Created Time  : 2026/06/29
# Author        : William GoGo
# Pipeline 第 12 步: Seurat 组对差异基因 (FindMarkers) + GO/KEGG 富集 + 火山图 + 热图。
# 复用 utils/seurat_helpers.R, utils/volcano_seurat.R, utils/enrich_plot_helpers.R。

suppressPackageStartupMessages({
    library(optparse)
    library(Seurat)
    library(dplyr)
    library(stringr)
    library(tibble)
    library(clusterProfiler)
    library(org.Hs.eg.db)
})

# ---- 解析参数 ----
option_list <- list(
    make_option(c("-i", "--input_rds"), type = "character",
                help = "11 步输出的 Seurat rds"),
    make_option(c("-o", "--output_dir"), type = "character",
                help = "输出根目录, 每个组对一个子目录"),
    make_option(c("--kegg_hierarchy"), type = "character",
                help = "KEGG 层级表 rds, 需含列 Gene_Symbol/KO_ID/Level_A/B/C"),
    make_option(c("--group_col"), type = "character", default = "group",
                help = "meta 中的分组列名 [default %default]"),
    make_option(c("--pairs"), type = "character", default = NULL,
                help = "指定组对, 形如 'HD,Class_III;HD,NoLN'; 默认两两组合"),
    make_option(c("--orgdb"), type = "character", default = "org.Hs.eg.db",
                help = "GO 富集用的 OrgDb 包名 [default %default]"),
    make_option(c("--padj_cutoff"), type = "double", default = 0.05),
    make_option(c("--logfc_cutoff"), type = "double", default = 0.5),
    make_option(c("--top_n_enrich"), type = "integer", default = 10),
    make_option(c("--utils_dir"), type = "character", default = NULL,
                help = "utils 目录, 默认相对脚本位置自动定位")
)
opt <- parse_args(OptionParser(option_list = option_list))
for (k in c("input_rds", "output_dir", "kegg_hierarchy")) {
    if (is.null(opt[[k]])) {
        print_help(OptionParser(option_list = option_list))
        stop("必须提供 --", k)
    }
}

# 解析 utils_dir
script_dir <- tryCatch({
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- args[grep("^--file=", args)]
    if (length(file_arg)) {
        normalizePath(dirname(sub("^--file=", "", file_arg[1])))
    } else {
        getwd()
    }
}, error = function(e) getwd())
utils_dir <- if (!is.null(opt$utils_dir)) {
    opt$utils_dir
} else {
    normalizePath(file.path(script_dir, "..", "utils"), mustWork = FALSE)
}
source(file.path(utils_dir, "seurat_helpers.R"))
source(file.path(utils_dir, "volcano_seurat.R"))
source(file.path(utils_dir, "enrich_plot_helpers.R"))

# 加载 OrgDb
suppressPackageStartupMessages(library(opt$orgdb, character.only = TRUE))
orgdb <- get(opt$orgdb)

# ---- 主流程 ----
message("读取 ", opt$input_rds)
seurat_obj <- readRDS(opt$input_rds)
Idents(seurat_obj) <- seurat_obj[[opt$group_col]][, 1]

seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

cleaned_kegg <- readRDS(opt$kegg_hierarchy)

# 组对
if (!is.null(opt$pairs)) {
    pair_list <- lapply(strsplit(opt$pairs, ";")[[1]],
                         function(p) trimws(strsplit(p, ",")[[1]]))
} else {
    unique_groups <- unique(seurat_obj@meta.data[[opt$group_col]])
    grp_pairs <- combn(unique_groups, 2)
    pair_list <- lapply(seq_len(ncol(grp_pairs)), function(i) grp_pairs[, i])
}

timestamp <- format(Sys.time(), "%Y%m%d")
directions <- c("Up", "Down", "Total")

dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

for (pair in pair_list) {
    comb <- paste(str_replace_all(pair[1], " ", ""),
                  str_replace_all(pair[2], " ", ""), sep = " vs ")
    pairDir <- file.path(opt$output_dir,
                         paste(str_replace_all(pair[1], " ", ""),
                               str_replace_all(pair[2], " ", ""), sep = "_vs_"))
    dir.create(pairDir, recursive = TRUE, showWarnings = FALSE)
    message("处理组对: ", comb)

    targetObj <- subset_groups(seurat_obj, pair, group_col = opt$group_col)
    combRes <- find_markers_pair(targetObj, pair[1], pair[2])

    degs_csv <- file.path(pairDir, paste0("diffExpre_", timestamp, ".csv"))
    write.csv(rownames_to_column(combRes$degs, "Gene_Symbol"), degs_csv, row.names = FALSE)

    annotated <- annotate_deg_direction(combRes$degs,
                                        padj_cutoff = opt$padj_cutoff,
                                        logfc_cutoff = opt$logfc_cutoff)

    top20_up <- annotated %>% filter(direction == "Up") %>% arrange(p_val_adj) %>% head(20)
    top20_down <- annotated %>% filter(direction == "Down") %>% arrange(p_val_adj) %>% head(20)

    comb_seuratObj <- ScaleData(combRes$seuratObj, features = rownames(combRes$seuratObj))
    heatmap_p <- deg_heatmap(comb_seuratObj, comb, pair,
                              c(rownames(top20_up), rownames(top20_down)))
    ggsave(file.path(pairDir,
                     paste0(str_replace_all(comb, " ", "-"), "_heatmap_", timestamp, ".png")),
           plot = heatmap_p, width = 8, height = 6, dpi = 300)

    deg_for_volcano <- rownames_to_column(annotated, "Gene_Symbol")
    volcano_plot(deg_for_volcano, comb, pairDir,
                 padj_cutoff = opt$padj_cutoff, logfc_cutoff = opt$logfc_cutoff)

    # GO/KEGG: 对 Up / Down / Total 各跑一次
    for (direction in directions) {
        genes <- if (direction == "Total") annotated
                 else annotated %>% filter(direction == !!direction)
        if (nrow(genes) == 0) next

        gene_ids <- tryCatch(
            mapIds(orgdb, keys = rownames(genes),
                   column = "ENTREZID", keytype = "SYMBOL", multiVals = "first"),
            error = function(e) { message("mapIds 失败: ", e$message); NULL }
        )
        if (is.null(gene_ids)) next

        goEnrich <- enrichGO(gene = gene_ids, OrgDb = orgdb,
                             keyType = "ENTREZID", ont = "ALL",
                             pAdjustMethod = "BH", qvalueCutoff = 0.05,
                             readable = TRUE)
        GO_dir <- file.path(pairDir, "GO_enrichment")
        dir.create(GO_dir, showWarnings = FALSE)
        combName <- paste(comb, direction, sep = "_")
        write.csv(as.data.frame(goEnrich),
                  file.path(GO_dir, paste0(combName, "_GO_", timestamp, ".csv")),
                  row.names = FALSE)
        if (nrow(as.data.frame(goEnrich)) > 0) {
            top <- enrichGO_topN_picking(goEnrich, opt$top_n_enrich)
            if (nrow(top) > 0) enrichGO_barplot(top, comb, direction, GO_dir)
        }

        KEGG_dir <- file.path(pairDir, "KEGG_enrichment")
        dir.create(KEGG_dir, showWarnings = FALSE)
        kegg_enrich_and_plot(genes, cleaned_kegg, comb, direction, KEGG_dir)
    }
}

message("Seurat DEG + 富集完成 → ", opt$output_dir)
