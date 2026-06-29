#!/usr/bin/env Rscript
# Created Time  : 2026/06/29
# Author        : William GoGo
# Pipeline 第 15 步: monocle3 拟时序分析
#   - 全细胞一次 (overall)
#   - 每个 group 各跑一次 (per-group)
# 输出 UMAP 轨迹图 + cluster→pseudotime CSV + 各组 CDS 列表 Rdata。

suppressPackageStartupMessages({
    library(optparse)
    library(monocle3)
    library(Seurat)
    library(ggplot2)
    library(dplyr)
    library(stringr)
})

option_list <- list(
    make_option(c("-i", "--input_rds"), type = "character",
                help = "11 步输出的 Seurat rds"),
    make_option(c("-o", "--output_dir"), type = "character",
                help = "输出目录"),
    make_option(c("--cluster_col"), type = "character", default = "leiden_r025"),
    make_option(c("--group_col"), type = "character", default = "group"),
    make_option(c("--target_clusters"), type = "character", default = NULL,
                help = "保留的 cluster ID 列表 (逗号分隔), 默认全部"),
    make_option(c("--root_cluster"), type = "character", default = NULL,
                help = "用作 root 的 cluster ID (例: '1'); 不传则不指定 root"),
    make_option(c("--n_dim"), type = "integer", default = 50,
                help = "preprocess_cds num_dim [default %default]"),
    make_option(c("--per_group"), action = "store_true", default = FALSE,
                help = "对每个 group 也单独跑一次")
)
opt <- parse_args(OptionParser(option_list = option_list))
for (k in c("input_rds", "output_dir")) {
    if (is.null(opt[[k]])) {
        print_help(OptionParser(option_list = option_list))
        stop("必须提供 --", k)
    }
}

dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 工具: 单次 monocle3 + 出图 ----
run_monocle3 <- function(seurat_subset, opt, label) {
    expression_matrix <- GetAssayData(seurat_subset, assay = "RNA", layer = "counts")
    cell_metadata <- seurat_subset@meta.data
    gene_annotation <- data.frame(gene_short_name = rownames(expression_matrix))
    rownames(gene_annotation) <- rownames(expression_matrix)

    cds <- new_cell_data_set(expression_matrix,
                             cell_metadata = cell_metadata,
                             gene_metadata = gene_annotation)
    cds <- preprocess_cds(cds, num_dim = opt$n_dim)
    cds <- align_cds(cds)
    cds <- reduce_dimension(cds)
    cds <- cluster_cells(cds)
    cds <- learn_graph(cds)

    if (!is.null(opt$root_cluster)) {
        root_cells <- colnames(cds)[
            as.character(colData(cds)[[opt$cluster_col]]) == opt$root_cluster]
        if (length(root_cells) > 0) {
            cds <- order_cells(cds, root_cells = root_cells)
        }
    }

    p <- plot_cells(cds, color_cells_by = "pseudotime",
                    group_cells_by = opt$cluster_col,
                    label_cell_groups = TRUE, label_leaves = FALSE,
                    label_groups_by_cluster = TRUE,
                    label_branch_points = TRUE, graph_label_size = 1.5)
    ggsave(file.path(opt$output_dir, paste0(label, "_umap_trajectoryPseudotime.png")),
           p, width = 14, height = 10, dpi = 300)

    p1 <- plot_cells(cds, color_cells_by = opt$cluster_col,
                     group_cells_by = "cluster",
                     label_cell_groups = FALSE, label_leaves = TRUE,
                     label_branch_points = TRUE, graph_label_size = 3)
    ggsave(file.path(opt$output_dir, paste0(label, "_umap_trajectoryCluster.png")),
           p1, width = 14, height = 10, dpi = 300)

    cluster2pseudotime <- data.frame(
        Cell_Barcode = colnames(seurat_subset),
        sample_ID = seurat_subset$sampleID,
        seurat_ClusterID = seurat_subset[[opt$cluster_col]][, 1],
        monocle3_ClusterID = clusters(cds, reduction_method = "UMAP"),
        Pseudotime = pseudotime(cds))
    write.csv(cluster2pseudotime,
              file.path(opt$output_dir, paste0(label, "_cluster2pseudotime.csv")),
              row.names = FALSE, fileEncoding = "UTF-8")

    return(cds)
}


# ---- 主流程 ----
seurat_obj <- readRDS(opt$input_rds)

if (!is.null(opt$target_clusters)) {
    target_clusters <- trimws(strsplit(opt$target_clusters, ",")[[1]])
    seurat_obj <- subset(seurat_obj,
                          subset = !!sym(opt$cluster_col) %in% target_clusters)
}

# overall
message("==== 全细胞 monocle3 ====")
cds_overall <- run_monocle3(seurat_obj, opt, label = "overall")

# per-group
if (opt$per_group) {
    grp2cds <- list()
    for (grp in unique(as.character(seurat_obj[[opt$group_col]][, 1]))) {
        message("==== group: ", grp, " ====")
        grp_cells <- rownames(seurat_obj@meta.data)[
            as.character(seurat_obj[[opt$group_col]][, 1]) == grp]
        grp_obj <- subset(seurat_obj, cells = grp_cells)
        label <- str_replace_all(grp, " ", "")
        grp2cds[[label]] <- run_monocle3(grp_obj, opt, label = label)
    }
    save(grp2cds,
         file = file.path(opt$output_dir, "grp2monocleCDS_lst.Rdata"))
}

message("monocle3 拟时序完成 → ", opt$output_dir)
