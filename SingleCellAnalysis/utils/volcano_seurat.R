#!/usr/bin/env Rscript
# Created Time  : 2026/06/29
# Author        : William GoGo
# Seurat DEG 火山图 + DoHeatmap 工具函数。
# source('utils/volcano_seurat.R') 调用。

suppressPackageStartupMessages({
    library(Seurat)
    library(ggplot2)
    library(ggrepel)
    library(dplyr)
    library(stringr)
})


#' Seurat FindMarkers 结果绘制火山图
#' 输入 degs_data 必须包含: Gene_Symbol, avg_log2FC, p_val_adj
#' @param degs_data data.frame
#' @param comb 主标题 (例: "HD vs Class_III")
#' @param workDir 输出目录
#' @param padj_cutoff p_val_adj 阈值
#' @param logfc_cutoff |avg_log2FC| 阈值
volcano_plot <- function(degs_data, comb, workDir,
                         padj_cutoff = 0.05, logfc_cutoff = 0.5) {
    degs_data$regulation <- ifelse(
        degs_data$p_val_adj < padj_cutoff & abs(degs_data$avg_log2FC) > logfc_cutoff,
        ifelse(degs_data$avg_log2FC > logfc_cutoff, "Up", "Down"),
        "NoSignificant")

    all_up   <- degs_data[degs_data$regulation == "Up", ]
    all_down <- degs_data[degs_data$regulation == "Down", ]
    up_total   <- nrow(all_up)
    down_total <- nrow(all_down)

    top10_up   <- all_up[order(all_up$p_val_adj), ][1:min(10, nrow(all_up)), ]
    top10_down <- all_down[order(all_down$p_val_adj), ][1:min(10, nrow(all_down)), ]
    label_gene <- rbind(top10_up, top10_down)

    # 修正极小 padj, 避免 Inf
    degs_data$p_val_adj <- ifelse(degs_data$p_val_adj < 1e-300, 1e-300, degs_data$p_val_adj)
    FC_max <- ceiling(max(abs(degs_data$avg_log2FC), na.rm = TRUE))

    p <- ggplot(degs_data, aes(x = avg_log2FC, y = -log10(p_val_adj), colour = regulation)) +
        geom_point(alpha = 0.7, size = 1.5) +
        geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff), lty = 4, col = "black", lwd = 0.8) +
        geom_hline(yintercept = -log10(padj_cutoff), lty = 4, col = "black", lwd = 0.8) +
        scale_color_manual(values = c("green", "gray60", "red"),
                           labels = c(paste0("Down(", down_total, ")"),
                                      "Not significant",
                                      paste0("Up(", up_total, ")"))) +
        geom_label_repel(data = label_gene,
                         aes(label = Gene_Symbol),
                         size = 2.5, max.overlaps = 100,
                         box.padding = 0.3, label.padding = 0.1,
                         segment.color = "black", show.legend = FALSE) +
        labs(x = "log2(Fold Change)", y = "-log10(Adjusted P-value)",
             title = comb, color = "Gene expression status") +
        scale_x_continuous(limits = c(-FC_max, FC_max)) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              plot.title = element_text(hjust = 0.5),
              legend.position = "right",
              legend.key.size = unit(1, "cm"),
              legend.key.height = unit(0.5, "cm"),
              legend.spacing.y = unit(0, "cm"))

    timestamp <- format(Sys.time(), "%Y%m%d")
    combName <- comb %>% str_replace(": ", "_") %>% str_replace_all(" ", "-")

    out_png <- file.path(workDir, paste0(paste(combName, "volcano", timestamp, sep = "_"), ".png"))
    ggsave(filename = out_png, plot = p, width = 8, height = 6, units = "in", dpi = 300)
    invisible(out_png)
}


#' 在合并后的 Seurat 对象上用 DoHeatmap 画 topN 基因
#' @param seurat_obj 已经 ScaleData 过的 Seurat 对象, group 列存在
#' @param comb 标题
#' @param pair 两个 group 名, 用作 factor levels
#' @param topN_genes 要画的基因向量
deg_heatmap <- function(seurat_obj, comb, pair, topN_genes) {
    seurat_obj@meta.data$group <- factor(seurat_obj@meta.data$group, levels = unlist(pair))
    cell_keep <- c(
        WhichCells(seurat_obj, expression = group == pair[1]),
        WhichCells(seurat_obj, expression = group == pair[2])
    )
    DoHeatmap(seurat_obj, group.by = "group", slot = "scale.data",
              cells = cell_keep, size = 3, angle = 0, features = topN_genes) +
        scale_fill_gradientn(colors = c("lightblue", "firebrick3")) +
        theme(plot.title = element_text(size = 14)) +
        labs(title = comb)
}
