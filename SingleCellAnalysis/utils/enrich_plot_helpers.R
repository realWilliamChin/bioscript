#!/usr/bin/env Rscript
# Created Time  : 2026/06/29
# Author        : William GoGo
# GO/KEGG 富集结果的 topN 选取与可视化通用函数。
# 通过 source('utils/enrich_plot_helpers.R') 在其他 R 脚本中使用。

suppressPackageStartupMessages({
    library(dplyr)
    library(stringr)
    library(ggplot2)
    library(RColorBrewer)
    library(clusterProfiler)
    library(ggforce)
})


#' 从字符串中清掉前导数字、空格和方括号/圆括号包裹的后缀
#' 例: "00010 Glycolysis [PATH:ko00010]" -> "Glycolysis"
extract_pathway_name <- function(pathway_string) {
    cleaned <- gsub("^\\d+\\s+|\\s*\\[[^]]*\\]$|\\s*\\([^)]*\\)$", "", pathway_string)
    return(trimws(cleaned))
}


#' 在 enrichGO 结果中按 BP/CC/MF 各取 top_n, 添加 -log10(p.adjust) 列
#' @param goEnrich enrichGO() 输出对象
#' @param top_n 每个 ontology 类别保留的数目
enrichGO_topN_picking <- function(goEnrich, top_n) {
    pick <- function(ont) {
        goEnrich@result %>%
            filter(ONTOLOGY == ont) %>%
            arrange(p.adjust) %>%
            slice_head(n = top_n)
    }
    bind_rows(pick("BP"), pick("CC"), pick("MF")) %>%
        mutate(log_pval = -log10(p.adjust))
}


#' 绘制 GO 条形图 (BP/CC/MF 分面)
#' 同时输出 pdf 与 png
#' @param enrichGO_topN enrichGO_topN_picking 的结果
#' @param comb 主标题前缀 (如 "HD vs Class_III")
#' @param direction "Up" / "Down" / "Total"
#' @param workDir 输出目录
enrichGO_barplot <- function(enrichGO_topN, comb, direction, workDir) {
    new_comb <- str_replace(comb, ": ", "_") %>% str_replace_all(" ", "-")
    combName <- paste(new_comb, direction, sep = "_")
    color_map <- c("BP" = "#33a02c", "CC" = "#e31a1c", "MF" = "#1f78b4")

    p <- ggplot(enrichGO_topN,
                aes(x = log_pval, y = reorder(Description, log_pval), fill = ONTOLOGY)) +
        geom_col(width = 0.7) +
        facet_wrap(~ONTOLOGY, ncol = 1, scales = "free_y") +
        scale_fill_manual(values = color_map) +
        labs(x = "-log10(p-value)", y = "GO Term",
             title = paste0(comb, " (", direction, ")\nTop GO Term Enrichment")) +
        theme_bw() +
        theme(panel.grid = element_blank(),
              axis.text.y = element_text(size = 10),
              axis.text.x = element_text(size = 10),
              strip.text = element_text(size = 12, face = "bold"),
              legend.position = "none",
              plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

    ggsave(file.path(workDir, paste0(combName, "_GO.pdf")), p,
           width = 14, height = 10, dpi = 300)
    ggsave(file.path(workDir, paste0(combName, "_GO.png")), p,
           width = 14, height = 10, dpi = 300)
}


#' 基于 KEGG 层级表 (Gene_Symbol/KO_ID/Level_A/Level_B/Level_C) 执行通用基因→KO→通路富集。
#' kegg_hierarchy 至少需要列: Gene_Symbol, KO_ID, Level_A, Level_B, Level_C
perform_kegg_enrichment <- function(genes, gene2ko, ko2pathway,
                                    universe_genes = NULL,
                                    pvalue_cutoff = 0.05,
                                    qvalue_cutoff = 0.2) {
    gene_ko_mapping <- gene2ko %>% filter(Gene_Symbol %in% genes)
    ko_list <- unique(gene_ko_mapping$KO_ID)
    message("输入基因对应的 KO 数量: ", length(ko_list))

    universe_ko <- NULL
    if (!is.null(universe_genes)) {
        universe_ko <- gene2ko %>%
            filter(Gene_Symbol %in% universe_genes) %>%
            pull(KO_ID) %>% unique()
    }

    enricher(gene = ko_list,
             pvalueCutoff = pvalue_cutoff,
             pAdjustMethod = "BH",
             qvalueCutoff = qvalue_cutoff,
             universe = universe_ko,
             minGSSize = 10,
             maxGSSize = 2000,
             TERM2GENE = ko2pathway %>% dplyr::select(Pathway_ID, KO_ID),
             TERM2NAME = ko2pathway %>% dplyr::select(Pathway_ID, Pathway_Description))
}


#' 一站式 KEGG 富集 + 气泡图导出 (top20 by p.adjust)
#' @param genes 数据框, 行名为基因符号 (例: FindMarkers 输出)
#' @param kegg_hierarchy KEGG 层级表
#' @param comb 主标题前缀
#' @param direc "Up"/"Down"/"Total"
#' @param workDir 输出目录
kegg_enrich_and_plot <- function(genes, kegg_hierarchy, comb, direc, workDir,
                                 universe_genes = NULL,
                                 pvalue_cutoff = 0.05, qvalue_cutoff = 0.2,
                                 top_n = 20) {
    gene2ko <- kegg_hierarchy %>%
        filter(!is.na(Gene_Symbol) & !is.na(KO_ID)) %>%
        dplyr::select(Gene_Symbol, KO_ID) %>% distinct()

    ko2pathway <- kegg_hierarchy %>%
        dplyr::select(KO_ID, Level_C) %>% distinct() %>%
        mutate(Pathway_ID = paste0("ko", str_extract(Level_C, "\\d+")),
               Pathway_Description = Level_C) %>%
        filter(!is.na(Pathway_ID))

    pathway_info <- kegg_hierarchy %>%
        mutate(Pathway_Num = str_extract(Level_C, "\\d+")) %>%
        filter(!is.na(Pathway_Num)) %>%
        dplyr::select(Pathway_Num, Level_A, Level_B, Level_C) %>% distinct()

    kegg_enrich <- perform_kegg_enrichment(
        genes = rownames(genes),
        gene2ko = gene2ko, ko2pathway = ko2pathway,
        universe_genes = universe_genes,
        pvalue_cutoff = pvalue_cutoff, qvalue_cutoff = qvalue_cutoff)
    if (is.null(kegg_enrich)) return(invisible(NULL))

    kegg_enrich_df <- as.data.frame(kegg_enrich) %>%
        mutate(Pathway_Num = str_extract(ID, "\\d+"), Pathway_ID = ID)

    timestamp <- format(Sys.time(), "%Y%m%d")
    new_string <- str_replace(comb, ": ", "_") %>% str_replace_all(" ", "-")
    combName <- paste(new_string, direc, sep = "_")

    if (nrow(kegg_enrich_df) == 0) {
        write.csv(kegg_enrich_df,
                  file.path(workDir, paste0(paste(combName, "KEGG", timestamp, sep = "_"), ".csv")),
                  row.names = FALSE)
        return(invisible(NULL))
    }

    kegg_enrich_df <- kegg_enrich_df %>%
        left_join(pathway_info, by = "Pathway_Num") %>%
        mutate(
            GeneRatio_num = as.numeric(str_split(GeneRatio, "/", simplify = TRUE)[, 1]) /
                            as.numeric(str_split(GeneRatio, "/", simplify = TRUE)[, 2]),
            BgRatio_num   = as.numeric(str_split(BgRatio, "/", simplify = TRUE)[, 1]) /
                            as.numeric(str_split(BgRatio, "/", simplify = TRUE)[, 2]),
            Fold_Enrichment = GeneRatio_num / BgRatio_num,
            Rich_Factor = Count / as.numeric(str_split(BgRatio, "/", simplify = TRUE)[, 2])
        ) %>%
        arrange(p.adjust) %>% mutate(Rank = row_number()) %>%
        mutate(Pathway_Name = sapply(Description, extract_pathway_name)) %>%
        mutate(LevelA = sapply(Level_A, extract_pathway_name)) %>%
        arrange(pvalue)

    write.csv(kegg_enrich_df[, !names(kegg_enrich_df) %in% colnames(kegg_enrich_df)[13:22]],
              file.path(workDir, paste0(paste(combName, "KEGG", timestamp, sep = "_"), ".csv")),
              row.names = FALSE)

    top_enrich_df <- kegg_enrich_df %>% arrange(p.adjust) %>% head(top_n)

    p <- ggplot(top_enrich_df,
                aes(x = GeneRatio_num, y = reorder(Pathway_Name, GeneRatio_num))) +
        geom_point(aes(size = Count, color = pvalue), alpha = 0.8) +
        scale_color_gradientn(colors = rev(brewer.pal(11, "Spectral")),
                              name = "P-value",
                              guide = guide_colorbar(reverse = TRUE, order = 1)) +
        scale_size_continuous(range = c(3, 8), name = "Gene Count",
                              breaks = scales::pretty_breaks(n = 4)(range(top_enrich_df$Count)),
                              guide = guide_legend(order = 2)) +
        theme_bw(base_size = 12) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
              axis.text.y = element_text(size = 10),
              axis.title = element_text(size = 12, face = "bold"),
              strip.text = element_text(size = 11, face = "bold", color = "white"),
              strip.background = element_rect(fill = "#2E86AB", color = "#2E86AB"),
              legend.position = "right",
              legend.box = "vertical",
              legend.title = element_text(face = "bold"),
              plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
              panel.grid.major = element_line(color = "gray90", size = 0.3),
              panel.grid.minor = element_blank()) +
        labs(title = paste("Top", top_n, "KEGG Pathways by Significance"),
             x = "Gene Ratio", y = "Pathway Description") +
        ggforce::facet_col(~LevelA, scales = "free_y", space = "free")

    fileName <- paste0("top", top_n, "_KEGG_", timestamp)
    ggsave(file.path(workDir, paste0(paste(combName, fileName, sep = "_"), ".pdf")),
           p, width = 14, height = 10, dpi = 300)
    ggsave(file.path(workDir, paste0(paste(combName, fileName, sep = "_"), ".png")),
           p, width = 14, height = 10, dpi = 300)
    invisible(kegg_enrich_df)
}
