packages <- c("optparse", "tidyverse", "ggpicrust2", "ggprism", "patchwork", "ALDEx2", "GGally", "openxlsx", "dplyr", "purrr")
suppressPackageStartupMessages(lapply(packages, library, character.only = TRUE))
rm(list = ls())

# 创建参数解析器
parser <- OptionParser(option_list = list(
  make_option(c("-w","--work_dir"), type = "character", default = NULL, help = "工作目录", metavar = "character"),
  make_option(c("-o","--output_dir"), type = "character", default = NULL, help = "输出目录", metavar = "character"),
  make_option(c("-s","--samples"), type = "character", default = NULL, help = "样本描述文件", metavar = "character"),
  make_option(c("-c","--comp"), type = "character", default = NULL, help = "比较信息文件，第一列为Treat，第二列为Control", metavar = "character")
))
opt <- parse_args(parser)

work_dir <- opt$work_dir
output_dir <- opt$output_dir
samples_file <- opt$samples
compare_file <- opt$comp

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 错误记录 error.txt -------------------------
log_plot_error <- function(plot_type, pair_name, extra_info = "") {
  msg <- paste0(
    "[", Sys.time(), "] ",
    "plot_type=", plot_type, "; ",
    "pair_name=", pair_name, "; ",
    extra_info, "\n\n"
  )
  cat(msg, file = file.path(output_dir, "error.txt"), append = TRUE)
}

# 通用KO差异分析处理函数（消除nonAdjust/Adjust重复代码）
process_ko_daa <- function(abundance_mat, sample_metadata, daa_method, p_adjust_method, adjust_label, pair_name, output_dir) {
  # 差异分析
  daa_df <- pathway_daa(
    abundance = abundance_mat,
    metadata = sample_metadata,
    group = "group",
    daa_method = daa_method,
    include_abundance_stats = TRUE,
    p.adjust = p_adjust_method
  )
  # 注释排序
  annotated_df <- pathway_annotation(
    pathway = "KO",
    daa_results_df = daa_df,
    ko_to_kegg = TRUE
  ) %>% arrange(p_values) %>% distinct(feature, .keep_all = TRUE)# 按p_values排序
  # 提取top10上下调
  top10_pos <- annotated_df %>% filter(log2_fold_change > 0) %>% slice_max(order_by = log2_fold_change, n = 10)
  top10_neg <- annotated_df %>% filter(log2_fold_change < 0) %>% slice_min(order_by = log2_fold_change, n = 10)
  top10_errorbar_df <- bind_rows(top10_pos, top10_neg)
  top10_errorbar_df$pathway_name <- paste0(top10_errorbar_df$feature, '_', top10_errorbar_df$pathway_name, '     ')
  # 提取top20上下调(显著)
  top20_pos <- annotated_df %>% filter(log2_fold_change > 0, p_adjust < 0.05) %>% slice_max(order_by = log2_fold_change, n = 20)
  top20_neg <- annotated_df %>% filter(log2_fold_change < 0, p_adjust < 0.05) %>% slice_min(order_by = log2_fold_change, n = 20)
  top20_heatmap_df <- bind_rows(top20_pos, top20_neg)

  # 检查是否有显著通路
  if (nrow(top20_heatmap_df) == 0) {
    warning(paste("No significant KEGG pathways found for", pair_name, adjust_label))
    # 保存结果时top20 sheet为空
    wb <- createWorkbook()
    addWorksheet(wb, "All")
    addWorksheet(wb, "Top10")
    addWorksheet(wb, "Top20")
    writeData(wb, "All", annotated_df)
    writeData(wb, "Top10", top10_errorbar_df)
    writeData(wb, "Top20", data.frame())
    saveWorkbook(wb, file.path(output_dir, paste0("KEGG_daa_", adjust_label, "_annotated_", pair_name, ".xlsx")), overwrite = TRUE)
    return(list(annotated = annotated_df, top10 = top10_errorbar_df, top20 = data.frame()))
  }

  # 合并丰度矩阵和top20信息，生成热图用矩阵
  top20_heatmap_df <- abundance_mat %>% rownames_to_column("feature") %>%
    inner_join(top20_heatmap_df %>% select(feature, pathway_name), by = "feature") %>%
    mutate(feature = paste(feature, pathway_name, sep = "_")) %>%
    .[, union("feature", sample_metadata$sample)] %>%
    column_to_rownames("feature")

  # 保存结果（分sheet保存All和Top10）
  wb <- createWorkbook()
  addWorksheet(wb, "All")
  addWorksheet(wb, "Top10")
  addWorksheet(wb, "Top20")
  writeData(wb, "All", annotated_df)
  writeData(wb, "Top10", top10_errorbar_df)
  top20_save_df <- data.frame(pathway = rownames(top20_heatmap_df), top20_heatmap_df, check.names = FALSE)
  writeData(wb, "Top20", top20_save_df, rowNames = FALSE)
  saveWorkbook(wb, file.path(output_dir, paste0("KEGG_daa_", adjust_label, "_annotated_", pair_name, ".xlsx")), overwrite = TRUE)
  return(list(annotated = annotated_df, top10 = top10_errorbar_df, top20 = top20_heatmap_df))
}

# 通用KO误差条形图绘制函数
plot_ko_errorbar <- function(abundance_mat, daa_results_df, group_info, adjust_label, pair_name, output_dir) {
  tryCatch(
    {
      plot_obj <- pathway_errorbar(
        abundance = abundance_mat,
        daa_results_df = daa_results_df,
        Group = group_info,
        p_values_threshold = 0.05,
        legend_text_size = 12,
        pvalue_size = 4,
        x_lab = "pathway_name"
      ) +
        scale_x_discrete(labels = scales :: label_wrap(50))
      ggsave(
        filename = file.path(output_dir, paste0("KEGG_", adjust_label, "_pathwayErrbar_", pair_name, ".jpg")),
        plot = plot_obj,
        dpi = 320,
        width = 12,
        height = 8
      )
    },
    error = function(e) {
      log_plot_error(
        plot_type = paste0("KEGG_pathway_errorbar_", adjust_label),
        pair_name = pair_name,
        extra_info = paste("error:", conditionMessage(e))
      )
    }
  )
}

# EC差异分析处理函数
process_ec_daa <- function(abundance_mat, sample_metadata, p_adjust_method, adjust_label, pair_name, output_dir) {
  # 差异分析
  daa_df <- pathway_daa(
    abundance = abundance_mat,
    metadata = sample_metadata,
    group = "group",
    daa_method = "LinDA",
    include_abundance_stats = TRUE,
    p.adjust = p_adjust_method
  )
  # 注释排序
  annotated_df <- pathway_annotation(
    pathway = "EC",
    daa_results_df = daa_df,
    ko_to_kegg = FALSE
  ) %>% arrange(p_adjust)
  # 提取top10上下调
  top10_pos <- annotated_df %>% filter(log2FoldChange > 0) %>% slice_max(order_by = log2FoldChange, n = 10)
  top10_neg <- annotated_df %>% filter(log2FoldChange < 0) %>% slice_min(order_by = log2FoldChange, n = 10)
  top10_errorbar_df <- bind_rows(top10_pos, top10_neg)
  top10_errorbar_df$description <- paste0(top10_errorbar_df$feature, '_', top10_errorbar_df$description, '     ')
  # 提取top20上下调(显著)，画图用
  top20_pos <- annotated_df %>% filter(log2FoldChange > 0, p_adjust < 0.05) %>% slice_max(order_by = log2FoldChange, n = 20)
  top20_neg <- annotated_df %>% filter(log2FoldChange < 0, p_adjust < 0.05) %>% slice_min(order_by = log2FoldChange, n = 20)
  top20_heatmap_df <- bind_rows(top20_pos, top20_neg)
  
  # 检查是否有显著通路
  if (nrow(top20_heatmap_df) == 0) {
    warning(paste("No significant EC pathways found for", pair_name, adjust_label))
    # 保存结果时top20 sheet为空
    wb <- createWorkbook()
    addWorksheet(wb, "All")
    addWorksheet(wb, "Top10")
    addWorksheet(wb, "Top20")
    writeData(wb, "All", annotated_df)
    writeData(wb, "Top10", top10_errorbar_df)
    writeData(wb, "Top20", data.frame())
    saveWorkbook(wb, file.path(output_dir, paste0("EC_daa_", adjust_label, "_annotated_", pair_name, ".xlsx")), overwrite = TRUE)
    return(list(annotated = annotated_df, top10 = top10_errorbar_df, top20 = data.frame()))
  }
  
  # 合并丰度矩阵和top20信息，生成热图用矩阵
  top20_heatmap_df <- abundance_mat %>% rownames_to_column("feature") %>%
    inner_join(top20_heatmap_df %>% select(feature, description), by = "feature") %>%
    mutate(feature = paste(feature, description, sep = "_")) %>%
    .[, union("feature", sample_metadata$sample)] %>%
    column_to_rownames("feature")
  # 保存结果
  wb <- createWorkbook()
  addWorksheet(wb, "All")
  addWorksheet(wb, "Top10")
  addWorksheet(wb, "Top20")
  writeData(wb, "All", annotated_df)
  writeData(wb, "Top10", top10_errorbar_df)
  top20_save_df <- data.frame(pathway = rownames(top20_heatmap_df), top20_heatmap_df, check.names = FALSE)
  writeData(wb, "Top20", top20_save_df, rowNames = FALSE)
  saveWorkbook(wb, file.path(output_dir, paste0("EC_daa_", adjust_label, "_annotated_", pair_name, ".xlsx")), overwrite = TRUE)
  return(list(annotated = annotated_df, top10 = top10_errorbar_df, top20 = top20_heatmap_df))
}

# 通用EC误差条形图绘制函数
plot_ec_errorbar <- function(abundance_mat, top10_df, group_info, adjust_label, pair_name, output_dir) {
  tryCatch(
    {
      top10_df$description <- paste0(top10_df$description, '   ')
      plot_obj <- pathway_errorbar(
        abundance = abundance_mat,
        daa_results_df = top10_df,
        Group = group_info,
        p_values_threshold = 0.01,
        legend_text_size = 12,
        pvalue_size = 4,
        max_features = 5,
        x_lab = "description"
      ) +
        scale_x_discrete(labels = scales :: label_wrap(50))
      ggsave(
        filename = file.path(output_dir, paste0("EC_", adjust_label, "_pathwayErrbar_", pair_name, ".jpg")),
        plot = plot_obj,
        dpi = 320,
        width = 12,
        height = 8
      )
    },
    error = function(e) {
      log_plot_error(
        plot_type = paste0("EC_pathway_errorbar_", adjust_label),
        pair_name = pair_name,
        extra_info = paste("error:", conditionMessage(e))
      )
    }
  )
}

# 通用通路热图绘制函数（支持EC/KO/MetaCyc等所有类型）
plot_pathway_heatmap <- function(heatmap_data, sample_metadata, group_col, data_type, adjust_label, pair_name, output_dir) {
  # 检查数据行数，小于2则无法绘制热图
  if (nrow(heatmap_data) < 2) {
    log_plot_error(
      plot_type = paste0(data_type, "_pathway_heatmap_", adjust_label),
      pair_name = pair_name,
      extra_info = "Insufficient data: less than 2 features to plot heatmap"
    )
    return(NULL)
  }
  
  tryCatch(
    {
      plot_obj <- pathway_heatmap(heatmap_data, sample_metadata, group_col) +
        ggtitle(pair_name) +
        scale_y_discrete(labels = scales :: label_wrap(50))
      ggsave(
        filename = file.path(output_dir, paste0(data_type, "_", adjust_label, "_heatmap_", pair_name, ".jpg")),
        plot = plot_obj,
        dpi = 320,
        width = 12,
        height = 14
      )
    },
    error = function(e) {
      log_plot_error(
        plot_type = paste0(data_type, "_pathway_heatmap_", adjust_label),
        pair_name = pair_name,
        extra_info = paste("error:", conditionMessage(e))
      )
      return(NULL)
    }
  )
  return(plot_obj)
}


global_sample_info <- read_delim(samples_file, delim = "\t", escape_double = FALSE, trim_ws = TRUE)
##合并和用 ‘KEGG_daa_df_test_20260122.xlsx’ 中的 p.adjust 列替代 ‘KEGG_daa_noAdjust_annotated_df_test_20260122.xlsx’ 中的 p.adjust 列并生成新文件，对 EC 和 metaCyc 也同样处理！！！

## 组间
# 创建 pairs list
metadata_pairs_list <- list()
if (!is.null(compare_file)) {
  # 读取比较信息文件，文件自带列名Treat和Control
  comp_df <- read_delim(compare_file, delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  for (i in 1:nrow(comp_df)) {
    treat <- comp_df$Treat[i]
    control <- comp_df$Control[i]
    # 按照control在前，treat在后的顺序排序
    df_sub <- global_sample_info %>% 
      filter(group %in% c(control, treat)) %>%
      mutate(group = factor(group, levels = c(control, treat))) %>%
      arrange(group)
    pair_name <- paste(control, treat, sep = "-vs-")
    metadata_pairs_list[[pair_name]] <- df_sub
  }
  comb_list <- lapply(1:nrow(comp_df), function(i) c(comp_df$Control[i], comp_df$Treat[i]))
} else {
  # 获取所有分组标签
  group_levels <- unique(global_sample_info[['group']])
  # 生成所有两两组合的分组
  comb_list <- combn(group_levels, 2, simplify = FALSE)
  # 生成metadata_pairs_list
  for (i in seq_along(comb_list)) {
    groups <- comb_list[[i]]
    df_sub <- global_sample_info %>% filter(group %in% groups)
    pair_name <- paste(groups, collapse = "-vs-")
    metadata_pairs_list[[pair_name]] <- df_sub
  }
}

## KO
KO_abundance_file <- file.path(work_dir, "KO_metagenome_out/pred_metagenome_unstrat.tsv")
ko_abundance_df <- read_delim(KO_abundance_file, delim = "\t", escape_double = FALSE, trim_ws = TRUE)
for (i in seq_along(comb_list)) {
  groups        <- comb_list[[i]]
  pair_name     <- paste(groups, collapse = "-vs-")
  pair_metadata <- metadata_pairs_list[[pair_name]]
  message("[", Sys.time(), "] 开始分析KO通路: ", pair_name, " (", i, "/", length(comb_list), ")")

  ## 1. 为当前分组对提取丰度子矩阵（保留第一列 function）
  sample_cols_idx   <- which(colnames(ko_abundance_df) %in% pair_metadata$sample)
  ko_abundance_pair <- ko_abundance_df[, c(1, sample_cols_idx)]

  ## 2. KO -> KEGG 丰度表
  ko2kegg_abundance_pair <- ko2kegg_abundance(data = ko_abundance_pair)

  ## 3. 构建 DAA 所需矩阵（行名为 KO ID / function 列）好像没用，加了个列名
  # ko_abundance_mat <- ko2kegg_abundance_pair %>% column_to_rownames("function")

  ## 4. 处理未校正和校正p值两种情况
  ko_res_nonAdjust <- process_ko_daa(
    abundance_mat = ko2kegg_abundance_pair,
    sample_metadata = pair_metadata,
    daa_method = "LinDA",
    p_adjust_method = "none",
    adjust_label = "nonAdjust",
    pair_name = pair_name,
    output_dir = output_dir
  )
  ko_res_Adjust <- process_ko_daa(
    abundance_mat = ko2kegg_abundance_pair,
    sample_metadata = pair_metadata,
    daa_method = "LinDA",
    p_adjust_method = "BH",
    adjust_label = "Adjust",
    pair_name = pair_name,
    output_dir = output_dir
  )

  ## 5. 绘制校正/未校正两套KEGG误差条形图
  plot_ko_errorbar(ko2kegg_abundance_pair, ko_res_Adjust$top10, pair_metadata$group, "Adjust", pair_name, output_dir)
  plot_ko_errorbar(ko2kegg_abundance_pair, ko_res_nonAdjust$top10, pair_metadata$group, "nonAdjust", pair_name, output_dir)

  ## 6. 绘制校正/未校正两套KEGG热图
  ko_pathway_Top20_heatmap_Adjust <- plot_pathway_heatmap(
    ko_res_nonAdjust$top20,
    pair_metadata,
    "group",
    "KEGG",
    "Adjust",
    pair_name,
    output_dir
  )
  ko_pathway_Top20_heatmap_nonAdjust <- plot_pathway_heatmap(
    ko_res_nonAdjust$top20,
    pair_metadata,
    "group",
    "KEGG",
    "nonAdjust",
    pair_name,
    output_dir
  )
}

## EC
EC_abundance_file <- file.path(work_dir, "EC_metagenome_out/pred_metagenome_unstrat.tsv")
EC_abundance <- read_tsv(EC_abundance_file, trim_ws = TRUE)

for (i in seq_along(comb_list)) {
  ## 1. 准备当前分组对的数据
  groups    <- comb_list[[i]]
  pair_name <- paste(groups, collapse = "-vs-")
  pair_metadata <- metadata_pairs_list[[pair_name]]
  message("[", Sys.time(), "] 开始分析EC酶: ", pair_name, " (", i, "/", length(comb_list), ")")

  sample_cols_idx    <- which(colnames(EC_abundance) %in% pair_metadata$sample)
  EC_abundance_pair  <- EC_abundance[, c(1, sample_cols_idx)]
  EC_abundance_mat   <- EC_abundance_pair %>% column_to_rownames("function")
  
  # pathway_heatmap 使用
  EC_example <- data.frame(
    sample_name = pair_metadata$sample,
    group       = pair_metadata$group,
    batch       = factor(rep("Batch1", length(pair_metadata$sample)))
  )

  ## 2. 批量处理未校正和校正p值两种情况（调用通用函数消除重复代码）
  res_nonAdjust <- process_ec_daa(EC_abundance_mat, pair_metadata, "none", "nonAdjust", pair_name, output_dir)
  res_Adjust <- process_ec_daa(EC_abundance_mat, pair_metadata, "BH", "Adjust", pair_name, output_dir)

  ## 3. 绘制校正/未校正两套误差条形图（调用通用绘图函数）
  plot_ec_errorbar(EC_abundance_mat, res_Adjust$top10, pair_metadata$group, "Adjust", pair_name, output_dir)
  plot_ec_errorbar(EC_abundance_mat, res_nonAdjust$top10, pair_metadata$group, "nonAdjust", pair_name, output_dir)
  
  ## 4. 绘制校正/未校正两套EC热图（调用通用热图函数）
  EC_pathway_Top20_heatmap_Adjust <- plot_pathway_heatmap(
    res_Adjust$top20,
    EC_example,
    "group",
    "EC",
    "Adjust",
    pair_name,
    output_dir
  )
  EC_pathway_Top20_heatmap_nonAdjust <- plot_pathway_heatmap(
    res_nonAdjust$top20,
    EC_example,
    "group",
    "EC",
    "nonAdjust",
    pair_name,
    output_dir
  )
}


## PCA（使用与热图一致的分组信息）
EC_example <- data.frame(
  sample_name = global_sample_info$sample,
  group       = global_sample_info$group,
  batch       = factor(rep("Batch1", length(global_sample_info$sample)))
)
EC_abundance_mat <- EC_abundance %>% column_to_rownames("function")
EC_desc_PCA <- pathway_pca(EC_abundance_mat, EC_example, "group")
ggsave(
  filename = file.path(output_dir, "EC_pca.jpg"),
  plot     = EC_desc_PCA,
  dpi      = 320,
  width    = 12,
  height   = 8
)

## 以下为 GSEA 分析代码（对所有分组对循环） ----------------------------------

## 1. 准备 GSEA 所需的 metadata 与 KO 丰度矩阵 -------------------------------
metadata_gsea_global <- global_sample_info %>% as.data.frame() %>% mutate(sample_bak = sample) %>% column_to_rownames("sample_bak")
ko_abundance_df <- ko_abundance_df %>% as.data.frame() %>% mutate(function_bak = make.unique(as.character(`function`))) %>% column_to_rownames("function_bak")
## 为 GSEA 准备 EC 丰度矩阵（行为 EC，列为样本）
EC_abundance_gsea <- EC_abundance %>% as.data.frame() %>% column_to_rownames("function")
ko_ref_df <- read.delim('/home/colddata/qinqiang/script/lib/KO_id_pathway_name.txt', sep='\t', header = FALSE, col.names = c("pathway_id", "pathway_name"), stringsAsFactors = FALSE)
## 2. 对所有分组对进行 KEGG / MetaCyc GSEA -----------------------------------
## MetaCyc
MetaCyc_abundance_file <- file.path(work_dir, "pathways_out/path_abun_unstrat.tsv")
MetaCyc_abundance <- read_delim(file =  MetaCyc_abundance_file, delim = "\t", escape_double = FALSE, trim_ws = TRUE)
data("ko_to_go_reference")
colnames(ko_to_go_reference)[1] <- 'pathway_id'
colnames(ko_to_go_reference)[2] <- 'pathway_name'
for (i in seq_along(comb_list)) {
  groups    <- comb_list[[i]]
  pair_name <- paste(groups, collapse = "-vs-")
  message("[", Sys.time(), "] 开始GSEA分析: ", pair_name, " (", i, "/", length(comb_list), ")")

  ## 当前分组对的 metadata
  metadata_gsea <- metadata_gsea_global %>% filter(group %in% groups)

  ## 当前分组对的 KO 丰度子矩阵
  sample_cols_idx_ko   <- which(colnames(ko_abundance_df) %in% metadata_gsea$sample)
  ko_abundance_pair    <- ko_abundance_df[, sample_cols_idx_ko, drop = FALSE]

  ## 2.1 KEGG GSEA（错误时记录并跳过）
  tryCatch(
    {
      kegg_gseaResults <- pathway_gsea(
        abundance    = ko_abundance_pair,
        metadata     = metadata_gsea,
        group        = "group",
        method       = "camera",
        pathway_type = "KEGG",
        rank_method = "signal2noise"
      )
      
      annotated_df <- kegg_gseaResults %>% select(-pathway_name) %>%
        left_join(ko_ref_df %>% select(pathway_id, pathway_name), by = "pathway_id")
      top10_pos <- annotated_df %>% filter(NES > 0, pvalue < 0.05) %>% slice_max(order_by = NES, n = 10)
      top10_neg <- annotated_df %>% filter(NES < 0, pvalue < 0.05) %>% slice_min(order_by = NES, n = 10)
      top10_df <- bind_rows(top10_pos, top10_neg)
      top10_df$pathway_name <- paste0(top10_df$pathway_id, '_', top10_df$pathway_name, '     ')
      
      plot_kegg_gsea <- visualize_gsea(
        gsea_results = top10_df, 
        plot_type = "barplot",
        metadata = metadata_gsea,
        group = 'group'
      ) +
        ggtitle(pair_name) +
        scale_x_discrete(labels = scales :: label_wrap(50))
      ggsave(
        filename = file.path(output_dir, paste0("GSEA_KEGG_Barplot_", pair_name, ".jpg")),
        plot     = plot_kegg_gsea,
        dpi      = 320,
        width    = 12,
        height   = 8
      )
      
      wb <- createWorkbook()
      addWorksheet(wb, "All")
      addWorksheet(wb, "Top10")
      writeData(wb, "All", annotated_df)
      writeData(wb, "Top10", top10_df)
      saveWorkbook(wb, file.path(output_dir, paste0("GSEA_KEGG_", pair_name, ".xlsx")), overwrite = TRUE)
    },
    error = function(e) {
      log_plot_error(
        plot_type = "KEGG_GSEA",
        pair_name = pair_name,
        extra_info = paste("error:", conditionMessage(e))
      )
    }
  )
  
  tryCatch(
    {
      GO_gsea_results <- pathway_gsea(
        abundance = ko_abundance_pair,
        metadata = metadata_gsea,
        group = "group",
        method = "camera",
        pathway_type = "GO",
        go_category = "MF",
        rank_method = "signal2noise"
      )
      GO_gsea_annotated_df <- GO_gsea_results %>% select(-pathway_name) %>%
        left_join(ko_to_go_reference %>% select(pathway_id, pathway_name), by = "pathway_id")

      GO_gsea_top10_pos <- GO_gsea_annotated_df %>% filter(NES > 0, pvalue < 0.05) %>% slice_max(order_by = NES, n = 10)
      GO_gsea_top10_neg <- GO_gsea_annotated_df %>% filter(NES < 0, pvalue < 0.05) %>% slice_min(order_by = NES, n = 10)
      GO_gsea_top10_df <- bind_rows(GO_gsea_top10_pos, GO_gsea_top10_neg)
      GO_gsea_top10_df$pathway_name <- paste0(GO_gsea_top10_df$pathway_id, '_', GO_gsea_top10_df$pathway_name, '     ')
      
      plot_go_gsea <- visualize_gsea(
        gsea_results = GO_gsea_top10_df, 
        plot_type = "barplot",
        metadata = metadata_gsea,
        group = 'group'
      ) +
        ggtitle(pair_name) +
        scale_x_discrete(labels = scales :: label_wrap(50))
      ggsave(
        filename = file.path(output_dir, paste0("GSEA_GO_Barplot_", pair_name, ".jpg")),
        plot     = plot_go_gsea,
        dpi      = 320,
        width    = 12,
        height   = 8
      )
      
      wb <- createWorkbook()
      addWorksheet(wb, "All")
      addWorksheet(wb, "Top10")
      writeData(wb, "All", GO_gsea_annotated_df)
      writeData(wb, "Top10", GO_gsea_top10_df)
      saveWorkbook(wb, file.path(output_dir, paste0("GSEA_GO_", pair_name, ".xlsx")), overwrite = TRUE)
    },
    error = function(e) {
      log_plot_error(
        plot_type = "GO_GSEA",
        pair_name = pair_name,
        extra_info = paste("error:", conditionMessage(e))
      )
    }
  )
  

  ## 当前分组对的 EC 丰度子矩阵（用于 MetaCyc GSEA）
  sample_cols_idx_ec   <- which(colnames(EC_abundance_gsea) %in% metadata_gsea$sample)
  EC_abundance_pair_gs <- EC_abundance_gsea[, sample_cols_idx_ec, drop = FALSE]

  ## 2.2 MetaCyc GSEA（错误时记录并跳过）
  tryCatch(
    {
      MetaCyc_gseaResults <- pathway_gsea(
        abundance    = EC_abundance_pair_gs,
        metadata     = metadata_gsea,
        group        = "group",
        method       = "camera",
        pathway_type = "MetaCyc"
      )
      colnames(MetaCyc_gseaResults)[1] <- "feature"
      annotated_df <- pathway_annotation(
        pathway = "MetaCyc",
        daa_results_df = MetaCyc_gseaResults,
      )
      colnames(annotated_df)[1] <- "pathway_id"
      
      top10_pos <- annotated_df %>% filter(NES > 0, pvalue < 0.05) %>% slice_max(order_by = NES, n = 10)
      top10_neg <- annotated_df %>% filter(NES < 0, pvalue < 0.05) %>% slice_min(order_by = NES, n = 10)
      top10_df <- bind_rows(top10_pos, top10_neg)
      top10_df$pathway_name <- paste0(top10_df$pathway_name, '_', top10_df$description, '     ')
      
      plot_metacyc_gsea <- visualize_gsea(
        gsea_results = top10_df,
        plot_type = "barplot"
      ) +
        ggtitle(pair_name) +
        scale_x_discrete(labels = scales :: label_wrap(50))

      ggsave(
        filename = file.path(output_dir, paste0("GSEA_MeffBarplot_", pair_name, ".jpg")),
        plot     = plot_metacyc_gsea,
        dpi      = 320,
        width    = 12,
        height   = 8
      )
      
      wb <- createWorkbook()
      addWorksheet(wb, "All")
      addWorksheet(wb, "Top10")
      writeData(wb, "All", annotated_df)
      writeData(wb, "Top10", top10_df)
      saveWorkbook(wb, file.path(output_dir, paste0("gsea_MetaCyc_", pair_name, ".xlsx")), overwrite = TRUE)
    },
    error = function(e) {
      log_plot_error(
        plot_type = "MetaCyc_GSEA",
        pair_name = pair_name,
        extra_info = paste("error:", conditionMessage(e))
      )
    }
  )
}

# extdata_dir <- file.path(find.package("ggpicrust2"), "extdata")
# if (!dir.exists(extdata_dir)) dir.create(extdata_dir, recursive = TRUE)
# extdata_dir
# # [1] "/home/data/opt/biosoft/R-4.4.3/build/lib64/R/library/ggpicrust2/extdata"
# download.file("https://github.com/cafferychen777/ggpicrust2/blob/main/data/metacyc_to_ec_reference.rda", file.path(extdata_dir, "metacyc_to_ec_reference.RData"))