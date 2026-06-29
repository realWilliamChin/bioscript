library(microeco)
library(openxlsx)
library(ggplot2)
library(optparse)
library(patchwork)
library(ggpubr)
rm(list=ls())

my_set_colors <- c(
  "#e91e63", "#3f51b5", "#dccc39", "#4caf50", "#607d8b",
  "#00bcd4", "#009688", "#8bc34a", "#f44336", "#9c27b0",
  "#ffeb3b", "#ffc107", "#ff5722", "#795548", "#9e9e9e",
  "#c8bfe7", "#b97a57", "#ffaec9", "#1ee1c4", "#2196f3"
)

option_list <- list(
  # 输入文件
  make_option(c("-o", "--otu"), type = "character", default = "Species_aboun.txt",
    help = "OTU/feature 矩阵文件路径，默认 asv_table.txt"),
  make_option(c("-s", "--sample"), type = "character", default = "samples_described.txt",
    help = "样本描述表文件路径，默认 samples_described.txt"),
  make_option(c("-t", "--tax"), type = "character", default = "Species_taxon.txt",
    help = "物种注释表文件路径，默认 asv_taxon.txt"),

  # 过滤相关
  make_option(c("-p", "--pollution"), type = "character", default = "Eukaryota,chloroplast,mitochondria",
    help = "污染物分类（逗号分隔），传递给 filter_pollution，默认值为 Eukaryota,chloroplast,mitochondria"),
  make_option(c("-r", "--rel_abund"), type = "double", default = 0.0001,
    help = "filter_taxa 的 rel_abund 阈值，默认 0.0001"),
  make_option(c("-f", "--freq"), type = "double", default = 0.1,
    help = "filter_taxa 的 freq 阈值，默认 0.1"),
  make_option(c("--outdir"), type = "character", default = "./",
    help = "输出目录，默认当前目录，会自动创建子目录")
)

# trans_abund 画图函数，包含多个分类水平的丰度柱状图与分组均值图
abundance_barplots <- function(dataset, tax_level, ntaxa = 10, output_dir) {
  abund_res <- trans_abund$new(dataset = dataset, taxrank = tax_level, ntaxa = ntaxa)
  abund_res_groupmean <- trans_abund$new(dataset = dataset, taxrank = tax_level, ntaxa = ntaxa, groupmean = "group")

  bar_plot <- abund_res$plot_bar(
    legend_text_italic = TRUE, color_values = my_set_colors, others_color = "darkgrey") +
    theme_pubr(
      base_size = 12, base_family = "Arial", border = FALSE, margin = TRUE,
      legend = 'right', x.text.angle = 45
    )
  #ggsave(filename = file.path(output_dir, paste0("abundance_bar_", tax_level, ".jpg")), plot = bar_plot, dpi = 320, width = 12, height = 8)
  # p_arial <- bar_plot + theme_pubr(
  #  text = element_text(family = "Ravie"), # 设置所有文本（包括坐标轴、图例）为Arial
  #  axis.title = element_text(family = "Ravie"), # 坐标轴标题
  #  axis.text = element_text(family = "Ravie"), # 坐标轴刻度标签
  #  legend.title = element_text(family = "Ravie"), # 图例标题
  #  legend.text = element_text(family = "Ravie") # 图例标签
  # )
  ggsave(filename = file.path(output_dir, paste0("abundance_bar_Ravie_", tax_level, ".jpg")), plot = bar_plot, dpi = 320, width = 14, height = 8)
  
  bar_plot_group <- abund_res_groupmean$plot_bar(
    legend_text_italic = TRUE, color_values = my_set_colors, others_color = "darkgrey") +
    theme_pubr(
      base_size = 12, base_family = "Arial", border = FALSE, margin = TRUE,
      legend = "right", x.text.angle = 45
    )
  ggsave(filename = file.path(output_dir, paste0("abundance_bar_groupmean_", tax_level, ".jpg")), plot = bar_plot_group, dpi = 320, width = 6, height = 8)
  write.xlsx(bar_plot_group$data, file = file.path(output_dir, paste0("abundance_bar_groupmean_", tax_level, ".xlsx")), sheetName = "Sheet1", rowNames = TRUE)
  
  abund_res_donut <- trans_abund$new(dataset = dataset, taxrank = tax_level, ntaxa = 5, groupmean = "group")
  #! donut graph by taxonomy 
  donut_plot <- abund_res_donut$plot_donut(
    xtext_keep = FALSE, legend_text_italic = TRUE, color_values = my_set_colors, others_color = "darkgrey") +
    theme_pubr(
      base_size = 12, base_family = "", border = TRUE, margin = TRUE,
      legend = "bottom") +
    geom_text(
      stat = 'identity',
      position = position_stack(vjust = 0.5),
      check_overlap = TRUE,
      size = 3
    )
  #donut_plot <- donut_plot + theme_pubr()
  ggsave(filename = file.path(output_dir, paste0("abundance_donut_", tax_level, ".jpg")), plot = donut_plot, dpi = 320, width = 12, height = 8)
  write.xlsx(bar_plot$plot_env$new_data, file = file.path(output_dir, paste0("abundance_bar_", tax_level, ".xlsx")), sheetName = "Sheet1", rowNames = TRUE)
  
  # 丰度组间数理分析
  # 判断组数, 如果样本组数大于 2 用 KW_dunn，2组 用 wilcox
  group_number <- length(unique(dataset$sample_table$group))
  if (group_number > 2) {method = "KW_dunn"} else {method = "wilcox"}
  abund_diff_obj <- trans_diff$new(
    dataset = dataset,
    method = "wilcox",
    group = "group",
    alpha = 0.5,
    lefse_subgroup = NULL,
    p_adjust_method = "none",
    taxa_level = 'all'
  )
  write.xlsx(abund_diff_obj$res_diff, file = file.path(output_dir, paste0("abundance_stats_", tax_level, ".xlsx")), sheetName = "Sheet1", rowNames = FALSE)
}

taxonomy_venn_plot <- function(dataset, tax_level, output_dir) {
  tmp <- trans_norm$new(dataset)
  mt_rarefied <- tmp$norm(method = "rarefy")
  mt_rarefied <- mt_rarefied$merge_samples("group")

  # Venn 构建和绘图过程中如果出现数据结构问题（如内部有 NULL），捕获错误并跳过该 level
  tryCatch({
    venn_project <- trans_venn$new(mt_rarefied, ratio = "seqratio")
    venn_project$plot_venn(
      petal_center_size = 50,
      petal_r = 1.5,
      petal_a = 3,
      petal_move_xy = 3.8,
      petal_color_center = "#BEBADA",
      color_circle = my_set_colors
    )
    #: TODO！不同组需要调整大小
    ggsave(file.path(output_dir, paste0("venn_", tax_level, ".jpg")), dpi = 320, width = 12, height = 10)

    venn_data_details <- t(venn_project$data_details)
    write.xlsx(venn_data_details, file = file.path(output_dir, paste0("venn_", tax_level, "_data_details.xlsx")), sheetName = "Sheet1", rowNames = TRUE, colNames = FALSE)
  }, error = function(e) {
    message(paste0(tax_level, " 当前分类水平的 Venn 分析出错，已跳过：", conditionMessage(e)))
    return(invisible(NULL))
  })
}

# 画样品分类进化树，只做 Genus 结果
tree_plot <- function(dataset, output_dir) {
  dataset$tidy_dataset()
  beta_res <- trans_beta$new(dataset, group = "group", measure = "bray")
  # use replace_name to set the label name, group parameter used to set the color
  tree_plot <- beta_res$plot_clustering(group = "group", replace_name = c("sample"), color_values = my_set_colors)
  ggsave(file.path(output_dir, "Genus_tree.jpg"), tree_plot, dpi = 320, width = 12, height = 10)
}

# 树和丰度柱状图结合，只画 Phylum 结果
tree_bar_combine_plot <- function(dataset, tax_level, ntaxa, output_dir) {
  abund_res <- trans_abund$new(dataset = dataset, taxrank = tax_level, ntaxa = ntaxa)
  beta_res <- trans_beta$new(dataset, group = "group", measure = "bray")
  tree_plot <- beta_res$plot_clustering(group = "group", replace_name = c("sample"), color_values = my_set_colors)

  tree_order <- tree_plot$plot_env$data2$sample
  bar_plot <- abund_res$plot_bar(
    xtext_keep = FALSE, legend_text_italic = TRUE, color_values = my_set_colors,
    coord_flip = TRUE, order_x = tree_order, others_color = "darkgrey")
  
  p_final = tree_plot + bar_plot + plot_layout(widths = c(2, 2))
  ggsave(file.path(output_dir, paste0(tax_level, "_tree_bar_combine.jpg")), p_final, dpi = 320, width = 15, height = 10)
}

# alpha-diversity 分析
alaha_analysis <- function(dataset, output_dir) {
  # alpha_stats 输出表
  # t1_single <- trans_alpha$new(dataset = dataset, group = "sample")
  # t1_group <- trans_alpha$new(dataset = dataset, group = "group")
  # write.xlsx(t1_single$data_stat, file = file.path(output_dir, paste0("alpha_stats_single.xlsx")), sheetName = "Sheet1", rowNames = FALSE)
  # write.xlsx(t1_group$data_stat, file = file.path(output_dir, paste0("alpha_stats_group.xlsx")), sheetName = "Sheet1", rowNames = FALSE)
  alpha_res <- trans_alpha$new(dataset = dataset, group = "group")
  alpha_res$cal_diff(method = "KW")
  # 在某些情况下（如结果结构缺少 Method 列）KW_dunn 可能会报错，这里做容错处理
  tryCatch({
    alpha_res$cal_diff(method = "KW_dunn", KW_dunn_letter = FALSE)
  }, error = function(e) {
    message("KW_dunn 多重比较计算出错，已跳过：", conditionMessage(e))
  })
  
  alpha_measures <- c("Chao1", "ACE", "Simpson", "Shannon", "InvSimpson", "Pielou", "Fisher", "Coverage")
  for (m in alpha_measures) {
    p <- tryCatch(
      alpha_res$plot_alpha(measure = m, add_sig_text_size = 6, y_increase = 0.25, add = "jitter", order_x_mean = TRUE) +
        theme_pubr(
          base_size = 12, base_family = "Arial", border = TRUE, margin = TRUE,
          legend = 'right', x.text.angle = 45
        ),
      error = function(e) {
        message(paste0(" alpha 指标 ", m, " 绘图失败，已跳过：", conditionMessage(e)))
      }
    )
    if (!is.null(p)) {
      ggsave(file.path(output_dir, paste0(m, '_bar.jpg')), p, dpi = 320, width = 8, height = 10)
    }
  }

  #输出summary结果
  write.xlsx(alpha_res$res_diff, file = file.path(output_dir, "alpha_summary_result.xlsx"), sheetName = "Sheet1", rowNames = FALSE)
  p <- tryCatch(
    alpha_res$plot_alpha(plot_type = "ggviolin", alpha = 0.2, y_increase = 0.4, add = "mean_se", add_sig_text_size = 6),
    error = function(e) {
      message("alpha 小提琴图绘图失败，已跳过：", conditionMessage(e))
    }
  )
  if (!is.null(p)) {
    ggsave(file.path(output_dir, "alpha_violin.jpg"), p, dpi = 320, width = 8, height = 10)
  }
}

beta_analysis <- function(dataset, output_dir) {
  beta_res <- trans_beta$new(dataset = dataset, group = "group", measure = "bray")
  write.xlsx(beta_res[["use_matrix"]], file = file.path(output_dir, "Beta_analysis_distance_matrix.xlsx"), sheetName = "Sheet1", rowNames = FALSE)
  
  beta_res$cal_anosim()
  beta_res$cal_manova()
  write.xlsx(beta_res$res_anosim, file = file.path(output_dir, "Beta_analysis_anosim_stat.xlsx"), sheetName = "Sheet1", rowNames = FALSE)
  write.xlsx(beta_res$res_manova, file = file.path(output_dir, "Beta_analysis_manova_stat.xlsx"), sheetName = "Sheet1", rowNames = FALSE)
  
  
  # plot compare the group distances
  beta_res$cal_group_distance(within_group = TRUE)
  beta_res$cal_group_distance_diff(method = 'wilcox')
  beta_res_group_distance_plot <- beta_res$plot_group_distance(add = "mean") +
    theme_pubr(
      base_size = 12,
      base_family = "",
      border = TRUE,
      margin = TRUE,
      legend = 'right',
      x.text.angle = 0
    )
  ggsave(file.path(output_dir, "beta_group_distance.jpg"), beta_res_group_distance_plot, dpi = 320, width = 8, height = 10)


  tryCatch({
    beta_res$cal_ordination(method = "PCoA")
    head(beta_res$res_ordination)
  }, error = function(e) {
    message(" 初始 PCoA 计算出错，已跳过：", conditionMessage(e))
  })

  beta_res$cal_manova(manova_all = FALSE, group = "group")
  write.xlsx(beta_res$res_manova, file = file.path(output_dir, "beta_analysis_bray_stat.xlsx"), sheetName = "Sheet1", rowNames = FALSE)

  beta_res$cal_anosim()
  # 可选的 PLS-DA 排序代码：默认不执行，只保留代码以便以后需要时开启
  if (isTRUE(run_plsda)) {
    tryCatch({
      beta_res$cal_ordination(method = "PLS-DA")
      # ‘method’ 默认 "PCoA"; 可选 "PCoA", "NMDS", "PCA", "DCA", "PLS-DA"(2组以上) 或 "OPLS-DA"(2组)
      # PCoA: principal coordinates analysis;
      # NMDS: non-metric multidimensional scaling;
      # PCA: principal component analysis;
      # DCA: detrended correspondence analysis;
      # PLS-DA: partial least squares discriminant analysis;
      # OPLS-DA: orthogonal partial least squares discriminant analysis.
      beta_ordination_plsda_plot <- beta_res$plot_ordination(
        plot_color = "group",
        plot_shape = "group",
        plot_type = c("point", "ellipse")
      )
      ggsave(file.path(output_dir, "beta_anosim.jpg"), beta_ordination_plsda_plot, dpi = 320, width = 8, height = 10)
      # write.xlsx(beta_res$res_ordination$scores, file = file.path(output_dir, "beta_anosim_result.xlsx"), sheetName = "Sheet1", rowNames = FALSE)
    }, error = function(e) {
      message(" PLS-DA/beta anosim ordination 出错，已跳过：", conditionMessage(e))
    })
  }

  tryCatch({
    beta_res$cal_ordination(method = "NMDS")
    # head(beta_res$res_ordination)
    beta_ordination_nmds_plot <- beta_res$plot_ordination(
      plot_color = "group", plot_shape = "group", plot_type = c("point", "ellipse")
      ) +
      theme_pubr(
        base_size = 12,
        base_family = "",
        border = TRUE,
        margin = TRUE,
        legend = 'right',
        x.text.angle = 0
      )
    ggsave(file.path(output_dir, "beta_NMDS.jpg"), beta_ordination_nmds_plot, dpi=320,width=8,height=10)
    #beta_res$res_ordination$scores
    # write.xlsx(beta_res$res_ordination$scores, file = file.path(output_dir, "beta_NMDS_result.xlsx"), sheetName = "Sheet1", rowNames = FALSE)
  }, error = function(e) {
    message(" NMDS/beta NMDS ordination 出错，已跳过：", conditionMessage(e))
  })

  #计算PCoA
  tryCatch({
    beta_res$cal_ordination(method = "PCoA")
    # head(beta_res$res_ordination)
    beta_ordination_pcoa_plot <- beta_res$plot_ordination(plot_color = "group", plot_shape = "group", plot_type = c("point", "ellipse")) +
      theme_pubr(
        base_size = 12,
        base_family = "",
        border = TRUE,
        margin = TRUE,
        legend = 'right',
        x.text.angle = 0
      )
    ggsave(file.path(output_dir, "beta_PCoA.jpg"), beta_ordination_pcoa_plot, dpi=320,width=8,height=10)
    # write.xlsx(beta_res$res_ordination$scores, file = file.path(output_dir, "beta_PCoA_result.xlsx"), sheetName = "Sheet1", rowNames = FALSE)
  }, error = function(e) {
    message(" PCoA/beta PCoA ordination 出错，已跳过：", conditionMessage(e))
  })

  #计算PCA
  # tryCatch({
  #   beta_res$cal_ordination(method = "PCA")
  #   #head(beta_res$res_ordination)
  #   beta_ordination_pca_plot <- beta_res$plot_ordination(plot_color = "group", plot_shape = "group", plot_type = c("point", "ellipse"))
  #   beta_ordination_pca_plot <- beta_ordination_pca_plot + theme_pubr()
  #   ggsave(file.path(output_dir, "beta_PCA.jpg"), beta_ordination_pca_plot, dpi=320,width=8,height=10)
  #   write.xlsx(beta_res$res_ordination$scores, file = file.path(output_dir, "beta_PCA_result.xlsx"), sheetName = "Sheet1", rowNames = FALSE)
  # }, error = function(e) {
  #   message(" PCA/beta PCA ordination 出错，已跳过：", conditionMessage(e))
  # })
}

net_plot <- function(dataset, tax_level, each_group, output_dir) {
  tryCatch({
    # 1. 构建网络：如果样本量少，建议显式设置较宽松的阈值
    net <- trans_network$new(dataset = dataset, 
                             cor_method = "spearman", 
                             filter_thres = 0.001,
                             COR_p_thres = 0.8) # 画不出来调大 COR_p_thres 这个值
    
    # 2. 计算网络
    net$cal_network(COR_cut = 0.6, COR_p_adjust = 'fdr', COR_p_thres = 0.8) # 画不出来调大 COR_p_thres 这个值
    
    if (is.null(net$res_network)) {
      message(paste0(tax_level, " 生成的网络没有任何连接（Edge），无法绘图"))
      return(invisible(NULL))
    }

    net$cal_module()
    net$get_node_table()
    
    # 3. 动态确定上色方案
    plot_color <- "module" # 默认按模块上色最稳妥
    if (tax_level %in% colnames(dataset$tax_table)) {
        plot_color <- tax_level
    }

    net_network_plot <- net$plot_network(method = "ggraph", node_color = plot_color) +
      theme_pubr(
        base_size = 12,
        base_family = "",
        border = TRUE,
        margin = TRUE,
        legend = 'right',
        x.text.angle = 0
      )
    
    ggsave(file.path(output_dir, paste0("net_plot_", tax_level, "_", each_group, ".jpg")), 
           net_network_plot, dpi = 320, width = 8, height = 10)
    
    plot_taxa_roles <- net$plot_taxa_roles(use_type=2, color_values=my_set_colors)
    ggsave(file.path(output_dir, paste0("module_analysis_", tax_level, "_", each_group,  ".jpg")), 
           plot_taxa_roles, dpi = 320, width = 8, height = 10)
    write.xlsx(plot_taxa_roles$data, file.path(output_dir, paste0("module_analysis_", tax_level, "_", each_group,  ".xlsx")), rowNames=FALSE)

    net$cal_network_attr()
    
    # 使用 openxlsx 导出
    write.xlsx(net_network_plot$data, file.path(output_dir, paste0("network_attribute_all_", tax_level, "_", each_group, ".xlsx")), rowNames=FALSE, colNames=TRUE)
    write.xlsx(net$res_network_attr, file.path(output_dir, paste0("network_attribute_", tax_level, "_", each_group, ".xlsx")), rowNames=TRUE, colNames=FALSE)
    write.xlsx(net$res_node_table, file.path(output_dir, paste0("network_node_module_", tax_level, "_", each_group, ".xlsx")), rowNames=FALSE)
    
  }, error = function(e) {
    message(paste0(each_group, " ", tax_level, " 网络分析 net_plot 出错：", conditionMessage(e)))
  })
}

# 如果组数量 = 3，画三元图
ternary_plot <- function(dataset, tax_level, output_dir, ntaxa = 10) {
  tryCatch({
    ternary_obj <- trans_abund$new(dataset = dataset, taxrank = tax_level, ntaxa = ntaxa, groupmean = "group")
    # 根据 microeco 包，方法名是 plot_tern 而不是 plot_ternary
    ternary_plot_fig <- ternary_obj$plot_tern(color_values = my_set_colors)
    ggsave(file.path(output_dir, paste0("ternary_plot_", tax_level, ".jpg")), ternary_plot_fig, dpi = 320, width = 12, height = 8)
  }, error = function(e) {
    message(paste0(tax_level, " 三元图绘制失败，已跳过：", conditionMessage(e)))
  })
}

# LDA 条形图
diff_bar_plot <- function(lefse_obj, output_dir) {
  top_n_bar <- min(20, nrow(lefse_obj$res_diff))
  bar_plot <- tryCatch(
    lefse_obj$plot_diff_bar(
      add_sig = TRUE,
      use_number = seq_len(top_n_bar),
      width = 0.8,
      # group_order = group_order,
      color_values = my_set_colors,
      group_two_sep = FALSE,
      add_sig_text_size = 3
    ) +
      theme_pubr(
        base_size = 12, base_family = "Arial", border = TRUE, margin = TRUE,
        legend = 'right'
      ),
    error = function(e) {
      message("LDA 条形图绘图失败，已跳过：", conditionMessage(e))
      NULL
    }
  )
  ggsave(file.path(output_dir, paste0("LDA_bar_", lefse_obj$taxa_level, ".png")), bar_plot, dpi = 320, width = 8, height = 10)
}

diff_abund_plot <- function(lefse_obj, output_dir) {
  # 丰度柱状图：关闭 add_sig，避免在空结果上标注导致报错I 
  top_n_abund <- min(20, nrow(lefse_obj$res_diff))
  abund_plot <- tryCatch(
    lefse_obj$plot_diff_abund(
      add_sig = FALSE,
      color_values = my_set_colors,
      use_number = seq_len(top_n_abund),
      # group_order = group_order,
      add_sig_text_size = 3
    ),
    error = function(e) {
      message("丰度柱状图绘图失败，已跳过：", conditionMessage(e))
      NULL
    }
  )
  ggsave(file.path(output_dir, paste0("Abundance_bar_", lefse_obj$taxa_level, ".png")), abund_plot, dpi = 320, width = 8, height = 13)
  # write.xlsx(lefse_obj[["res_diff"]], file = file.path(output_dir, "Lefse_diff_result.xlsx"), sheetName = "Sheet1", rowNames = FALSE)
}

diff_cladogram_plot <- function(lefse_obj, output_dir) {
  # Cladogram，只能画 tax_level = default 的结果
  clado_plot <- tryCatch(
    lefse_obj$plot_diff_cladogram(
      color = my_set_colors,
      use_taxa_num = 500,
      use_feature_num = 20,
      clade_label_level = 5,
      node_size_offset = 1
      # group_order = group_order
    ),
    error = function(e) {
      message("Cladogram 绘图失败，已跳过：", conditionMessage(e))
      NULL
    }
  )
  ggsave(file.path(output_dir, paste0("Clade_tree_", lefse_obj$taxa_level, ".png")), clado_plot, dpi = 320, width = 12, height = 12)
}



# ============================================= 流程开始 =============================================
## 是否运行 beta 分析中的 PLS-DA 排序（默认 FALSE，代码保留但不执行）
run_plsda <- FALSE
opt <- parse_args(OptionParser(option_list = option_list))
otu_path <- opt$otu
sample_path <- opt$sample
tax_path <- opt$tax
pollution_taxa <- trimws(unlist(strsplit(opt$pollution, ",")))
rel_abund <- opt$rel_abund
freq <- opt$freq
outdir <- opt$outdir

# 调试使用
run_plsda <- FALSE
otu_path <- "table_resampling.txt"
sample_path <- "samples_described.txt"
tax_path <- "table_resampling_tax.txt"
pollution_taxa <- "Eukaryota,chloroplast,mitochondria"
rel_abund <- 0.0001
freq <- 0.1
outdir <- "./"

# 创建输出目录（如果不存在则自动创建）
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
setwd(outdir)

# 读取输入表
feature_table <- read.table(otu_path, row.names = 1, header = TRUE, check.names = FALSE)
sample_table <- read.table(sample_path, header = TRUE, check.names = FALSE)
rownames(sample_table)<-sample_table$sample
tax_table<-read.table(tax_path, row.names = 1, header = TRUE, check.names = FALSE)

# 计算 group 和 sample 的 number
group_list <- unique(sample_table$group)
group_number <- length(group_list)
sample_list <- unique(sample_table$sample)
sample_number <- length(sample_list)
print(paste0("group 数量：", group_number))
print(paste0("sample 数量：", sample_number))
print(paste0("group 列表：", paste(group_list, collapse = ", ")))
print(paste0("sample 列表：", paste(sample_list, collapse = ", ")))

# 根据 tax_table 中实际存在的列，自动确定可用的分类水平，避免不存在的 tax_level 报错
all_tax_levels <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
tax_levels <- intersect(all_tax_levels, colnames(tax_table))
print("检测到的 tax_levels：")
print(tax_levels)
if (length(tax_levels) == 0) {
  stop("在 tax_table 中未找到任何标准分类列：", paste(all_tax_levels, collapse = ", "))
}

# 构建 dataset
dataset <- microtable$new(
  sample_table = sample_table,
  otu_table = feature_table,
  tax_table = tax_table
)
# 过滤 pollution / abundance
dataset <- dataset$filter_pollution(taxa = pollution_taxa)
# 过滤 dataset
dataset <- dataset$filter_taxa(rel_abund = rel_abund, freq = freq)
dataset <- dataset$cal_abund()
dataset <- dataset$cal_alphadiv()
dataset <- dataset$cal_betadiv()
dataset_list <- list()
for (tax_level in tax_levels) {
  dataset_list[[tax_level]] <- clone(dataset$merge_taxa(tax_level))
  dataset_list[[tax_level]] <- dataset_list[[tax_level]]$cal_abund()
  dataset_list[[tax_level]] <- dataset_list[[tax_level]]$cal_alphadiv()
  dataset_list[[tax_level]] <- dataset_list[[tax_level]]$cal_betadiv()  # beta 分析需要的
}

dir.create(file.path(outdir, 'Alpha'), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, 'Venn'), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, 'Barplot'), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, 'Beta'), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, 'Tree'), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, 'Lefse'), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, 'Tree_bar_combine'), recursive = TRUE, showWarnings = FALSE)

# 画丰度柱状图
for (tax_level in tax_levels) {
  abundance_barplots(dataset_list[[tax_level]], tax_level, 10, file.path(outdir, 'Barplot'))
}

# 画 Venn 图
for (tax_level in c("Phylum", "Genus", "Species")) {
  taxonomy_venn_plot(dataset_list[[tax_level]], tax_level, file.path(outdir, 'Venn'))
}

# 如果组数量 = 3，画三元图
if (group_number == 3) {
  dir.create(file.path(outdir, 'Ternary'), recursive = TRUE, showWarnings = FALSE)
  for (tax_level in tax_levels) {
    ternary_plot(dataset_list[[tax_level]], tax_level, file.path(outdir, 'Ternary'), 10)
  }
}

# 画树图
tree_plot(dataset_list[["Genus"]], file.path(outdir, 'Tree'))
tree_bar_combine_plot(dataset_list[["Phylum"]], "Phylum", 10, file.path(outdir, 'Tree_bar_combine'))


# 画网络图
for (tax_level in c('Genus', 'Species')) {
  for (each_group in group_list) {
    dir.create(file.path(outdir, 'Net', each_group), recursive = TRUE, showWarnings = FALSE)
    dataset_subgroup <- clone(dataset_list[[tax_level]])
    dataset_subgroup$sample_table <- subset(dataset_subgroup$sample_table, group == each_group)
    dataset_subgroup$tidy_dataset()
    net_plot(dataset_subgroup, tax_level, each_group, file.path(outdir, 'Net', each_group))
  }
}

# alpha 分析
alaha_analysis(dataset, file.path(outdir, 'Alpha'))

# beta 分析
beta_analysis(dataset, file.path(outdir, 'Beta'))

# Genus 和 all 做 lefse 分析
lefse_obj_all <- trans_diff$new(dataset = dataset, method = "lefse", group = "group", alpha = 0.5,
  lefse_subgroup = NULL, p_adjust_method = "none", taxa_level = 'all')
lefse_obj_genus <- trans_diff$new(dataset = dataset, method = "lefse", group = "group", alpha = 0.5,
  lefse_subgroup = NULL, p_adjust_method = "none", taxa_level = 'Genus')
lefse_table_all <- lefse_obj_all$res_diff
lefse_table_genus <- lefse_obj_genus$res_diff
write.xlsx(lefse_table_all, file = file.path(outdir, 'Lefse', "LDA_result_all.xlsx"), sheetName = "Sheet1", rowNames = FALSE)
write.xlsx(lefse_table_genus, file = file.path(outdir, 'Lefse', "LDA_result_genus.xlsx"), sheetName = "Sheet1", rowNames = FALSE)
write.xlsx(lefse_obj_all$res_abund, file = file.path(outdir, 'Lefse', "Abund_result_all.xlsx"), sheetName = "Sheet1", rowNames = FALSE)
write.xlsx(lefse_obj_genus$res_abund, file = file.path(outdir, 'Lefse', "Abund_result_genus.xlsx"), sheetName = "Sheet1", rowNames = FALSE)

for (each_lefse_obj in list(lefse_obj_all, lefse_obj_genus)) {
  diff_bar_plot(each_lefse_obj, file.path(outdir, 'Lefse'))
  diff_abund_plot(each_lefse_obj, file.path(outdir, 'Lefse'))
}
diff_cladogram_plot(lefse_obj_all, file.path(outdir, 'Lefse'))
