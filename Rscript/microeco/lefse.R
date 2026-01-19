library(microeco)
library(openxlsx)
library(ggplot2)
library(optparse)
library(patchwork)
rm(list=ls())

option_list <- list(
  # 输入文件
  make_option(c("-o", "--otu"), type = "character", default = "Species_aboun.txt",
    help = "OTU/feature 矩阵文件路径，默认 Species_aboun.txt"),
  make_option(c("-s", "--sample"), type = "character", default = "samples_described.txt",
    help = "样本描述表文件路径，默认 samples_described.txt"),
  make_option(c("-t", "--tax"), type = "character", default = "Species_taxon.txt",
    help = "物种注释表文件路径，默认 Species_taxon.txt"),

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
otu_path <- "Species_aboun.txt"
sample_path <- "samples_described.txt"
tax_path <- "Species_taxon.txt"
pollution_taxa <- "Eukaryota,chloroplast,mitochondria"
rel_abund <- 0.0001
freq <- 0.1
outdir <- "./"

# 创建输出目录（如果不存在则自动创建）
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
setwd(outdir)

my_set_colors <- c(
  "#e91e63", "#3f51b5", "#dccc39", "#4caf50", "#607d8b",
  "#00bcd4", "#009688", "#8bc34a", "#f44336", "#9c27b0",
  "#ffeb3b", "#ffc107", "#ff5722", "#795548", "#9e9e9e",
  "#c8bfe7", "#b97a57", "#ffaec9", "#1ee1c4", "#2196f3"
)

# 读取输入表
feature_table <- read.table(otu_path, row.names = 1, header = TRUE, check.names = FALSE)
sample_table <- read.table(sample_path, header = TRUE, check.names = FALSE)
rownames(sample_table)<-sample_table$sample
tax_table<-read.table(tax_path, row.names = 1, header = TRUE, check.names = FALSE)

# 计算 group 和 sample 的 number
group_number <- length(unique(sample_table$group))
sample_number <- nrow(sample_table)
print(paste0("group 数量：", group_number))
print(paste0("sample 数量：", sample_number))

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
dataset_list <- list()
for (tax_level in tax_levels) {
  dataset_list[[tax_level]] <- clone(dataset$merge_taxa(tax_level))
  dataset_list[[tax_level]] <- dataset_list[[tax_level]]$cal_abund()
  dataset_list[[tax_level]] <- dataset_list[[tax_level]]$cal_alphadiv()
  dataset_list[[tax_level]] <- dataset_list[[tax_level]]$cal_betadiv()  # beta 分析需要的
}

# trans_abund 画图函数，包含多个分类水平的丰度柱状图与分组均值图
plot_abundance_barplots <- function(dataset, tax_level, ntaxa = 10, output_dir) {
  abund_res <- trans_abund$new(dataset = dataset, taxrank = tax_level, ntaxa = ntaxa)
  abund_res_groupmean <- trans_abund$new(dataset = dataset, taxrank = tax_level, ntaxa = ntaxa, groupmean = "group")

  bar_plot <- abund_res$plot_bar(legend_text_italic = TRUE, color_values = my_set_colors)
  ggsave(filename = file.path(output_dir, paste0("abundance_bar_", tax_level, ".jpg")), plot = bar_plot, dpi = 320, width = 12, height = 8)
  
  bar_plot_group <- abund_res_groupmean$plot_bar(legend_text_italic = TRUE, color_values = my_set_colors)
  ggsave(filename = file.path(output_dir, paste0("abundance_bar_groupmean_", tax_level, ".jpg")), plot = bar_plot_group, dpi = 320, width = 4, height = 8)

  #! donut graph by taxonomy 
  donut_plot <- abund_res_groupmean$plot_donut(xtext_keep = FALSE, legend_text_italic = TRUE, color_values = my_set_colors)
  ggsave(filename = file.path(output_dir, paste0("abundance_donut_", tax_level, ".jpg")), plot = donut_plot, dpi = 320, width = 12, height = 8)

  write.xlsx(bar_plot$plot_env$new_data, file = file.path(output_dir, paste0("abundance_bar_", tax_level, ".xlsx")), sheetName = "Sheet1", rowNames = TRUE)
}
for (tax_level in tax_levels) {
  dir.create(file.path(outdir, 'Barplot'), recursive = TRUE, showWarnings = FALSE)
  plot_abundance_barplots(dataset_list[[tax_level]], tax_level, 10, file.path(outdir, 'Barplot'))
}


plot_taxonomy_venn <- function(dataset, output_dir) {
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
for (tax_level in c("Phylum", "Genus", "Species")) {
  dir.create(file.path(outdir, 'Venn'), recursive = TRUE, showWarnings = FALSE)
  plot_taxonomy_venn(dataset_list[[tax_level]], file.path(outdir, 'Venn'))
}


# 画样品分类进化树，只做 Genus 结果
plot_tree <- function(dataset, output_dir) {
  dataset$tidy_dataset()
  beta_res <- trans_beta$new(dataset, group = "group", measure = "bray")
  # use replace_name to set the label name, group parameter used to set the color
  tree_plot <- beta_res$plot_clustering(group = "group", replace_name = c("sample"), color_values = my_set_colors)
  ggsave(file.path(output_dir, "Genus_tree.jpg"), tree_plot, dpi = 320, width = 12, height = 10)
}
dir.create(file.path(outdir, 'Tree'), recursive = TRUE, showWarnings = FALSE)
plot_tree(dataset_list[["Genus"]], file.path(outdir, 'Tree'))

# 树和丰度柱状图结合，只画 Phylum 结果
plot_tree_bar_combine <- function(dataset, tax_level, ntaxa, output_dir) {
  abund_res <- trans_abund$new(dataset = dataset, taxrank = tax_level, ntaxa = ntaxa)
  beta_res <- trans_beta$new(dataset, group = "group", measure = "bray")
  tree_plot <- beta_res$plot_clustering(group = "group", replace_name = c("sample"), color_values = my_set_colors)

  tree_order <- tree_plot$plot_env$data2$sample
  bar_plot <- abund_res$plot_bar(xtext_keep = TRUE, legend_text_italic = TRUE, color_values = my_set_colors, coord_flip = TRUE, order_x = tree_order)
  
  
  p_final = tree_plot + bar_plot + plot_layout(widths = c(2, 2))
  ggsave(file.path(output_dir, paste0(tax_level, "_tree_bar_combine.jpg")), p_final, dpi = 320, width = 12, height = 10)
}
dir.create(file.path(outdir, 'Tree_bar_combine'), recursive = TRUE, showWarnings = FALSE)
plot_tree_bar_combine(dataset_list[["Phylum"]], "Phylum", 10, file.path(outdir, 'Tree_bar_combine'))


# alpha-diversity 分析
alaha_analysis <- function(dataset, tax_level, output_dir) {
  alpha_res <- trans_alpha$new(dataset = dataset, group = "group")
  alpha_res$cal_diff(method = "KW")
  # 在某些情况下（如结果结构缺少 Method 列）KW_dunn 可能会报错，这里做容错处理
  tryCatch({
    alpha_res$cal_diff(method = "KW_dunn", KW_dunn_letter = FALSE)
  }, error = function(e) {
    message(paste0(tax_level, " KW_dunn 多重比较计算出错，已跳过：", conditionMessage(e)))
  })
  
  alpha_measures <- c("Chao1", "ACE", "Simpson", "Shannon", "InvSimpson", "Pielou", "Fisher", "Coverage")
  for (m in alpha_measures) {
    p <- tryCatch(
      alpha_res$plot_alpha(measure = m, add_sig_text_size = 6, y_increase = 0.25, add = "jitter", order_x_mean = TRUE),
      error = function(e) {
        message(paste0(tax_level, " alpha 指标 ", m, " 绘图失败，已跳过：", conditionMessage(e)))
        return(NULL)
      }
    )
    if (!is.null(p)) {
      ggsave(file.path(output_dir, paste0(m, '_', tax_level, '_bar.jpg')), p, dpi = 320, width = 8, height = 10)
    }
  }

  #输出summary结果
  write.xlsx(alpha_res$res_diff, file = file.path(output_dir, paste0("alpha_summary_result_", tax_level, ".xlsx")), sheetName = "Sheet1", rowNames = FALSE)
  p <- tryCatch(
    alpha_res$plot_alpha(plot_type = "ggviolin", alpha = 0.2, y_increase = 0.4, add = "mean_se", add_sig_text_size = 6),
    error = function(e) {
      message(paste0(tax_level, " alpha 小提琴图绘图失败，已跳过：", conditionMessage(e)))
      return(NULL)
    }
  )
  if (!is.null(p)) {
    ggsave(file.path(output_dir, paste0("alpha_violin_", tax_level, ".jpg")), p, dpi = 320, width = 8, height = 10)
  }
}
for (tax_level in tax_levels) {
  dir.create(file.path(outdir, 'Alpha'), recursive = TRUE, showWarnings = FALSE)
  alaha_analysis(dataset_list[[tax_level]], tax_level, file.path(outdir, 'Alpha'))
}



beta_analysis <- function(dataset, tax_level, output_dir) {
  beta_res <- trans_beta$new(dataset = dataset, group = "group", measure = "bray")
  tryCatch({
    beta_res$cal_ordination(method = "PCoA")
    head(beta_res$res_ordination)
  }, error = function(e) {
    message(paste0(tax_level, " 初始 PCoA 计算出错，已跳过：", conditionMessage(e)))
  })

  beta_res$cal_manova(manova_all = FALSE, group = "group")
  write.xlsx(beta_res$res_manova, file = file.path(output_dir, paste0("beta_stat_by_group_", tax_level, ".xlsx")), sheetName = "Sheet1", rowNames = FALSE)

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
      ggsave(file.path(output_dir, paste0("beta_anosim_", tax_level, ".jpg")),
             beta_ordination_plsda_plot, dpi = 320, width = 8, height = 10)
      write.xlsx(beta_res$res_ordination$scores,
                 file = file.path(output_dir, paste0("beta_anosim_result_", tax_level, ".xlsx")),
                 sheetName = "Sheet1", rowNames = FALSE)
    }, error = function(e) {
      message(paste0(tax_level, " PLS-DA/beta anosim ordination 出错，已跳过：", conditionMessage(e)))
    })
  }

  tryCatch({
    beta_res$cal_ordination(method = "NMDS")
    head(beta_res$res_ordination)
    beta_ordination_nmds_plot <- beta_res$plot_ordination(plot_color = "group", plot_shape = "group", plot_type = c("point", "ellipse"))
    ggsave(file.path(output_dir, paste0("beta_NMDS_", tax_level, ".jpg")), beta_ordination_nmds_plot, dpi=320,width=8,height=10)
    #beta_res$res_ordination$scores
    write.xlsx(beta_res$res_ordination$scores, file = file.path(output_dir, paste0("beta_NMDS_result_", tax_level, ".xlsx")), sheetName = "Sheet1", rowNames = FALSE)
  }, error = function(e) {
    message(paste0(tax_level, " NMDS/beta NMDS ordination 出错，已跳过：", conditionMessage(e)))
  })

  #计算PCoA
  tryCatch({
    beta_res$cal_ordination(method = "PCoA")
    #head(beta_res$res_ordination)
    beta_ordination_pcoa_plot <- beta_res$plot_ordination(plot_color = "group", plot_shape = "group", plot_type = c("point", "ellipse"))
    ggsave(file.path(output_dir, paste0("beta_PCoA_", tax_level, ".jpg")), beta_ordination_pcoa_plot, dpi=320,width=8,height=10)
    write.xlsx(beta_res$res_ordination$scores, file = file.path(output_dir, paste0("beta_PCoA_result_", tax_level, ".xlsx")), sheetName = "Sheet1", rowNames = FALSE)
  }, error = function(e) {
    message(paste0(tax_level, " PCoA/beta PCoA ordination 出错，已跳过：", conditionMessage(e)))
  })

  #计算PCA
  tryCatch({
    beta_res$cal_ordination(method = "PCA")
    #head(beta_res$res_ordination)
    beta_ordination_pca_plot <- beta_res$plot_ordination(plot_color = "group", plot_shape = "group", plot_type = c("point", "ellipse"))
    ggsave(file.path(output_dir, paste0("beta_PCA_", tax_level, ".jpg")), beta_ordination_pca_plot, dpi=320,width=8,height=10)
    write.xlsx(beta_res$res_ordination$scores, file = file.path(output_dir, paste0("beta_PCA_result_", tax_level, ".xlsx")), sheetName = "Sheet1", rowNames = FALSE)
  }, error = function(e) {
    message(paste0(tax_level, " PCA/beta PCA ordination 出错，已跳过：", conditionMessage(e)))
  })
}
for (tax_level in tax_levels) {
  dir.create(file.path(outdir, 'Beta'), recursive = TRUE, showWarnings = FALSE)
  tryCatch({
    beta_analysis(dataset_list[[tax_level]], tax_level, file.path(outdir, 'Beta'))
  }, error = function(e) {
    message(paste0(tax_level, " 分类水平 的 beta 分析整体出错，已跳过：", conditionMessage(e)))
  })
}


net_plot <- function(dataset, tax_level, output_dir) {
  tryCatch({
    # 1. 构建网络：如果样本量少，建议显式设置较宽松的阈值
    net <- trans_network$new(dataset = dataset, 
                             cor_method = "spearman", 
                             filter_thres = 0.001,
                             COR_p_thres = 0.05) # 显式设置P值阈值
    
    # 2. 计算网络，捕获可能的空网络错误
    net$cal_network(COR_cut = 0.5) # 适当降低相关系数要求
    
    if (is.null(net$res_network)) {
      message(paste0(tax_level, " 生成的网络没有任何连接（Edge），无法绘图"))
      return(invisible(NULL))
    }

    net$cal_module()
    
    # 3. 动态确定上色方案
    plot_color <- "module" # 默认按模块上色最稳妥
    if (tax_level %in% colnames(dataset$tax_table)) {
        plot_color <- tax_level
    }

    net_network_plot <- net$plot_network(method = "ggraph", node_color = plot_color)
    
    ggsave(file.path(output_dir, paste0("net_plot_", tax_level, ".jpg")), 
           net_network_plot, dpi = 320, width = 8, height = 10)

    net$cal_network_attr()
    
    # 使用 openxlsx 导出
    write.xlsx(net$res_network_attr, file.path(output_dir, paste0("network_attribute_", tax_level, ".xlsx")))
    write.xlsx(net$res_node_table, file.path(output_dir, paste0("network_node_module_", tax_level, ".xlsx")))
    
  }, error = function(e) {
    message(paste0(tax_level, " 网络分析 net_plot 出错：", conditionMessage(e)))
  })
}

for (tax_level in tax_levels) {
  dir.create(file.path(outdir, 'Net'), recursive = TRUE, showWarnings = FALSE)
  net_plot(dataset_list[[tax_level]], tax_level, file.path(outdir, 'Net'))
}


# 如果组数量 = 3，画三元图
plot_ternary_plot <- function(dataset, tax_level, output_dir, ntaxa = 10) {
  ternary_plot <- trans_abund$new(dataset = dataset, taxrank = tax_level, ntaxa = ntaxa, groupmean = "group")
  ternary_plot$plot_ternary(color_values = my_set_colors)
  ggsave(file.path(output_dir, paste0("ternary_plot_", tax_level, ".jpg")), ternary_plot, dpi = 320, width = 8, height = 10)
}
if (group_number == 3) {
  for (tax_level in tax_levels) {
    plot_ternary_plot(dataset_list[[tax_level]], tax_level, file.path(outdir, 'Ternary'), 10)
  }
}





library('microeco')
library(openxlsx)
library(ggplot2)
library(optparse)
rm(list=ls())

#TODO: 加入参数解析功能

my_set_colors <- c(
  "#e91e63", "#3f51b5", "#dccc39", "#4caf50", "#607d8b",
  "#00bcd4", "#009688", "#8bc34a", "#f44336", "#9c27b0",
  "#ffeb3b", "#ffc107", "#ff5722", "#795548", "#9e9e9e",
  "#c8bfe7", "#b97a57", "#ffaec9", "#1ee1c4", "#2196f3",
)

feature_table<-read.table("Species_aboun.txt",row.names=1,header = T,check.names=F)
sample_table<-read.table("samples_described.txt",header = T,check.names=F)
rownames(sample_table)<-sample_table$sample
tax_table<-read.table("Species_taxon.txt",row.names=1,header = T,check.names=F)
dataset <- microtable$new(  
  sample_table = sample_table,  
  otu_table = feature_table,  
  tax_table = tax_table
)

lefse_res <- trans_diff$new(dataset=dataset,
                          method = "lefse",
                          group = "group",
                          alpha=0.5,
                          lefse_subgroup=NULL,
                          p_adjust_method = "none"
)
res.dt<-lefse_res$res_diff
write.xlsx(res.dt, file = "LDA_result.xlsx", sheetName = "Sheet1", rowNames = FALSE)

p1<-lefse_res$plot_diff_bar(add_sig = T,
                        use_number = 1:20,
                        width=0.8,
                        group_order = unique(sample_table$group),
                        color_values = my_set_colors,
                        )
ggsave("LDA_bar.png", p1, dpi=320,width=8,height=10)

# 画丰富度柱状图：关闭 add_sig，避免在空结果上标注导致报错
use_n_abund <- min(20, nrow(res.dt))
p2 <- lefse_res$plot_diff_abund(add_sig = F,
                          color_values = my_set_colors,
                          use_number = 1:use_n_abund,
                        #width=0.8,
                          group_order = unique(sample_table$group),
)
ggsave("Abundance_bar.png",p2, dpi=320,width=8,height=13)

p3<-lefse_res$plot_diff_cladogram(color = my_set_colors,
                              use_taxa_num = 500,
                              use_feature_num = 20,
                              clade_label_level = 5,
                              #clade_label_size = 100,
                              #select_show_labels = "c__Bacilli",
                              node_size_offset = 1,
                              #annotation_shape = 21,
                              #annotation_shape_size = 4.5,
                              group_order = unique(sample_table$group)
                        #clade_label_size =1
                      )

ggsave("Clade_tree.png",p3,dpi=320,width=12,height=12)
