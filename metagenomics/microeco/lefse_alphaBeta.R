library('microeco')
library(openxlsx)
library(ggplot2)
library(optparse)
rm(list=ls())

# alpha 多样性(Chao1/ACE/Simpson/Shannon/InvSimpson/Pielou/Fisher/Coverage 等多指标 + KW/Dunn 检验)
# 与 beta 多样性(Bray/Jaccard 距离 + PCoA/NMDS/PCA/PLS-DA 排序 + ANOSIM)联合分析脚本
# 输入沿用 microeco 体系下 lefse.R 的三个标准文件: OTU 矩阵 / 样本表 / 物种注释表

option_list <- list(
  make_option(c("-o", "--otu"), type = "character", default = "Species_aboun.txt",
    help = "OTU/feature 矩阵文件路径，默认 Species_aboun.txt"),
  make_option(c("-s", "--sample"), type = "character", default = "samples_described.txt",
    help = "样本描述表文件路径，默认 samples_described.txt"),
  make_option(c("-t", "--tax"), type = "character", default = "Species_taxon.txt",
    help = "物种注释表文件路径，默认 Species_taxon.txt"),
  make_option(c("--outdir"), type = "character", default = "./",
    help = "输出目录，默认当前目录，会自动创建")
)

opt <- parse_args(OptionParser(option_list = option_list))
otu_path <- opt$otu
sample_path <- opt$sample
tax_path <- opt$tax
outdir <- opt$outdir

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

my_set_colors <- c(
  "#e91e63", "#3f51b5", "#dccc39", "#4caf50", "#607d8b",
  "#00bcd4", "#009688", "#8bc34a", "#f44336", "#9c27b0",
  "#ffeb3b", "#ffc107", "#ff5722", "#795548", "#9e9e9e",
  "#c8bfe7", "#b97a57", "#ffaec9", "#1ee1c4", "#2196f3"
)

feature_table <- read.table(otu_path, row.names = 1, header = TRUE, check.names = FALSE)
sample_table <- read.table(sample_path, header = TRUE, check.names = FALSE)
rownames(sample_table) <- sample_table$sample
tax_table <- read.table(tax_path, row.names = 1, header = TRUE, check.names = FALSE)
dataset <- microtable$new(
  sample_table = sample_table,
  otu_table = feature_table,
  tax_table = tax_table
)

# ============== alpha 多样性 ==============
alpha_res <- trans_alpha$new(dataset = dataset, group = "group")

alpha_res$cal_diff(method = "KW")
alpha_res$cal_diff(method = "KW_dunn", KW_dunn_letter = FALSE)

head(alpha_res$res_diff)

# y_increase can adjust the distance from the letters to the highest point
# add_sig_text_size: letter size adjustment
alpha_res$plot_alpha(measure = "Chao1", add_sig_text_size = 6, y_increase = 0.25, add = "jitter", order_x_mean = TRUE)
ggsave(file.path(outdir, "Chao1_bar.png"), dpi=320, width=8, height=10)

alpha_res$plot_alpha(measure = "ACE", add_sig_text_size = 6, y_increase = 0.25, add = "jitter", order_x_mean = TRUE)
ggsave(file.path(outdir, "ACE_bar.png"), dpi=320, width=8, height=10)

alpha_res$plot_alpha(measure = "Simpson", add_sig_text_size = 6, y_increase = 0.25, add = "jitter", order_x_mean = TRUE)
ggsave(file.path(outdir, "Simpson_bar.png"), dpi=320, width=8, height=10)

alpha_res$plot_alpha(measure = "Shannon", add_sig_text_size = 6, y_increase = 0.25, add = "jitter", order_x_mean = TRUE)
ggsave(file.path(outdir, "Shannon_bar.png"), dpi=320, width=8, height=10)

alpha_res$plot_alpha(measure = "InvSimpson", add_sig_text_size = 6, y_increase = 0.25, add = "jitter", order_x_mean = TRUE)
ggsave(file.path(outdir, "InvSimpson_bar.png"), dpi=320, width=8, height=10)

alpha_res$plot_alpha(measure = "Pielou", add_sig_text_size = 6, y_increase = 0.25, add = "jitter", order_x_mean = TRUE)
ggsave(file.path(outdir, "Pielou_bar.png"), dpi=320, width=8, height=10)

alpha_res$plot_alpha(measure = "Fisher", add_sig_text_size = 6, y_increase = 0.25, add = "jitter", order_x_mean = TRUE)
ggsave(file.path(outdir, "Fisher_bar.png"), dpi=320, width=8, height=10)

alpha_res$plot_alpha(measure = "Coverage", add_sig_text_size = 6, y_increase = 0.25, add = "jitter", order_x_mean = TRUE)
ggsave(file.path(outdir, "Coverage_bar.png"), dpi=320, width=8, height=10)

# 输出 summary 结果
write.xlsx(alpha_res$res_diff, file = file.path(outdir, "alpha_summary_result.xlsx"), sheetName = "Sheet1", rowNames = FALSE)


# Note that from v1.9.0, the parameter plot_type control which type of plot is employed. All the options starting with "gg" (e.g., "ggboxplot", "ggdotplot", "ggviolin", "ggstripchart", "ggerrorplot") means they are the functions coming from the ggpubr package.
alpha_res$plot_alpha(plot_type = "ggviolin", alpha = 0.2, y_increase = 0.4, add = "mean_se", add_sig_text_size = 6)

# ============== beta 多样性 ==============
dataset$cal_betadiv()
beta_res <- trans_beta$new(dataset = dataset, group = "group", measure = "bray") # The input measure should only have one element!
beta_res <- trans_beta$new(dataset = dataset, group = "group", measure = "jaccard")

beta_res$cal_ordination(method = "PCoA")
head(beta_res$res_ordination)

# 计算 anosim
beta_res$cal_anosim()
beta_res$cal_ordination(method = "PLS-DA")
# 'method' default "PCoA"; "PCoA", "NMDS", "PCA", "DCA", "PLS-DA"(2组以上) or "OPLS-DA"(2组). PCoA: principal coordinates analysis;
# NMDS: non-metric multidimensional scaling, PCA: principal component analysis; DCA: detrended correspondence analysis; PLS-DA: partial least squares
# discriminant analysis; OPLS-DA: orthogonal partial least squares discriminant analysis.
beta_res$plot_ordination(plot_color = "group", plot_shape = "group", plot_type = c("point", "ellipse"))
ggsave(file.path(outdir, "beta_anosim.png"), dpi=320, width=8, height=10)
write.xlsx(beta_res$res_ordination$scores, file = file.path(outdir, "beta_anosim_result.xlsx"), sheetName = "Sheet1", rowNames = FALSE)

beta_res$cal_ordination(method = "NMDS")
head(beta_res$res_ordination)
beta_res$plot_ordination(plot_color = "group", plot_shape = "group", plot_type = c("point", "ellipse"))
ggsave(file.path(outdir, "beta_NMDS.png"), dpi=320, width=8, height=10)
write.xlsx(beta_res$res_ordination$scores, file = file.path(outdir, "beta_NMDS_result.xlsx"), sheetName = "Sheet1", rowNames = FALSE)

# PCoA
beta_res$cal_ordination(method = "PCoA")
beta_res$plot_ordination(plot_color = "group", plot_shape = "group", plot_type = c("point", "ellipse"))
ggsave(file.path(outdir, "beta_PCoA.png"), dpi=320, width=8, height=10)
write.xlsx(beta_res$res_ordination$scores, file = file.path(outdir, "beta_PCoA_result.xlsx"), sheetName = "Sheet1", rowNames = FALSE)

# PCA
beta_res$cal_ordination(method = "PCA")
beta_res$plot_ordination(plot_color = "group", plot_shape = "group", plot_type = c("point", "ellipse"))
ggsave(file.path(outdir, "beta_PCA.png"), dpi=320, width=8, height=10)
write.xlsx(beta_res$res_ordination$scores, file = file.path(outdir, "beta_PCA_result.xlsx"), sheetName = "Sheet1", rowNames = FALSE)

# ============== Phylum 三元图 ==============
abund_res <- trans_abund$new(dataset = dataset, taxrank = "Phylum", ntaxa = 8, groupmean = "group")
abund_res$plot_tern()
