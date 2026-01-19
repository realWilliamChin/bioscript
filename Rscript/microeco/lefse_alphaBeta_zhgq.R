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
  "#c8bfe7", "#b97a57", "#ffaec9", "#1ee1c4", "#2196f3"
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

#计算alpha-diversity
alpha_res <- trans_alpha$new(dataset = dataset, group = "group")

alpha_res$cal_diff(method = "KW")
alpha_res$cal_diff(method = "KW_dunn", KW_dunn_letter = FALSE)

head(alpha_res$res_diff)

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

alpha_res$plot_alpha(measure = "ACE", add_sig_text_size = 6, y_increase = 0.25, add = "jitter", order_x_mean = TRUE)
ggsave("ACE_bar.png", dpi=320,width=8,height=10)

alpha_res$plot_alpha(measure = "Simpson", add_sig_text_size = 6, y_increase = 0.25, add = "jitter", order_x_mean = TRUE)
ggsave("Simpson_bar.png", dpi=320,width=8,height=10)

alpha_res$plot_alpha(measure = "Shannon", add_sig_text_size = 6, y_increase = 0.25, add = "jitter", order_x_mean = TRUE)
ggsave("Shannon_bar.png", dpi=320,width=8, height=10)

alpha_res$plot_alpha(measure = "InvSimpson", add_sig_text_size = 6, y_increase = 0.25, add = "jitter", order_x_mean = TRUE)
ggsave("InvSimpson_bar.png", dpi=320,width=8,height=10)

alpha_res$plot_alpha(measure = "Pielou", add_sig_text_size = 6, y_increase = 0.25, add = "jitter", order_x_mean = TRUE)
ggsave("Pielou_bar.png", dpi=320,width=8,height=10)

alpha_res$plot_alpha(measure = "Fisher", add_sig_text_size = 6, y_increase = 0.25, add = "jitter", order_x_mean = TRUE)
ggsave("Fisher_bar.png", dpi=320,width=8,height=10)

alpha_res$plot_alpha(measure = "Coverage", add_sig_text_size = 6, y_increase = 0.25, add = "jitter", order_x_mean = TRUE)
ggsave("Coverage_bar.png", dpi=320,width=8,height=10)

#输出summary结果
write.xlsx(alpha_res$res_diff, file = "alpha_summary_result.xlsx", sheetName = "Sheet1", rowNames = FALSE)


#Note that from v1.9.0, the parameter plot_type control which type of plot is employed. All the options starting with “gg” (e.g., “ggboxplot”, “ggdotplot”, “ggviolin”, “ggstripchart”, “ggerrorplot”) means they are the functions coming from the ggpubr package.
alpha_res$plot_alpha(plot_type = "ggviolin", alpha = 0.2, y_increase = 0.4, add = "mean_se", add_sig_text_size = 6)

dataset$cal_betadiv()
beta_res <- trans_beta$new(dataset = dataset, group = "group", measure = "bray") #The input measure should only have one element!
beta_res <- trans_beta$new(dataset = dataset, group = "group", measure = "jaccard")

beta_res$cal_ordination(method = "PCoA")
head(beta_res$res_ordination)

#计算anosim
beta_res$cal_anosim()
beta_res$cal_ordination(method = "PLS-DA")
#‘method’ default "PCoA"; "PCoA", "NMDS", "PCA", "DCA", "PLS-DA"(2组以上) or "OPLS-DA"(2组). PCoA: principal coordinates analysis; 
# NMDS: non-metric multidimensional scaling, PCA: principal component analysis; DCA: detrended correspondence analysis; PLS-DA: partial least squares
# discriminant analysis; OPLS-DA: orthogonal partial least squares discriminant analysis. For the methods details,
beta_res$plot_ordination(plot_color = "group", plot_shape = "group", plot_type = c("point", "ellipse"))
ggsave("beta_anosim.png", dpi=320,width=8,height=10)
write.xlsx(beta_res$res_ordination$scores, file = "beta_anosim_result.xlsx", sheetName = "Sheet1", rowNames = FALSE)

beta_res$cal_ordination(method = "NMDS")
head(beta_res$res_ordination)
beta_res$plot_ordination(plot_color = "group", plot_shape = "group", plot_type = c("point", "ellipse"))
ggsave("beta_NMDS.png", dpi=320,width=8,height=10)
#beta_res$res_ordination$scores
write.xlsx(beta_res$res_ordination$scores, file = "beta_NMDS_result.xlsx", sheetName = "Sheet1", rowNames = FALSE)

#计算PCoA
beta_res$cal_ordination(method = "PCoA")
#head(beta_res$res_ordination)
beta_res$plot_ordination(plot_color = "group", plot_shape = "group", plot_type = c("point", "ellipse"))
ggsave("beta_PCoA.png", dpi=320,width=8,height=10)
write.xlsx(beta_res$res_ordination$scores, file = "beta_PCoA_result.xlsx", sheetName = "Sheet1", rowNames = FALSE)

#计算PCA
beta_res$cal_ordination(method = "PCA")
#head(beta_res$res_ordination)
beta_res$plot_ordination(plot_color = "group", plot_shape = "group", plot_type = c("point", "ellipse"))
ggsave("beta_PCA.png", dpi=320,width=8,height=10)
write.xlsx(beta_res$res_ordination$scores, file = "beta_PCA_result.xlsx", sheetName = "Sheet1", rowNames = FALSE)

#abund_res <- trans_abund$new(dataset = dataset, taxrank = "Phylum", ntaxa = 8)
abund_res <- trans_abund$new(dataset = dataset, taxrank = "Phylum", ntaxa = 8, groupmean = "group")
abund_res$plot_tern()
