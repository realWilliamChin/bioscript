library(openxlsx)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(optparse)
rm(list = ls())

my_set_colors <- c(
  "#f44336", "#e91e63", "#9c27b0", "#3f51b5", "#2196f3",
  "#00bcd4", "#009688", "#4caf50", "#8bc34a", "#cddc39",
  "#ffeb3b", "#ffc107", "#ff5722", "#795548", "#9e9e9e",
  "#607d8b", "#c8bfe7", "#b97a57", "#ffaec9", "#1ee1c4"
)

option_list <- list(
  make_option(c("-f", "--datafile"), type="character", default=NULL, 
              help="两组之间比对的 fpkm 表达量 excel file", metavar="character"),
  make_option(c("-o", "--outputfile"), type="character", default=NULL, 
              help="输出的图片名称", metavar="character"),
  # 画图参数
  make_option("--no-cluster-rows", action="store_false", dest="cluster_rows",
              help='指定画图时不对行聚类'),
  make_option("--cluster-rows", action="store_true", dest="cluster_rows",
              help='指定画图时不对行聚类')
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check that all required arguments are provided
if (is.null(opt$datafile)) {
  print_help(opt_parser)
  stop("请输入 datafile 文件", call. = FALSE)
} 

datafile <- opt$datafile
output <- opt$outputfile
cluster_rows <- opt$cluster_rows


wb <- loadWorkbook(datafile)
sheet_names <- getSheetNames(datafile)

mat <- read.xlsx(datafile, rowNames = TRUE, sheet = 1)
sample_annotation <- read.xlsx(datafile, rowNames = TRUE, sheet = 2)
sample_group_number <- length(unique(sample_annotation$group))
my_sample_colors <- sample(my_set_colors, sample_group_number)
names(my_sample_colors) <- unique(sample_annotation$group)

# group 是 sample_annotation 的 group 列名，Ontology 和 gene_annotation 一致
ann_colors <- list(group = my_sample_colors)
annotation_row <- NA
if (length(sheet_names) == 3) {
  gene_annotation <- read.xlsx(datafile, rowNames = TRUE, sheet = 3)
  gene_group_number <- length(unique(gene_annotation$Ontology))
  my_gene_colors <- sample(my_set_colors, gene_group_number)
  names(my_gene_colors) <- unique(gene_annotation$Ontology)
  ann_colors <- list(group = my_sample_colors, Ontology = my_gene_colors)
  annotation_row <- gene_annotation
}

show_rownames <- if (nrow(mat) > 200) FALSE else TRUE
# 使用pheatmap绘制热图
p1 <- pheatmap(mat,
  cluster_rows = cluster_rows, # 不对行（基因）进行聚类
  cluster_cols = FALSE, # 不对列（样品）进行聚类
  annotation_row = annotation_row, # 添加基因注释
  annotation_col = sample_annotation, # 添加样品注释
  colorRampPalette(c("green", "black", "red"))(50),
  scale = "row",
  annotation_colors = ann_colors, # 设置注释颜色
  show_rownames = show_rownames, # 如果你不想显示基因名，可以设置为FALSE
  show_colnames = TRUE,
  annotation_names_row = FALSE,
  annotation_names_col = FALSE
) # 如果你不想显示样品名，可以设置为FALSE

p1_height <- 3 + (nrow(mat) / 8)
p1_width <- 3 + ncol(mat) * 0.8
if ((p1_width - p1_height) > 5 * p1_height) {
  p1_width <- p1_height * 5
} else if ((p1_height - p1_width) > 3 * p1_width) {
  p1_width <- p1_height / 3
}
ggsave(output, p1, dpi = 320, width = p1_width, height = p1_height, limitsize = FALSE)


# ---------------- 以下是测试代码 ----------------
# 示例数据
# set.seed(123) # 设置随机种子以便结果可复现
# mat <- matrix(rnorm(200), 20, 10) # 20个基因，10个样品的随机数据
# rownames(mat)=paste0("gene",1:20)
# colnames(mat)=paste0("sample",1:10)
# mat
# # 假设的分组信息
# gene_groups <- rep(c("GroupA", "GroupB"), each = 10) # 基因分为两组A和B
# gene_groups
# sample_groups <- rep(c("GroupC", "GroupD"), each = 5) # 样品分为两组C和D
# sample_groups
# 创建基因和样品的注释数据框
# gene_annotation <- data.frame(Group_gene = gene_groups)
# rownames(gene_annotation)<-rownames(mat)
# gene_annotation
# sample_annotation <- data.frame(Group_sample = sample_groups)
# sample_annotation
# rownames(sample_annotation)<-colnames(mat)
# ------------------------------------------------
