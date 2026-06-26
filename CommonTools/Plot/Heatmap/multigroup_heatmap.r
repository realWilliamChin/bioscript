suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(ggplot2))

option_list <- list(
  make_option(c("--fpkm"),
    type = "character", default = "fpkm_matrix_filtered.txt",
    help = "fpkm_matrix_filtered.txt", metavar = "character"
  ),
  make_option(c("--output"),
    type = "character", default = "all_gene_heatmap.jpeg",
    help = "输出目录，默认当前目录", metavar = "character"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check that all required arguments are provided
if (!file.exists(opt$fpkm)) {
  print_help(opt_parser)
  stop("请提供 fpkm_matrix_filtered.txt 文件", call. = FALSE)
}

# Assign the first argument to prefix
fpkm_file <- opt$fpkm
output_file <- opt$output
read.table(fpkm_file, sep = "\t", header = T, row.names = 1, check.names = F) -> fpkm
all_fpkm <- fpkm[rowSums(fpkm) > 0, ]
if (nrow(all_fpkm) > 65535) {
  # 随机挑选 65535 行的索引
  sample_indices <- sample(1:nrow(all_fpkm), size = 65535, replace = FALSE)

  # 按照原始顺序重新排列数据
  subset_data <- all_fpkm[sample_indices, ]
  subset_data <- subset_data[order(sample_indices), ]

  # 创建热图
  all.heatmap <- pheatmap(subset_data, scale = "row", cluster_cols = FALSE, show_rownames = FALSE)
} else {
  # 如果数据不超过 65535 行，直接创建热图
  all.heatmap <- pheatmap(all_fpkm, scale = "row", cluster_cols = FALSE, show_rownames = FALSE)
}

ggsave(output_file, all.heatmap, dpi = 300, width = 10, height = 10)