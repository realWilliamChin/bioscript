suppressPackageStartupMessages(library(ggthemes))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggcorrplot))
suppressPackageStartupMessages(library(corrplot))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(FactoMineR))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(c("--fpkm"),
    type = "character", default = "fpkm_matrix_filtered.txt",
    help = "fpkm_matrix_filtered.txt", metavar = "character"
  ),
  make_option(c("--output"),
    type = "character", default = "density.jpeg",
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

melt(fpkm, variable.name = "sample", value.name = "fpkm") -> data.m
data.m[data.m[, 2] > 0, ] -> data.m
d <- ggplot(data.m, aes(x = log10(fpkm), col = sample)) +
  geom_density(aes(fill = sample), colour = NA, alpha = .2) +
  geom_line(stat = "density", size = 1.5) +
  xlab("log2(FPKM)") +
  theme_base() +
  theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), axis.text.x = element_text(size = rel(3)), axis.text.y = element_text(size = rel(3)))
ggsave(output_file, plot = d, height = 10, width = 16.8, dpi = 300)