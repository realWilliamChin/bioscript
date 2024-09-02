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
    type = "character", default = "correlation.jpeg",
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

fpkm.m <- as.matrix(fpkm)
fpkm.cor <- cor(fpkm.m)
# p.mat <- cor_pmat(fpkm.m)
# p.mat
# cor.p<-ggcorrplot(fpkm.cor,lab=TRUE)
# cor.p<-corrplot(res, method = "number", col = RColorBrewer::brewer.pal(n=11, name = "RdYlGn"),title = "")
# ?corrplot
# ?ggcorrplot
write.table(fpkm.cor, "sample_cor.txt", sep = "\t", quote = F)
# ggsave("correlation.jpeg",cor.p,dpi=300)
min(fpkm.cor)
# fpkm.m<-as.data.frame(fpkm)
# fpkm.cor<-cor(fpkm.m)
resfactor <- 5
png(output_file, res = 72 * resfactor, height = 1200 * resfactor, width = 1200 * resfactor)
corrplot(fpkm.cor, is.corr = F, col = rev(COL2("PiYG")), method = "color", addCoef.col = "black", tl.col = "black", col.lim = c(min(fpkm.cor) - 0.01, max(fpkm.cor)), cl.ratio = 0.1)
dev.off()

