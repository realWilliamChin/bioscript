#######################################################
# version: 8.0.0
#######################################################
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
  make_option(c("--reads"),
    type = "character", default = "reads_matrix_filtered.txt",
    help = "提供 reads_matrix_filtered.txt 文件", metavar = "character"
  ),
  make_option(c("--samples"),
    type = "character", default = "samples_described.txt",
    help = "提供 samples_described.txt 文件", metavar = "integer"
  ),
  make_option(c("--compare"),
    type = "character", default = "compare_info.txt",
    help = "提供 compare_info.txt 文件", metavar = "character"
  ),
  make_option(c("--filtertype"),
    type = "character", default = "padj",
    help = "过滤类型，选择 padj 或者 pvalue，默认是 padj，如果差异太少，则可以尝试重新选择 padj", metavar = "character"
  ),
  make_option(c("--filtervalue"),
    type = "double", default = 0.05,
    help = "通常是 0.05 如果值太少，或者太多可以进行调整", metavar = "double"
  ),
  make_option(c("--degvalue"),
    type = "double", default = 2.00,
    help = "deg 值", metavar = "double"
  ),
  make_option(c("--outputdir"),
    type = "character", default = "./",
    help = "输出目录，默认当前目录", metavar = "character"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check that all required arguments are provided
if (!file.exists(opt$fpkm)) {
  print_help(opt_parser)
  stop("请提供 fpkm_matrix_filtered.txt 文件", call. = FALSE)
} else if (!file.exists(opt$reads)) {
  print_help(opt_parser)
  stop("请提供 reads_matrix_filtered.txt 文件", call. = FALSE)
} else if (!file.exists(opt$samples)) {
  print_help(opt_parser)
  stop("请提供 samples_described.txt 文件", call. = FALSE)
} else if (!file.exists(opt$compare)) {
  print_help(opt_parser)
  stop("请提供 compare_info.txt 文件", call. = FALSE)
}

# 处理输出文件夹
if (!dir.exists(opt$outputdir)) {
  dir.create(opt$outputdir)
  cat("文件夹已创建：", opt$outputdir, "\n")
} else {
  cat("文件夹已存在，自动覆盖已有的文件：", opt$outputdir, "\n")
}
if (substr(opt$outputdir, nchar(opt$outputdir), nchar(opt$outputdir)) != "/") {
  opt$outputdir <- paste0(opt$outputdir, "/")
}

# 单独运行数据
# fpkm_file <- 'fpkm_matrix_filtered.txt'
# reads_file <- 'reads_matrix_filtered.txt'
# samples_file <- 'samples_described.txt'
# compare_file <- 'compare_info.txt'
# filter_type <- 'padj'
# filter_value <- '0.05'
# deg_value <- 1.5
# bs_pos <- log2(deg_value)
# output_dir <- './'
# bs_neg <- -bs_pos

# Assign the first argument to prefix
fpkm_file <- opt$fpkm
reads_file <- opt$reads
samples_file <- opt$samples
compare_file <- opt$compare
filter_type <- opt$filtertype
filter_value <- opt$filtervalue
deg_value <- opt$degvalue
bs_pos <- log2(deg_value)
output_dir <- opt$outputdir
bs_neg <- -bs_pos

deg_dir <- paste0(output_dir, "DEG_analysis_results/")
deg_exp_data_dir <- paste0(output_dir, "DEG_analysis_results/Expression_data/")
deg_exp_graph_dir <- paste0(output_dir, "DEG_analysis_results/Expression_data_graphs/")
exp_evaluation_dir <- paste0(output_dir, "Expression_data_evaluation/")

dir.create(deg_dir)
dir.create(deg_exp_data_dir)
dir.create(deg_exp_graph_dir)
dir.create(exp_evaluation_dir)

read.table(fpkm_file, sep = "\t", header = T, row.names = 1, check.names = F) -> fpkm

all_fpkm <- fpkm[rowSums(fpkm) > 0, ]
# all.heatmap<-pheatmap(all_fpkm,scale="row",cluster_cols=F,show_rownames=F)
# 检查数据行数是否超过 65535
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

ggsave(paste0(exp_evaluation_dir, "all_gene_heatmap.jpeg"), all.heatmap, dpi = 300, width = 10, height = 10)

sample_info <- read.table(samples_file, sep = "\t", header = T, check.names = F, stringsAsFactors = F)
# sample_info[sample_info$group=="CK",]
reads_data <- read.table(reads_file, sep = "\t", row.names = 1, header = T, check.names = F, stringsAsFactors = F)
reads_data <- na.omit(reads_data)
comp_info <- read.table(compare_file, sep = "\t", header = T, check.names = F, stringsAsFactors = F)
# head(reads_data)
# head(reads_data[,sample_info[sample_info$group %in% comp_info[1,],]$sample])
# tmp.data<-reads_data[,sample_info[sample_info$group %in% comp_info[1,],]$sample]
# rownames(comp_info)
i <- ""
total.deg <- ""
stat.deg <- data.frame()
for (i in seq_along(1:nrow(comp_info))) {
  group_vs_group_name <- paste(comp_info[i, 1], comp_info[i, 2], sep = "_vs_")
  # print(i)
  data.treat <- reads_data[, sample_info[sample_info$group == comp_info[i, 1], ]$sample]
  data.control <- reads_data[, sample_info[sample_info$group == comp_info[i, 2], ]$sample]
  rnaseqMatrix <- cbind(data.treat, data.control)
  fpkm.treat <- fpkm[, sample_info[sample_info$group == comp_info[i, 1], ]$sample]
  fpkm.control <- fpkm[, sample_info[sample_info$group == comp_info[i, 2], ]$sample]
  fpkm.deg <- cbind(fpkm.treat, fpkm.control)
  rnaseqMatrix <- round(rnaseqMatrix)
  rnaseqMatrix <- rnaseqMatrix[rowSums(cpm(rnaseqMatrix) > 1) >= 2, ]
  conditions <- data.frame(conditions = factor(c(rep(comp_info[i, 1], ncol(data.treat)), rep(comp_info[i, 2], ncol(data.control)))))
  # conditions
  rownames(conditions) <- colnames(rnaseqMatrix)
  ddsFullCountTable <- DESeqDataSetFromMatrix(
    countData = rnaseqMatrix,
    colData = conditions,
    design = ~conditions
  )
  dds <- DESeq(ddsFullCountTable)
  contrast <- c("conditions", comp_info[i, 1], comp_info[i, 2])
  res <- results(dds, contrast)
  baseMeanA <- rowMeans(counts(dds, normalized = TRUE)[, colData(dds)$conditions == comp_info[i, 1]])
  baseMeanB <- rowMeans(counts(dds, normalized = TRUE)[, colData(dds)$conditions == comp_info[i, 2]])
  res <- cbind(baseMeanA, baseMeanB, as.data.frame(res))
  res <- cbind(sampleA = comp_info[i, 1], sampleB = comp_info[i, 2], as.data.frame(res))
  res$padj[is.na(res$padj)] <- 1
  res <- as.data.frame(res[order(res$pvalue), ])
  outfile <- paste0(group_vs_group_name, "_DE_results")
  outfile.ma <- paste0(group_vs_group_name, "_DE_results_readCounts.matrix")
  rnaseqMatrix <- cbind(as.data.frame(rownames(rnaseqMatrix)), rnaseqMatrix)
  colnames(rnaseqMatrix)[1] <- "GeneID"
  write.table(rnaseqMatrix, file = outfile.ma, sep = "	", quote = FALSE, row.names = F)
  volcano <- res
  volcano$padj <- ifelse(volcano$padj < 0.000000000000001, 0.000000000000001, volcano$padj)
  volcano$pvalue <- ifelse(is.na(volcano$pvalue), 1, volcano$pvalue)

  if (filter_type == "padj") {
    volcano_filter_col <- volcano$padj
  } else {
    volcano_filter_col <- volcano$pvalue
  }

  volcano$regulation <- as.factor(ifelse(volcano_filter_col < filter_value & abs(volcano$log2FoldChange) >= bs_pos, ifelse(volcano$log2FoldChange >= bs_pos, "Up", "Down"), "NoSignificant"))
  volcano$FC <- 2^volcano$log2FoldChange
  fpkm.tmp <- fpkm.deg[rownames(volcano), ]
  colnames(fpkm.tmp) <- paste(colnames(fpkm.tmp), "_FPKM", sep = "")
  volcano <- cbind(volcano, fpkm.tmp)
  volcano <- cbind(as.data.frame(rownames(volcano)), volcano)
  colnames(volcano)[1] <- "GeneID"
  write.table(volcano, file = outfile, sep = "	", quote = FALSE, row.names = F)
  total.deg <- c(total.deg, rownames(volcano)[volcano$regulation == "Up" | volcano$regulation == "Down"])
  total_deg_num <- nrow(volcano[volcano$regulation == "Up" | volcano$regulation == "Down", ])
  up_deg_num <- nrow(volcano[volcano$regulation == "Up", ])
  down_deg_num <- nrow(volcano[volcano$regulation == "Down", ])
  deg_stat_group <- paste(comp_info[i, 1], comp_info[i, 2], sep = "_vs_")
  stat.deg <- rbind(stat.deg, rbind(c(deg_stat_group, total_deg_num, up_deg_num, down_deg_num)))
  if (up_deg_num == 0 | down_deg_num == 0) {
    next
  }
  p.volcano <- ggplot(data = volcano, aes(x = log2FoldChange, y = -log10(padj), colour = regulation)) +
    geom_point() +
    scale_color_manual(values = c("green", "grey", "red")) +
    geom_vline(xintercept = c(bs_neg, bs_pos), lty = 4, col = "black", lwd = 0.8) +
    geom_hline(yintercept = -log10(0.05), lty = 4, col = "black", lwd = 0.8) +
    theme_base()
  outfile.volcano <- paste0(deg_exp_graph_dir, group_vs_group_name, "_volcano.jpeg")
  ggsave(outfile.volcano, p.volcano, dpi = 300, width = 10, height = 10)
  tmp.deg <- as.character(rownames(volcano)[volcano$regulation == "Up" | volcano$regulation == "Down"])
  fpkm_tmp <- na.omit(fpkm.deg[tmp.deg, ])
  fpkm_tmp <- log2(fpkm_tmp[rowSums(fpkm_tmp) > 0, ] + 1)
  p.tmpdeg.heatmap <- pheatmap(fpkm_tmp, scale = "row", cluster_cols = F, show_rownames = F)
  outfile.tmpdeg.heatmap <- paste0(deg_exp_graph_dir, group_vs_group_name, "_heatmap.jpeg")
  ggsave(outfile.tmpdeg.heatmap, p.tmpdeg.heatmap, dpi = 300, width = 10, height = 10)
  outfile.Up_ID <- paste0(deg_dir, group_vs_group_name, "_Up_ID.txt")
  UP_ID.dt <- data.frame(GeneID = as.character(rownames(volcano)[volcano$regulation == "Up"]))
  write.table(UP_ID.dt, file = outfile.Up_ID, sep = "	", quote = FALSE, row.names = F, col.names = FALSE)
  outfile.Down_ID <- paste0(deg_dir, group_vs_group_name, "_Down_ID.txt")
  Down_ID.dt <- data.frame(GeneID = as.character(rownames(volcano)[volcano$regulation == "Down"]))
  write.table(Down_ID.dt, file = outfile.Down_ID, sep = "	", quote = FALSE, row.names = F, col.names = FALSE)
}
rm(i)
############ draw heatmap #############
total.deg <- as.character(unique(total.deg))
# length(total.deg)

deg_fpkm <- na.omit(fpkm[total.deg, ])

deg_fpkm <- log2(deg_fpkm[rowSums(deg_fpkm) > 0, ] + 1)
p.heatmap <- pheatmap(deg_fpkm, scale = "row", cluster_cols = F, show_rownames = F)
ggsave("deg_heatmap.jpeg", p.heatmap, dpi = 300, width = 10, height = 10)
############ draw boxplot density########

melt(fpkm, variable.name = "sample", value.name = "fpkm") -> data.m
data.m[data.m[, 2] > 0, ] -> data.m
p <- ggplot(data.m, aes(x = sample, y = log2(fpkm), fill = sample)) +
  geom_boxplot() +
  theme_base() +
  ylab("log2(FPKM)") +
  xlab("Sample") +
  theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), axis.text.x = element_text(size = rel(1.5), angle = 90, hjus = 1, vjust = .5))
ggsave(filename = "Expression_data_evaluation/fpkm_boxplot.jpeg", plot = p, height = 10, width = 16.8, dpi = 300)
d <- ggplot(data.m, aes(x = log10(fpkm), col = sample)) +
  geom_density(aes(fill = sample), colour = NA, alpha = .2) +
  geom_line(stat = "density", size = 1.5) +
  xlab("log2(FPKM)") +
  theme_base() +
  theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), axis.text.x = element_text(size = rel(3)), axis.text.y = element_text(size = rel(3)))
ggsave(filename = "Expression_data_evaluation/fpkm_density.jpeg", plot = d, height = 10, width = 16.8, dpi = 300)

# install.packages("ggcorrplot")
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

# 假设fpkm.cor是一个100行150列的矩阵
# nrow(fpkm.cor)  # 100
# ncol(fpkm.cor)  # 150

cell_size <- 260

png("Expression_data_evaluation/correlation.jpeg", res = 72 * 5, 
    height = (nrow(fpkm.cor) * cell_size) + 1000, width = (ncol(fpkm.cor) * cell_size) + 1000) 
corrplot(fpkm.cor, is.corr = F, col = rev(COL2("PiYG")), method = "color", 
         addCoef.col = "black", tl.col = "black", col.lim = c(min(fpkm.cor) - 0.01, max(fpkm.cor)), 
         cl.ratio = 0.1)
dev.off()

# pdf("correlation.pdf", res = 72 * 5, 
#     height = (nrow(fpkm.cor) * cell_size) + 1000, width = (ncol(fpkm.cor) * cell_size) + 1000) 
# corrplot(fpkm.cor, is.corr = F, col = rev(COL2("PiYG")), method = "color", addCoef.col = "black", tl.col = "black", col.lim = c(min(fpkm.cor) - 0.01, max(fpkm.cor)), cl.ratio = 0.1)
# dev.off()


# head(deg_fpkm)
# stat.deg
# as.data.frame(matrix(stat.deg,nrow=nrow(comp_info),ncol=4,byrow=T),colnames=c("Groups","Total_DEGs","UP"))
# stat.deg
colnames(stat.deg) <- c("Comparisons", "Total DEGs", "Up regulated", "Down regulated")

deg_first_line <- paste0("# 筛选条件：",filter_type, " < ",filter_value,"; FoldChange > ",deg_value,"\n")
# write.lines(deg_first_line, paste0(deg_dir, "DEG_summary.txt"))
cat(deg_first_line, file = paste0(deg_dir, "DEG_summary.txt"))

write.table(stat.deg, file = paste0(deg_dir, "DEG_summary.txt"), sep = "\t", quote = F, row.names = F, append = TRUE)


##############################################################################################################
# setwd("d:/pll/R_work/anova")
# getwd()
data <- read.table(fpkm_file, sep = "\t", row.names = 1, header = T, check.names = F)
group <- read.table(samples_file, sep = "\t", header = T, check.names = F)
head(group)

colnames(data)
# match(colnames(data),group$sample)
# order(match(colnames(data),group$sample))
# group[match(colnames(data),group$sample),]
tmp.group <- group$group[match(colnames(data), group$sample)]
# tmp.group
# class(tmp.group)
# head(data)
p <- NULL
for (i in rownames(data)) {
  # print(i);
  # print(data[i,])
  tmp.count <- as.numeric(data[i, ])
  # tmp.dt<-data.frame(tmp.count,group=rep(c("S1","S3","S2"),each=3))
  tmp.dt <- data.frame(tmp.count, group = tmp.group)
  colnames(tmp.dt)[1] <- i
  tmp.dt$group <- as.factor(tmp.dt$group)
  mod <- aov(tmp.dt[, 1] ~ group, data = tmp.dt)
  p <- c(p, summary(mod)[[1]][5]$Pr[1])
}
length(p)
head(p)
data$p_value <- p
data <- data[order(data$p_value), ]
p.adjust(data$p_value, method = "BH") -> p.aj
data$BH_p_value <- p.aj
df.rname <- data.frame(GeneID = rownames(data))
data <- cbind(df.rname, data)
write.table(data, "anova_analysis_p.txt", sep = "\t", quote = F, row.names = F)
# dim(data)
# rep(c("S1","S2","S3"),each=3)
# as.numeric(data[1,])
# tmp.count<-as.numeric(data[rownames(data)[1],])
# tmp.dt<-data.frame(tmp.count,group=rep(c("S1","S2","S3"),each=3))
# colnames(tmp.dt)[1]<-rownames(data)[1]
# tmp.dt
# tmp.dt$group<-as.factor(tmp.dt$group)
# mod <- aov( ENSMUSG00000029661~group , data = tmp.dt)
# haha<-summary(mod)
# summary(mod)[[1]][5]$Pr[1]
#####################################################################

# get the log2 transformed FPKM values
gene <- read.delim(fpkm_file, row.names = 1, sep = "\t", check.names = FALSE, header = T)
gene <- log2(gene + 1)
# transform the values to the same scale
gene <- t(gene)

# PCA plot
gene.pca <- PCA(gene, ncp = 2, scale.unit = TRUE, graph = FALSE)

# get the PCA coordinates
pca_sample <- data.frame(gene.pca$ind$coord[, 1:2])
head(pca_sample)

pca_eig1 <- round(gene.pca$eig[1, 2], 2)
pca_eig2 <- round(gene.pca$eig[2, 2], 2)

# get the group information
group <- read.delim(samples_file, row.names = 2, sep = "\t", check.names = FALSE, header = T)
group <- group[rownames(pca_sample), ]

p <- ggplot(data = pca_sample, aes(x = Dim.1, y = Dim.2)) +
  geom_point(aes(color = group), size = 5) +
  theme(
    panel.grid = element_blank(), panel.background = element_rect(color = "black", fill = "transparent"),
    legend.key = element_rect(fill = "transparent")
  ) +
  labs(x = paste("PCA1:", pca_eig1, "%"), y = paste("PCA2:", pca_eig2, "%"), color = "")
p <- p + geom_text_repel(data = pca_sample, aes(Dim.1, Dim.2, label = rownames(pca_sample)))
cluster_border <- ddply(pca_sample, "group", function(df) df[chull(df[[1]], df[[2]]), ])
p <- p + geom_polygon(data = cluster_border, aes(color = group), fill = NA, show.legend = FALSE)
ggsave("Expression_data_evaluation/PCA.jpeg", p, dpi = 300, width = 10, height = 10)

