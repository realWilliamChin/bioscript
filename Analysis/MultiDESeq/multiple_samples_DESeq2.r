pkgs <- c(
  'ggthemes', 'ggplot2', 'pheatmap', 'reshape2', 'ggcorrplot',
  'corrplot', 'DESeq2', 'edgeR', 'FactoMineR', 'ggrepel',
  'plyr', 'optparse'
)
suppressPackageStartupMessages(invisible(lapply(pkgs, require, character.only = TRUE)))

source('/home/colddata/qinqiang/script/Plot/Heatmap/heatmap_1.r', echo = TRUE, encoding = 'UTF-8')
source('/home/colddata/qinqiang/script/Plot/PCA/pca_1.r', echo = TRUE, encoding = 'UTF-8')
source('/home/colddata/qinqiang/script/Plot/Corrplot/correlation_1.r', echo = TRUE, encoding = 'UTF-8')


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
    help = "提供 samples_described.txt 文件", metavar = "character"
  ),
  make_option(c("--compare"),
    type = "character", default = "compare_info.txt",
    help = "提供 compare_info.txt 文件", metavar = "character"
  ),
  make_option(c("--filtercol"),
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
# filter_col <- 'padj'
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
filter_col <- opt$filtercol
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
  png(file.path(exp_evaluation_dir, 'All_gene_heatmap.png'), width = 10, height = 10, units = "in", res = 300)
  all.heatmap <- pheatmap(subset_data, scale = "row", cluster_cols = FALSE, show_rownames = FALSE)
  dev.off()
} else {
  # 如果数据不超过 65535 行，直接创建热图
  png(file.path(exp_evaluation_dir, 'All_gene_heatmap.png'), width = 10, height = 10, units = "in", res = 300)
  all.heatmap <- pheatmap(all_fpkm, scale = "row", cluster_cols = FALSE, show_rownames = FALSE)
  dev.off()
}

sample_info <- read.table(samples_file, sep = "\t", header = T, check.names = F, stringsAsFactors = F)
# sample_info[sample_info$group=="CK",]
reads_data <- read.table(reads_file, sep = "\t", row.names = 1, header = T, check.names = F, stringsAsFactors = F)
reads_data <- na.omit(reads_data)
comp_info <- read.table(compare_file, sep = "\t", header = T, check.names = F, stringsAsFactors = F) 


total.deg <- ""
stat.deg <- data.frame()
for (i in seq_along(1:nrow(comp_info))) {
  group_vs_group_name <- paste(comp_info[i, 1], comp_info[i, 2], sep = "-vs-")  
  # 从reads_data中提取处理组和对照组的样本数据
  data.treat <- reads_data[, sample_info[sample_info$group == comp_info[i, 1], ]$sample]
  data.control <- reads_data[, sample_info[sample_info$group == comp_info[i, 2], ]$sample]
  
  # 合并处理组和对照组数据，构建RNA-seq矩阵
  rnaseqMatrix <- cbind(data.treat, data.control)
  
  # 从fpkm数据中提取对应的处理组和对照组样本数据（用于后续可视化）
  fpkm.treat <- fpkm[, sample_info[sample_info$group == comp_info[i, 1], ]$sample]
  fpkm.control <- fpkm[, sample_info[sample_info$group == comp_info[i, 2], ]$sample]
  fpkm.deg <- cbind(fpkm.treat, fpkm.control)
  
  # 将reads数据四舍五入为整数（DESeq2要求count数据为整数）
  rnaseqMatrix <- round(rnaseqMatrix)
  
  # 过滤低表达基因：保留在至少2个样本中CPM>1的基因
  rnaseqMatrix <- rnaseqMatrix[rowSums(cpm(rnaseqMatrix) > 1) >= 2, ]
  
  # 构建样本分组信息，用于DESeq2分析
  conditions <- data.frame(conditions = factor(c(rep(comp_info[i, 1], ncol(data.treat)), rep(comp_info[i, 2], ncol(data.control)))))
  # conditions
  rownames(conditions) <- colnames(rnaseqMatrix)
  
  # 创建DESeq2数据集对象
  ddsFullCountTable <- DESeqDataSetFromMatrix(
    countData = rnaseqMatrix,      # 输入reads count数据
    colData = conditions,          # 样本分组信息
    design = ~conditions           # 实验设计公式
  )
  
  # 运行DESeq2差异表达分析
  dds <- DESeq(ddsFullCountTable)
  
  # 设置比较对比：处理组 vs 对照组
  contrast <- c("conditions", comp_info[i, 1], comp_info[i, 2])
  
  # 获取差异表达分析结果
  res <- results(dds, contrast)
  
  # 计算处理组和对照组的标准化表达量均值
  baseMeanA <- rowMeans(counts(dds, normalized = TRUE)[, colData(dds)$conditions == comp_info[i, 1]])
  baseMeanB <- rowMeans(counts(dds, normalized = TRUE)[, colData(dds)$conditions == comp_info[i, 2]])
  
  # 将均值信息添加到结果中
  res <- cbind(baseMeanA, baseMeanB, as.data.frame(res))
  
  # 添加样本组别信息到结果中
  res <- cbind(sampleA = comp_info[i, 1], sampleB = comp_info[i, 2], as.data.frame(res))
  
  # 将NA的padj值设为1（避免后续分析出错）
  res$padj[is.na(res$padj)] <- 1
  
  # 按p值排序结果
  res <- as.data.frame(res[order(res$pvalue), ])
  
  # 设置输出文件名
  outfile <- paste0(group_vs_group_name, "_DE_results")
  outfile.ma <- paste0(group_vs_group_name, "_DE_results_readCounts.matrix")
  
  # 为reads count矩阵添加基因ID列，准备输出
  rnaseqMatrix <- cbind(as.data.frame(rownames(rnaseqMatrix)), rnaseqMatrix)
  colnames(rnaseqMatrix)[1] <- "GeneID"
  
  # 输出reads count矩阵到文件
  write.table(rnaseqMatrix, file = outfile.ma, sep = "	", quote = FALSE, row.names = F)
  
  # 准备火山图数据
  volcano <- res
  
  # 处理极小的padj值，避免火山图显示问题
  volcano$padj <- ifelse(volcano$padj < 0.000000000000001, 0.000000000000001, volcano$padj)
  volcano$pvalue <- ifelse(is.na(volcano$pvalue), 1, volcano$pvalue)

  # 为火山图选择统一的显著性列，并处理NA/非正值，设置阈值alpha
  volcano[[filter_col]][is.na(volcano[[filter_col]])] <- 1
  volcano[[filter_col]][volcano[[filter_col]] <= 0] <- .Machine$double.xmin
  alpha <- filter_value

  # 根据过滤标准和log2FoldChange阈值，标记基因的调控状态
  # Up: 上调基因，Down: 下调基因，NoSignificant: 无显著差异
  volcano$regulation <- as.factor(ifelse(volcano[[filter_col]] < filter_value & abs(volcano$log2FoldChange) >= bs_pos, ifelse(volcano$log2FoldChange >= bs_pos, "Up", "Down"), "NoSignificant"))
  
  # 计算倍数变化（FC = 2^log2FoldChange）
  volcano$FC <- 2^volcano$log2FoldChange
  
  # 将FPKM表达量数据添加到结果中（用于结果展示）
  fpkm.tmp <- fpkm.deg[rownames(volcano), ]
  colnames(fpkm.tmp) <- paste(colnames(fpkm.tmp), "_FPKM", sep = "")
  volcano <- cbind(volcano, fpkm.tmp)
  
  # 添加基因ID列到结果中
  volcano <- cbind(as.data.frame(rownames(volcano)), volcano)
  colnames(volcano)[1] <- "GeneID"
  
  # 输出完整的差异表达分析结果到文件
  write.table(volcano, file = outfile, sep = "	", quote = FALSE, row.names = F)
  
  # 收集所有差异表达基因的ID
  total.deg <- c(total.deg, rownames(volcano)[volcano$regulation == "Up" | volcano$regulation == "Down"])
  
  # 统计差异基因数量
  total_deg_num <- nrow(volcano[volcano$regulation == "Up" | volcano$regulation == "Down", ])
  up_deg_num <- nrow(volcano[volcano$regulation == "Up", ])
  down_deg_num <- nrow(volcano[volcano$regulation == "Down", ])
  
  # 生成统计信息并添加到统计表中
  deg_stat_group <- paste(comp_info[i, 1], comp_info[i, 2], sep = "-vs-")
  stat.deg <- rbind(stat.deg, rbind(c(deg_stat_group, total_deg_num, up_deg_num, down_deg_num)))
  
  # 如果上调或下调基因数量为0，跳过后续的可视化步骤
  if (up_deg_num == 0 | down_deg_num == 0) {
    next
  }
  
  # 绘制火山图：展示差异表达基因的分布
  p.volcano <- ggplot(data = volcano, aes(x = log2FoldChange, y = -log10(.data[[filter_col]]), colour = regulation)) +
    geom_point() +
    scale_color_manual(values = c("green", "grey", "red")) +  # 绿色=下调，红色=上调，灰色=无显著差异
    geom_vline(xintercept = c(bs_neg, bs_pos), lty = 4, col = "black", linewidth = 0.8) +  # 添加log2FoldChange阈值线
    geom_hline(yintercept = -log10(alpha), lty = 4, col = "black", linewidth = 0.8) +  # 添加显著性阈值线
    theme_base()
  
  # 保存火山图
  outfile.volcano <- paste0(deg_exp_graph_dir, group_vs_group_name, "_volcano.jpeg")
  ggsave(outfile.volcano, p.volcano, dpi = 300, width = 10, height = 10)
  
  # 提取差异表达基因的FPKM数据用于热图绘制
  tmp.deg <- as.character(rownames(volcano)[volcano$regulation == "Up" | volcano$regulation == "Down"])
  fpkm_tmp <- na.omit(fpkm.deg[tmp.deg, ])
  fpkm_tmp <- log2(fpkm_tmp[rowSums(fpkm_tmp) > 0, ] + 1)  # log2(FPKM+1)转换
  
  # 绘制差异基因表达热图
  p.tmpdeg.heatmap <- pheatmap(fpkm_tmp, scale = "row", cluster_cols = F, show_rownames = F)
  
  # 保存热图
  outfile.tmpdeg.heatmap <- paste0(deg_exp_graph_dir, group_vs_group_name, "_heatmap.jpeg")
  ggsave(outfile.tmpdeg.heatmap, p.tmpdeg.heatmap, dpi = 300, width = 10, height = 10)
  
  # 输出上调基因ID列表
  outfile.Up_ID <- paste0(deg_dir, group_vs_group_name, "_Up_ID.txt")
  UP_ID.dt <- data.frame(GeneID = as.character(rownames(volcano)[volcano$regulation == "Up"]))
  write.table(UP_ID.dt, file = outfile.Up_ID, sep = "	", quote = FALSE, row.names = F, col.names = FALSE)
  
  # 输出下调基因ID列表
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
  geom_line(stat = "density", linewidth = 1.5) +
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

deg_first_line <- paste0("# 筛选条件：",filter_col, " < ",filter_value,"; FoldChange > ",deg_value,"\n")
cat(deg_first_line, file = paste0(deg_dir, "DEG_summary.txt"))
write.table(stat.deg, file = paste0(deg_dir, "DEG_summary.txt"), sep = "\t", quote = F, row.names = F, append = TRUE, col.names = TRUE)


##############################################################################################################
# setwd("d:/pll/R_work/anova")
# getwd()
data <- read.table(fpkm_file, sep = "\t", row.names = 1, header = T, check.names = F)
group <- read.table(samples_file, sep = "\t", header = T, check.names = F, stringsAsFactors = F)
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
group_info <- read.delim(samples_file, sep = "\t", check.names = FALSE, header = T)
group_info <- group_info[match(rownames(pca_sample), group_info$sample), ]

# 添加group列到pca_sample
pca_sample$group <- group_info$group

p <- ggplot(data = pca_sample, aes(x = Dim.1, y = Dim.2)) +
  geom_point(aes(color = group), size = 5) +
  theme(
    panel.grid = element_blank(), panel.background = element_rect(color = "black", fill = "transparent"),
    legend.key = element_rect(fill = "transparent")
  ) +
  labs(x = paste("PC1:", pca_eig1, "%"), y = paste("PC2:", pca_eig2, "%"), color = "")
p <- p + geom_text_repel(data = pca_sample, aes(Dim.1, Dim.2, label = rownames(pca_sample)))
cluster_border <- ddply(pca_sample, "group", function(df) df[chull(df[[1]], df[[2]]), ])
p <- p + geom_polygon(data = cluster_border, aes(color = group), fill = NA, show.legend = FALSE)
ggsave("Expression_data_evaluation/PCA.jpeg", p, dpi = 300, width = 10, height = 10)

