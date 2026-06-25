library(QTLseqr)
library(ggplot2)
rm(list=ls())

prefix = 'yumi'
HighBulk <- "Mut"
LowBulk <- "WT"
file <- "BSA_merged_snp.table"

# 染色体信息列表
Chroms <- paste0(rep("", 10), 1:10)
Chroms <- read.table('chroms_list.txt', header=FALSE)[[1]]

df <-
  importFromGATK(
    file = file,
    highBulk = HighBulk,
    lowBulk = LowBulk,
    chromList = Chroms)

# SNP 质检，决定后面 df filt 的参数
dp_plot <- ggplot(data = df) +
  geom_histogram(aes(x = DP.HIGH + DP.LOW )) +
  xlim(0, 100)
ggsave(paste0(prefix, "_DP.jpeg"), dp_plot, width = 10, height = 8, units = "in", dpi = 300)

gq_high_plot <- ggplot(data = df) +
  geom_histogram(aes(x = GQ.HIGH)) +
  xlim(0,150)
ggsave(paste0(prefix, "_GQ_HIGH.jpeg"), gq_high_plot, width = 10, height = 8, units = "in", dpi = 300)

ref_frq_plot <- ggplot(data = df) +
  geom_histogram(aes(x = REF_FRQ)) +
  xlim(0,0.1) + 
  ylim(0, 100)
ggsave(paste0(prefix, "_REF_FRQ.jpeg"), ref_frq_plot, width = 10, height = 8, units = "in", dpi = 300)

snp_index_plot <- ggplot(data = df) +
  geom_histogram((aes(x = SNPindex.HIGH))) +
  xlim(0,1) +
  ylim(0,1)
ggsave(paste0(prefix, "_SNPindex_HIGH.jpeg"), snp_index_plot, width = 10, height = 8, units = "in", dpi = 300)


df_filt <-
  filterSNPs(
    SNPset = df,
    refAlleleFreq = 0.001,  # 认为可变位点的比例
    minTotalDepth = 5,      # 看图
    maxTotalDepth = 100,    # 看图
    depthDifference = 24,   # 关注 HighBulk 的测序深度
    minSampleDepth = 10,    # vcf 文件读出来的
    minGQ = 10)             # SNP 位点的质量数


df_filt_G <- runGprimeAnalysis(
  SNPset = df_filt,
  windowSize = 4e6,             # 如果程序报错，就增大窗口大小
  outlierFilter = "deltaSNP")   # Run G' analysis

#write.csv(df_filt_G, file="ZZW-MAX-MIN_G.csv", row.names = FALSE)

df_filt_Q <- runQTLseqAnalysis(
  SNPset = df_filt_G,
  windowSize = 4e6,         # 如果程序报错，就增大窗口大小
  popStruc = "F2",
  bulkSize = c(24, 1),      # 依据混池样本数量重新填写
  replications = 10000,
  intervals = c(95, 99))    # Run QTLseq analysis

# 画图使用的所有 SNP 位点
write.csv(df_filt_Q, file=paste0(prefix,"_QTL_snp_all.csv"), row.names = FALSE)

# Gprime Plot
# plotGprimeDist(SNPset = df_filt_G, outlierFilter = "Hampel")
q_value <- 0.55  # 根据丰图调整 q，达到预期的 Gprime 值
gprime_plot <- plotQTLStats(
  SNPset = df_filt_G,
  var = "Gprime",
  plotThreshold = TRUE,
  q = q_value
)
gprime_plot
ggsave(paste0(prefix, "_Gprime_q", q_value, ".jpeg"), gprime_plot, width = 10, height = 8, units = "in", dpi = 300)

# deltaSNP Plot
deltasnp_plot <- plotQTLStats(SNPset = df_filt_Q, var = "deltaSNP", plotIntervals = TRUE)
deltasnp_plot
ggsave(paste0(prefix, "_deltaSNP.jpeg"), deltasnp_plot, width = 10, height = 8, units = "in", dpi = 300)

# 显著 qtl 区域表单
getQTLTable(SNPset = df_filt_Q, method="Gprime", alpha = q_value, export = TRUE, fileName = paste0(prefix,"_G_sig_qtl_region.csv"))

