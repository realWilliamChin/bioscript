library(QTLseqr)
library(ggplot2)
rm(list=ls())

# 【项目模板脚本】下面的 prefix/HighBulk/LowBulk/file 等参数是按具体项目硬编码的，
# 每次新项目使用时需要修改这些参数为对应混池名称与文件路径，不是通用 CLI 工具。
prefix = 'Nacl'
HighBulk <- "32Nacl"
LowBulk <- "33Nacl"
file <- "snp.table"

# 染色体信息列表
Chroms <- paste0(rep("", 10), 1:10)
Chroms <- read.table('chroms_list.txt', header=FALSE)[[1]]

df <-
  importFromGATK(
    file = file,
    highBulk = HighBulk,
    lowBulk = LowBulk,
    chromList = Chroms)

gp_high_plot <- ggplot(data = df) +
  geom_histogram(aes(x = DP.HIGH)) +
  xlim(0,200)
ggsave(paste0(prefix, "_GP_HIGH.png"), gp_high_plot, width = 10, height = 8, units = "in", dpi = 300)

gp_low_plot <- ggplot(data = df) +
  geom_histogram(aes(x = DP.LOW)) +
  xlim(0,200)
ggsave(paste0(prefix, "_GP_LOW.png"), gp_low_plot, width = 10, height = 8, units = "in", dpi = 300)

# SNP 质检，决定后面 df filt 的参数
dp_plot <- ggplot(data = df) +
  geom_histogram(aes(x = DP.HIGH + DP.LOW )) +
  xlim(0, 100)
ggsave(paste0(prefix, "_DP.png"), dp_plot, width = 10, height = 8, units = "in", dpi = 300)

gq_high_plot <- ggplot(data = df) +
  geom_histogram(aes(x = GQ.HIGH)) +
  xlim(0,150)
ggsave(paste0(prefix, "_GQ_HIGH.png"), gq_high_plot, width = 10, height = 8, units = "in", dpi = 300)
gq_low_plot <- ggplot(data = df) +
  geom_histogram(aes(x = GQ.LOW))
ggsave(paste0(prefix, "_GQ_LOW.png"), gq_low_plot, width = 10, height = 8, units = "in", dpi = 300)

ref_frq_plot <- ggplot(data = df) +
  geom_histogram(aes(x = REF_FRQ))
  # xlim(0,0.1) +
  # ylim(0, 100)
ggsave(paste0(prefix, "_REF_FRQ.png"), ref_frq_plot, width = 10, height = 8, units = "in", dpi = 300)

snp_index_hight_plot <- ggplot(data = df) +
  geom_histogram((aes(x = SNPindex.HIGH)))
  # xlim(0,1) +
  # ylim(0,1)
ggsave(paste0(prefix, "_SNPindex_HIGH.png"), snp_index_hight_plot, width = 10, height = 8, units = "in", dpi = 300)

snp_index_low_plot <- ggplot(data = df) +
  geom_histogram((aes(x = SNPindex.LOW)))
ggsave(paste0(prefix, "_SNPindex_LOW.png"), snp_index_low_plot, width = 10, height = 8, units = "in", dpi = 300)


df_filt <-
  filterSNPs(
    SNPset = df,
    refAlleleFreq = 0.05,  # 认为可变位点的比例
    minTotalDepth = 40,      # 看图
    maxTotalDepth = 85,    # 看图
    depthDifference = 30,   # 关注 HighBulk 的测序深度
    minSampleDepth = 30,    # vcf 文件读出来的
    minGQ = 99)             # SNP 位点的质量数


df_filt_G <- runGprimeAnalysis(
  SNPset = df_filt,
  windowSize = 4e6,             # 如果程序报错，就增大窗口大小
  outlierFilter = "deltaSNP")   # Run G' analysis

#write.csv(df_filt_G, file="ZZW-MAX-MIN_G.csv", row.names = FALSE)

df_filt_Q <- runQTLseqAnalysis(
  SNPset = df_filt_G,
  windowSize = 4e6,         # 如果程序报错，就增大窗口大小
  popStruc = "F2",
  bulkSize = c(30, 1),      # 依据混池样本数量重新填写
  replications = 10000,
  intervals = c(95, 99))    # Run QTLseq analysis

# 画图使用的所有 SNP 位点
write.csv(df_filt_Q, file=paste0(prefix,"_QTL_snp_all.csv"), row.names = FALSE)

# Gprime Plot
# plotGprimeDist(SNPset = df_filt_G, outlierFilter = "Hampel")
q_value <- 0.01  # 根据丰图调整 q，达到预期的 Gprime 值
gprime_plot <- plotQTLStats(
  SNPset = df_filt_G,
  var = "Gprime",
  plotThreshold = TRUE,
  q = q_value
)
gprime_plot
ggsave(paste0(prefix, "_Gprime_q", q_value, ".png"), gprime_plot, width = 50, height = 20, units = "in", dpi = 300, limitsize = FALSE)

# deltaSNP Plot
deltasnp_plot <- plotQTLStats(SNPset = df_filt_Q, var = "deltaSNP", plotIntervals = TRUE)
deltasnp_plot
ggsave(paste0(prefix, "_deltaSNP.png"), deltasnp_plot, width = 50, height = 20, units = "in", dpi = 300, limitsize = FALSE)

# 显著 qtl 区域表单
getQTLTable(SNPset = df_filt_Q, method="Gprime", alpha = q_value, export = TRUE, fileName = paste0(prefix,"_G_sig_qtl_region.csv"))

