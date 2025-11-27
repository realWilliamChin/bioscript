library(QTLseqr)
library(ggplot2)
rm(list=ls())

prefix = 'yumi'

HighBulk <- "Mut"
LowBulk <- "WT"
file <- "BSA_merged_snp.table" ##set sample and file name

Chroms <- paste0(rep("", 10), 1:10) ##choose which chromosomes will be included in the analysis 

Chroms <- read.table('chroms_list.txt', header=FALSE)[[1]]

df <-
  importFromGATK(
    file = file,
    highBulk = HighBulk,
    lowBulk = LowBulk,
    chromList = Chroms) ##Import SNP data from file


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
    refAlleleFreq = 0.001,
    minTotalDepth = 5,
    maxTotalDepth = 100,
    depthDifference = 24,
    minSampleDepth = 10,
    minGQ = 10) ##Filter SNPs based on some criteria


df_filt_G <- runGprimeAnalysis(
  SNPset = df_filt,
  windowSize = 4e6,   ##如果程序报错，就增大窗口大小
  outlierFilter = "deltaSNP") ##Run G' analysis

#write.csv(df_filt_G, file="ZZW-MAX-MIN_G.csv", row.names = FALSE)

df_filt_Q <- runQTLseqAnalysis(
  SNPset = df_filt_G,
  windowSize = 4e6,   ##如果程序报错，就增大窗口大小
  popStruc = "F2",
  bulkSize = c(24, 1),
  replications = 10000,
  intervals = c(95, 99)) ##Run QTLseq analysis

write.csv(df_filt_Q, file=paste0(prefix,"_QTL.csv"), row.names = FALSE)


# plotGprimeDist(SNPset = df_filt_G, outlierFilter = "Hampel")

gprime_plot <- plotQTLStats(SNPset = df_filt_G, var = "Gprime", plotThreshold = TRUE, q = 0.55)  ##如果没有阈值线出现，就调节q的数值
gprime_plot
ggsave(paste0(prefix, "_Gprime_q0.55.jpeg"), gprime_plot, width = 10, height = 8, units = "in", dpi = 300)
# results <- getQTLTable(SNPset = df_filt_G, method = "Gprime",alpha = 0.05, export = FALSE)


deltasnp_plot <- plotQTLStats(SNPset = df_filt_Q, var = "deltaSNP", plotIntervals = TRUE) ##Plot
ggsave(paste0(prefix, "_deltaSNP.jpeg"), deltasnp_plot, width = 10, height = 8, units = "in", dpi = 300)

getQTLTable(SNPset = df_filt_Q, method="Gprime", alpha = 0.55, export = TRUE, fileName = paste0(prefix,"_G_sig_chr_region.csv")) ##export summary CSV   ##如果没有显著QTL生成，就调节alpha的数值

