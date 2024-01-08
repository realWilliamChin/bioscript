library(ggthemes)
library(ggplot2)
library(pheatmap)
library(reshape2)
library(ggcorrplot)
library(corrplot)
library(DESeq2)
library(edgeR)
library(FactoMineR)
library(ggrepel)
library(plyr)
library(dplyr)
library(tidyr)
rm(list=ls())

# samples_file <- 'samples_described.txt'
# comp_file <- 'compare_info.txt'
# reads_file <- 'reads_matrix_filtered.txt'
# fpkm_file <- 'fpkm_matrix_filtered.txt'
raw_all_deg_summary <- function(samples_file, comp_file, reads_file, fpkm_file) {
  reads_data<-read.table(reads_file,sep="\t",row.names=1,header=T,check.names=F,stringsAsFactors = F)
  reads_data<-na.omit(reads_data)
  comp_info<-read.table(comp_file,sep="\t",header=T,check.names=F,stringsAsFactors = F)
  sample_info<-read.table(samples_file,sep="\t",header=T,check.names=F,stringsAsFactors = F)
  read.table(fpkm_file,sep="\t",header=T,row.names=1,check.names=F) -> fpkm
  deg_numbers <- c(1, 1.1, 1.2, 1.3, 1.4, 1.5, 2, 3, 4)
  deg_count_summary <- data.frame("group"=character(0), "deg1"=numeric(0), "deg1.1"=numeric(0), "deg1.2"=numeric(0), "deg1.3"=numeric(0),
                                  "deg1.4"=numeric(0), "deg1.5"=numeric(0), "deg2"=numeric(0), "deg3"=numeric(0), "deg4"=numeric(0))
  for(i in seq_along(1:nrow(comp_info))){
    group_vs_group_name <- paste(comp_info[i,1],comp_info[i,2],sep="_vs_")
    
    # 分出每组比较的 reads 和 fpkm 表
    data.treat <- reads_data[,sample_info[sample_info$group == comp_info[i,1],]$sample]
    data.control <- reads_data[,sample_info[sample_info$group == comp_info[i,2],]$sample]
    rnaseqMatrix <- cbind(data.treat,data.control)
    
    fpkm.treat <- fpkm[,sample_info[sample_info$group == comp_info[i,1],]$sample]
    fpkm.control <- fpkm[,sample_info[sample_info$group == comp_info[i,2],]$sample]
    fpkm.deg <- cbind(fpkm.treat,fpkm.control)
    
    rnaseqMatrix = round(rnaseqMatrix)
    rnaseqMatrix = rnaseqMatrix[rowSums(cpm(rnaseqMatrix) > 1) >= 2,] # 过滤掉表达水平较低的基因
    
    # 差异分析
    conditions = data.frame(conditions=factor(c(rep(comp_info[i,1], ncol(data.treat)), rep(comp_info[i,2], ncol(data.control)))))
    rownames(conditions) = colnames(rnaseqMatrix)
    ddsFullCountTable <- DESeqDataSetFromMatrix(
      countData = rnaseqMatrix,
      colData = conditions,
      design = ~ conditions)
    dds = DESeq(ddsFullCountTable)
    contrast = c("conditions",comp_info[i,1],comp_info[i,2])
    res = results(dds, contrast)
    
    # DE_results
    baseMeanA <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$conditions == comp_info[i,1]])
    baseMeanB <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$conditions == comp_info[i,2]])
    res = cbind(baseMeanA, baseMeanB, as.data.frame(res))
    res = cbind(sampleA=comp_info[i,1], sampleB=comp_info[i,2], as.data.frame(res))
    res$padj[is.na(res$padj)]  <- 1
    res = as.data.frame(res[order(res$pvalue),])
    rnaseqMatrix <- cbind(as.data.frame(rownames(rnaseqMatrix)), rnaseqMatrix)
    colnames(rnaseqMatrix)[1]<-"GeneID"
    
    volcano <- res
    volcano$padj <- ifelse(volcano$padj<0.000000000000001,0.000000000000001,volcano$padj)
    #volcano$regulation = as.factor(ifelse(volcano$padj < 0.05 & abs(volcano$log2FoldChange) >= bs_pos, ifelse(volcano$log2FoldChange >=bs_pos ,'Up','Down'),'NoSignificant'))
    volcano$FC = 2^volcano$log2FoldChange
    fpkm.tmp <- fpkm.deg[rownames(volcano),]
    colnames(fpkm.tmp) <- paste(colnames(fpkm.tmp),'_FPKM',sep='')
    volcano <- cbind(volcano,fpkm.tmp)
    volcano <- cbind(as.data.frame(rownames(volcano)),volcano)
    colnames(volcano)[1] <- "GeneID"
    
    result_numbers <- c()
    for (fc_num in deg_numbers) {
      # 计算大于 fc 的基因数，padj < 0.05
      gtfc_num <- nrow(volcano[abs(volcano$log2FoldChange) >= log2(fc_num) & volcano$padj < 0.05, ])
      result_numbers <- c(result_numbers, gtfc_num)
    }
    new_row <- data.frame("group"=group_vs_group_name, "deg1"=result_numbers[1],
                          "deg1.1"=result_numbers[2], "deg1.2"=result_numbers[3],
                          "deg1.3"=result_numbers[4], "deg1.4"=result_numbers[5],
                          "deg1.5"=result_numbers[6], "deg2"=result_numbers[7],
                          "deg3"=result_numbers[8], "deg4"=result_numbers[9])
    deg_count_summary <- rbind(deg_count_summary, new_row)
  }

  # 画 deg_count_summary 的点线图

  row.names(deg_count_summary) <- deg_count_summary$group
  deg_count_summary$group <- NULL

  # 这里我们将行名转换为一个普通的列，假设行名代表样本名
  deg_count_summary <- data.frame(group = rownames(deg_count_summary), deg_count_summary)

  # 将数据从宽格式转换为长格式
  long_data <- pivot_longer(
    deg_count_summary, 
    cols = -group, 
    names_to = "condition", 
    values_to = "value"
  )
  # 获取每个样本第一个数据点的值
  first_points <- long_data %>%
    group_by(group) %>%
    slice(1) %>%
    ungroup() %>%
    arrange(value)

  # 获取排序后样本名的向量
  ordered_sample_names <- first_points$group

  # 找到每个样本最后一个数据点的位置
  last_points <- long_data %>%
    group_by(group) %>%
    summarize(Condition = last(condition), Value = last(value))

  # 绘制点线图
  p <- ggplot(long_data, aes(x = condition, y = value, group = group)) +
    geom_line(aes(color = group)) +
    geom_point(aes(color = group)) +
    theme_minimal() +
    labs(x = "Condition", y = "Value") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_discrete(limits = unique(long_data$condition)) # 设置X轴的刻度

  # 确定图像的尺寸
  width <- 11
  height <- 8.5

  # 保存图片
  ggsave("deg_line_plot.png", plot = p, width = width, height = height, dpi = 300)
}
draw_all_deg_summary('samples_described.txt', 'compare_info.txt', 'reads_matrix_filtered.txt', 'fpkm_matrix_filtered.txt')
