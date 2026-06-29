#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# 功能：从VCF文件直接绘制1Mb窗口SNP密度图
# 使用方法：Rscript snp_density_from_vcf.r -i input.vcf -o output_prefix

library(optparse)
library(CMplot)
library(vcfR)

# 解析命令行参数
option_list <- list(
  make_option(c("-i", "--input"), type="character", help="输入VCF文件路径", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="snp_density", help="输出文件前缀 [default: %default]", metavar="character"),
  make_option(c("-w", "--window"), type="integer", default=1e6, help="窗口大小(bp) [default: %default]", metavar="integer"),
  make_option(c("--dpi"), type="integer", default=300, help="图片分辨率 [default: %default]", metavar="integer"),
  make_option(c("--width"), type="integer", default=14, help="图片宽度(inch) [default: %default]", metavar="integer"),
  make_option(c("--height"), type="integer", default=10, help="图片高度(inch) [default: %default]", metavar="integer")
)

parser <- OptionParser(option_list=option_list)
args <- parse_args(parser)

# 检查输入参数
if (is.null(args$input)) {
  print_help(parser)
  stop("必须提供输入VCF文件路径")
}

# 读取VCF文件
message("正在读取VCF文件: ", args$input)
vcf <- read.vcfR(args$input, verbose = FALSE)

# 提取SNP信息
message("正在处理SNP信息...")
snp_data <- as.data.frame(getFIX(vcf))

# 格式化数据，满足CMplot要求
# CMplot需要三列：SNP名称, 染色体, 位置
cm_data <- data.frame(
  SNP = ifelse(!is.na(snp_data$ID) & snp_data$ID != ".", snp_data$ID, paste0(snp_data$CHROM, ":", snp_data$POS)),
  Chrom = snp_data$CHROM,
  Pos = as.numeric(snp_data$POS)
)

# 过滤掉非正常染色体（可选，根据需要调整）
# 支持人类染色体格式chr1-chr22,chrX,chrY,chrMT或者1-22,X,Y,MT
valid_chr <- c(paste0("chr", c(1:22, "X", "Y", "MT")), as.character(c(1:22, "X", "Y", "MT")))
cm_data <- cm_data[cm_data$Chrom %in% valid_chr, ]

# 转换染色体编号为数字格式（保证绘图顺序正确）
cm_data$Chrom <- gsub("chr", "", cm_data$Chrom, ignore.case = TRUE)
cm_data$Chrom <- factor(cm_data$Chrom, levels = c(as.character(1:22), "X", "Y", "MT"))

message("总共有 ", nrow(cm_data), " 个有效SNP用于绘图")

# 绘制SNP密度图
message("正在绘制SNP密度图...")
CMplot(cm_data,
       plot.type = "d",
       bin.size = args$window,
       col = c("darkgreen", "yellow", "red"),  # 和原图配色一致
       file = "jpg",
       file.name = args$output,
       dpi = args$dpi,
       width = args$width,
       height = args$height,
       file.output = TRUE,
       verbose = TRUE,
       main = paste0("The number of SNPs within ", args$window/1e6, "Mb window size")
)

message("绘图完成！输出文件: ", args$output, ".jpg")