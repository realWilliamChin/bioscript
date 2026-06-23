library(ggplot2)  
# setwd("d:/pll/R_work/tmp/")
rm(list=ls())
library(vegan)  
library(ggplot2)  
library(optparse)

# 创建参数解析器
option_list <- list(
  make_option(c("-i", "--input"), 
              type = "character", 
              default = "Species.txt",
              help = "输入文件路径（支持 txt 和 csv 格式），默认为 Species.txt",
              metavar = "FILE"),
  make_option(c("-o", "--output"),
              type = "character",
              default = ".",
              help = "输出目录路径，默认为当前目录",
              metavar = "DIR")
)

parser <- OptionParser(option_list = option_list,
                       description = "物种累积曲线分析脚本")

# 解析参数
args <- parse_args(parser)
input_file <- args$input
# --output 为输出目录，图片文件名固定为 species_accumulation_curve.png
output_dir <- args$output
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
output_jpg <- file.path(output_dir, "species_accumulation_curve.png")

# 根据文件扩展名选择读取方法
file_ext <- tolower(tools::file_ext(input_file))
if(file_ext == "csv") {
  otu <- read.csv(input_file, row.names = 1, stringsAsFactors = FALSE, check.names = FALSE, header = TRUE)
} else {
  # 默认为 txt (tab分隔)
  otu <- read.delim(input_file, row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE, header = TRUE)
}
otu=t(otu)
nrow(otu)
sp <- specaccum(otu, method = 'random')
sp.min<-200*trunc(min((sp$richness)/200))
sp.max<-200*ceiling(max((sp$richness)/200))
seq(sp.min,sp.max,by=200)
summary(sp)
png(output_jpg, width = 2000, height = 1600, res=320)
plot(
  sp, ci.type = 'poly', col = 'blue', lwd = 2, ci.lty = 0, ci.col = 'white',
  bty="o",
  xlim = c(0, max(nrow(otu))),
  ylim=c(sp.min,sp.max),
  xaxt = "n",
  yaxt="n",
  xlab="Sample Number",
  ylab="ASV Count")
axis(1, at = 1:nrow(otu), labels = 1:nrow(otu))
axis(2, at = seq(sp.min,sp.max,by=200), labels = seq(sp.min,sp.max,by=200))
boxplot(sp, col = '#3f51b5', add = TRUE, pch = '+')
box()
dev.off()

