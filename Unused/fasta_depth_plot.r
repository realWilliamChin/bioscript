setwd("/home/colddata/qinqiang/Project/2024_01_19_xijun_unicycler_renjing/minimap2")
library(ggplot2)
library(ggthemes)
rm(list=ls())
reads_data <- read.table("S20624_depth_average.txt",sep="\t",header=T,check.names=F,stringsAsFactors = F)

# 绘制点线图
p <- ggplot(reads_data, aes(x = position, y = depth)) +
  geom_line() +
  theme_minimal() +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15)) +
  labs(x = "Position（kb）", y = "Sequenceing depth（X）")
p

# 确定图像的尺寸
width <- 11
height <- 8.5

# 保存图片
ggsave("S20624_depth_plot.jpeg", plot = p, width = width, height = height, dpi = 300)
