# setwd("/home/colddata/qinqiang/Project/2024_06_07_renminyiyuan_lixuewei_ji/联合分析")
library(ggplot2)
library(ggthemes)
library(ggrepel)

args=commandArgs(T)
input_file <- args[1]
output_file <- args[2]

data.xx<-read.table(input_file,header = T,check.names = F,quote="",sep='\t')
head(data.xx)
colnames(data.xx)

# 检查是否存在 GeneSymbol 列
if ("GeneSymbol" %in% colnames(data.xx)) {
  text_repel_layer <- geom_text_repel(
    data = data.xx,
    aes(
      x = log2FC_RNA,
      y = log2FC_protein,
      label = GeneSymbol
    ),
    size = 3,  # 减小字体大小
    color = 'black',
    box.padding = unit(1, "lines"),  # 增加标签之间的间距
    point.padding = unit(1, "lines"),  # 调整点与标签之间的距离
    segment.color = 'black',
    show.legend = FALSE,
    max.overlaps = 20,  # 限制最大重叠数
    nudge_x = 0.1,  # 微调标签的x位置
    nudge_y = 0.1,  # 微调标签的y位置
    min.segment.length = 0.2  # 控制标签与点之间的连线长度
  )
} else {
  text_repel_layer <- NULL  # 如果 GeneSymbol 列不存在，设置为 NULL
}


ggplot(
  data = data.xx,
  aes(
    x=log2FC_RNA,
    y=log2FC_protein,
    colour=Group)
  ) +
  geom_point() +
  theme_classic() +
  geom_vline(xintercept = 0.26, color = "grey", linetype = "dashed") +
  geom_vline(xintercept = -0.26, color = "grey", linetype = "dashed") +
  geom_hline(yintercept = 0.26, color = "grey", linetype = "dashed") +
  geom_hline(yintercept = -0.26, color = "grey", linetype = "dashed") +
  text_repel_layer
  

ggsave(output_file, width = 6, height = 6, dpi = 300)
#length(data.xx$Group[data.xx$Group=="group5"])
