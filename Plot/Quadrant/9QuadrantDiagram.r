# setwd("/home/colddata/qinqiang/Project/2024_06_07_renminyiyuan_lixuewei_ji/联合分析")
library(ggplot2)
library(ggthemes)

args=commandArgs(T)
input_file <- args[1]
output_file <- args[2]

data.xx<-read.table(input_file,header = T,check.names = F)
head(data.xx)
colnames(data.xx)
ggplot(data.xx,aes(x=log2FC_RNA,y=log2FC_protein,colour=Group))+geom_point()+
  theme_classic()+
  geom_vline(xintercept = 0.26, color = "grey", linetype = "dashed")+
  geom_vline(xintercept = -0.26, color = "grey", linetype = "dashed")+
  geom_hline(yintercept = 0.26, color = "grey", linetype = "dashed") +
  geom_hline(yintercept = -0.26, color = "grey", linetype = "dashed") 
ggsave(output_file, width = 6, height = 6, dpi = 300)
#length(data.xx$Group[data.xx$Group=="group5"])
