setwd("/home/colddata/qinqiang/Project/2024_06_07_renminyiyuan_lixuewei_ji/联合分析")
library(ggplot2)
library(ggthemes)
data.xx<-read.table("result_group.txt",header = T,check.names = F)
head(data.xx)
colnames(data.xx)
ggplot(data.xx,aes(x=log2FC_RNA,y=log2FC_protein,colour=Group))+geom_point()+
  theme_classic()+
  geom_vline(xintercept = 0.26, color = "grey", linetype = "dashed")+
  geom_vline(xintercept = -0.26, color = "grey", linetype = "dashed")+
  geom_hline(yintercept = 0.26, color = "grey", linetype = "dashed") +
  geom_hline(yintercept = -0.26, color = "grey", linetype = "dashed") 

#length(data.xx$Group[data.xx$Group=="group5"])
