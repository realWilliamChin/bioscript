setwd("C:/Users/Analysis/OneDrive/Work/zhiwusuo/01_Work/2023_03_30_马铃薯/wego")
dir()
library("openxlsx")
# go<-read.xlsx("Book4.xlsx")
file <- "target_gene_for_potato_snp_ZK_specific_min4_Go.txt.Level3.count.xls"
output_file <- sub("\\.xls$", ".jpeg", file)
go<-read.table(file, sep = "\t", header = T)
colnames(go)[2]<-"GoTerms"
colnames(go)
library(ggplot2)
library(ggthemes)
go <- go[go[,4] > 20,]
factor(go$GoTerms,levels=go$GoTerms)->go$GoTerms
bar.p<-ggplot(go,aes(x=GoTerms,y=Counts,fill=Class))+
  geom_bar(stat="identity")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90,
                                   vjust=0.2,hjust=1,size=15),
        legend.title=element_blank(),
        legend.position = c(0.8,0.8),
        panel.grid = element_blank())
ggsave(output_file,bar.p,dpi=320,height = 12,width =15)
