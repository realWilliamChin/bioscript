getwd()
setwd("C:/Users/Analysis/OneDrive/Work/zhiwusuo/01_Work/2023_05_26_丝瓜/交付/Prep_files/report")
##########################################################
rm(list=ls())
library(ggplot2)
library(ggthemes)
library(gdata)
library("svglite")
getwd()
library("openxlsx")
#library("xlsx")

infiles<-grep('_report.xlsx',dir(),value = T)
infiles
for (i in infiles) {
  out.name<-gsub("_report.xlsx","",i)
  out.name
  dt.bp<-read.xlsx(i,startRow = 9,check.names = F,sheet=1)
  class(dt.bp)
  head(dt.bp)
  colnames(dt.bp)[1]<-"Category"
  colnames(dt.bp)[2]<-"Count"
  colnames(dt.bp)[5]<-"Fold_enrichment"
  colnames(dt.bp)[6]<-"P_value"
  colnames(dt.bp)[9]<-"Q_value"
  dt.bp<-dt.bp[order(dt.bp[,6]),]
  dt.bp<-dt.bp[dt.bp[,2]>=3 & dt.bp[,6]<0.05,]
  dt.bp<-head(dt.bp[order(dt.bp[,6]),],n=15)
  dt.bp<-dt.bp[order(-dt.bp[,6]),]
  dt.bp[,c(1,2,6)]
  if (nrow(dt.bp) > 0) {
    p<-ggplot(dt.bp,aes(y=factor(dt.bp[,1],levels=dt.bp[,1]),x=Fold_enrichment,size=Count,colour=P_value))+
      geom_point()+scale_y_discrete(limits=dt.bp[,1])+
      scale_colour_gradient(low="red",high="green")+
      ggtitle(paste0(out.name,"_GOBP",sep=""))+ylab(colnames(data)[1])+theme_base()
    p
    ggsave(paste0(out.name,"_enrich_GOBP_Bubble.png",sep=""),p,dpi=320,width = 20,height = 10)
  }
  
  dt.cc<-read.xlsx(i,startRow = 9,check.names = F,sheet=2)
  class(dt.cc)
  head(dt.cc)
  colnames(dt.cc)[1]<-"Category"
  colnames(dt.cc)[2]<-"Count"
  colnames(dt.cc)[5]<-"Fold_enrichment"
  colnames(dt.cc)[6]<-"P_value"
  colnames(dt.cc)[9]<-"Q_value"
  dt.cc<-dt.cc[order(dt.cc[,6]),]
  dt.cc<-dt.cc[dt.cc[,2]>=3 & dt.cc[,6]<0.05,]
  dt.cc<-head(dt.cc[order(dt.cc[,6]),],n=15)
  dt.cc<-dt.cc[order(-dt.cc[,6]),]
  dt.cc[,c(1,2,6)]
  if (nrow(dt.cc) > 0) {
    p<-ggplot(dt.cc,aes(y=factor(dt.cc[,1],levels=dt.cc[,1]),x=Fold_enrichment,size=Count,colour=P_value))+
      geom_point()+scale_y_discrete(limits=dt.cc[,1])+
      scale_colour_gradient(low="red",high="green")+
      ggtitle(paste0(out.name,"_GOCC",sep=""))+ylab(colnames(data)[1])+theme_base()
    p
    ggsave(paste0(out.name,"_enrich_GOCC_Bubble.png",sep=""),p,dpi=320,width = 20,height = 10)
  }
  
  
  dt.mf<-read.xlsx(i,startRow = 9,check.names = F,sheet=3)
  class(dt.mf)
  head(dt.mf)
  colnames(dt.mf)[1]<-"Category"
  colnames(dt.mf)[2]<-"Count"
  colnames(dt.mf)[5]<-"Fold_enrichment"
  colnames(dt.mf)[6]<-"P_value"
  colnames(dt.mf)[9]<-"Q_value"
  dt.mf<-dt.mf[order(dt.mf[,6]),]
  dt.mf<-dt.mf[dt.mf[,2]>=3 & dt.mf[,6]<0.05,]
  dt.mf<-head(dt.mf[order(dt.mf[,6]),],n=15)
  dt.mf<-dt.mf[order(-dt.mf[,6]),]
  dt.mf[,c(1,2,6)]
  if (nrow(dt.mf) > 0) {
    p<-ggplot(dt.mf,aes(y=factor(dt.mf[,1],levels=dt.mf[,1]),x=Fold_enrichment,size=Count,colour=P_value))+
      geom_point()+scale_y_discrete(limits=dt.mf[,1])+
      scale_colour_gradient(low="red",high="green")+
      ggtitle(paste0(out.name,"_GOMF",sep=""))+ylab(colnames(data)[1])+theme_base()
    p
    ggsave(paste0(out.name,"_enrich_GOMF_Bubble.png",sep=""),p,dpi=320,width = 20,height = 10)
  }
  
  
  dt.kegg<-read.xlsx(i,startRow = 9,check.names = F,sheet=4)
  class(dt.kegg)
  head(dt.kegg)
  colnames(dt.kegg)[1]<-"Category"
  colnames(dt.kegg)[2]<-"Count"
  colnames(dt.kegg)[5]<-"Fold_enrichment"
  colnames(dt.kegg)[6]<-"P_value"
  colnames(dt.kegg)[9]<-"Q_value"
  dt.kegg<-dt.kegg[order(dt.kegg[,6]),]
  dt.kegg<-dt.kegg[dt.kegg[,2]>=3 & dt.kegg[,6]<0.05,]
  dt.kegg<-head(dt.kegg[order(dt.kegg[,6]),],n=10)
  dt.kegg<-dt.kegg[order(-dt.kegg[,6]),]
  dt.kegg[,c(1,2,6)]
  if (nrow(dt.kegg) > 0) {
    p<-ggplot(dt.kegg,aes(y=factor(dt.kegg[,1],levels=dt.kegg[,1]),x=Fold_enrichment,size=Count,colour=P_value))+
      geom_point()+scale_y_discrete(limits=dt.kegg[,1])+
      scale_colour_gradient(low="red",high="green")+
      ggtitle(paste0(out.name,"_KEGG",sep=""))+ylab(colnames(data)[1])+theme_base()
    p
    ggsave(paste0(out.name,"_enrich_KEGG_Bubble.png",sep=""),p,dpi=320,width = 20,height = 10)
  }
  
  dt.bp<-dt.bp[order(-dt.bp[,2]),]
  dt.bp$Class<-rep("BP",nrow(dt.bp))
  dt.cc<-dt.cc[order(-dt.cc[,2]),]
  dt.cc$Class<-rep("CC",nrow(dt.cc))
  dt.mf<-dt.mf[order(-dt.mf[,2]),]
  dt.mf$Class<-rep("MF",nrow(dt.mf))
  dt.total<-rbind(dt.bp,dt.cc,dt.mf)
  bar.p<-ggplot(dt.total,aes(x=factor(dt.total[,1],levels=dt.total[,1]),y=Count,fill=Class))+
    geom_bar(stat="identity")+
    theme_bw()+
    theme(axis.text.x = element_text(angle=90,
                                     vjust=0.2,hjust=1,size=15),
          legend.title=element_blank(),
          legend.position = c(0.8,0.8),
          panel.grid = element_blank())
  ggsave(paste0(out.name,"_enrich_GObarplot.png",sep=""),bar.p,dpi=320,width = 20,height = 30)
}
#i<-"B8-10_vs_B8-18_deg_Up_report.xlsx"


#xl_list<-list()
#xl_list$BP<-dt.tib
#xl_list$MF<-dt.tib
#write.xlsx(xl_list,"test1.xlsx")
#write.xlsx(dt.tib,"test1.xlsx",sheetName="MF")

