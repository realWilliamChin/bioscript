getwd()
setwd("d:/pll/R_work/metabolism/") #change work directory
rm(list=ls())
library("openxlsx")
reads_data<-read.table("reads_matrix_filtered.txt",sep="\t",row.names=1,header=T,check.names = F) # change input file name
sample_info<-read.table("samples_described.txt",sep="\t",header=T,check.names=F,stringsAsFactors = F)
comp_info<-read.table("compare_info.txt",sep="\t",header=T,check.names=F,stringsAsFactors = F)
i<-NULL
total.df<-data.frame(Count=1:nrow(reads_data))
#total_title<-NULL
for(i in seq_along(1:nrow(comp_info))){
  data.treat<-reads_data[,sample_info[sample_info$group == comp_info[i,1],]$sample]
  data.control<-reads_data[,sample_info[sample_info$group == comp_info[i,2],]$sample]
  n<-NULL
  out.df<-data.frame()
  for (n in (1:nrow(reads_data))) {
    tmp.p<-t.test(x=as.numeric(data.treat[n,]),y=as.numeric(data.control[n,]),paired=T)$p.value
    tmp.fc<-mean(as.numeric(data.treat[n,]))/mean(as.numeric(data.control[n,]))
    #out.df<-rbind(out.df,cbind(data.treat[n,],data.control[n,],data.frame(p_value=tmp.p,FC=tmp.fc)))
    out.df<-rbind(out.df,data.frame(meanA=mean(as.numeric(data.treat[n,])),meanB=mean(as.numeric(data.control[n,])),p_value=tmp.p,FC=tmp.fc))
  }
  total_p<-paste0(comp_info[i,1],"_vs_",comp_info[i,2],"_p")
  total_fc<-paste0(comp_info[i,1],"_vs_",comp_info[i,2],"_fc")
  #total_title<-c(total_title,total_p,total_fc)
  out.df<-cbind(data.frame(ID=rownames(reads_data),groupA=rep(comp_info[i,1],nrow(reads_data)),groupB=rep(comp_info[i,2],nrow(reads_data))),out.df)
  out.file=paste0(comp_info[i,1],"_vs_",comp_info[i,2],".xlsx",sep="")
  total.df<-cbind(total.df,out.df[,(ncol(out.df)-1):ncol(out.df)])
  colnames(total.df)[(ncol(total.df)-1):ncol(total.df)] = c(total_p,total_fc)
  #write.table(out.df,file=out.file,quote=F,sep="\t",row.names = F)
  write.xlsx(out.df,out.file)
}
#head(total.df)
total.df<-total.df[,-1]
head(total.df)
total.df<-cbind(data.frame(KEGG=rownames(reads_data)),total.df)
#dim(reads_data)
write.xlsx(total.df,"all_result.xlsx")
write.table(total.df,file="all_result.txt",quote=F,sep="\t",row.names = F)

