library(ggplot2)
library(ggthemes)
getwd()
#setwd("d:/pll/R_work")
#library(ggpubr)
library(agricolae)
library(vegan)
data<-read.csv("Alpha_stat.csv",header=T,row.names=1)
group=data["group"]
rownames(group)
data<-data[,-ncol(data)]
data<-data[match(rownames(group),rownames(data)),]
p.k<-NULL
for (i in seq_along(colnames(data))){cbind(data[,i],group)->data.new
	colnames(data.new)[1]<-colnames(data)[i]
	out<-colnames(data)[i]
	p.kruskal<-kruskal.test(data.new[,1]~group,data=data.new)
	p.k<-c(p.k,round(p.kruskal$p.value,digits=2))
	p<-ggplot(data.new,aes(x=group,y=data.new[,1],fill=group))+geom_boxplot()+ylab(out)+
	theme_base()+ggtitle(paste("P=",round(p.kruskal$p.value,digits=2),sep=""))+
	theme(plot.title=element_text(hjust=0.5,vjust=-6))
	outfile<-paste(out,".jpeg",sep='')
	ggsave(outfile,p,dpi=300)
}
#data.new<-t(data)
#head(data.new)
#head(data)
data[nrow(data)+1,]<-p.k
rownames(data)[nrow(data)]<-"P_value"
#p.tmp<-data.frame(Alpha=colnames(data),P=p.k)
#write.table(p.tmp,"alpha_index_p.txt",sep="\t",quote=F,row.names = F)
write.table(data,"alpha_index_p.txt",sep="\t",quote=F,row.names = T)
