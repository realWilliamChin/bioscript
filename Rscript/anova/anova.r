##############################################################################################################
setwd("d:/pll/R_work/RNA-seq/Niubang/anova/")
############Niubang germinated
library(readxl)
data<-read_xlsx("gene_id_seed_vs_germinated_fpkm.xlsx",1)
data<-as.data.frame(data)
head(data)
rownames(data)<-data$GeneID
data<-data[,-1]
data<-subset(data,rowSums(data)>0)
#data<-read.table("ABA_FPKM.txt",sep="\t",row.names=1,header=T,check.names=F)
group<-read.table("gene_id_seed_vs_germinated_group.txt",sep="\t",header=T,check.names = F)
#head(group)
group
colnames(data)
#match(colnames(data),group$sample)
#order(match(colnames(data),group$sample))
#group[match(colnames(data),group$sample),]
tmp.group<-group$group[match(colnames(data),group$sample)]
tmp.group
#class(tmp.group)
#head(data)
p<-NULL
for (i in rownames(data)) {
  #print(i);
  #print(data[i,])
  tmp.count<-as.numeric(data[i,])
  #tmp.dt<-data.frame(tmp.count,group=rep(c("S1","S3","S2"),each=3))
  tmp.dt<-data.frame(tmp.count,group=tmp.group)
  colnames(tmp.dt)[1]<-i
  tmp.dt$group<-as.factor(tmp.dt$group)
  mod <- aov( tmp.dt[,1]~group , data = tmp.dt)
  p<-c(p,summary(mod)[[1]][5]$Pr[1])
}
length(p)
#head(p)
data$p_value<-p
data<-data[order(data$p_value),]
#head(data)
#?p.adjust
p.adjust(data$p_value,method = "BH")->p.aj
#head(p.aj)
data$BH_p_value<-p.aj
#head(data)
df.rname<-data.frame(GeneID=rownames(data))
data<-cbind(df.rname,data)
head(data)
write.table(data,"anova_analysis_seed_vs_germinated_p.txt",sep="\t",quote=F,row.names = F)
################################从excel文件里读取后进行操作(ABA)
library(readxl)
data<-read_xlsx("ABA_FPKM.xlsx",1)
data<-as.data.frame(data)
head(data)
colnames(data)
rownames(data)<-data[,1]
data<-data[,-1]
head(data)
dim(data)
group<-read.table("ABA_group.txt",sep="\t",header=T,check.names = F)
head(group)
group
colnames(data)
#match(colnames(data),group$sample)
#order(match(colnames(data),group$sample))
#group[match(colnames(data),group$sample),]
tmp.group<-group$group[match(colnames(data),group$sample)]
tmp.group
#class(tmp.group)
#head(data)
p<-NULL
for (i in rownames(data)) {
  #print(i);
  #print(data[i,])
  tmp.count<-as.numeric(data[i,])
  #tmp.dt<-data.frame(tmp.count,group=rep(c("S1","S3","S2"),each=3))
  tmp.dt<-data.frame(tmp.count,group=tmp.group)
  colnames(tmp.dt)[1]<-i
  tmp.dt$group<-as.factor(tmp.dt$group)
  mod <- aov( tmp.dt[,1]~group , data = tmp.dt)
  p<-c(p,summary(mod)[[1]][5]$Pr[1])
}
length(p)
#head(p)
data$p_value<-p
data<-data[order(data$p_value),]
#head(data)
#?p.adjust
p.adjust(data$p_value,method = "BH")->p.aj
#head(p.aj)
data$BH_p_value<-p.aj
#head(data)
df.rname<-data.frame(GeneID=rownames(data))
data<-cbind(df.rname,data)
head(data)
write.table(data,"anova_analysis_ABA_p.txt",sep="\t",quote=F,row.names = F)
############################(PUT)############################################
data<-read_xlsx("Put_FPKM.xlsx",1)
data<-as.data.frame(data)
head(data)
colnames(data)
rownames(data)<-data[,1]
data<-data[,-1]
head(data)
dim(data)
group<-read.table("Put_group.txt",sep="\t",header=T,check.names = F)
head(group)
group
colnames(data)
#match(colnames(data),group$sample)
#order(match(colnames(data),group$sample))
#group[match(colnames(data),group$sample),]
tmp.group<-group$group[match(colnames(data),group$sample)]
tmp.group
#class(tmp.group)
#head(data)
p<-NULL
for (i in rownames(data)) {
  #print(i);
  #print(data[i,])
  tmp.count<-as.numeric(data[i,])
  #tmp.dt<-data.frame(tmp.count,group=rep(c("S1","S3","S2"),each=3))
  tmp.dt<-data.frame(tmp.count,group=tmp.group)
  colnames(tmp.dt)[1]<-i
  tmp.dt$group<-as.factor(tmp.dt$group)
  mod <- aov( tmp.dt[,1]~group , data = tmp.dt)
  p<-c(p,summary(mod)[[1]][5]$Pr[1])
}
length(p)
#head(p)
data$p_value<-p
data<-data[order(data$p_value),]
#head(data)
#?p.adjust
p.adjust(data$p_value,method = "BH")->p.aj
#head(p.aj)
data$BH_p_value<-p.aj
#head(data)
df.rname<-data.frame(GeneID=rownames(data))
data<-cbind(df.rname,data)
head(data)
write.table(data,"anova_analysis_Put_p.txt",sep="\t",quote=F,row.names = F)
############################(CK)############################################
data<-read_xlsx("CK_FPKM.xlsx",1)
data<-as.data.frame(data)
head(data)
colnames(data)
rownames(data)<-data[,1]
data<-data[,-1]
head(data)
dim(data)
group<-read.table("CK_group.txt",sep="\t",header=T,check.names = F)
head(group)
group
colnames(data)
#match(colnames(data),group$sample)
#order(match(colnames(data),group$sample))
#group[match(colnames(data),group$sample),]
tmp.group<-group$group[match(colnames(data),group$sample)]
tmp.group
#class(tmp.group)
#head(data)
p<-NULL
for (i in rownames(data)) {
  #print(i);
  #print(data[i,])
  tmp.count<-as.numeric(data[i,])
  #tmp.dt<-data.frame(tmp.count,group=rep(c("S1","S3","S2"),each=3))
  tmp.dt<-data.frame(tmp.count,group=tmp.group)
  colnames(tmp.dt)[1]<-i
  tmp.dt$group<-as.factor(tmp.dt$group)
  mod <- aov( tmp.dt[,1]~group , data = tmp.dt)
  p<-c(p,summary(mod)[[1]][5]$Pr[1])
}
length(p)
#head(p)
data$p_value<-p
data<-data[order(data$p_value),]
#head(data)
#?p.adjust
p.adjust(data$p_value,method = "BH")->p.aj
#head(p.aj)
data$BH_p_value<-p.aj
#head(data)
df.rname<-data.frame(GeneID=rownames(data))
data<-cbind(df.rname,data)
head(data)
write.table(data,"anova_analysis_CK_p.txt",sep="\t",quote=F,row.names = F)

########################PCA
gene <- read.delim('T_ALL_all_RNAseq_data_RPKM.txt', row.names = 1, sep = '\t',check.names = FALSE,header=T)
gene<-log2(gene+1)
#??????????ֵ????????ת?ã?ʹ??Ϊ????????Ϊ????
gene <- t(gene)

#????ʹ?? FactoMineR ???еķ?????ʵ?? PCA ?????;???????
library(FactoMineR)

#?????л???????ֵ?? PCA ????
gene.pca <- PCA(gene, ncp = 2, scale.unit = TRUE, graph = FALSE)
#plot(gene.pca)  #PCA ??ͼ

#??ȡ?????? PCA ǰ��???е?????
pca_sample <- data.frame(gene.pca$ind$coord[ ,1:2])
head(pca_sample)

#??ȡ PCA ǰ��???Ĺ??׶?
pca_eig1 <- round(gene.pca$eig[1,2], 2)
pca_eig2 <- round(gene.pca$eig[2,2],2 )

#??ȡ???ϲ???????????Ϣ
group <- read.delim('baixuebing_group.txt', row.names = 2, sep = '\t', check.names = FALSE,header=T)
group <- group[rownames(pca_sample), ]

#ggplot2 ???ƶ?άɢ??ͼ
library(ggplot2)
p <- ggplot(data = pca_sample, aes(x = Dim.1, y = Dim.2)) +
  geom_point(aes(color = group), size = 5) +  #???????????????????????????????????????
  #scale_color_manual(values = c('orange', 'purple','blue','black')) +  #???????????????
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent')) +  #????????????????????????
  labs(x =  paste('PCA1:', pca_eig1, '%'), y = paste('PCA2:', pca_eig2, '%'), color = '')  #??? PCA ???????????????????????????????????????
p
#p
#????ͼ1
#??ʶ???????ƣ?ʹ?? ggplot2 ????չ?? ggrepel ��????
#?????????Ʊ?ǩ
library(ggrepel)
#pca_sample
#group
p<-p+geom_text_repel(data=pca_sample,aes(Dim.1, Dim.2, label=rownames(pca_sample)))

#p

#????ͼ2
#???? 95% ??????Բ???????ڱ?ʾ???????࣬??ֻ???????ڸ????????????? 5 ????????
#????????????
#p + stat_ellipse(aes(color = group), level = 0.95, show.legend = FALSE)

#????ͼ3
#??????????Ӱ
#p + stat_ellipse(aes(fill = group), geom = 'polygon', level = 0.95, alpha = 0.1, show.legend = FALSE) +
# scale_fill_manual(values = c('orange', 'purple','blue','black'))


#????ͼ4
#??????��??ͬ?????????߽?????ʽ???????ڸ????????????? 3 ????????
#??????????
library(plyr)
cluster_border <- ddply(pca_sample, 'group', function(df) df[chull(df[[1]], df[[2]]), ])
p<-p + geom_polygon(data = cluster_border, aes(color = group),fill=NA, show.legend = FALSE)
ggsave("baixuebing_PCA.jpeg",p,dpi=300,width=10,height=10)

#????ͼ5
#????????Ӱ
#p<-p + geom_polygon(data = cluster_border, aes(fill = group), alpha = 0.2, show.legend = FALSE) #+
# #scale_fill_manual(values = c('orange', 'purple','blue','black'))





#getwd()
data<-read.table("fpkm_matrix_filtered.txt",sep="\t",row.names=1,header=T,check.names=F)
group<-read.table("samples_described.txt",sep="\t",header=T,check.names = F)
head(group)
colnames(data)
#match(colnames(data),group$sample)
#order(match(colnames(data),group$sample))
#group[match(colnames(data),group$sample),]
tmp.group<-group$group[match(colnames(data),group$sample)]
#tmp.group
#class(tmp.group)
#head(data)
p<-NULL
for (i in rownames(data)) {
  #print(i);
  #print(data[i,])
  tmp.count<-as.numeric(data[i,])
  #tmp.dt<-data.frame(tmp.count,group=rep(c("S1","S3","S2"),each=3))
  tmp.dt<-data.frame(tmp.count,group=tmp.group)
  colnames(tmp.dt)[1]<-i
  tmp.dt$group<-as.factor(tmp.dt$group)
  mod <- aov( tmp.dt[,1]~group , data = tmp.dt)
  p<-c(p,summary(mod)[[1]][5]$Pr[1])
}
length(p)
head(p)
data$p_value<-p
data<-data[order(data$p_value),]
write.table(data,"anova_analysis_p.txt",sep="\t",quote=F)
#dim(data)
#rep(c("S1","S2","S3"),each=3)
#as.numeric(data[1,])
#tmp.count<-as.numeric(data[rownames(data)[1],])
#tmp.dt<-data.frame(tmp.count,group=rep(c("S1","S2","S3"),each=3))
#colnames(tmp.dt)[1]<-rownames(data)[1]           
#tmp.dt
#tmp.dt$group<-as.factor(tmp.dt$group)
#mod <- aov( ENSMUSG00000029661~group , data = tmp.dt)
#haha<-summary(mod)
#summary(mod)[[1]][5]$Pr[1]
#####################################################################
#setwd("D:/pll/R_work/baisanye")
#??ȡ????????ֵ????
#?Ƽ?ʹ?? log ת?????Ļ???????ֵ?????Ͳ?ͬ????????ˮƽ??��????????????????
gene <- read.delim('fpkm_matrix_filtered.txt', row.names = 1, sep = '\t',check.names = FALSE,header=T)
gene<-log2(gene+1)
#??????????ֵ????????ת?ã?ʹ??Ϊ????????Ϊ????
gene <- t(gene)

#????ʹ?? FactoMineR ???еķ?????ʵ?? PCA ?????;???????
library(FactoMineR)

#?????л???????ֵ?? PCA ????
gene.pca <- PCA(gene, ncp = 2, scale.unit = TRUE, graph = FALSE)
#plot(gene.pca)  #PCA ??ͼ

#??ȡ?????? PCA ǰ��???е?????
pca_sample <- data.frame(gene.pca$ind$coord[ ,1:2])
head(pca_sample)

#??ȡ PCA ǰ��???Ĺ??׶?
pca_eig1 <- round(gene.pca$eig[1,2], 2)
pca_eig2 <- round(gene.pca$eig[2,2],2 )

#??ȡ???ϲ???????????Ϣ
group <- read.delim('samples_described.txt', row.names = 2, sep = '\t', check.names = FALSE,header=T)
group <- group[rownames(pca_sample), ]

#ggplot2 ???ƶ?άɢ??ͼ
library(ggplot2)
p <- ggplot(data = pca_sample, aes(x = Dim.1, y = Dim.2)) +
  geom_point(aes(color = group), size = 5) +  #???????????????????????????????????????
  #scale_color_manual(values = c('orange', 'purple','blue','black')) +  #???????????????
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent')) +  #????????????????????????
  labs(x =  paste('PCA1:', pca_eig1, '%'), y = paste('PCA2:', pca_eig2, '%'), color = '')  #??? PCA ???????????????????????????????????????
#p
#p
#????ͼ1
#??ʶ???????ƣ?ʹ?? ggplot2 ????չ?? ggrepel ��????
#?????????Ʊ?ǩ
library(ggrepel)
#pca_sample
#group
p<-p+geom_text_repel(data=pca_sample,aes(Dim.1, Dim.2, label=rownames(pca_sample)))

#p

#????ͼ2
#???? 95% ??????Բ???????ڱ?ʾ???????࣬??ֻ???????ڸ????????????? 5 ????????
#????????????
#p + stat_ellipse(aes(color = group), level = 0.95, show.legend = FALSE)

#????ͼ3
#??????????Ӱ
#p + stat_ellipse(aes(fill = group), geom = 'polygon', level = 0.95, alpha = 0.1, show.legend = FALSE) +
 # scale_fill_manual(values = c('orange', 'purple','blue','black'))


#????ͼ4
#??????��??ͬ?????????߽?????ʽ???????ڸ????????????? 3 ????????
#??????????
library(plyr)
cluster_border <- ddply(pca_sample, 'group', function(df) df[chull(df[[1]], df[[2]]), ])
p<-p + geom_polygon(data = cluster_border, aes(color = group),fill=NA, show.legend = FALSE)
ggsave("PCA.jpeg",p,dpi=300,width=10,height=10)

#????ͼ5
#????????Ӱ
#p<-p + geom_polygon(data = cluster_border, aes(fill = group), alpha = 0.2, show.legend = FALSE) #+
 # #scale_fill_manual(values = c('orange', 'purple','blue','black'))

