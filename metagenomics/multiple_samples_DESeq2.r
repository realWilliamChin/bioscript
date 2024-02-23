suppressPackageStartupMessages(library(ggthemes))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggcorrplot))
suppressPackageStartupMessages(library(corrplot))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(FactoMineR))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(openxlsx))

rm(list=ls())
getwd()

# Define the command-line options
option_list = list(
  make_option(c("-t", "--tableType"), type="character", default=NULL, 
              help="table type Species ...", metavar="character"),
  make_option(c("-p", "--bsPos"), type="numeric", default=2, 
              help="positive log2FoldChange value", metavar="numeric")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Check that all required arguments are provided
if (is.null(opt$tableType)){
  print_help(opt_parser)
  stop("Please provide the table type", call.=FALSE)
}else if (is.null(opt$bsPos)){
  print_help(opt_parser)
  stop("Please provide the positive log2FoldChange value", call.=FALSE)
}

# args=commandArgs(T)
bs_pos <- opt$bsPos
bs_neg <- -bs_pos
bs_pos

dir.create("Analysis")
dir.create("Analysis/Expression_data")
dir.create("Analysis/Expression_data_graphs")
dir.create("Expression_data_evaluation")


table_type <- opt$tableType
fpkm_file <- paste0(table_type, "_Summary_count.txt")
reads_file <- paste0(table_type, ".txt")

read.table(fpkm_file,sep="\t",header=T,row.names=1,check.names=F)->fpkm

all_fpkm<-fpkm[rowSums(fpkm)>0,]
#all.heatmap<-pheatmap(all_fpkm,scale="row",cluster_cols=F,show_rownames=F)
# 检查数据行数是否超过 65535
if (nrow(all_fpkm) > 65535) {
  # 随机挑选 65535 行的索引
  sample_indices <- sample(1:nrow(all_fpkm), size = 65535, replace = FALSE)
  
  # 按照原始顺序重新排列数据
  subset_data <- all_fpkm[sample_indices, ]
  subset_data <- subset_data[order(sample_indices), ]
  
  # 创建热图
  all.heatmap <- pheatmap(subset_data, scale = "row", cluster_cols = FALSE, show_rownames = FALSE)
} else {
  # 如果数据不超过 65535 行，直接创建热图
  all.heatmap <- pheatmap(all_fpkm, scale = "row", cluster_cols = FALSE, show_rownames = FALSE)
}

all_gene_heatmap_filename <- paste0("Expression_data_evaluation/all_", table_type, "_heatmap.jpeg")
ggsave(all_gene_heatmap_filename,all.heatmap,dpi=300,width=10,height=10)

sample_info<-read.table("samples_described.txt",sep="\t",header=T,check.names=F,stringsAsFactors = F)
#sample_info[sample_info$group=="CK",]
reads_data<-read.table(reads_file,sep="\t",row.names=1,header=T,check.names=F,stringsAsFactors = F)
reads_data<-na.omit(reads_data)
comp_info<-read.table("compare_info.txt",sep="\t",header=T,check.names=F,stringsAsFactors = F)
#head(reads_data)
#head(reads_data[,sample_info[sample_info$group %in% comp_info[1,],]$sample])
#tmp.data<-reads_data[,sample_info[sample_info$group %in% comp_info[1,],]$sample]
#rownames(comp_info)
i<-''
total.deg<-''
stat.deg<-data.frame()
for(i in seq_along(1:nrow(comp_info))){
  group_vs_group_name<-paste(comp_info[i,1],comp_info[i,2],sep="_vs_")
  #print(i)
  data.treat<-reads_data[,sample_info[sample_info$group == comp_info[i,1],]$sample]
  data.control<-reads_data[,sample_info[sample_info$group == comp_info[i,2],]$sample]
  rnaseqMatrix<-cbind(data.treat,data.control)
  fpkm.treat<-fpkm[,sample_info[sample_info$group == comp_info[i,1],]$sample]
  fpkm.control<-fpkm[,sample_info[sample_info$group == comp_info[i,2],]$sample]
  fpkm.deg<-cbind(fpkm.treat,fpkm.control)
  rnaseqMatrix = round(rnaseqMatrix)
  rnaseqMatrix = rnaseqMatrix[rowSums(cpm(rnaseqMatrix) > 1) >= 2,]
  conditions = data.frame(conditions=factor(c(rep(comp_info[i,1], ncol(data.treat)), rep(comp_info[i,2], ncol(data.control)))))
  #conditions
  rownames(conditions) = colnames(rnaseqMatrix)
  ddsFullCountTable <- DESeqDataSetFromMatrix(
    countData = rnaseqMatrix,
    colData = conditions,
    design = ~ conditions)
  dds = DESeq(ddsFullCountTable)
  contrast=c("conditions",comp_info[i,1],comp_info[i,2])
  res = results(dds, contrast)
  baseMeanA <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$conditions == comp_info[i,1]])
  baseMeanB <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$conditions == comp_info[i,2]])
  res = cbind(baseMeanA, baseMeanB, as.data.frame(res))
  res = cbind(sampleA=comp_info[i,1], sampleB=comp_info[i,2], as.data.frame(res))
  res$padj[is.na(res$padj)]  <- 1
  res = as.data.frame(res[order(res$pvalue),])
  outfile=paste0(group_vs_group_name,"_DE_results")
  outfile.ma=paste0(group_vs_group_name,"_DE_results_readCounts.matrix")
  rnaseqMatrix<-cbind(as.data.frame(rownames(rnaseqMatrix)),rnaseqMatrix)
  colnames(rnaseqMatrix)[1]<-table_type
  write.table(rnaseqMatrix, file=outfile.ma, sep='	', quote=FALSE,row.names=F)
  
  volcano<-res
  volcano$padj<-ifelse(volcano$padj<0.000000000000001,0.000000000000001,volcano$padj)
  volcano$regulation = as.factor(ifelse(volcano$padj < 0.05 & abs(volcano$log2FoldChange) >= bs_pos, ifelse(volcano$log2FoldChange >=bs_pos ,'Up','Down'),'NoSignificant'))
  volcano$FC = 2^volcano$log2FoldChange
  fpkm.tmp<-fpkm.deg[rownames(volcano),]
  colnames(fpkm.tmp)<-paste(colnames(fpkm.tmp),'_RA',sep='')
  volcano<-cbind(volcano,fpkm.tmp)
  volcano<-cbind(as.data.frame(rownames(volcano)),volcano)
  colnames(volcano)[1]<-table_type
  # 修改保存为 excel 格式 (2024_02_23)
  # write.table(volcano, file=outfile, sep='	', quote=FALSE,row.names=F)
  write.xlsx(volcano, file = paste0(outfile, ".xlsx"), row.names = FALSE)
  total.deg<-c(total.deg,rownames(volcano)[volcano$regulation=="Up" | volcano$regulation=="Down" ])
  total_deg_num<-nrow(volcano[volcano$regulation=="Up" | volcano$regulation=="Down", ])
  up_deg_num<- nrow(volcano[volcano$regulation=="Up",])
  down_deg_num<-nrow(volcano[volcano$regulation=="Down",])
  deg_stat_group<-paste(comp_info[i,1],comp_info[i,2],sep="_vs_")
  stat.deg<-rbind(stat.deg,rbind(c(deg_stat_group,total_deg_num,up_deg_num,down_deg_num)))
  if(up_deg_num==0 | down_deg_num==0){next}
  p.volcano<-ggplot(data = volcano, aes(x = log2FoldChange, y = -log10(padj), colour=regulation)) +
    geom_point() +
    scale_color_manual(values=c("green", "grey","red")) +
    geom_vline(xintercept=c(bs_neg,bs_pos),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8)+theme_base()
  outfile.volcano=paste0("Analysis/Expression_data_graphs/",group_vs_group_name,"_volcano.jpeg")
  ggsave(outfile.volcano,p.volcano,dpi=300,width=10,height=10)
  tmp.deg<-as.character(rownames(volcano)[volcano$regulation=="Up" | volcano$regulation=="Down" ])
  fpkm_tmp<-na.omit(fpkm.deg[tmp.deg,])
  fpkm_tmp<-log2(fpkm_tmp[rowSums(fpkm_tmp)>0,]+1)
  p.tmpdeg.heatmap<-pheatmap(fpkm_tmp,scale="row",cluster_cols=F,show_rownames=F)
  outfile.tmpdeg.heatmap=paste0("Analysis/Expression_data_graphs/",group_vs_group_name,"_heatmap.jpeg")
  ggsave(outfile.tmpdeg.heatmap,p.tmpdeg.heatmap,dpi=300,width=10,height=10)
  outfile.Up_ID<-paste0("Analysis/", group_vs_group_name,"_Up_ID.txt")
  UP_ID.dt<-data.frame(GeneID=as.character(rownames(volcano)[volcano$regulation=="Up"]))
  write.table(UP_ID.dt, file=outfile.Up_ID, sep='	', quote=FALSE,row.names=F, col.names=FALSE)
  outfile.Down_ID<-paste0("Analysis/",group_vs_group_name,"_Down_ID.txt")
  Down_ID.dt<-data.frame(GeneID=as.character(rownames(volcano)[volcano$regulation=="Down"]))
  write.table(Down_ID.dt, file=outfile.Down_ID, sep='	', quote=FALSE,row.names = F, col.names=FALSE)
}
rm(i)
############draw heatmap #############
total.deg<-as.character(unique(total.deg))
#length(total.deg)

deg_fpkm<-na.omit(fpkm[total.deg,])

deg_fpkm<-log2(deg_fpkm[rowSums(deg_fpkm)>0,]+1)
p.heatmap<-pheatmap(deg_fpkm,scale="row",cluster_cols=F,show_rownames=F)
ggsave("deg_heatmap.jpeg",p.heatmap,dpi=300,width=10,height=10)
############draw boxplot density########

melt(fpkm,variable.name="sample",value.name="fpkm")->data.m
data.m[data.m[,2]>0,]->data.m
p<-ggplot(data.m,aes(x=sample,y=log2(fpkm),fill=sample))+
  geom_boxplot()+theme_base()+ylab("log2(Relative_Abundance)")+xlab("Sample")+
  theme(axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),axis.text.x=element_text(size=rel(1.5),angle=90,hjus=1,vjust=.5))

# boxplot 根据样本数量调整宽度
boxplot_width <- (length(unique(data.m$sample)) - 1) / 2
ggsave(filename="Expression_data_evaluation/RA_boxplot.jpeg",plot=p,height=10,width=boxplot_width,dpi=300)
d<-ggplot(data.m,aes(x=log10(fpkm),col=sample))+geom_density(aes(fill=sample),colour=NA,alpha=.2)+geom_line(stat="density",size=1.5)+xlab("log2(Relative_Abundance)")+theme_base()+theme(axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),axis.text.x=element_text(size=rel(3)),axis.text.y=element_text(size=rel(3)))
ggsave(filename="Expression_data_evaluation/RA_density.jpeg",plot=d,height=10,width=16.8,dpi=300)

#install.packages("ggcorrplot")
fpkm.m<-as.matrix(fpkm)
fpkm.cor<-cor(fpkm.m)
#p.mat <- cor_pmat(fpkm.m)
#p.mat
#cor.p<-ggcorrplot(fpkm.cor,lab=TRUE)
#cor.p<-corrplot(res, method = "number", col = RColorBrewer::brewer.pal(n=11, name = "RdYlGn"),title = "")
#?corrplot
#?ggcorrplot
write.table(fpkm.cor,"sample_cor.txt",sep="\t",quote=F)
#ggsave("correlation.jpeg",cor.p,dpi=300)
min(fpkm.cor)
#fpkm.m<-as.data.frame(fpkm)
#fpkm.cor<-cor(fpkm.m)
resfactor = 5
png("Expression_data_evaluation/correlation.png",res=72*resfactor,height=1800*resfactor,width=1800*resfactor)
corrplot(fpkm.cor,is.corr = F,col= rev(COL2('PiYG')),method="color",addCoef.col = 'black',tl.col="black",col.lim=c(min(fpkm.cor)-0.01,max(fpkm.cor)),cl.ratio=0.1)
dev.off()

pdf("correlation.pdf",width = 15,height = 15)
corrplot(fpkm.cor,is.corr = F,col= rev(COL2('PiYG')),method="color",addCoef.col = 'black',tl.col="black",col.lim=c(min(fpkm.cor)-0.01,max(fpkm.cor)),cl.ratio=0.1)
dev.off()


#head(deg_fpkm)
#stat.deg
#as.data.frame(matrix(stat.deg,nrow=nrow(comp_info),ncol=4,byrow=T),colnames=c("Groups","Total_DEGs","UP"))
#stat.deg
colnames(stat.deg)=c("Groups","Total DEGs","Up regulated","Down regulated")
write.table(stat.deg,file="Analysis/DEG_summary.txt",sep="\t",quote=F,row.names=F)


##############################################################################################################
#setwd("d:/pll/R_work/anova")
#getwd()
data<-read.table(fpkm_file,sep="\t",row.names=1,header=T,check.names=F)
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
p.adjust(data$p_value,method = "BH")->p.aj
data$BH_p_value<-p.aj
df.rname<-data.frame(GeneID=rownames(data))
data<-cbind(df.rname,data)
write.table(data,"anova_analysis_p.txt",sep="\t",quote=F,row.names = F)
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
gene <- read.delim(fpkm_file, row.names = 1, sep = '\t',check.names = FALSE,header=T)
gene<-log2(gene+1)
gene <- t(gene)

# PCA
gene.pca <- PCA(gene, ncp = 2, scale.unit = TRUE, graph = FALSE)

pca_sample <- data.frame(gene.pca$ind$coord[ ,1:2])
head(pca_sample)

#??ȡ PCA ǰ��???Ĺ??׶?
pca_eig1 <- round(gene.pca$eig[1,2], 2)
pca_eig2 <- round(gene.pca$eig[2,2],2 )

group <- read.delim('samples_described.txt', row.names = 2, sep = '\t', check.names = FALSE,header=T)
group <- group[rownames(pca_sample), ]

#ggplot2
p <- ggplot(data = pca_sample, aes(x = Dim.1, y = Dim.2)) +
  geom_point(aes(color = group), size = 5) +
  #scale_color_manual(values = c('orange', 'purple','blue','black')) + 
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent')) +  
  labs(x =  paste('PCA1:', pca_eig1, '%'), y = paste('PCA2:', pca_eig2, '%'), color = '')

#pca_sample
#group
p<-p+geom_text_repel(data=pca_sample,aes(Dim.1, Dim.2, label=rownames(pca_sample)))


cluster_border <- ddply(pca_sample, 'group', function(df) df[chull(df[[1]], df[[2]]), ])
p<-p + geom_polygon(data = cluster_border, aes(color = group),fill=NA, show.legend = FALSE)
ggsave("Expression_data_evaluation/PCA.jpeg",p,dpi=300,width=10,height=10)

#????ͼ5
#????????Ӱ
#p<-p + geom_polygon(data = cluster_border, aes(fill = group), alpha = 0.2, show.legend = FALSE) #+
# #scale_fill_manual(values = c('orange', 'purple','blue','black'))


