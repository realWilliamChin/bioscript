#setwd("/home/colddata/qinqiang/Project/2024_01_24_songqian_hongjiyinzu/2024_03_11/group1/R")
library(openxlsx)
library(vegan)
rm(list=ls())

Alpha_analysis_dir = '01_Alpha_analysis/'
Beta_analysis_dir = '02_Beta_analysis/'
Classification_dir = '03_Classification_Distribution/'
dir.create(Alpha_analysis_dir)
dir.create(Beta_analysis_dir)
dir.create(Classification_dir)

convert_to_numeric <- function(df) {
  # 遍历DataFrame的每一列
  for (col in names(df)) {
    # 尝试将列转换为numeric
    converted_col <- as.numeric(as.character(df[[col]]))
    # 检查是否存在非数值型的值
    if (any(is.na(converted_col))) {
      # 输出警告信息
      warning(paste("Column", col, "contains non-numeric values and will not be converted."))
    } else {
      # 更新原始列为数值型
      df[[col]] <- converted_col
    }
  }
  # 返回转换后的DataFrame
  return(df)
}

#data <- read.table("Species.txt",header=T, sep="\t",check.names = F,row.names=1,fileEncoding = "UTF-8")
data<-read.xlsx("Species.xlsx")
data <- convert_to_numeric(data)
head(data)
dim(data)
rownames(data)
data_transposed <- t(data[, -1])
head(data_transposed)
dim(data_transposed)
# Shannon指数
shannon <- diversity(data_transposed, index="shannon")
class(shannon)
shannon
# Simpson指数
simpson <- diversity(data_transposed, index="simpson")
simpson
# Chao1指数
chao1 <- estimateR(data_transposed, method="Chao")[2,]
chao1
# ACE指数
# ace <- estimateR(data_transposed, method="ACE")[4,]
# ace

# 创建alpha指数数据框
diversity_indices <- data.frame(Shannon=shannon, Simpson=simpson, Chao1=chao1)
rownames(diversity_indices)
#write.xlsx(diversity_indices,"alpha_index.xlsx",rowNames=T)
diversity_indices
group <- read.delim('samples_described.txt', row.names = 2, sep = '\t', check.names = FALSE,header=T)
head(group)
#diversity_indices$group<-group$group
#diversity_indices
#data.tmp<-as.data.frame(cbind(diversity_indices[,1],group$group))
#class(data.tmp)
#data.tmp
library(ggplot2)
library(ggthemes)
#########################alpha指数总体变化比较##########################
p.k<-NULL
for (i in seq_along(colnames(diversity_indices))){as.data.frame(cbind(diversity_indices[,i],group$group))->data.new
  colnames(data.new)[1]<-colnames(diversity_indices)[i]
  colnames(data.new)[2]<-"group"
  data.new$group<-as.factor(data.new$group)
  data.new[,1]<-as.numeric(data.new[,1])
  head(data.new)
  out<-colnames(diversity_indices)[i]
  p.kruskal<-kruskal.test(data.new[,1]~group,data=data.new)
  p.k<-c(p.k,round(p.kruskal$p.value,digits=2))
  p.box<-ggplot(data.new,aes(x=group,y=data.new[,1],fill=group))+geom_boxplot()+ylab(out)+
    theme_base()+ggtitle(paste("P=",round(p.kruskal$p.value,digits=2),sep=""))+
    theme(plot.title=element_text(hjust=0.5,vjust=-6))
  outfile<-paste(out,".jpeg",sep='')
  ggsave(outfile,p.box,dpi=300)
}
class(diversity_indices$Shannon)
p.k
diversity_indices[nrow(diversity_indices)+1,]<-p.k
rownames(diversity_indices)[nrow(diversity_indices)]<-"P_value"
#write.table(diversity_indices,"alpha_index_p.txt",sep="\t",quote=F,row.names = T)
#?write.xlsx
#write.xlsx(diversity_indices,"alpha_index_p.xlsx",rowNames=T)
#################################alpha指数两组之间比较#################
diversity_paired<-diversity_indices[1:(nrow(diversity_indices)-1),]
diversity_paired
diversity_paired$group<-group$group[match(rownames(diversity_paired),rownames(group))]
comp_info<-read.table("compare_info.txt",sep="\t",header=T,check.names=F,stringsAsFactors = F)
comp_info
wilcox.test(Shannon~group,subset(diversity_paired,group %in% as.character(comp_info[1,]))[c(1,4)])
#ggplot(subset(diversity_paired,group %in% as.character(comp_info[1,]))[c(1,4)],aes(x=group,y=Shannon,fill=group))+geom_boxplot()
p.paired<-NULL
col_names<-NULL
as.character(comp_info[1,])[1]
for (i in 1:nrow(comp_info)) {
  alpha_paired<-subset(diversity_paired,group %in% as.character(comp_info[i,]))
  col_names<-c(col_names,paste0(as.character(comp_info[i,])[1],"_vs_",as.character(comp_info[i,])[2]))
  for (t in 1:(ncol(alpha_paired)-1)) {
    alpha_box<-alpha_paired[c(t,ncol(alpha_paired))]
    p.wilcox<-wilcox.test(alpha_box[,1]~group,alpha_box)
    paired_boxplot<-ggplot(alpha_box,aes(x=group,y=alpha_box[,1],fill=group))+geom_boxplot()+ylab(colnames(alpha_box)[1])+
      theme_base()+ggtitle(paste("P=",round(p.wilcox$p.value,digits=2),sep=""))+
      theme(plot.title=element_text(hjust=0.5,vjust=-6))
    paired_outfile<-paste0(Alpha_analysis_dir, comp_info[i,1],"_vs_",comp_info[i,2],"_",colnames(alpha_box)[1],".jpeg")
    ggsave(paired_outfile,paired_boxplot,dpi=300)
    p.paired<-c(p.paired,round(p.wilcox$p.value,digits=2))
    
  }
}
p.paired
col_names
p_paired_df<-as.data.frame(matrix(p.paired,nrow=ncol(diversity_paired)-1))
p_paired_df
colnames(p_paired_df)<-col_names
rownames(p_paired_df)<-colnames(diversity_paired)[1:(ncol(diversity_paired)-1)]
#write.xlsx(p_paired_df,"alpha_paired_p.xlsx",rowNames=T)
#######将alpha指数及组与组之间的alpha指数比较结果输出到同一个excel文件
alpha_list<-list("Sheet1"=diversity_indices,"alpha_paired_p"=p_paired_df)
write.xlsx(alpha_list,paste0(Alpha_analysis_dir,"alpha_index_p.xlsx"),rowNames=T)
#########################PCA 分析############################
head(data)
dim(data)
rownames(data)
pca_data<-t(data[,-1])
rownames(pca_data)
head(pca_data)
dim(pca_data)
rownames(pca_data)
library(FactoMineR)
data.pca <- PCA(pca_data, ncp = 2, scale.unit = TRUE, graph = FALSE)
summary(data.pca$var)
head(data.pca$var$cor)
head(data.pca$var$contrib)
pca_sample <- data.frame(data.pca$ind$coord[ ,1:2])
head(pca_sample)
pca_eig1 <- round(data.pca$eig[1,2], 2)
pca_eig2 <- round(data.pca$eig[2,2],2 )
group <- read.delim('samples_described.txt', row.names = 2, sep = '\t', check.names = FALSE,header=T)
#group <- group[rownames(pca_sample), ]
pca_sample$group<-group[rownames(pca_sample), ]
pca_sample
library(ggplot2)
p <- ggplot(data = pca_sample, aes(x = Dim.1, y = Dim.2)) +
  geom_point(aes(color = group), size = 5) +  
  #scale_color_manual(values = c('orange', 'purple','blue','black')) +  #???????????????
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent')) +  
  labs(x =  paste('PCA1:', pca_eig1, '%'), y = paste('PCA2:', pca_eig2, '%'), color = '')  
p
library(ggrepel)
#pca_sample
#group
p<-p+geom_text_repel(data=pca_sample,aes(Dim.1, Dim.2, label=rownames(pca_sample)))
p
ggsave(paste0(Beta_analysis_dir, "PCA.jpeg"),p,dpi=300,width=10,height=10)
#############################beta 多样性分析###################################
library(FactoMineR)
library(ggrepel)
library(plyr)
library(ggthemes)
library(vegan)
distance <- vegdist(pca_data, method = 'bray')
#head(distance)
dist.df<-as.data.frame(as.matrix(distance))
head(dist.df)
write.xlsx(dist.df,file = paste0(Beta_analysis_dir, "bray_curtis_dist.xlsx"),rowNames=T)
head(dist.df)
pcoa <- cmdscale(distance, k = (nrow(pca_data) - 1), eig = TRUE)
plot_data <- data.frame({pcoa$point})[1:2]
head(plot_data)
dim(plot_data)
rownames(plot_data)
plot_data$Sample_name <- rownames(plot_data)
names(plot_data)[1:2] <- c('PCoA1', 'PCoA2')
head(plot_data)
plot_data
rownames(plot_data)
#tmp.group<-as.data.frame(group[rownames(plot_data),],row.names=rownames(plot_data))
plot_data$group<-group[rownames(plot_data),]
plot_data
eig = pcoa$eig
eig
library(tidyverse)
p.PCoA <- ggplot(data = plot_data, aes(x=PCoA1, y=PCoA2)) +
  geom_point(alpha=.7, size=2,aes(color=group)) +
  geom_text_repel(data=plot_data,aes(PCoA1, PCoA2,label=Sample_name)) +
  theme_base()+
  #stat_chull() +
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""))
p.PCoA
ggsave(paste0(Beta_analysis_dir,"Species_PCoA.jpeg"),p.PCoA,dpi=320,height = 10,width=10)
#######################################
###############################NMDS 分析###############################
df_nmds <- metaMDS(distance, k = 2)
summary(df_nmds)
#应力函数值（<=0.2合理）
df_nmds_stress <- df_nmds$stress
df_nmds_stress
#检查观测值非相似性与排序距离之间的关系——没有点分布在线段较远位置表示该数据可以使用NMDS分析
stressplot(df_nmds)
#提取作图数据
df_points <- as.data.frame(df_nmds$points)
class(df_points)
df_points
#添加samp1es变量
df_points$samples <- rownames(df_points)
#修改列名
names(df_points)[1:2] <- c('NMDS1', 'NMDS2')
df_points$group<-group[rownames(df_points),]
df_points

p.NMDS<-ggplot(data=df_points,aes(x=NMDS1,y=NMDS2))+#指定数据、X轴、Y轴，颜色
  theme_bw()+#主题设置
  geom_point(aes(color = group), shape = 19, size=3)+#绘制点图并设定大小
  theme(panel.grid = element_blank())+
  geom_vline(xintercept = 0,lty="dashed", size = 1, color = 'grey50')+
  geom_hline(yintercept = 0,lty="dashed", size = 1, color = 'grey50')+#图中虚线
  geom_text(aes(label=samples, y=NMDS2+0.01,x=NMDS1+0.01,
                vjust=0, color = group),size=3.5, show.legend = F)+#添加数据点的标签
  stat_ellipse(data=df_points,
               geom = "polygon",level=0.95,
               linetype = 2,size=0.5,
               aes(fill=group),
               alpha=0.2)+
  #scale_color_manual(values = color) +#点的颜色设置
  #scale_fill_manual(values = color)+#椭圆颜色
  theme(axis.title.x=element_text(size=12),#修改X轴标题文本
        axis.title.y=element_text(size=12,angle=90),#修改y轴标题文本
        axis.text.y=element_text(size=10),#修改x轴刻度标签文本
        axis.text.x=element_text(size=10),#修改y轴刻度标签文本
        panel.grid=element_blank())+#隐藏网格线
  ggtitle(paste('Stress=',round(df_nmds_stress, 3)))#添加应力函数值

p.NMDS

ggsave(paste0(Beta_analysis_dir, "Species_NMDS.jpeg"),p.NMDS,dpi=320,height = 10,width=10)

##########################Anosim分析
veg.dist<-distance
group
class(veg.dist)
#head(veg.dist)
species.ano<-anosim(veg.dist,grouping = group$group)
print(species.ano)
summary(species.ano)
plot(species.ano, xlab='Group', ylab='Rank of Distances')
png(file=paste0(Beta_analysis_dir,"Anosim.png"),width = 800,height = 800,res=100)
plot(species.ano)
dev.off()
#ggsave(paste0(Beta_analysis_dir,"Anosim.jpeg"),p.anosim)

#################################Adonis 非参数多元方差分析（PERMANOVA）

veg.dist.ma<-as.matrix(veg.dist)
adonis.group<-group
adonis.group$group<-adonis.group[rownames(veg.dist.ma),]
adonis.group
Adonis.result <- adonis2(veg.dist.ma ~ group, data = adonis.group)
print(Adonis.result)

p.PCoA.Adonis<-p.PCoA+ggtitle(paste("PCoA based on Bray-Curtis Distance", "\nAdonis: R2 =",
                     round(Adonis.result$R2[1], 2), "p =", Adonis.result$'Pr(>F)'[1]))

ggsave(paste0(Beta_analysis_dir, "Species_PCoA_Adonis.jpeg"),p.PCoA.Adonis,dpi=320,height = 10,width=10)
########################堆积条形图######################
######################Phylum#################
library(reshape2)
phylum_prop<-read.table("Phylum_Top15.txt",header=T,check.names=F,stringsAsFactors = F)
phylum_prop<-phylum_prop[,-ncol(phylum_prop)]
colnames(phylum_prop)[1]<-"Phylum"
colnames(phylum_prop)
#phylum_prop$Phylum<-('k__Bacteria\\|','',phylum_prop$Phylum)
head(phylum_prop)
phylum_prop$Phylum <- factor(phylum_prop$Phylum, levels = rev(phylum_prop$Phylum))
phylum_prop[2]
phylum_prop.melt<-melt(phylum_prop)
head(phylum_prop.melt)
phylum_prop.plot<-ggplot(phylum_prop.melt, aes(x = variable, y = value)) +
  geom_bar(aes(fill = Phylum), stat = "identity", width = 0.3) +
  labs(x = "Samples", y = "Relative Abundance (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 75, hjust = 1, vjust = 1, size = 18),
        axis.text.y = element_text(size = 18),
        axis.title = element_text(size = 25),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20)) +
  scale_fill_manual(values = c("#FF9999", "#66B2FF", "#99FF99", "#FFCC99", "#FF66CC",
                               "#CC99FF", "#FFFF66", "#CCE5FF", "#FF6666", "#FFCC66",
                               "#66FFCC", "#CC99CC", "#FF9966", "#99CCFF", "#CCCCFF",
                               "#B3E6B3"))
ggsave(paste0(Classification_dir,"phylum_distribution.jpeg"),phylum_prop.plot,dpi=320,width=30,height=20)
#############################Class#######################
class_prop<-read.table("Class_Top15.txt",header=T,check.names=F,stringsAsFactors = F)
class_prop<-class_prop[,-ncol(class_prop)]
colnames(class_prop)[1]<-"Class"
colnames(class_prop)
class_prop$Class<-gsub("k__Bacteria.*\\|c__",'c__',class_prop$Class,perl=T)
head(class_prop)
class_prop$Class <- factor(class_prop$Class, levels = rev(class_prop$Class))
class_prop[2]
class_prop.melt<-melt(class_prop)
head(class_prop.melt)
class_prop.plot<-ggplot(class_prop.melt, aes(x = variable, y = value)) +
  geom_bar(aes(fill = Class), stat = "identity", width = 0.3) +
  labs(x = "Samples", y = "Relative Abundance (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 75, hjust = 1, vjust = 1, size = 18),
        axis.text.y = element_text(size = 18),
        axis.title = element_text(size = 25),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20)) +
  scale_fill_manual(values = c("#FF9999", "#66B2FF", "#99FF99", "#FFCC99", "#FF66CC",
                               "#CC99FF", "#FFFF66", "#CCE5FF", "#FF6666", "#FFCC66",
                               "#66FFCC", "#CC99CC", "#FF9966", "#99CCFF", "#CCCCFF",
                               "#B3E6B3"))
ggsave(paste0(Classification_dir, "class_distribution.jpeg"),class_prop.plot,dpi=320,width=30,height=20)
##########################################order##################################
order_prop<-read.table("Order_Top15.txt",header=T,check.names=F,stringsAsFactors = F)
order_prop<-order_prop[,-ncol(order_prop)]
colnames(order_prop)[1]<-"Order"
colnames(order_prop)
order_prop$Order<-gsub("k__Bacteria.*\\|o__",'o__',order_prop$Order,perl=T)
head(order_prop)
order_prop$Order <- factor(order_prop$Order, levels = rev(order_prop$Order))
order_prop[2]
order_prop.melt<-melt(order_prop)
head(order_prop.melt)
order_prop.plot<-ggplot(order_prop.melt, aes(x = variable, y = value)) +
  geom_bar(aes(fill = Order), stat = "identity", width = 0.3) +
  labs(x = "Samples", y = "Relative Abundance (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 75, hjust = 1, vjust = 1, size = 18),
        axis.text.y = element_text(size = 18),
        axis.title = element_text(size = 25),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20)) +
  scale_fill_manual(values = c("#FF9999", "#66B2FF", "#99FF99", "#FFCC99", "#FF66CC",
                               "#CC99FF", "#FFFF66", "#CCE5FF", "#FF6666", "#FFCC66",
                               "#66FFCC", "#CC99CC", "#FF9966", "#99CCFF", "#CCCCFF",
                               "#B3E6B3"))
ggsave(paste0(Classification_dir, "order_distribution.jpeg"),order_prop.plot,dpi=320,width=30,height=20)
########################################family################################
family_prop<-read.table("Family_Top15.txt",header=T,check.names=F,stringsAsFactors = F)
family_prop<-family_prop[,-ncol(family_prop)]
colnames(family_prop)[1]<-"Family"
colnames(family_prop)
family_prop$Family<-gsub("k__Bacteria.*\\|f__",'f__',family_prop$Family,perl=T)
head(family_prop)
family_prop$Family <- factor(family_prop$Family, levels = rev(family_prop$Family))
family_prop[2]
family_prop.melt<-melt(family_prop)
head(family_prop.melt)
family_prop.plot<-ggplot(family_prop.melt, aes(x = variable, y = value)) +
  geom_bar(aes(fill = Family), stat = "identity", width = 0.3) +
  labs(x = "Samples", y = "Relative Abundance (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 75, hjust = 1, vjust = 1, size = 18),
        axis.text.y = element_text(size = 18),
        axis.title = element_text(size = 25),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20)) +
  scale_fill_manual(values = c("#FF9999", "#66B2FF", "#99FF99", "#FFCC99", "#FF66CC",
                               "#CC99FF", "#FFFF66", "#CCE5FF", "#FF6666", "#FFCC66",
                               "#66FFCC", "#CC99CC", "#FF9966", "#99CCFF", "#CCCCFF",
                               "#B3E6B3"))
ggsave(paste0(Classification_dir, "family_distribution.jpeg"),family_prop.plot,dpi=320,width=30,height=20)

#######################################genus#################################
genus_prop<-read.table("Genus_Top15.txt",header=T,check.names=F,stringsAsFactors = F)
genus_prop<-genus_prop[,-ncol(genus_prop)]
colnames(genus_prop)[1]<-"Genus"
colnames(genus_prop)
genus_prop$Genus<-gsub("k__Bacteria.*\\|g__",'g__',genus_prop$Genus,perl=T)
head(genus_prop)
genus_prop$Genus <- factor(genus_prop$Genus, levels = rev(genus_prop$Genus))
genus_prop[2]
genus_prop.melt<-melt(genus_prop)
head(genus_prop.melt)
genus_prop.plot<-ggplot(genus_prop.melt, aes(x = variable, y = value)) +
  geom_bar(aes(fill = Genus), stat = "identity", width = 0.3) +
  labs(x = "Samples", y = "Relative Abundance (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 75, hjust = 1, vjust = 1, size = 18),
        axis.text.y = element_text(size = 18),
        axis.title = element_text(size = 25),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20)) +
  scale_fill_manual(values = c("#FF9999", "#66B2FF", "#99FF99", "#FFCC99", "#FF66CC",
                               "#CC99FF", "#FFFF66", "#CCE5FF", "#FF6666", "#FFCC66",
                               "#66FFCC", "#CC99CC", "#FF9966", "#99CCFF", "#CCCCFF",
                               "#B3E6B3"))
ggsave(paste0(Classification_dir, "genus_distribution.jpeg"),genus_prop.plot,dpi=320,width=30,height=20)
###################################species#################################
species_prop<-read.table("Species_Top15.txt",header=T,check.names=F,stringsAsFactors = F)
species_prop<-species_prop[,-ncol(species_prop)]
colnames(species_prop)[1]<-"Species"
colnames(species_prop)
species_prop$Species<-gsub("k__Bacteria.*\\|s__",'s__',species_prop$Species,perl=T)
head(species_prop)
species_prop$Species <- factor(species_prop$Species, levels = rev(species_prop$Species))
species_prop[2]
species_prop.melt<-melt(species_prop)
head(species_prop.melt)
species_prop.plot<-ggplot(species_prop.melt, aes(x = variable, y = value)) +
  geom_bar(aes(fill = Species), stat = "identity", width = 0.3) +
  labs(x = "Samples", y = "Relative Abundance (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 75, hjust = 1, vjust = 1, size = 18),
        axis.text.y = element_text(size = 18),
        axis.title = element_text(size = 25),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20)) +
  scale_fill_manual(values = c("#FF9999", "#66B2FF", "#99FF99", "#FFCC99", "#FF66CC",
                               "#CC99FF", "#FFFF66", "#CCE5FF", "#FF6666", "#FFCC66",
                               "#66FFCC", "#CC99CC", "#FF9966", "#99CCFF", "#CCCCFF",
                               "#B3E6B3"))
ggsave(paste0(Classification_dir, "species_distribution.jpeg"),species_prop.plot,dpi=320,width=30,height=20)
