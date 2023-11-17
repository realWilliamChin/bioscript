setwd("C:/Users/Analysis/Desktop/work/20221202_医学文章数据准备/抗氧化SHSY5Y实验_转录组整理/wgcna")
dir.create("Cytoscape")
library(DESeq2)
library(WGCNA)
library(genefilter)
library(tidyverse)
allowWGCNAThreads(32)
options(stringsAsFactors = FALSE)
getwd()
rm(list=ls())
###############读入raw data，DESeq normalize之后保存备份#####################
sp_name<-"kangai"
infile1<-paste0(sp_name,"_reads_matrix_filtered.txt",sep="")
infile2<-paste0(sp_name,"_samples_described.txt",sep="")
data<-read.table(infile1,sep="\t",header=T,stringsAsFactors = F,check.names = F,row.names=1,quote='')
head(data)
colnames(data)
meta<-read.table(infile2,sep="\t",header=T,check.names = F)
rownames(meta)<-meta$sample
#meta<-meta[1]
meta<-meta[match(colnames(data),meta$sample),]
meta$group<-factor(meta$group)
meta
all(colnames(data) %in% rownames(meta))
all(colnames(data) == rownames(meta))

meta
dds <- DESeqDataSetFromMatrix(countData = round(data), colData = meta, design = ~ group)
dds <- DESeq(dds)
vsd <- varianceStabilizingTransformation(dds)
head(vsd)
class(vsd)
#wpn_vsd<- getVarianceStabilizedData(dds)
#rv_wpn <- rowVars(wpn_vsd)
#q75_wpn <- quantile( rowVars(wpn_vsd), .75)  # <= original
#q95_wpn <- quantile( rowVars(wpn_vsd), .95)  # <= changed to 95 quantile to reduce dataset
#expr_normalized <- wpn_vsd[ rv_wpn > q95_wpn, ]
normalized_counts <- getVarianceStabilizedData(dds)
head(normalized_counts)
range(normalized_counts)
class(normalized_counts)
dim(normalized_counts)
rv_normalized_counts <- rowVars(normalized_counts)
head(rv_normalized_counts)
summary(rv_normalized_counts)
rownames(normalized_counts)
#dim(data)
#View(counts(dds))
#dds <- estimateSizeFactors(dds)
#sizeFactors(dds)
#normalized_counts <- counts(dds, normalized=TRUE)
head(normalized_counts)
normalized_counts<-cbind(as.data.frame(rownames(normalized_counts)),normalized_counts)
colnames(normalized_counts)[1]<-"GeneID"
outfile1<-paste0(sp_name,"_DESeq_counts.txt",sep="")
write.table(normalized_counts, file=outfile1, sep="\t", quote=F, col.names=T,row.names = F)
dim(normalized_counts)
################################anova分析筛选变化大的基因##########################
tmp.group<-meta$group
rownames(normalized_counts)
colnames(normalized_counts)
normalized_counts<-normalized_counts[,-1]
p<-NULL
for (i in rownames(normalized_counts)) {
  #print(i);
  #print(data[i,])
  tmp.count<-as.numeric(normalized_counts[i,])
  #tmp.dt<-data.frame(tmp.count,group=rep(c("S1","S3","S2"),each=3))
  tmp.dt<-data.frame(tmp.count,group=tmp.group)
  colnames(tmp.dt)[1]<-i
  tmp.dt$group<-as.factor(tmp.dt$group)
  mod <- aov( tmp.dt[,1]~group , data = tmp.dt)
  p<-c(p,summary(mod)[[1]][5]$Pr[1])
}
length(p)
normalized_counts$p_value<-p
normalized_counts<-normalized_counts[order(normalized_counts$p_value),]
p.adjust(normalized_counts$p_value,method = "BH")->p.aj
normalized_counts$BH_p_value<-p.aj
df.rname<-data.frame(GeneID=rownames(normalized_counts))
normalized_counts<-cbind(df.rname,normalized_counts)
head(normalized_counts)
outfile2<-paste0(sp_name,"_DESeq_anova_p.txt",sep="")
write.table(normalized_counts,outfile2,sep="\t",quote=F,row.names = F)
###################################获得wgcna的输入数据############################
dim(normalized_counts)
head(normalized_counts)
dt_wgcna<-normalized_counts[normalized_counts$BH_p_value<=0.01,]
dt_wgcna<-dt_wgcna[1:10000,1:(ncol(dt_wgcna)-2)]
dim(dt_wgcna)
class(dt_wgcna)
colnames(dt_wgcna)
rownames(dt_wgcna)
meta$sample
expr_normalized_df<-pivot_longer(dt_wgcna,-GeneID)
#expr_normalized_df$name
expr_normalized_df %>% ggplot(., aes(x = factor(name,levels=meta$sample), y = value)) +
  geom_violin() +
  geom_point() +
  theme_bw() +
  theme(
    axis.text.x = element_text( angle = 90)
  ) +
  ylim(0, NA) +
  labs(
    title = "Normalized Gene Expression",
    x = "Sample",
    y = "normalized expression"
  )->exp_violin.plot
ggsave("normalized_exp_violin.pdf",exp_violin.plot,dpi=320)
##########################start wgcna analysis dt_wgcna ##################
#########################dt_wgcna 为数据框，第一列为GeneID，行名也是GeneID
datExpr0 = as.data.frame(t(dt_wgcna[, -c(1)]))
colnames(datExpr0)
rownames(datExpr0)
#输入数据的预处理
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK
#如果数据集中含有过多的缺失值，则对数据集执行下列代码
if (!gsg$allOK){
  
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
###################################################
pdf(file = "wgcna-01hclust.pdf", width = 12, height = 9);
sampleTree = hclust(dist(datExpr0), method = "average");
par(cex = 0.6); 
par(mar = c(0,4,2,0))  
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,  cex.axis = 1.5, cex.main = 2)
abline(h = 50000,col="red")
dev.off()
clust = cutreeStatic(sampleTree, cutHeight = 50000, minSize = 10) 
table(clust)
########如果要删除异常样本，则运行下面脚本，如果要保留全部样本则忽略下面4行命令
keepSamples = (clust==1) 
datExpr = datExpr0[keepSamples, ] 
nGenes = ncol(datExpr) 
nSamples = nrow(datExpr)
#################上面4行命令可以选择性忽略###############################
#################如保留所有样品，则跳过上面4行命令，直接运行下面命令#############
datExpr = datExpr0 
save(datExpr,file="wgcna-01-dataInput.RData")
lnames = load(file = "wgcna-01-dataInput.RData"); 
lnames
#确定最佳BETA值
#一系列数值，从中选择最优的
powers = c(c(1:10), seq(from = 12, to=20, by=2)) 
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5) 
#结果作图
pdf(file="wgcna-02-1.pdf",width=9,height=5) 
par(mfrow = c(1,2));
cex1 = 0.8;  
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],  xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",  main = paste("Scale independence"));  
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=cex1,col="red"); 
#在图中高度为0.8的位置画一条红线，这个值对应的是相关系数的平方R^2，此值必须大于0.8，一般选择0.9


abline(h=0.8,col="red")  
plot(sft$fitIndices[,1], sft$fitIndices[,5],  xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))  
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
sft$powerEstimate
#如果上一步最佳估计的结果是[1] NA，自己根据绘制的图，选择一个最佳的power，选择wgcna-02-1.pdf右边图里下降曲线变平滑的那个值

# !!!!这里一个具体的例子，选了10.但是不同的数据，肯定是不同的。
sft$powerEstimate<-14
#一步法构建共表达矩阵,注意这个power值的选择，此power值应选择wgcna-02-1.pdf右边图里下降曲线变平滑处的值
net = blockwiseModules(datExpr, power = sft$powerEstimate, TOMType = "unsigned", minModuleSize = 30, reassignThreshold =0, mergeCutHeight = 0.25, numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = TRUE,  saveTOMFileBase = "fpkmTOM", verbose = 3)
table(net$colors)
table(labels2colors(net$colors))

pdf(file="wgcna-02-2.pdf",width=12,height=9) 
mergedColors = labels2colors(net$colors) 
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors",  dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)

dev.off()
#将计算的得到的变量保存，以便下次使用。
moduleLabels = net$colors  
moduleColors = labels2colors(net$colors)
MEs = net$MEs;  
geneTree = net$dendrograms[[1]];  
save(net, MEs, moduleLabels, moduleColors, geneTree,net, file = "wgcna-02-networkConstruction-auto.RData")
rownames(MEs)

#载入上步保存的数据
lnames = load(file = "wgcna-01-dataInput.RData")
#查看载入的变量名
lnames
lnames = load(file = "wgcna-02-networkConstruction-auto.RData");
lnames
names(datExpr)
#TOM矩阵计算，层次聚类和动态切割
############这个项目里面，power=9.不同的项目肯定不同！！！！################
TOM = TOMsimilarityFromExpr(datExpr, power = sft$powerEstimate)

save(net, MEs, moduleLabels, moduleColors, geneTree,net,TOM, file = "wgcna-03-TOM2networkConstruction-auto.RData")
# 计算基因之间的相异度
dissTOM = 1-TOM
hierTOM = hclust(as.dist(dissTOM),method="average")

geneTree = hclust(as.dist(dissTOM), method = "average")
# Plot the resulting clustering tree (dendrogram)
#windows()
#sizeGrWindow(12,9)
pdf(file="2ndwgcna-02-2.pdf",width=12,height=9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",     labels = FALSE, hang = 0.04)
dev.off()

#画TOM图，非常耗时
dissTOM = 1-TOM;
plotTOM = dissTOM^7;
diag(plotTOM) = NA;
tiff(file="fpkm-heatmap.tif")
#it takes 10 min
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes");
dev.off()

# 使用动态剪切树构建模块：
minModuleSize = 100;
# 动态切割树:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize);
table(dynamicMods)
dynamicMods

dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

png("dynamic_tree_cut.png",width = 800,height = 600)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05,main = "Gene dendrogram and module colors")
dev.off()

MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);#计算模块和模块之间的相关性和相异性
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
#sizeGrWindow(7, 6)
#png("Clustering_of_module_eigengenes1.png",width = 800,height = 600)

png("Clustering_of_module_eigengenes1_0.2.png",width = 800,height = 600)

plot(METree, main = "Clustering of module eigengenes",     xlab = "", sub = "")
#######################这个地方需要手动选择理想的cutoff,这个值越大，模块数越少##################################
##################MEDissThres值的选择根据Clustering_of_module_eigengenes1.png图来决定，值越大，模块数越少
MEDissThres = 0.2

# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
dev.off()

png("Clustering_of_module_eigengenes1_0.3.png",width = 800,height = 600)

plot(METree, main = "Clustering of module eigengenes",     xlab = "", sub = "")
#######################这个地方需要手动选择理想的cutoff,这个值越大，模块数越少##################################
##################MEDissThres值的选择根据Clustering_of_module_eigengenes1.png图来决定，值越大，模块数越少
MEDissThres = 0.3

# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
dev.off()
png("Clustering_of_module_eigengenes1_0.1.png",width = 800,height = 600)

plot(METree, main = "Clustering of module eigengenes",     xlab = "", sub = "")
#######################这个地方需要手动选择理想的cutoff,这个值越大，模块数越少##################################
##################MEDissThres值的选择根据Clustering_of_module_eigengenes1.png图来决定，值越大，模块数越少
# kanzhege Clustering_of_module_eigengenes1_0.1.png
MEDissThres = 0.1

# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
dev.off()

# Call an automatic merging function MEDissThres=0.1
# !!!
MEDissThres = 0.1
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
str(merge)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs
png("dynamicColors_mergedColors.png",width = 800,height = 600)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),c("Dynamic Tree Cut", "Merged dynamic"),dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)
dev.off()
#绘制两两模块间的邻接矩阵
png("wgcna.adjacency.heatmap.png",height = 1000,width = 900)
plotEigengeneNetworks(MEs, "Eigengene adjacency heatmap",plotDendrograms = F,marDendro = c(4,4,2,4))
dev.off()

moduleColors = merge$colors
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0); 
MEs
colnames(MEs)
module_order = names(MEs) %>% gsub("ME","", .)
MEs$treatment = row.names(MEs)
head(MEs)
class(MEs)
row.names(MEs)
colnames(MEs)
MEs.out<-MEs[c(ncol(MEs),1:(ncol(MEs)-1))]
colnames(MEs.out)<-gsub("ME", "", colnames(MEs.out))
write.table(MEs.out,"sample_module_cor.txt",sep="\t",row.names=F,quote=F)

mME = MEs %>%
  pivot_longer(-treatment) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )

head(mME)
class(mME)
mME$treatment<-factor(mME$treatment,levels=unique(mME$treatment))

mME %>% ggplot(., aes(x=treatment, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-trait Relationships", y = "Modules", fill="corr")->mME.plot
ggsave("sample_module_cor.pdf",mME.plot,dpi=320)


MEs<-MEs[1:(ncol(MEs)-1)]
colnames(MEs)
for(module in substring(colnames(MEs),3)){
  #if(module == "grey") next
  ME=as.data.frame(MEs[,paste("ME",module,sep="")])
  colnames(ME)=module
  datModExpr=datExpr[,moduleColors==module]
  datKME = signedKME(datModExpr, ME)
  
  datKME=cbind(datKME,rep(module,length(datKME)))
  write.table(datKME,quote = F,row.names = T,append = T,file = "All_Gene_KME.txt",col.names = F,sep="\t")
}

for(module in substring(colnames(MEs),3)){
  #if(module == "grey") next
  ME=MEs[,paste("ME",module,sep="")]
  png(paste("wgcna.", module, ".express.barplot.png", sep=""),height = 700,width = 900)
  par(mfrow=c(2,1),mar=c(0.3,5.5,3,2))
  plotMat(t(scale(datExpr[,moduleColors==module])),rlabels=F,main=module,cex.main=2,clabels=F)
  names(ME)<-rownames(MEs)
  par(mar=c(5,4.2,0,0.7))
  barplot(ME,col=module,main="",cex.main=2,ylab="eigengene expression",xlab="",las=2)
  dev.off()
}

# !!! threshold = 0.2 看一下这个文件edge_weighted_value.txt
for(module in substring(colnames(MEs),3)){
  #if(module == "grey") next
  probes = colnames(datExpr)
  inModule = is.finite(match(moduleColors, module))
  modProbes = probes[inModule]
  modTOM = TOM[inModule, inModule]
  dimnames(modTOM) = list(modProbes, modProbes)
  cyt = exportNetworkToCytoscape(modTOM,edgeFile = paste("CytoscapeInput-edges-", module, ".txt", sep=""),nodeFile = paste("CytoscapeInput-nodes-", module, ".txt", sep=""),weighted = TRUE,threshold = 0.2,nodeNames = modProbes,nodeAttr = moduleColors[inModule])
}
##############################################################################################################
dir()
edge_files<-grep("CytoscapeInput-edges-",dir(),value = T)
class(edge_files)
weight_value<-NULL
weight_value_df<-as.data.frame(matrix(1:8,nrow=1))
i<-NULL
for (i in 1:length(edge_files)) {
  temp_df <-read.table(edge_files[i],header=T,sep="\t",check.names = F,quote='',stringsAsFactors = F)
  temp_color<-gsub("CytoscapeInput-edges-",'',edge_files[i])
  temp_color<-gsub(".txt",'',temp_color)
  if (nrow(temp_df)==0) {
    weight_value_df[i,]<-c(temp_color,nrow(temp_df),as.numeric(quantile(0)),0)
    next
    
  }
  weight_value_df[i,]<-c(temp_color,nrow(temp_df),as.numeric(quantile(temp_df[,3])),mean(temp_df[,3]))
  weight_value<-c(weight_value,temp_df[,3])
}
weight_value<-as.numeric(weight_value)
length(weight_value)
colnames(weight_value_df)<-c("module","edge_count","q0","q25","q50","q75","q100","mean_value")
write.table(weight_value_df,"edge_weighted_value.txt",sep="\t",quote=F,row.names=F)
weight_value_df
#WeightValue_60<-quantile(weight_value,0.6)
par(mar=c(5.1,4.1,4.1,2.1))
pdf("histogram_edge_weight_value.pdf",width=9,height=5)
hist(weight_value,breaks=100,col = "light blue")
dev.off()
png("histogram_edge_weight_value.png",height = 1000,width = 900)
hist(weight_value,breaks=100,col = "light blue")
dev.off()
for (sub_edge_file in edge_files) {
  edge_tmp_df<-read.table(sub_edge_file,header=T,sep="\t",check.names = F,quote='',stringsAsFactors = F)
  if(nrow(edge_tmp_df) == 0){next}
  #edge_tmp_df<-edge_tmp_df[edge_tmp_df[,3]>=WeightValue_60,]
  node_freq<-as.data.frame(table(edge_tmp_df[,1]))
  node_freq<-node_freq[order(-node_freq[,2]),]
  if (nrow(node_freq)>100) {
    node_freq<-node_freq[1:100,]
    
  }
  sub_edge_file_out<-paste0(sub_edge_file,"_Top100.txt",sep="")
  edge_out_df<-edge_tmp_df[edge_tmp_df[,1] %in% node_freq[,1],]
  write.table(edge_out_df,file = sub_edge_file_out,row.names=F,quote=F,sep="\t")
  sub_node_file <- gsub("edges","nodes",sub_edge_file)
  node_tmp_df<-read.table(sub_node_file,header=T,sep="\t",check.names = F,quote='',stringsAsFactors = F)
  node_out_name<-unique(c(edge_out_df[,1],edge_out_df[,2]))
  node_out_df<-node_tmp_df[node_tmp_df[,1] %in% node_out_name,]
  sub_node_file_out<-paste0(sub_node_file,"_Top100.txt",sep="")
  write.table(node_out_df,file = sub_node_file_out,row.names=F,quote=F,sep="\t")
}

#head(node_freq)
