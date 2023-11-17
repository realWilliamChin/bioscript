library(Mfuzz)
library(factoextra)
library(cluster)
setwd("d:/pll/R_work/RNA-seq/time_series/")
rm(list=ls())
protein <- read.delim('lappa_seed_annovafiltered_fpkm.txt', row.names = 1, check.names = FALSE,header=T)
colnames(protein)

####protein<-protein[,-9]
dim(protein)
class(protein)
protein <- as.matrix(protein)
df <- protein
df <- na.omit(df)
df <- scale(df)


dim(df)
############# choose the number of clusters ############
wss <- (nrow(df)-1)*sum(apply(df,2,var))
for (i in 2:15) 
  wss[i] <- sum(kmeans(df,centers=i)$withinss)
png("kmeans_number.jpeg")
plot(1:15, wss, type="b", xlab="Number of Clusters",ylab="Within groups sum of squares")
dev.off()
############# select the cluster numbers according to the picture: kmeans_number.jpeg############
cluster_num <- 8 ###########check the number and change it everytime 
dim(df)
mfuzz_class <- new('ExpressionSet',exprs = protein)
mfuzz_class <- filter.NA(mfuzz_class, thres = 0.25)
mfuzz_class <- fill.NA(mfuzz_class, mode = 'mean')
mfuzz_class <- filter.std(mfuzz_class, min.std = 0)
mfuzz_class <- standardise(mfuzz_class)
set.seed(123)
mfuzz_cluster <- mfuzz(mfuzz_class, c = cluster_num, m = mestimate(mfuzz_class))
mfuzz.plot2(mfuzz_class, cl = mfuzz_cluster, mfrow = c(3, 3), time.labels = colnames(protein))
cluster_size <- mfuzz_cluster$size
names(cluster_size) <- 1:cluster_num
cluster_size
head(mfuzz_cluster$membership)
protein_cluster <- mfuzz_cluster$cluster
class(mfuzz_cluster$cluster)
protein_cluster <- cbind(protein[names(protein_cluster), ], protein_cluster)
protein_cluster<-as.data.frame(protein_cluster)
rownames(protein_cluster)
colnames(protein_cluster)
protein_cluster<-cbind(rownames(protein_cluster),protein_cluster)
colnames(protein_cluster)[1]<-"GeneID"
write.table(protein_cluster, 'all_gene_cluster.txt', sep = '\t', quote = FALSE,row.names = F)
class(protein_cluster)
for(i in 1:cluster_num){
  outfile<-paste0("cluster",i,".txt")
  write.table(subset(protein_cluster,protein_cluster==i),outfile,sep="\t",quote=F,row.names = F)
}




