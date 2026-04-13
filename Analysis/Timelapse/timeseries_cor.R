library(Mfuzz)
library(corrplot)
library(WGCNA)
rm(list=ls())
expression_data<-read.table("fpkm_matrix_filter_timelapse.txt",sep="\t",row.names=1,header = T,check.names = F)
colnames(expression_data)
expression_data <- log2(expression_data + 1)
expression_data <- scale(expression_data)
expression_data <- na.omit(expression_data)
############# choose the number of clusters ############
############# choose the number of clusters ############
wss <- (nrow(expression_data)-1)*sum(apply(expression_data,2,var))
for (i in 2:15) { 
  wss[i] <- sum(kmeans(expression_data,centers=i)$withinss)}
png("kmeans_number.jpeg")
plot(1:15, wss, type="b", xlab="Number of Clusters",ylab="Within groups sum of squares")
dev.off()
############# select the cluster numbers according to the picture: kmeans_number.jpeg############
cluster_num <- 7 ######每次看图决定聚类cluster个数
###########################################################################################
eset <- new("ExpressionSet", exprs = as.matrix(expression_data))
eset<- standardise(eset)
head(eset)
length(colnames(expression_data))

mfuzz_cluster <- mfuzz(eset, c = cluster_num, m = mestimate(eset))
mfrow_row <- floor(sqrt(cluster_num))
mfrow_col <- ceiling(cluster_num / mfrow_row)
# 根据 cluster_num 设置图片的大小
timelapse_width <- 5 * mfrow_col
timelapse_height <- 5 * mfrow_row
png("timelapse.png", width = timelapse_width, height = timelapse_height, units = 'in', res = 300)
mfuzz.plot2(eset, cl = mfuzz_cluster, mfrow = c(mfrow_row, mfrow_col), time.labels = colnames(expression_data),x11 = FALSE)
dev.off()
####################提取每一个cluster单独的图
for(i in 1:cluster_num){
  out_pic<-paste0("cluster",i,".png")
  png(out_pic)
  mfuzz.plot2(eset, cl = mfuzz_cluster, mfrow = c(1, 1),time.labels = colnames(expression_data),single = i,x11 = FALSE)
  dev.off()}


# 提取每个聚类的基因
mfuzz_cluster$size

cluster_genes <- list()
sample_pc<-matrix(rep(1,length(colnames(expression_data))*cluster_num),nrow =length(colnames(expression_data)))
sample_pc
colnames(sample_pc)<-paste0("cluster",1:cluster_num)
rownames(sample_pc)<-colnames(expression_data)
sample_pc
sample_pc[,1]
for (i in 1:cluster_num) {
  cluster_genes[[i]] <- names(mfuzz_cluster$cluster[mfuzz_cluster$cluster == i])
  expression_df<-expression_data[cluster_genes[[i]],]
  #pca_result<-prcomp(t(expression_df))
  #sample_pc[,i]<-pca_result$x[,1]
  sample_pc[,i] <- colMeans(expression_df)
  out_file<-paste0("cluster",i,"_gene.txt")
  write.table(as.data.frame(cluster_genes[[i]]),out_file,quote=F,row.names = F,col.names = F)
}
#expression_df<-expression_data[cluster_genes[[1]],]
#pca_result<-prcomp(t(expression_df))
#pca_result$x[,1]
#sample_pc[,1]<-pca_result$x[,1]
#pca_result$x

sample_pc
group<-factor(rownames(sample_pc),levels = colnames(expression_data))
group_dummy <- model.matrix(~ group - 1)
group_dummy
cluster_group_cor <- cor(sample_pc, group_dummy, use = "p")
cluster_group_cor
cluster_group_p <- corPvalueStudent(cluster_group_cor, nrow(sample_pc))
colnames(cluster_group_p)<-gsub("group",'',colnames(cluster_group_p))
cluster_group_p
colnames(cluster_group_cor)<-colnames(cluster_group_p)
cluster_group_cor
my_colors <- colorRampPalette(c("green", "white", "red"))(200)
cor_matrix<-cluster_group_cor
p_matrix<-cluster_group_p
plot_corr_with_p <- function(cor_matrix, p_matrix) {
  # 绘制相关性热图
  corrplot(
    cor_matrix,
    method = "color",          # 使用颜色表示相关性
    type = "full",             # 显示完整的矩阵
    order = "original",        # 保持原始顺序
    diag = TRUE,              # 显示对角线
    tl.col = "black",          # 标签颜色
    tl.srt = 45,               # 标签旋转角度
    addCoef.col = "black",     # 添加相关系数值，颜色为黑色
    number.cex = 0.8,          # 相关系数字体大小
    col = my_colors,           # 使用自定义颜色调色板
    cl.pos = "r"               # 图例在右侧
  )
  
  # 在图中添加 P 值
  n <- nrow(cor_matrix)
  m <- ncol(cor_matrix)
  for (i in 1:n) {
    for (j in 1:m) {
      if (!is.na(p_matrix[i, j])) {
        # 在相关系数下方显示 P 值
        text(j, n-i+1, 
             labels = sprintf("(%.3f)", p_matrix[i, j]), 
             cex = 0.6, col = "black", pos = 1, offset = 1.2)
      }
    }
  }
}
png("corrplot_with_p_values.png", width = 800, height = 800, res = 150)  # 设置文件名称和尺寸
plot_corr_with_p(cluster_group_cor, cluster_group_p)  # 调用自定义函数绘制图形
dev.off()  # 关闭图形设备
cluster_group_cor
cluster_group_p
############the end
#################################################################################################################

# 将分组信息转换为虚拟变量矩阵
#group_dummy <- model.matrix(~ group - 1)
#group_dummy 
# 进行典型相关分析
#library(CCA)
#install.packages("CCA")
#cca_result <- cc(t(gene_expression), group_dummy)

# 提取第一对典型变量的相关系数
#correlation <- cca_result$cor[1]

# 输出结果
#print(paste("基因整体表达与分组的典型相关系数:", correlation))
#?cc

#####################################
# 生成示例数据
#set.seed(123)
#gene_expression <- matrix(rnorm(100 * 5), nrow = 100)
#rownames(gene_expression) <- paste0("Gene", 1:100)
#colnames(gene_expression) <- paste0("Sample", 1:5)
#group <- factor(1:5)

# 对基因表达数据进行主成分分析
#pca_result <- prcomp(t(gene_expression))

# 提取所有主成分
#pc_scores <- pca_result$x
#tmp.df<-as.data.frame(pc_scores[,1])
#rownames(tmp.df)
#tmp.dt<-as.data.frame(pc_scores[,1:5])
#colnames(tmp.dt)<-paste0("cluster",1:5)
#rownames(tmp.dt)
#tmp.dt
#group_dummy <- model.matrix(~ group - 1)
#group_dummy
#module_group_cor <- cor(tmp.dt, group_dummy, use = "p")
#tmp.dt
#module_group_cor
#library(WGCNA)
#module_group_p <- corPvalueStudent(module_group_cor, nrow(tmp.dt))
#module_group_p
#library(corrplot)
#getwd()
#my_colors <- colorRampPalette(c("green", "white", "red"))(200)
#cor_matrix<-module_group_cor
#p_matrix<-module_group_p
#plot_corr_with_p <- function(cor_matrix, p_matrix) {
  # 绘制相关性热图
 # corrplot(
  #  cor_matrix,
   # method = "color",          # 使用颜色表示相关性
    #type = "full",             # 显示完整的矩阵
    #order = "original",        # 保持原始顺序
    #diag = TRUE,              # 显示对角线
    #tl.col = "black",          # 标签颜色
    #tl.srt = 45,               # 标签旋转角度
    #addCoef.col = "black",     # 添加相关系数值，颜色为黑色
    #number.cex = 0.8,          # 相关系数字体大小
    #col = my_colors,           # 使用自定义颜色调色板
    #cl.pos = "r"               # 图例在右侧
  #)
  
  # 在图中添加 P 值
 # n <- nrow(cor_matrix)
#  for (i in 1:n) {
 #   for (j in 1:n) {
  #    if (!is.na(p_matrix[i, j])) {
   #     # 在相关系数下方显示 P 值
    #    text(j, i, 
     #        labels = sprintf("(%.3f)", p_matrix[i, j]), 
      #       cex = 0.8, col = "gray40", pos = 1, offset = 1.2)
      #}
    #}
  #}
#}
#png("corrplot_with_p_values.png", width = 800, height = 800, res = 150)  # 设置文件名称和尺寸
#plot_corr_with_p(module_group_cor, module_group_p)  # 调用自定义函数绘制图形
#dev.off()  # 关闭图形设备
