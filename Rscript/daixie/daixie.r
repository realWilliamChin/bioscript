setwd("C:/Users/Analysis/OneDrive/Work/zhiwusuo/01_Work/04_daixie/2023_01_27_大鼠/2023_11_28_蛋白分析/shizhuofei")
library(ggthemes)
library(ggplot2)
library(pheatmap)
library(reshape2)
library(ggcorrplot)
library(corrplot)
library(FactoMineR)
library(plyr)
library(ropls)
library(ggrepel)
library(mixOmics)
rm(list=ls())

# 读取文件
reads_data<-read.table("All_sample_data.txt",sep="\t",row.names=1,header=T,check.names=F,stringsAsFactors = F)
sample_info<-read.table("samples_described.txt",sep="\t",header=T,check.names=F,stringsAsFactors = F)
# 如果需要合并定义则读取单独定义文件，没有则跳过
definition_df<-read.table("def.txt",sep="\t",row.names=1,header=T,check.names=F,stringsAsFactors = F,quote="")
definition_df$Metabolite <- rownames(definition_df)

fpkm<-as.data.frame(t(apply(reads_data,1,function(x){(x-mean(x))/(sd(x)**0.5)})))


# heatmap
# 化合物多于 100 的 ，show_rownames=F
all.heatmap<-pheatmap(fpkm,scale="row",cluster_cols=F,show_rownames=F)
ggsave("All_metabolites_heatmap.jpeg",all.heatmap,dpi=300,width=10,height=10,limitsize=FALSE)


# correlation
reads_data<-na.omit(reads_data)
reads_data<-reads_data + 0.000000001
fpkm.m<-as.matrix(fpkm)
fpkm.cor<-cor(fpkm.m)
correlation_df<-as.data.frame(fpkm.cor)
correlation_df$ID<-rownames(correlation_df)
correlation_df <- correlation_df[, c("ID", setdiff(names(correlation_df), "ID"))]
write.table(correlation_df,"Metabolite_correlation.txt",sep="\t",quote=F,row.names=F)
min(fpkm.cor)
resfactor = 5
png("Metabolite_correlation_graph.png",res=72*resfactor,height=1200*resfactor,width=1200*resfactor)
corrplot(fpkm.cor,is.corr = F,col= rev(COL2('PiYG')),method="color",addCoef.col = 'black',tl.col="black",col.lim=c(min(fpkm.cor)-0.01,max(fpkm.cor)),cl.ratio=0.1)
dev.off()
pdf("Metabolite_correlation_graph.pdf",width = 15,height = 15)
corrplot(fpkm.cor,is.corr = F,col= rev(COL2('PiYG')),method="color",addCoef.col = 'black',tl.col="black",col.lim=c(min(fpkm.cor)-0.01,max(fpkm.cor)),cl.ratio=0.1)
dev.off()


# pca 图
gene <- log2(fpkm+1)
gene <- t(gene)
gene.pca <- PCA(gene, ncp = 2, scale.unit = TRUE, graph = FALSE)
gene.pca$eig
summary(gene.pca)
pca_sample <- data.frame(gene.pca$ind$coord[ ,1:2])
head(pca_sample)
pca_eig1 <- round(gene.pca$eig[1,2], 2)
pca_eig2 <- round(gene.pca$eig[2,2],2 )
group <- read.delim('samples_described.txt', row.names = 2, sep = '\t', check.names = FALSE,header=T)
group <- group[rownames(pca_sample), ]
p <- ggplot(data = pca_sample, aes(x = Dim.1, y = Dim.2)) +
  geom_point(aes(color = group), size = 5) + 
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent')) + 
  labs(x =  paste('PCA1:', pca_eig1, '%'), y = paste('PCA2:', pca_eig2, '%'), color = '')
p<-p+geom_text_repel(data=pca_sample,aes(Dim.1, Dim.2, label=rownames(pca_sample)))
cluster_border <- ddply(pca_sample, 'group', function(df) df[chull(df[[1]], df[[2]]), ])
p<-p + geom_polygon(data = cluster_border, aes(color = group),fill=NA, show.legend = FALSE)
ggsave("Metabolite_PCA_analysis.jpeg",p,dpi=300,width=10,height=10)


# 代谢的图
dir.create("多组分析")
fpkm_t<-t(fpkm)
metabolites<-as.matrix(fpkm_t)
groups=sample_info$group

df_plsda <- plsda(fpkm_t, groups, ncomp = 4)

# scree plot
scree_df <- as.data.frame(df_plsda$prop_expl_var$X)
scree_df$comp=rownames(scree_df)
colnames(scree_df)[1]="value"

df<-data.frame(x=letters[1:6],y=seq(10,60,10))
# 创建 scree plot 图
scree_plot <- ggplot(scree_df, aes(x = comp, y = value)) +
#  geom_bar(stat = "identity", fill = "steelblue") +
  labs(x = "Principal Component", y = "Explained Variance") +
  ggtitle("Scree Plot") +
  theme_minimal()+
  geom_line(aes(group = 1), color = "red")+
  geom_point()

ggsave('多组分析/Metabolite_quantitation_scree_plot.jpeg', scree_plot)

df <- unclass(df_plsda)

df1 = as.data.frame(df$variates$X)
df1$group = sample_info$group
df1$samples = rownames(df1)

explain = df$prop_expl_var$X
x_lable <- round(explain[1],digits=3)
y_lable <- round(explain[2],digits=3)

num_samples <- nrow(df1)
plot_width <- 5 + num_samples * 0.1
plot_height <- 4 + num_samples * 0.1

p1 <- ggplot(df1, aes(x = comp1, y = comp2, color = group, shape = group)) +
  theme_bw() +
  geom_point(size = 1.8) +
  theme(panel.grid = element_blank()) +  # 设置高度
  geom_text(aes(label = samples, y = comp2 + 0.4, x = comp1 + 0.5, vjust = 0), size = 3.5) +
  labs(x = paste0("P1 (", x_lable * 100, "%)"), y = paste0("P2 (", y_lable * 100, "%)")) +
  stat_ellipse(data = df1, geom = "polygon", level = 0.95,
               linetype = 2, linewidth = 0.5, aes(fill = group),
               alpha = 0.2, show.legend = TRUE) +
  scale_color_discrete() +
  scale_fill_discrete() +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12, angle = 90),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        panel.grid = element_blank())

ggsave("多组分析/Multigroup_PCA_Distribution_Graph.jpeg", p1, width = plot_width, height = plot_height)


# vip 值使用 fpkm, 其他值计算都是用 reads
plsda_model <- opls(x = metabolites, y = groups, predI = 1)
vip_values <- plsda_model@vipVn
df.vip<-as.data.frame(vip_values)
reads_data_with_def <- cbind(reads_data, VIP = df.vip$vip_values)

reads_data_with_def$Metabolite <- rownames(reads_data_with_def)
reads_data_with_def <- reads_data_with_def[, c("Metabolite", setdiff(names(reads_data_with_def), "Metabolite"))]

if (exists("definition_df")) {
  reads_data_with_def<-merge(reads_data_with_def, definition_df, by = "Metabolite", all.x=TRUE)
}
reads_data_with_def <- reads_data_with_def[order(-reads_data_with_def$VIP, na.last = TRUE), ]

write.table(reads_data_with_def, file="多组分析/Metabolite_quantitation_VIP.txt", sep='\t', row.names = FALSE, col.names = TRUE, quote=FALSE)


# 自动生成组间分析
group_levels <- unique(sample_info$group)
comparisons <- combn(group_levels, 2, simplify = FALSE)

# 指定组间分析
comp_info<-read.table("compare_info.txt",sep="\t",header=T,check.names=F,stringsAsFactors = F)
comparisons <- list()
for(i in seq_along(1:nrow(comp_info))){
  comparisons <- append(comparisons, list(as.character(comp_info[i,])))
}

plot_types <- c('correlation', 'outlier', 'overview', 'permutation', 
                'predict-train', 'x-loading', 'x-score', 
                'x-variance', 'xy-score')



dir.create('组间分析')

deg_data <- data.frame(group = character(0), All=numeric(0), Up = numeric(0), Down = numeric(0))

# 循环中注意可能需要修改 corssvalI 值
for (i in seq_along(comparisons)) {
  
  key_name <- paste(comparisons[[i]], collapse = "_vs_")
  dir_create_path <- paste0('组间分析/', key_name)
  dir.create(dir_create_path, showWarnings = FALSE)
  print(paste("Processing:", key_name)) # 打印当前处理的组合名称
  groups <- comparisons[[i]]
  groupA_cols <- sample_info$sample[sample_info$group == groups[1]]
  groupB_cols <- sample_info$sample[sample_info$group == groups[2]]
  
  # 获取当前两组的样本的表达量数据
  current_samples <- sample_info$sample[sample_info$group %in% groups]
  
  current_expression_data <- reads_data[, current_samples, drop = FALSE]
  current_expression_fpkm_data <- fpkm[, current_samples, drop = FALSE]
  #p_value_current_expression_data <- reads_data[, current_samples, drop = FALSE]
  # 检查因变量的水平数量
  y_factor <- factor(sample_info$group[sample_info$sample %in% current_samples])
  if(length(unique(y_factor)) < 2) {
    print(paste("y_factor 小于 2 个", key_name))
    next # 跳过当前的迭代
  }
  
  # 检查是否有缺失值
  if(any(is.na(current_expression_data))) {
    print(paste("Missing values found in expression data for:", key_name))
    next # 跳过当前的迭代
  }
  
  # 进行OPLS分析
  transposed_expression_data <- t(current_expression_fpkm_data)
  
  opls_model <- try(
    # crossvalI 默认是 7, crossvalI 需要小于等于两组样本的数量
    opls(x = transposed_expression_data, y = y_factor, predI = 1, orthoI = 2, crossvalI = 5),
    silent = TRUE
    )
  
  # 检查模型是否成功建立
  if(class(opls_model) == "try-error") {
    print(paste("Failed to build model for:", key_name))
    next # 跳过当前的迭代
  }
  
  # 获取VIP值，如果模型失败则赋予 0
  vip_values <- tryCatch({
    opls_model@vipVn
  }, error = function(e) {
    rep(0, ncol(current_expression_data))
  })
  
  # 计算FoldChange
  baseMeanA <- rowMeans(current_expression_data[, groupA_cols])
  baseMeanB <- rowMeans(current_expression_data[, groupB_cols])
  FoldChange <- ifelse(baseMeanB > 0, abs(baseMeanA / baseMeanB), 0)
  
  # 定义一个函数进行成对t检验
  paired_t_test <- function(x, y) {
    # 如果x中所有值都相同
    if(length(unique(x)) == 1) {
      x[1] <- x[1] + 0.000000001
    }
    
    # 如果y中所有值都相同
    if(length(unique(y)) == 1) {
      y[1] <- y[1] + 0.000000001
    }
    
    test_result <- t.test(x, y)
    return(test_result$p.value)
  }
  
  # 使用apply函数
  p_values <- apply(current_expression_data, 1, function(row) paired_t_test(row[groupA_cols], row[groupB_cols]))
  #padj <- p.adjust(p_values, "BH")
  # 将 VIP 值和 FoldChange 保存到文件
  current_expression_data_def <- cbind(current_expression_data, baseMeanA, baseMeanB, FoldChange, pvalues = p_values)
  current_expression_data_def <- as.data.frame(current_expression_data_def)
  
  # 创建一个长度与current_expression_data行数相同的全零向量
  vip_values_adjusted <- rep(0, nrow(current_expression_data_def))
  # 更新vip_values_adjusted向量中相应的位置
  # 这假设vip_values的名称与current_expression_data的行名对应
  if (!is.null(names(vip_values))) {
    matching_rows <- match(names(vip_values), rownames(current_expression_data_def))
    vip_values_adjusted[matching_rows] <- vip_values
  }
  current_expression_data_def <- cbind(current_expression_data_def, VIP = vip_values_adjusted)
  
  current_expression_data_def <- current_expression_data_def[order(current_expression_data_def$pvalues, na.last = TRUE), ]
  current_expression_data_def$Metabolite <- rownames(current_expression_data_def)
  current_expression_data_def <- current_expression_data_def[, c("Metabolite", setdiff(names(current_expression_data_def), "Metabolite"))]
  current_expression_data_def$padj <- p.adjust(current_expression_data_def$pvalues, "BH")
  
  if (exists("definition_df")) {
    current_expression_data_def<-merge(current_expression_data_def, definition_df, by = "Metabolite", all.x=TRUE)
  }
  
  current_expression_data_def <- current_expression_data_def[order(-current_expression_data_def$VIP, na.last = TRUE), ]
  current_expression_data_def[is.na(current_expression_data_def)] <- 0
  
  # 将更新后的数据框保存为文本文件
  write.table(current_expression_data_def, file = paste0('组间分析/', key_name, '/', key_name, '_VIP.txt'), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  # class 分组计数
  class_count <- ""
  if (exists("definition_df")) {
    greater_than_one_data_def <- current_expression_data_def[current_expression_data_def$VIP > 1,]
    class_count <- aggregate(greater_than_one_data_def$Metabolite, by = list(greater_than_one_data_def$class), length)
    class_count <- class_count[class_count$Group.1 != '',]
    class_count <- class_count[order(-class_count$x),]
    write.table(class_count, file=paste0('组间分析/',key_name,'/Significant_compound_count_by_class.txt'),
                sep="\t", row.names=FALSE,col.names=FALSE, quote = FALSE)
  }
  
  # 计算 deg
  deg_df <- current_expression_data_def[current_expression_data_def$VIP > 1,]
  deg_up <- nrow(deg_df[deg_df$FoldChange>=1.2,])
  deg_down <- nrow(deg_df[deg_df$FoldChange<=0.8,])
  new_row <- data.frame(group = key_name, All= deg_up + deg_down, Up = deg_up, Down = deg_down)
  deg_data <- rbind(deg_data, new_row)
  
  # 生成图形
  for (type in plot_types) {
    # 创建文件名
    file_name <- paste0('组间分析/', key_name, '/', key_name, "_OPLS_DA_", type, ".png")
    
    # 检查模型是否有效
    if(inherits(opls_model, "opls")) {
      png(file_name)
      try({
        plot(opls_model, typeVc = type)
      }, silent = TRUE)
      dev.off() # 确保关闭图形设备
    }
  }
}
write.table(deg_data, file='DEG_summary.txt', sep='\t', row.names=FALSE, col.names=TRUE, quote=FALSE)


