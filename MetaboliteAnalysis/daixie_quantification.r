setwd("L:/work/05_daixie/2023_12_12_xuelaoshi_caomei")
setwd("/home/colddata/qinqiang/work/05_daixie/2024_03_25_张译心/sepal_44_group")
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
library(openxlsx)
rm(list=ls())

# 读取文件
#reads_data<-read.table("All_sample_data.txt",sep="\t",row.names=1,header=T,check.names=F,stringsAsFactors = F)
reads_data <- read.xlsx("All_sample_data.xlsx", sheet = 1, rowNames = TRUE)
sample_info<-read.table("samples_described.txt",sep="\t",header=T,check.names=F,stringsAsFactors = F)

reads_data<-reads_data[sample_info$sample]
# 数据预处理
reads_data <- reads_data[rowSums(reads_data != 0) > 0, ]
# 清理每个元素头尾的空格
#reads_data <- apply(reads_data, 2, function(x) trimws(x, which = c("both")))
# 检查是否油 NA，油则退出
if (any(is.na(reads_data))) {
  print("检查数据，有 NA")
  quit()
}
# 如果需要合并定义则读取单独定义文件，没有则跳过
#definition_df<-read.table("def.txt",sep="\t",row.names=1,header=T,check.names=F,stringsAsFactors = F,quote="")
if (file.exists("def.xlsx")) {
  definition_df <- read.xlsx("def.xlsx", sheet = 1, rowNames = TRUE)
  definition_df$Metabolite <- rownames(definition_df)
}



# heatmap
# 化合物多于 100 的 ，show_rownames=F
if (nrow(reads_data) > 100) {
  all.heatmap<-pheatmap(reads_data,scale="row",cluster_cols=F,show_rownames=F)
} else {
  all.heatmap<-pheatmap(reads_data,scale="row",cluster_cols=F,show_rownames=T)
}
ggsave("All_metabolites_heatmap.jpeg",all.heatmap,dpi=300,width=20,height=10,limitsize=FALSE)


# correlation
reads_data<-na.omit(reads_data)
reads_data<-reads_data + 0.000000001
fpkm.m<-as.matrix(reads_data)
fpkm.cor<-cor(fpkm.m)
correlation_df<-as.data.frame(fpkm.cor)
correlation_df$ID<-rownames(correlation_df)
correlation_df <- correlation_df[, c("ID", setdiff(names(correlation_df), "ID"))]
# write.table(correlation_df,"Metabolite_correlation.txt",sep="\t",quote=F,row.names=F)
write.xlsx(correlation_df, "Metabolite_correlation.xlsx", sheetName = "Sheet1", rowNames = FALSE)
min(fpkm.cor)
# 根据样本数量，设置图的大小
sample_num <- ncol(fpkm.cor)
correlation_plot_width <- 5 + sample_num * 0.5
correlation_plot_height <- 4 + sample_num * 0.5
png("Metabolite_correlation_graph.png", width = correlation_plot_width, height = correlation_plot_height, units = 'in', res = 300)
corrplot(fpkm.cor,is.corr = F,col= rev(COL2('PiYG')),method="color",addCoef.col = 'black',tl.col="black",col.lim=c(min(fpkm.cor)-0.01,max(fpkm.cor)),cl.ratio=0.1)
dev.off()
pdf("Metabolite_correlation_graph.pdf", width = correlation_plot_width, height = correlation_plot_height)
corrplot(fpkm.cor,is.corr = F,col= rev(COL2('PiYG')),method="color",addCoef.col = 'black',tl.col="black",col.lim=c(min(fpkm.cor)-0.01,max(fpkm.cor)),cl.ratio=0.1)
dev.off()


# pca 图
gene.pca <- PCA(t(reads_data), ncp = 2, scale.unit = TRUE, graph = FALSE)
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
multigroup_dir <- "多组分析"
dir.create(multigroup_dir)
# >>>>>> 多组分析中有多个组的
#select_sample_info <- read.table("multigroup4_samples_described.txt",sep="\t",header=T,check.names=F,stringsAsFactors = F)
#select_reads <- reads_data[, select_sample_info$sample, drop = FALSE]

#reads_data_t<-t(select_reads)
#groups=select_sample_info$group
#select_reads_data <- reads_data[, select_sample_info$sample, drop = FALSE]
# <<<<<<

reads_data_t<-t(reads_data)
groups=sample_info$group

metabolites<-as.matrix(reads_data_t)

# 如果代谢物的数量小于 10 用 4，10-20 用 6，20 以上用 10
if (ncol(metabolites) < 10) {
  ncomp = 4
} else if (ncol(metabolites) < 20) {
  ncomp = 6
} else {
  ncomp = 10
}

df_plsda <- plsda(reads_data_t, groups, ncomp = ncomp)

# scree plot
scree_df <- as.data.frame(df_plsda$prop_expl_var$X)
scree_df$comp=rownames(scree_df)
colnames(scree_df)[1]="value"

df<-data.frame(x=letters[1:6],y=seq(10,60,10))
# 创建 scree plot 图
scree_df$comp <- factor(scree_df$comp, levels=scree_df$comp)
scree_plot <- ggplot(scree_df, aes(x = comp, y = value)) +
#  geom_bar(stat = "identity", fill = "steelblue") +
  labs(x = "Principal Component", y = "Explained Variance") +
  ggtitle("Scree Plot") +
  theme_minimal()+
  geom_line(aes(group = 1), color = "red")+
  geom_point()

ggsave(paste0(multigroup_dir, '/Metabolite_quantitation_scree_plot.jpeg'), scree_plot, width = ncomp * 0.7, height = 4)

comp_load_df <- as.data.frame(df_plsda$loadings$X)
comp_load_df <- cbind(rownames(comp_load_df), comp_load_df)
colnames(comp_load_df)[1] = "compound_name"
# write.table(comp_load_df, file="多组分析/pc_loading_value.txt", sep='\t', row.names=FALSE,col.names = TRUE,quote = FALSE)
write.xlsx(comp_load_df, file=paste0(multigroup_dir, "/pc_loading_value.xlsx"), sheetName = "Sheet1", rowNames = FALSE)

df <- unclass(df_plsda)

df1 = as.data.frame(df$variates$X)
# df1$group = sample_info$group
df1$group = groups
df1$samples = rownames(df1)

explain = df$prop_expl_var$X
x_label <- round(explain[1],digits=3)
y_label <- round(explain[2],digits=3)

num_samples <- nrow(df1)
plot_width <- 5 + num_samples * 0.2
plot_height <- 4 + num_samples * 0.2

p1 <- ggplot(df1, aes(x = comp1, y = comp2, color = group)) +
  theme_bw() +
  geom_point(size = 1.8) +
  geom_text_repel(
    aes(label = samples),
    size = 2,
    box.padding = 0.35,
    point.padding = 0.5
  ) +
  labs(x = paste0("P1 (", x_label * 100, "%)"), y = paste0("P2 (", y_label * 100, "%)")) +
  stat_ellipse(data = df1, geom = "polygon", level = 0.95,
               linetype = 2, linewidth = 0.5, aes(fill = group),
               alpha = 0.2, show.legend = TRUE) +
  scale_color_discrete() +
  scale_fill_discrete() +
  theme(
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12, angle = 90),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    panel.grid = element_blank()
  )
# 提取当前坐标轴的极限
current_limits <- ggplot_build(p1)$layout$panel_params[[1]]
x_range <- current_limits$x.range
y_range <- current_limits$y.range

# 根据当前值每个增加2
x_min <- x_range[1] - 2
x_max <- x_range[2] + 2
y_min <- y_range[1]# - 1
y_max <- y_range[2]# + 1
p1 <- p1 + coord_cartesian(xlim = c(x_min, x_max), ylim = c(y_min, y_max))

ggsave(paste0(multigroup_dir, "/Multigroup_Plsda_Distribution_Graph.jpeg"), p1, width = plot_width, height = plot_height)



# vip 值使用 fpkm, 其他值计算都是用 reads
plsda_model <- opls(x = metabolites, y = groups, predI = 1)
vip_values <- plsda_model@vipVn
df.vip<-as.data.frame(vip_values)

# 整个多组
reads_data_with_def <- cbind(reads_data, VIP = df.vip$vip_values)
# 多个多组
# reads_data_with_def <- cbind(select_reads_data, VIP = df.vip$vip_values)

reads_data_with_def$Metabolite <- rownames(reads_data_with_def)
reads_data_with_def <- reads_data_with_def[, c("Metabolite", setdiff(names(reads_data_with_def), "Metabolite"))]

# 没有定义跳过
class_count <- ""
if (exists("definition_df")) {
  reads_data_with_def<-merge(reads_data_with_def, definition_df, by = "Metabolite", all.x=TRUE)
  greater_than_one_data_def <- reads_data_with_def[reads_data_with_def$VIP > 1,]
  # 如果 greater_than_one_data_def 没有数据，就是用 reads_data_with_def
  if (nrow(greater_than_one_data_def) == 0) {
    print("没有 VIP 大于 1，输出全部")
    greater_than_one_data_def <- reads_data_with_def
  }

  class_count <- aggregate(greater_than_one_data_def$Metabolite, by = list(greater_than_one_data_def$Class), length)
  class_count <- class_count[class_count$Group.1 != '',]
  class_count <- class_count[order(-class_count$x),]
  # write.table(class_count, file='多组分析/Significant_compound_count_by_class.txt',
  #             sep="\t", row.names=FALSE,col.names=FALSE, quote = FALSE)
  write.xlsx(class_count, file=paste0(multigroup_dir, "/Significant_compound_count_by_class.xlsx"),
             sheetName = "Sheet1", rowNames=FALSE,colNames=FALSE)
}

reads_data_with_def <- reads_data_with_def[order(-reads_data_with_def$VIP, na.last = TRUE), ]

#write.table(reads_data_with_def, file="多组分析/Metabolite_quantitation_VIP.txt", sep='\t', row.names = FALSE, col.names = TRUE, quote=FALSE)
write.xlsx(reads_data_with_def, file=paste0(multigroup_dir, "/Metabolite_quantitation_VIP.xlsx"), sheetName = "Sheet1", rowNames = FALSE, colNames = TRUE)


# 自动生成组间分析
# group_levels <- unique(sample_info$group)
# comparisons <- combn(group_levels, 2, simplify = FALSE)

# 指定组间分析
dir.create('组间分析')
comp_info<-read.table("compare_info.txt",sep="\t",header=T,check.names=F,stringsAsFactors = F)
comparisons <- list()
for(i in seq_along(1:nrow(comp_info))){
  comparisons <- append(comparisons, list(as.character(comp_info[i,])))
}

plot_types <- c('correlation', 'outlier', 'overview', 'permutation', 
                'predict-train', 'x-loading', 'x-score', 
                'x-variance', 'xy-score')

deg_data <- data.frame(group = character(0), All=numeric(0), Up = numeric(0), Down = numeric(0),
                       Up_list = character(0), Down_list = character(0))
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
  # current_expression_fpkm_data <- fpkm[, current_samples, drop = FALSE]
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
  transposed_expression_data <- t(current_expression_data)
  
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
  # write.table(current_expression_data_def, file = paste0('组间分析/', key_name, '/', key_name, '_VIP.txt'), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  write.xlsx(current_expression_data_def, file=paste0('组间分析/', key_name, '/', key_name, '_VIP.xlsx'), sheetName = "Sheet1", rowNames=FALSE, colNames=TRUE)
  
  # 计算 deg
  deg_df <- current_expression_data_def[current_expression_data_def$VIP > 1,]
  deg_up <- nrow(deg_df[deg_df$FoldChange>=1.2,])
  deg_up_idlist <- deg_df[deg_df$FoldChange>=1.2,]$Metabolite
  deg_down <- nrow(deg_df[deg_df$FoldChange<=0.8,])
  deg_down_idlist <- deg_df[deg_df$FoldChange<=0.8,]$Metabolite
  new_row <- data.frame(group = key_name, All= deg_up + deg_down, Up = deg_up, Down = deg_down, 
                        Up_list = paste0(deg_up_idlist, collapse = ","), Down_list = paste0(deg_down_idlist, collapse = ","))
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
# write.table(deg_data, file='组间分析/Differential_metabolite_count_summary.txt', sep='\t', row.names=FALSE, col.names=TRUE, quote=FALSE)
write.xlsx(deg_data, file='组间分析/Differential_metabolite_count_summary.xlsx', sheetName = "Sheet1", rowNames=FALSE, colNames=TRUE)


# class 分组计数
all_files <- list.files(path = "组间分析", recursive = TRUE)
vip_files <- all_files[grep("VIP", all_files)]  
class_count_list <- list()  

for (vip_file in vip_files) {
  print(vip_file)
  vip_df <- read.xlsx(paste0("组间分析/",vip_file), sheet = 1, rowNames = FALSE)
  vip_df_vipgt1 <- vip_df[vip_df$VIP > 1 & (vip_df$FoldChange > 1.2 | vip_df$FoldChange < 0.8),]
  # vip_df_vipgt1 <- vip_df
  class_count <- aggregate(Metabolite ~ Class, data=vip_df_vipgt1, FUN=length)
  class_count <- class_count[class_count$Class != '',]
  names(class_count)[2] <- strsplit(vip_file, "/")[[1]][1]  # 修改列名为对应的 count_name
  class_count_list <- append(class_count_list, list(class_count))  
}

# 去除重复的行?
#class_count_list <- lapply(class_count_list, function(df) df[!duplicated(df$class), ])

# 合并 class_count_list 中所有的 class_count，根据第一列的 class 合并，合并方式为并集  
class_count_result <- Reduce(function(x, y) merge(x, y, by = "Class", all = TRUE), class_count_list)  
class_count_result[is.na(class_count_result)] <- 0
# write.table(class_count_result, file=paste0('组间分析/Significant_compound_count_by_class.txt'),
#             sep="\t", row.names=FALSE,col.names=TRUE, quote = FALSE)
write.xlsx(class_count_result, file=paste0('组间分析/Significant_compound_count_by_class.xlsx'),
           sheetName = "Sheet1", rowNames=FALSE,colNames=TRUE)
