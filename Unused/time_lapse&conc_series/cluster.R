library(Mfuzz)

library(factoextra)
library(cluster)
setwd("D:/chenshen/R_work/cluster")

#读取矩阵表格，在我网盘中，示例数据为“mmc2.union_all_protein_exp.txt”
#该示例中，行为基因或蛋白名称，列为时间样本（按时间顺序提前排列好，若存在生物学重复需提前取均值）
protein <- read.delim('Relief_Group_1.txt', row.names = 1, check.names = FALSE)
#protein<-read_xlsx("gene_id_germinated.xlsx",1)
protein <- as.matrix(protein)


#load data
df <- protein

#remove rows with missing values
df <- na.omit(df)

#scale each variable to have a mean of 0 and sd of 1
df <- scale(df)


#view first six rows of dataset
head(df)


fviz_nbclust(df, kmeans, method = "wss")

#使用 Bioconductor 安装 Mfuzz 包
#install.packages('BiocManager')
#BiocManager::install('Mfuzz')

#加载 Mfuzz 包
library(Mfuzz)


#构建对象
mfuzz_class <- new('ExpressionSet',exprs = protein)

#预处理缺失值或者异常值
mfuzz_class <- filter.NA(mfuzz_class, thres = 0.25)
mfuzz_class <- fill.NA(mfuzz_class, mode = 'mean')
mfuzz_class <- filter.std(mfuzz_class, min.std = 0)

#标准化数据
mfuzz_class <- standardise(mfuzz_class)

#Mfuzz 基于 fuzzy c-means 的算法进行聚类，详情 ?mfuzz
#需手动定义目标聚类群的个数，例如这里我们为了重现原作者的结果，设定为 10，即期望获得 10 组聚类群
#需要设定随机数种子，以避免再次运行时获得不同的结果
set.seed(123)
cluster_num <- 7
mfuzz_cluster <- mfuzz(mfuzz_class, c = cluster_num, m = mestimate(mfuzz_class))

#作图，详情 ?mfuzz.plot2
#time.labels 参数设置时间轴，需要和原基因表达数据集中的列对应
#颜色、线宽、坐标轴、字体等细节也可以添加其他参数调整，此处略，详见函数帮助
mfuzz.plot2(mfuzz_class, cl = mfuzz_cluster, mfrow = c(4,2 ), time.labels = colnames(protein), las = 2)
#查看每个聚类群中各自包含的蛋白数量
cluster_size <- mfuzz_cluster$size
names(cluster_size) <- 1:cluster_num
cluster_size

#查看每个蛋白所属的聚类群
head(mfuzz_cluster$cluster)

#Mfuzz 通过计算一个叫 membership 的统计量判断蛋白质所属的聚类群，以最大的 membership 值为准
#查看各蛋白的 membership 值
head(mfuzz_cluster$membership)
protein_cluster <- mfuzz_cluster$cluster
protein_cluster <- cbind(protein[names(protein_cluster), ], protein_cluster)
head(protein_cluster)
write.table(protein_cluster, 'Relief_Group_1_culster.txt', sep = '\t', col.names = NA, quote = FALSE)

#如果您想提取数据分析过程中，标准化后的表达值（绘制曲线图用的那个值，而非原始蛋白表达值）
protein_cluster <- mfuzz_cluster$cluster
protein_standard <- mfuzz_class@assayData$exprs
protein_standard_cluster <- cbind(protein_standard[names(protein_cluster), ], protein_cluster)
head(protein_standard_cluster)
#write.table(protein_standard_cluster, 'protein_standard_cluster.txt', sep = '\t', col.names = NA, quote = FALSE)
