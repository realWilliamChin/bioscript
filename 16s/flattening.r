# 加载所需的R包
library(vegan)  # 用于微生物多样性分析
library(openxlsx)  # 用于读取和写入Excel文件
library(optparse)  # 用于解析命令行参数

# 设置命令行参数
option_list <- list(
  make_option(c("-i", "--input"), type="character", default="asv.txt",
              help="输入文件路径 [默认= %default]"),
  make_option(c("-o", "--output"), type="character", default="asv_Flattening.txt",
              help="输出文件路径 [默认= %default]")
)

# 解析命令行参数
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# 读取OTU表
asv <- read.table(
  opt$input,
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE
)

# 查看每个样本的OTU总数
print("原始样本OTU总数：")
print(colSums(asv))  # 输出每个样本的OTU总和

# 进行抽平分析
# 使用rrarefy函数对OTU表进行抽平，抽平到最小样本的OTU总数
otu_Flattening <- as.data.frame(t(rrarefy(t(asv), mean(colSums(asv)))))

# 查看抽平后每个样本的OTU总数
print("抽平后样本OTU总数：")
print(colSums(otu_Flattening))  # 输出抽平后每个样本的OTU总和

# 将行名转换为数据框的一列，并命名为"ASV_ID"
otu_Flattening$ASV_ID <- rownames(otu_Flattening)
# 将ASV_ID列移到第一列
otu_Flattening <- otu_Flattening[, c("ASV_ID", setdiff(names(otu_Flattening), "ASV_ID"))]

# 将抽平后的OTU表保存到工作目录，以便后续分析
write.table(
  otu_Flattening, 
  file = opt$output,  # 使用命令行参数指定的输出文件路径
  sep = "\t",  # 使用制表符作为分隔符
  col.names = TRUE,   # 保留列名
  row.names = FALSE,    # 保留行名
  quote = FALSE  # 避免输出双引号
)

print(paste("抽平后的数据已保存到：", opt$output))
