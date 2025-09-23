# 加载所需的R包（尽量安静，若缺失则给出明确错误）
suppressMessages({
  if (!requireNamespace("vegan", quietly = TRUE)) {
    stop("R package 'vegan' is not installed. Please install it: install.packages('vegan')")
  }
  library(vegan)
  if (!requireNamespace("optparse", quietly = TRUE)) {
    stop("R package 'optparse' is not installed. Please install it: install.packages('optparse')")
  }
  library(optparse)
})

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
if (!file.exists(opt$input)) {
  stop(paste0("Input file not found: ", opt$input))
}

asv <- tryCatch({
  read.table(
    opt$input,
    header = TRUE,
    row.names = 1,
    sep = "\t",
    check.names = FALSE,
    comment.char = "",
    quote = "",
    stringsAsFactors = FALSE
  )
}, error = function(e) {
  stop(paste0("Failed to read input table: ", conditionMessage(e)))
})

if (nrow(asv) == 0 || ncol(asv) == 0) {
  stop("Input ASV table has zero rows or zero columns.")
}

# 确保所有样本列为非负整数（rrarefy 需要计数数据）
non_numeric_cols <- names(asv)[!vapply(asv, is.numeric, logical(1))]
if (length(non_numeric_cols) > 0) {
  # 尝试强制转换
  asv[non_numeric_cols] <- lapply(asv[non_numeric_cols], function(x) suppressWarnings(as.numeric(x)))
}
if (any(!vapply(asv, is.numeric, logical(1)))) {
  stop("All sample columns must be numeric counts.")
}
if (any(is.na(asv))) {
  stop("Input ASV table contains NA values after parsing. Please clean the data.")
}
if (any(as.matrix(asv) < 0)) {
  stop("Input ASV table contains negative counts, which are invalid.")
}

# 查看每个样本的OTU总数
print("原始样本OTU总数：")
print(colSums(asv))  # 输出每个样本的OTU总和

# 进行抽平分析：抽平到最小测序深度（每列之和的最小值）
lib_sizes <- colSums(asv)
if (any(is.na(lib_sizes))) {
  stop("Column sums contain NA. Check input data.")
}
target_depth <- min(lib_sizes)
if (target_depth <= 0) {
  stop(paste0("Invalid target depth for rarefaction: ", target_depth, ". Each sample must have positive total counts."))
}
otu_Flattening <- as.data.frame(t(rrarefy(t(as.matrix(asv)), sample = target_depth)))

# 查看抽平后每个样本的OTU总数
print("抽平后样本OTU总数：")
print(colSums(otu_Flattening))  # 输出抽平后每个样本的OTU总和

# 将行名转换为数据框的一列，并命名为"ASV_ID"
otu_Flattening$ASV_ID <- rownames(otu_Flattening)
# 将ASV_ID列移到第一列
otu_Flattening <- otu_Flattening[, c("ASV_ID", setdiff(names(otu_Flattening), "ASV_ID"))]

# 将抽平后的OTU表保存到工作目录，以便后续分析
tryCatch({
  write.table(
    otu_Flattening, 
    file = opt$output,  # 使用命令行参数指定的输出文件路径
    sep = "\t",  # 使用制表符作为分隔符
    col.names = TRUE,   # 保留列名
    row.names = FALSE,    # 保留行名
    quote = FALSE  # 避免输出双引号
  )
}, error = function(e) {
  stop(paste0("Failed to write output: ", conditionMessage(e)))
})

print(paste("抽平后的数据已保存到：", opt$output))
