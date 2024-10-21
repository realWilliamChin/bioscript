library(ggplot2)
library(optparse)
library(openxlsx)

option_list <- list(
  make_option(c("-i", "--input"),
    type = "character", default = NULL,
    help = "输入画图所需文件", metavar = "character"
  ),
  make_option(c("-o", "--output"),
    type = "character", default = NULL,
    help = "输入输出图片的名字", metavar = "character"
  ),
  make_option(c("--keggclean"),
    type = "character", default = NULL,
    help = "提供 kegg 注释出来的 KEGG_clean.txt 文件", metavar = "integer"
  ),
  make_option(c("--datadir"),
    type = "character", default = "./",
    help = "提供 compareinfo up 和 down gene id 的目录，默认是当前目录", metavar = "character"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


# Check that all required arguments are provided
if (!file.exists(opt$input)) {
  print_help(opt_parser)
  stop("请输入输入文件", call. = FALSE)
} else if (!file.exists(opt$output)) {
  print_help(opt_parser)
  stop("请提供 compare_info.txt 文件", call. = FALSE)
}

# Assign the first argument to prefix
gene_go <- opt$genego
compare_info <- opt$compare
kegg_clean <- opt$keggclean
data_dir <- opt$datadir
output_dir <- opt$output_dir

if (tolower(file_ext) == "xlsx" || tolower(file_ext) == "xls") {
  # 如果文件是 Excel 文件（xlsx 或 xls），使用 readxl 包中的 read_excel 函数读取
  data <- readxl::read_excel(file_path)
} else if (tolower(file_ext) == "txt" || tolower(file_ext) == "tab") {
  # 如果文件是以 tab 分隔的文本文件，使用 read.table 函数读取
  data <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
} else {
  stop("Unsupported file format. Please provide an Excel file (xlsx or xls) or a tab-separated text file (txt or tab).")
}

plot_bubble <- function(
  dt,
  ID_col,
  num_col,
  value_color_col,
  value_count_col,
  plot_title = "Bubble Graph",
  ylab_title = ID_col
  ) {
  p <- ggplot(dt, aes(y = factor(ID_col, levels = ID_col), x = num_col, size = value_count_col, colour = value_color_col)) +
    geom_point() +  # 绘制散点图
    scale_size_continuous(range = c(2, 8))+  # 调整点的大小范围
    scale_y_discrete(limits = dt$ID_col) +
    scale_colour_gradient(low = "red", high = "green") +  # 设置颜色渐变
    ggtitle(plot_title) +  # 设置图表标题
    ylab(ylab_title) +  # 设置 y 轴标签
    theme_base()

  return p
}

p <- plot

ggsave(paste0(, "_enrich_GOMF_Bubble.png", sep = ""), p, dpi = 320, width = 20, height = 10)

