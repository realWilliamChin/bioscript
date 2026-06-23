#!/usr/bin/env Rscript

# =============================================================================
# 宏基因组分类学堆积条形图绘制脚本
# 
# 功能：为单个分类学水平绘制堆积条形图
# 输入：单个分类学水平的丰度数据文件（格式：第一列为分类单元，后续列为样本丰度）
# 输出：高质量JPEG格式的堆积条形图
# 
# 作者：基于原始代码优化
# 日期：2024
# =============================================================================

# 加载必要的库
suppressPackageStartupMessages({
  library(ggplot2)    # 绘图核心库
  library(reshape2)   # 数据重塑
  library(argparse)   # 命令行参数解析
  library(rlang)      # 提供 sym() 用于 aes 中按变量名引用列
})

# 创建参数解析器
parser <- ArgumentParser(description = "绘制宏基因组分类学堆积条形图")
parser$add_argument("--input_file", type = "character", required = TRUE,
                   help = "输入数据文件路径")
parser$add_argument("--output_file", type = "character", required = TRUE,
                   help = "输出图片文件路径")
parser$add_argument("--taxa_level", type = "character", required = TRUE,
                   help = "要绘制的分类学水平 (例如: Phylum, Class)")
parser$add_argument("--dpi", type = "integer", default = 320,
                   help = "图片分辨率")

# 解析命令行参数
args <- parser$parse_args()

# 定义颜色调色板（20种颜色，确保足够区分不同分类单元）
color_palette <- c("#FF9999", "#66B2FF", "#99FF99", "#FFCC99", "#FF66CC",
                   "#CC99FF", "#FFFF66", "#CCE5FF", "#FF6666", "#FFCC66",
                   "#66FFCC", "#CC99CC", "#FF9966", "#99CCFF", "#CCCCFF",
                   "#B3E6B3", "#FFB366", "#CC99FF", "#99FFCC", "#FF99CC")

# 定义分类学水平对应的正则表达式模式（用于清理分类名称）
taxa_patterns <- list(
  "Phylum" = list(pattern = "k__Bacteria\\|p__", replacement = "p__"),
  "Class" = list(pattern = "k__Bacteria.*\\|c__", replacement = "c__"),
  "Order" = list(pattern = "k__Bacteria.*\\|o__", replacement = "o__"),
  "Family" = list(pattern = "k__Bacteria.*\\|f__", replacement = "f__"),
  "Genus" = list(pattern = "k__Bacteria.*\\|g__", replacement = "g__"),
  "Species" = list(pattern = "k__Bacteria.*\\|s__", replacement = "s__")
)

# =============================================================================
# 计算自适应图片尺寸的函数
# 根据样本数量和分类单元数量自动调整图片尺寸
# =============================================================================
calculate_adaptive_size <- function(n_samples, n_taxa) {
  # 基础尺寸
  base_width <- 8
  base_height <- 6
  
  # 根据样本数量调整宽度（每个样本增加0.3英寸）
  width <- base_width + (n_samples * 0.3)
  
  # 根据分类单元数量调整高度（每个分类单元增加0.2英寸）
  height <- base_height + (n_taxa * 0.2)
  
  # 设置最小和最大限制
  width <- max(10, min(width, 40))   # 宽度范围：10-40英寸
  height <- max(8, min(height, 30))  # 高度范围：8-30英寸
  
  return(list(width = width, height = height))
}

# =============================================================================
# 绘制堆积条形图的主函数
# 
# 参数：
#   data_file: 输入数据文件路径
#   taxa_level: 分类学水平名称
#   output_file: 输出图片文件路径
#   args: 命令行参数对象
# =============================================================================
plot_stacked_bar <- function(data_file, taxa_level, output_file, args) {
  cat("正在处理", taxa_level, "水平的数据...\n")
  
  # 检查文件是否存在
  if (!file.exists(data_file)) {
    cat("警告: 文件", data_file, "不存在，跳过...\n")
    return(FALSE)
  }
  
  # 读取数据文件
  taxa_data <- read.table(
    data_file,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    sep = "\t"
  )
  
  # 设置第一列名称为分类学水平
  colnames(taxa_data)[1] <- taxa_level
  
  # 应用正则表达式清理分类名称（去除冗余的前缀）
  if (taxa_level %in% names(taxa_patterns)) {
    pattern_info <- taxa_patterns[[taxa_level]]
    taxa_data[[taxa_level]] <- gsub(pattern_info$pattern, pattern_info$replacement, 
                                   taxa_data[[taxa_level]], perl = TRUE)
  }
  
  # 设置因子水平（确保图例顺序与数据顺序一致）
  taxa_data[[taxa_level]] <- factor(taxa_data[[taxa_level]], 
                                   levels = rev(taxa_data[[taxa_level]]))
  
  # 数据重塑（从宽格式转换为长格式，便于ggplot绘图）
  taxa_data_melt <- melt(taxa_data)
  taxa_data_melt <- na.omit(taxa_data_melt)
  
  # 计算自适应图片尺寸
  n_samples <- length(unique(taxa_data_melt$variable))
  n_taxa <- length(unique(taxa_data_melt[[taxa_level]]))
  size_info <- calculate_adaptive_size(n_samples, n_taxa)
  
  # 计算字体大小（根据图片尺寸自适应）
  font_size <- max(10, min(18, size_info$width * 0.6))
  title_size <- max(14, min(25, size_info$width * 0.8))
  legend_size <- max(12, min(20, size_info$width * 0.7))
  
  # 创建堆积条形图
  p <- ggplot(taxa_data_melt, aes(x = variable, y = value)) +
    # 绘制堆积条形图
    geom_bar(aes(fill = !!sym(taxa_level)), stat = "identity", width = 0.3) +
    # 设置坐标轴标签和标题
    labs(x = "Samples", y = "Relative Abundance (%)", 
         title = paste(taxa_level, "Distribution")) +
    # 使用简洁的主题
    theme_minimal() +
    # 自定义主题设置
    theme(
      # X轴标签：75度倾斜，避免重叠
      axis.text.x = element_text(angle = 75, hjust = 1, vjust = 1, 
                                size = font_size),
      # Y轴标签
      axis.text.y = element_text(size = font_size),
      # 坐标轴标题
      axis.title = element_text(size = title_size),
      # 图例文字
      legend.text = element_text(size = legend_size),
      # 图例标题
      legend.title = element_text(size = legend_size),
      # 图片标题居中
      plot.title = element_text(size = title_size, hjust = 0.5)
    ) +
    # 设置颜色填充（根据分类单元数量选择颜色）
    scale_fill_manual(values = color_palette[1:length(unique(taxa_data[[taxa_level]]))])
  
  # 保存图片
  ggsave(output_file, p, dpi = args$dpi, 
         width = size_info$width, height = size_info$height)
  
  cat("已保存", output_file, 
      sprintf("(尺寸: %.1f x %.1f 英寸, 样本数: %d, 分类单元数: %d)\n", 
              size_info$width, size_info$height, n_samples, n_taxa))
  
  return(TRUE)
}

# =============================================================================
# 主函数：处理所有分类学水平
# =============================================================================
main <- function() {
  cat("开始绘制宏基因组分类学堆积条形图...\n")
  cat("输入文件:", args$input_file, "\n")
  cat("输出文件:", args$output_file, "\n")
  cat("分类学水平:", args$taxa_level, "\n")

  # 创建输出目录（如果不存在）
  output_dir <- dirname(args$output_file)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # 绘制图形
  success <- plot_stacked_bar(args$input_file, args$taxa_level, args$output_file, args)
  
  # 输出总结
  cat(strrep("=", 50), "\n")
  if (success) {
    cat("处理完成！图形已保存至:", args$output_file, "\n")
  } else {
    cat("处理失败。\n")
  }
}

# 运行主函数（仅在非交互模式下执行）
if (!interactive()) {
  main()
} 