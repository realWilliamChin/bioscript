#!/usr/bin/env Rscript
# CellChat 风格气泡图 (Bubble Plot)
# 用法: Rscript Plot/Bubble/cellchat_bubble.R -i input.txt -o output_prefix
#
# 输入文件必需列（可通过参数自定义列名）:
#   interaction_name : 配体-受体对 (y轴)
#   source           : 发送细胞
#   target           : 接收细胞
#   prob             : 通讯概率 (映射大小)
#   pval             : p-value (映射颜色)
#
# 说明: 默认按标准 CellChat 风格映射 prob=大小, pval=颜色。
#   颜色/大小区间自动按输入数据的最小最大值设置。

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(openxlsx))

option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "输入数据文件路径 (txt/tsv/xlsx)", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = "cellchat_bubble",
              help = "输出文件前缀 [default: %default]", metavar = "character"),
  make_option(c("--source_col"), type = "character", default = "source",
              help = "发送细胞(source)列名 [default: %default]", metavar = "character"),
  make_option(c("--target_col"), type = "character", default = "target",
              help = "接收细胞(target)列名 [default: %default]", metavar = "character"),
  make_option(c("--lr_col"), type = "character", default = "interaction_name",
              help = "配体-受体对列名 (y轴) [default: %default]", metavar = "character"),
  make_option(c("--prob_col"), type = "character", default = "prob",
              help = "通讯概率列名 (映射大小) [default: %default]", metavar = "character"),
  make_option(c("--pval_col"), type = "character", default = "pval",
              help = "p-value列名 (映射颜色) [default: %default]", metavar = "character"),
  make_option(c("--sep"), type = "character", default = "->",
              help = "source与target合并为x轴的分隔符 [default: %default]", metavar = "character"),
  make_option(c("--facet_by"), type = "character", default = NULL,
              help = "按此列分面，模拟顶部source分组条效果 (如 source) [default: %default]", metavar = "character"),
  make_option(c("--width"), type = "numeric", default = 16,
              help = "图片宽度 (inch) [default: %default]", metavar = "numeric"),
  make_option(c("--height"), type = "numeric", default = 10,
              help = "图片高度 (inch) [default: %default]", metavar = "numeric"),
  make_option(c("--dpi"), type = "numeric", default = 320,
              help = "PNG分辨率 [default: %default]", metavar = "numeric"),
  make_option(c("--title"), type = "character", default = NULL,
              help = "图标题 [default: %default]", metavar = "character"),
  make_option(c("--font_family"), type = "character", default = "Helvetica",
              help = "字体 (Helvetica/Arial/Times) [default: %default]", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("必须提供输入文件 (-i)", call. = FALSE)
}

if (!file.exists(opt$input)) {
  stop(paste0("输入文件不存在: ", opt$input), call. = FALSE)
}

# ---- 读取数据 ----
ext <- tolower(tools::file_ext(opt$input))
if (ext %in% c("xlsx", "xls")) {
  df <- read.xlsx(opt$input, sheet = 1)
} else {
  df <- read.table(opt$input, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                   check.names = FALSE, quote = "", comment.char = "")
}

required_cols <- c(opt$source_col, opt$target_col, opt$lr_col, opt$prob_col, opt$pval_col)
missing_cols <- setdiff(required_cols, colnames(df))
if (length(missing_cols) > 0) {
  stop(paste0("输入文件缺少必需列: ", paste(missing_cols, collapse = ", ")), call. = FALSE)
}

# ---- 构建 x 轴标签 ----
df$x_axis <- paste0(df[[opt$source_col]], opt$sep, df[[opt$target_col]])

# 保持 y 轴顺序：数据文件中越靠上的行在图中越靠上
y_levels <- rev(unique(df[[opt$lr_col]]))
x_levels <- unique(df$x_axis)
df[[opt$lr_col]] <- factor(df[[opt$lr_col]], levels = y_levels)
df$x_axis <- factor(df$x_axis, levels = x_levels)

# ---- 按数据范围设置区间 ----
prob_range <- range(df[[opt$prob_col]], na.rm = TRUE)
pval_range <- range(df[[opt$pval_col]], na.rm = TRUE)
# pval 最大值小于 0.05 时，颜色条上限固定为 0.05
if (pval_range[2] < 0.05) {
  pval_range[2] <- 0.05
}

# ---- 绘图 ----
p <- ggplot(df, aes_string(
    x = "x_axis",
    y = opt$lr_col,
    size = opt$prob_col,
    color = opt$pval_col
  )) +
  geom_point(alpha = 0.9, stroke = 0) +
  scale_size_continuous(
    range = c(1.5, 8),
    name = opt$prob_col,
    limits = prob_range
  ) +
  scale_color_gradientn(
    colours = c("#2166AC", "#67A9CF", "#D1E5F0", "#FDDBC7", "#EF8A62", "#B2182B"),
    name = opt$pval_col,
    limits = pval_range
  ) +
  labs(x = NULL, y = NULL, title = opt$title) +
  theme_bw(base_family = opt$font_family) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8, color = "black", face = "bold"),
    axis.text.y = element_text(size = 9, color = "black", face = "bold"),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5),
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold")
  )

# ---- 分面: 模拟顶部 source 分组条 ----
facet_col <- opt$facet_by
if (!is.null(facet_col) && facet_col %in% colnames(df)) {
  df[[facet_col]] <- factor(df[[facet_col]], levels = unique(df[[facet_col]]))
  p <- p + facet_grid(as.formula(paste(". ~", facet_col)),
                      scales = "free_x", space = "free_x") +
    theme(
      strip.background = element_rect(fill = "grey85", color = "black"),
      strip.text = element_text(size = 10, face = "bold"),
      panel.spacing.x = unit(0.5, "lines")
    )
}

# ---- 保存 ----
out_png <- paste0(opt$output, ".png")
out_pdf <- paste0(opt$output, ".pdf")

ggsave(out_png, p, dpi = opt$dpi, width = opt$width, height = opt$height, limitsize = FALSE)
ggsave(out_pdf, p, width = opt$width, height = opt$height, limitsize = FALSE)

message(paste0("图片已保存: ", out_png, ", ", out_pdf))
