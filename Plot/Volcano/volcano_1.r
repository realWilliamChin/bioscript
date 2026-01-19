# 代谢火山图，纵坐标 VIP 值，横坐标 FC

pkgs <- c('ggthemes', 'ggplot2', 'optparse', 'ggrepel')
suppressPackageStartupMessages(invisible(lapply(pkgs, require, character.only = TRUE)))

plot_volcano <- function(
  df,
  x_col = "FoldChange",
  y_col = "padj",
  reg_col = "regulation",
  alpha = 0.05,
  bs_neg = 0.8,
  bs_pos = 1.2,
  title = NULL,
  label_col = "Annotation",
  center_zero = FALSE,
  x_left = NA_real_,
  x_right = NA_real_
) {

  missing_cols <- setdiff(c(x_col, y_col, reg_col), colnames(df))
  if (length(missing_cols) > 0) {
    stop(sprintf("数据缺少必要列：%s", paste(missing_cols, collapse = ", ")))
  }

  # # 若没有提供分组列，且存在 VIP 与 FoldChange，则按 VIP 与 FoldChange 自动判定上下调
  # if (!(reg_col %in% colnames(df)) && all(c("VIP", "FoldChange") %in% colnames(df))) {
  #   reg <- rep("stable", nrow(df))
  #   reg[df$VIP > 1 & df$FoldChange > 1.2] <- "up"
  #   reg[df$VIP > 1 & df$FoldChange < 0.8] <- "down"
  #   df[[reg_col]] <- factor(reg, levels = c("down", "stable", "up"))
  # }

  # 按横轴类型确定竖线位置（默认：FoldChange 使用 0.8/1.2；log2FoldChange 通常用 -1/1，需在调用时指定）

  p <- ggplot(data = df, aes_string(x = x_col, y = y_col, colour = reg_col)) +
    geom_point() +
    scale_color_manual(values = c("green", "grey", "red")) +
    geom_vline(xintercept = c(bs_neg, bs_pos), lty = 4, col = "black", linewidth = 0.8) +
    geom_hline(yintercept = 1, lty = 4, col = "black", linewidth = 0.8) +
    theme_base()

  # 若存在标注列且非空，添加自动避让文本标注
  if (!is.null(label_col) && nzchar(label_col) && (label_col %in% colnames(df))) {
    lab_vals <- as.character(df[[label_col]])
    keep <- !is.na(lab_vals) & nzchar(lab_vals)
    if (any(keep)) {
      p <- p + ggrepel::geom_text_repel(
        data = df[keep, , drop = FALSE],
        aes_string(label = label_col),
        max.overlaps = Inf,
        size = 3,
        show.legend = FALSE
      )
    }
  }

  # 需要以 0 为中心对称时，调整 x 轴范围
  if (isTRUE(center_zero)) {
    x_vals <- df[[x_col]]
    if (is.numeric(x_vals)) {
      max_abs <- suppressWarnings(max(abs(range(x_vals, na.rm = TRUE)), na.rm = TRUE))
      if (is.finite(max_abs) && max_abs > 0) {
        p <- p + xlim(-max_abs, max_abs)
      }
    }
  }

  # 手动指定左右轴范围时优先使用，传入 NA 表示保留自动范围
  if (is.finite(x_left) || is.finite(x_right)) {
    limits <- c(x_left, x_right)
    p <- p + coord_cartesian(xlim = limits)
  }

  if (!is.null(title) && nzchar(title)) {
    p <- p + labs(title = title)
  }
  return(p)
}


is_main <- function() {
  # 检查是否通过 source() 调用
  if (sys.nframe() == 0) {
    return(TRUE) # 直接执行
  } else if (any(grepl("source", sapply(sys.calls(), function(x) deparse(x)[1])))) {
    return(FALSE) # 被 source() 调用
  } else {
    return(TRUE) # 其他情况视为直接执行
  }
}

if (is_main()) {
  option_list <- list(
    make_option(c("-i", "--input"), type = "character", help = "输入数据文件路径（制表/逗号分隔）", metavar = "FILE"),
    make_option(c("-o", "--output"), type = "character", default = "volcano.png", help = "输出图片文件，默认 volcano.png"),
    make_option(c("--x_col"), type = "character", default = "FoldChange", help = "横轴列名，默认 FoldChange（如用 log2FoldChange 请自行设置阈值）"),
    make_option(c("--y_col"), type = "character", default = "padj", help = "纵轴列名，默认 padj"),
    make_option(c("--reg_col"), type = "character", default = "regulation", help = "分组列名，默认 regulation"),
    make_option(c("--bs_neg"), type = "double", default = 0.8, help = "左侧竖线阈值（FoldChange），默认 0.8；若 x_col=log2FoldChange 可设为 np.log10(0.8)"),
    make_option(c("--bs_pos"), type = "double", default = 1.2, help = "右侧竖线阈值（FoldChange），默认 1.2；若 x_col=log2FoldChange 可设为 np.log10(1.2)"),
    make_option(c("--width"), type = "double", default = 10, help = "图片宽度(inches)，默认 10"),
    make_option(c("--height"), type = "double", default = 10, help = "图片高度(inches)，默认 10"),
    make_option(c("--dpi"), type = "integer", default = 300, help = "图片 DPI，默认 300"),
    make_option(c("--title"), type = "character", default = "", help = "图标题，可为空"),
    make_option(c("--label_col"), type = "character", default = "", help = "标注列名，默认 Annotation；为空则不标注"),
    make_option(c("--center_zero"), type = "logical", default = TRUE, help = "x 轴 0 居中且两侧距离相等，默认 FALSE"),
    make_option(c("--x_left"), type = "double", default = NA, help = "x 轴左侧最小值，默认自动"),
    make_option(c("--x_right"), type = "double", default = NA, help = "x 轴右侧最大值，默认自动")
  )
  opt_parser <- OptionParser(option_list = option_list)
  opt <- parse_args(opt_parser)

  df <- tryCatch({
    read.table(opt$input, header = TRUE, sep = '\t', quote = "\"", check.names = FALSE, comment.char = "")
  }, error = function(e) {
    stop(sprintf("读取输入失败：%s", e$message))
  })

  p <- plot_volcano(
    df = df,
    x_col = opt$x_col,
    y_col = opt$y_col,
    reg_col = opt$reg_col,
    alpha = opt$alpha,
    bs_neg = opt$bs_neg,
    bs_pos = opt$bs_pos,
    title = opt$title,
    label_col = opt$label_col,
    center_zero = opt$center_zero,
    x_left = opt$x_left,
    x_right = opt$x_right
  )

  ggsave(filename = opt$output, plot = p, width = opt$width, height = opt$height, dpi = opt$dpi, units = "in")

}