library(pheatmap)
library(openxlsx)
library(grid)
library(rlang)

#' 精确计算热图所需尺寸的高级函数
#'
#' @param matrix_data 数据矩阵
#' @param cellwidth 单个单元格宽度（点）
#' @param cellheight 单个单元格高度（点）
#' @param fontsize_row 行标签字体大小
#' @param fontsize_col 列标签字体大小
#' @param annotation_row 行注释数据框
#' @param annotation_col 列注释数据框
#' @param show_rownames 是否显示行名
#' @param show_colnames 是否显示列名
#' @param treeheight_row 行聚类树高度
#' @param treeheight_col 列聚类树高度
#' @param annotation_names_row 是否显示行注释名称
#' @param annotation_names_col 是否显示列注释名称
#' @param legend 是否显示图例
#' @param margin_left 左侧额外边距（英寸）
#' @param margin_right 右侧额外边距（英寸）
#' @param margin_top 顶部额外边距（英寸）
#' @param margin_bottom 底部额外边距（英寸）
#' @param silent 是否静默运行（不绘制图形）
#' @param ... 传递给pheatmap的其他参数
#'
#' @return 包含宽度和高度的列表（单位：英寸）
calculate_heatmap_dimensions <- function(matrix_data,
                                          cellwidth = 12,
                                          cellheight = 12,
                                          fontsize_row = 10,
                                          fontsize_col = 10,
                                          annotation_row = NA,
                                          annotation_col = NA,
                                          show_rownames = TRUE,
                                          show_colnames = TRUE,
                                          treeheight_row = 30,
                                          treeheight_col = 30,
                                          annotation_names_row = TRUE,
                                          annotation_names_col = TRUE,
                                          legend = TRUE,
                                          margin_left = 0.5,
                                          margin_right = 0.5,
                                          margin_top = 0.5,
                                          margin_bottom = 0.5,
                                          silent = TRUE,
                                          ...) {
  
  # 在独立的临时设备上绘制临时热图以获取布局信息，避免污染当前设备
  previous_device <- grDevices::dev.cur()
  tmp_file <- tempfile(fileext = ".png")
  grDevices::png(filename = tmp_file, width = 10, height = 10, units = "in", res = 10)
  on.exit({
    try(grDevices::dev.off(), silent = TRUE)
    # 恢复之前的设备（如果存在）
    if (!is.null(previous_device) && previous_device > 1) {
      try(grDevices::dev.set(previous_device), silent = TRUE)
    }
    # 清理临时文件
    try(unlink(tmp_file), silent = TRUE)
  }, add = TRUE)

  temp_heatmap <- pheatmap(matrix_data,
                          cellwidth = cellwidth,
                          cellheight = cellheight,
                          fontsize_row = fontsize_row,
                          fontsize_col = fontsize_col,
                          annotation_row = annotation_row,
                          annotation_col = annotation_col,
                          show_rownames = show_rownames,
                          show_colnames = show_colnames,
                          treeheight_row = treeheight_row,
                          treeheight_col = treeheight_col,
                          annotation_names_row = annotation_names_row,
                          annotation_names_col = annotation_names_col,
                          legend = legend,
                          silent = TRUE,
                          ...)
  
  # 获取gtable对象的尺寸
  gtable_widths <- temp_heatmap$gtable$widths
  gtable_heights <- temp_heatmap$gtable$heights
  
  # 转换为英寸
  total_width_inch <- convertUnit(sum(gtable_widths), "inch", valueOnly = TRUE)
  total_height_inch <- convertUnit(sum(gtable_heights), "inch", valueOnly = TRUE)
  
  # 添加额外边距
  total_width_inch <- total_width_inch + margin_left + margin_right
  total_height_inch <- total_height_inch + margin_top + margin_bottom
  
  return(list(width = total_width_inch, height = total_height_inch))
}

#' 自动设置热图参数函数
#'
#' @param matrix_data 数据矩阵
#' @param autoset_image_specification 是否使用自动规格设置
#' @param ... 用户提供的其他参数
#'
#' @return 包含所有参数的列表
auto_set_heatmap_parameters <- function(matrix_data, autoset_image_specification = TRUE, ...) {
  # 获取用户提供的参数
  user_params <- list(...)
  
  # 如果不需要自动设置，直接返回用户参数
  if (!autoset_image_specification) {
    return(user_params)
  }
  
  # 参数设置表
  n_values <- c(2, 4, 8, 16, 64, 128, 256, 512, 2000, 3000, 5000, 10000, 65535)
  cellwidth_values <- c(120, 60, 30, 20, 15, 7.5, 3.3, 1.2, 0.6, 0.35, 0.1, 0.1, 0.1)
  fontsize_values <- c(10, 10, 10, 10, 10, 8, 4, 1, 1, 1, 1, 1, 1)
  dpi_values <- c(300, 300, 300, 300, 300, 300, 512, 1024, 1024, 1024, 1024, 1024, 1024)
  
  # 获取数据矩阵的维度
  nrow_data <- nrow(matrix_data)
  ncol_data <- ncol(matrix_data)
  
  # 根据行数选择参数
  row_idx <- findInterval(nrow_data, n_values, rightmost.closed = TRUE)
  row_idx <- min(row_idx, length(n_values))
  
  # 根据列数选择参数
  col_idx <- findInterval(ncol_data, n_values, rightmost.closed = TRUE)
  col_idx <- min(col_idx, length(n_values))
  
  # 自动设置的参数
  auto_params <- list(
    cellwidth = cellwidth_values[col_idx],
    cellheight = cellwidth_values[row_idx],
    fontsize_row = fontsize_values[row_idx],
    fontsize_col = fontsize_values[col_idx],
    dpi = max(dpi_values[row_idx], dpi_values[col_idx]),
    color = colorRampPalette(c("navy", "white", "firebrick3"))(100)
  )
  
  # 当 n_values 超过 3000 时不添加 names
  if (nrow_data > 3000) {
    auto_params$show_rownames <- FALSE
  }
  if (ncol_data > 3000) {
    auto_params$show_colnames <- FALSE
  }
  
  # 合并参数：用户参数优先，自动参数作为默认值
  final_params <- user_params
  for (param_name in names(auto_params)) {
    if (!param_name %in% names(final_params)) {
      final_params[[param_name]] <- auto_params[[param_name]]
    }
  }
  
  return(final_params)
}

#' 自动检测文件后缀函数
#' @param filename 文件名
#' @return 检测到的文件格式（"pdf", "png", "tiff", "svg", "jpg", "jpeg"）
detect_file_extension <- function(filename) {
  if (is.null(filename) || filename == "") {
    return(NULL)
  }
  
  # 提取文件扩展名
  extension <- tolower(tools::file_ext(filename))
  
  # 支持的格式映射
  format_mapping <- list(
    "pdf" = "pdf",
    "png" = "png", 
    "tiff" = "tiff",
    "tif" = "tiff",
    "svg" = "svg",
    "jpg" = "jpg",
    "jpeg" = "jpg"
  )
  
  # 检查是否支持该格式
  if (extension %in% names(format_mapping)) {
    return(format_mapping[[extension]])
  } else {
    # 默认返回png格式
    warning(sprintf("不支持的文件格式 '%s'，将使用PNG格式", extension))
    return("png")
  }
}

#' 智能热图生成函数
#'
#' @param matrix_data 数据矩阵
#' @param filename 输出文件名（可选）
#' @param output_format 输出格式（"pdf", "png", "tiff", "svg", "jpg", "jpeg"），如果为NULL则自动检测
#' @param autoset_image_specification 是否使用自动规格设置
#' @param ... 传递给calculate_heatmap_dimensions_advanced和pheatmap的参数
#'
#' @return 生成的热图对象
smart_heatmap <- function(matrix_data, 
                         filename = NULL,
                         output_format = NULL,
                         autoset_image_specification = TRUE,
                         ...) {
  
  # 使用自动参数设置函数
  params <- auto_set_heatmap_parameters(matrix_data, autoset_image_specification, ...)
  
  # 获取dpi参数（用于输出设备设置）
  dpi <- if ("dpi" %in% names(params)) params$dpi else 300

  # 自动检测输出格式（如果未指定）
  if (is.null(output_format) && !is.null(filename)) {
    output_format <- detect_file_extension(filename)
  }

  # 计算所需尺寸
  dimensions <- do.call(calculate_heatmap_dimensions, c(list(matrix_data), params))
  
  # 针对位图设备限制进行像素尺寸自适应（避免超过设备最大尺寸）
  max_pixels_per_side <- 65535
  adjust_dpi_if_needed <- function(cur_dpi, width_in, height_in) {
    width_px <- width_in * cur_dpi
    height_px <- height_in * cur_dpi
    if (width_px > max_pixels_per_side || height_px > max_pixels_per_side) {
      scale_factor <- min(max_pixels_per_side / width_px, max_pixels_per_side / height_px)
      # 至少保持 72 DPI，防止过低导致设备报错或图形过于模糊
      new_dpi <- max(72, floor(cur_dpi * scale_factor))
      return(new_dpi)
    }
    return(cur_dpi)
  }
  
  # 如果指定了文件名，设置输出设备
  if (!is.null(filename)) {
    # 确保设备能被安全关闭
    device_opened <- FALSE
    if (output_format == "pdf") {
      pdf(file = filename, 
          width = dimensions$width, 
          height = dimensions$height)
      device_opened <- TRUE
    } else if (output_format == "png") {
      dpi <- adjust_dpi_if_needed(dpi, dimensions$width, dimensions$height)
      png(filename = filename,
          width = dimensions$width * dpi, 
          height = dimensions$height * dpi,
          res = dpi)
      device_opened <- TRUE
    } else if (output_format == "tiff") {
      dpi <- adjust_dpi_if_needed(dpi, dimensions$width, dimensions$height)
      tiff(filename = filename,
           width = dimensions$width * dpi, 
           height = dimensions$height * dpi,
           res = dpi)
      device_opened <- TRUE
    } else if (output_format == "svg") {
      svg(filename = filename,
          width = dimensions$width, 
          height = dimensions$height)
      device_opened <- TRUE
    } else if (output_format == "jpg" || output_format == "jpeg") {
      dpi <- adjust_dpi_if_needed(dpi, dimensions$width, dimensions$height)
      jpeg(filename = filename,
           width = dimensions$width * dpi, 
           height = dimensions$height * dpi,
           res = dpi)
      device_opened <- TRUE
    } else {
      stop("不支持的输出格式。请选择: 'pdf', 'png', 'tiff', 'svg', 'jpg', 或 'jpeg'")
    }
    
    # 绘制热图并保证关闭设备
    on.exit({ if (isTRUE(device_opened)) try(dev.off(), silent = TRUE) }, add = TRUE)
    heatmap_result <- do.call(pheatmap, c(list(matrix_data), params))
    
    message(sprintf("热图已保存到: %s (尺寸: %.2f × %.2f 英寸)", 
                   filename, dimensions$width, dimensions$height))
  } else {
    # 如果没有指定文件名，直接在设备上绘制
    heatmap_result <- do.call(pheatmap, c(list(matrix_data), params))
  }
  
  return(invisible(heatmap_result))
}

#' 估算最长标签所需的宽度
#'
#' @param labels 标签向量
#' @param fontsize 字体大小
#'
#' @return 估算的宽度（英寸）
estimate_label_width <- function(labels, fontsize = 10) {
  if (length(labels) == 0) return(0)
  
  # 找到最长的标签
  max_label <- labels[which.max(nchar(labels))]
  
  # 估算宽度（假设平均字符宽度为字体大小的0.6倍）
  # 这个系数可以根据实际字体调整
  estimated_width <- (nchar(max_label) * fontsize * 0.6) / 72 # 转换为英寸
  
  return(estimated_width)
}

#' 估算最长标签所需的高度
#'
#' @param labels 标签向量
#' @param fontsize 字体大小
#'
#' @return 估算的高度（英寸）
estimate_label_height <- function(labels, fontsize = 10) {
  if (length(labels) == 0) return(0)
  
  # 高度主要取决于字体大小和行数
  estimated_height <- fontsize / 72 # 转换为英寸
  
  return(estimated_height)
}


# ============= 示例数据 ================
# set.seed(2398)
# n_rows <- 3
# n_cols <- 1024

# data_matrix <- matrix(rnorm(n_rows * n_cols), nrow = n_rows, ncol = n_cols)
# rownames(data_matrix) <- paste0("Gene_", 1:n_rows, "_with_very_long_name_for_testing")
# colnames(data_matrix) <- paste0("Sample_", 1:n_cols, "_condition_timepoint")

# # 创建行注释
# row_annotation <- data.frame(
#   Pathway = rep(c("Apoptosis", "Cell Cycle", "Metabolism"), length.out = n_rows),
#   Expression = rnorm(n_rows)
# )
# rownames(row_annotation) <- rownames(data_matrix)

# # 创建列注释
# col_annotation <- data.frame(
#   Treatment = rep(c("Control", "Drug_A", "Drug_B"), length.out = n_cols),
#   Timepoint = rep(1:3, length.out = n_cols)
# )
# rownames(col_annotation) <- colnames(data_matrix)
# =======================================

#nrow or ncol 2   4   8   16  64  128 256 512
#cellwidth    120 60  30  20  15  7.5 3.3 1.2
#fontsize     10  10  10  10  10  8   4   1
#dpi          300 300 300 300 300 300 512 1024
# 实际的热图生成调用
# smart_heatmap(
#   matrix_data = reads_data,
#   filename = output_pic,
#   output_format = "png",
#   # 以下参数将根据数据矩阵的行数和列数自动设置：
#   # cellwidth, cellheight, fontsize_row, fontsize_col, dpi, color
#   show_rownames = TRUE,
#   show_colnames = TRUE,
#   main = "Heatmap"
# )

# 示例：如果需要自定义颜色，可以这样调用：
# smart_heatmap(
#   matrix_data = data_matrix,
#   filename = "advanced_heatmap.svg",
#   cellwidth = 1.2,
#   cellheight = 120,
#   fontsize_col = 1,
#   fontsize_row = 10,
#   annotation_row = row_annotation,
#   annotation_col = col_annotation,
#   show_rownames = TRUE,
#   show_colnames = TRUE,
#   #treeheight_row = 20,
#   #treeheight_col = 20,
#   main = "Class:",
#   dpi = 2048
# )