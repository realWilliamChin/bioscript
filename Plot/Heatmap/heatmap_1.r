library(pheatmap)
library(openxlsx)
library(grid)
library(optparse)

# 全局配置，集中管理魔法数字
.HEATMAP_CONFIG <- list(
  MAX_PIXELS_PER_SIDE = 32767L,
  MAX_DIMENSION_INCHES = 200,
  MIN_DPI = 72L,
  DEFAULT_DPI = 300L,
  AUTO_PARAMS = list(
    n_values = c(2, 4, 8, 16, 64, 128, 256, 512, 2000, 3000, 5000, 10000, 65535),
    cellwidth_values = c(120, 60, 30, 20, 15, 7.5, 3.3, 1.2, 0.6, 0.35, 0.1, 0.1, 0.1),
    fontsize_values = c(10, 10, 10, 10, 10, 8, 4, 1, 1, 1, 1, 1, 1),
    dpi_values = c(300, 300, 300, 300, 300, 300, 512, 1024, 1024, 1024, 1024, 1024, 1024)
  ),
  NAME_THRESHOLD = 3000L,
  FORMAT_MAPPING = list(
    pdf = "pdf",
    png = "png",
    tiff = "tiff",
    tif = "tiff",
    svg = "svg",
    jpg = "jpg",
    jpeg = "jpg"
  )
)

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
  on.exit(
    {
      try(grDevices::dev.off(), silent = TRUE)
      # 恢复之前的设备（如果存在）
      if (!is.null(previous_device) && previous_device > 1) {
        try(grDevices::dev.set(previous_device), silent = TRUE)
      }
      # 清理临时文件
      try(unlink(tmp_file), silent = TRUE)
    },
    add = TRUE
  )

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
    ...
  )

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

#' 读取以首列为索引的表格（txt/tsv/csv 或 Excel）
#'
#' @param path 文件路径（txt/csv/tsv/xlsx/xls）
#' @param sheet Excel 工作表名或序号；非 Excel 忽略
#' @param expect_numeric 是否期望主体为数值（用于数据矩阵）
#' @return data.frame，行名来自首列
read_indexed_table <- function(path, sheet = NULL, expect_numeric = FALSE) {
  stopifnot(!is.null(path) && nzchar(path))
  ext <- tolower(tools::file_ext(path))
  df <- NULL
  if (ext %in% c("xlsx", "xls")) {
    if (is.null(sheet)) {
      df <- openxlsx::read.xlsx(path, sheet = 1, rowNames = FALSE, check.names = FALSE)
    } else {
      df <- openxlsx::read.xlsx(path, sheet = sheet, rowNames = FALSE, check.names = FALSE)
    }
  } else {
    sep <- if (ext %in% c("tsv")) "\t" else if (ext %in% c("csv")) "," else "\t"
    df <- utils::read.table(path, header = TRUE, sep = sep, quote = "\"", comment.char = "", check.names = FALSE, stringsAsFactors = FALSE)
  }
  if (ncol(df) < 2) stop(sprintf("文件 %s 列数不足，首列需为索引，后续需包含数据列", path))
  rownames(df) <- as.character(df[[1]])
  df[[1]] <- NULL
  if (expect_numeric) {
    suppressWarnings({
      for (j in seq_len(ncol(df))) df[[j]] <- as.numeric(df[[j]])
    })
    if (any(is.na(as.matrix(df)))) {
      # 允许 NA，由上游 clean_matrix_data 处理；这里只输出提示
      message(sprintf("注意：%s 中存在无法转为数值的单元，已置为 NA", basename(path)))
    }
  }
  return(df)
}

#' 解析 sheet 参数（数字/字符均可）
parse_sheet_param <- function(sheet_value) {
  if (is.null(sheet_value) || !nzchar(sheet_value)) {
    return(NULL)
  }
  suppressWarnings({
    as_num <- as.numeric(sheet_value)
    if (!is.na(as_num)) {
      return(as_num)
    }
  })
  return(sheet_value)
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

  # 获取数据矩阵的维度
  nrow_data <- nrow(matrix_data)
  ncol_data <- ncol(matrix_data)
  auto_cfg <- .HEATMAP_CONFIG$AUTO_PARAMS

  # 根据行数选择参数
  row_idx <- findInterval(nrow_data, auto_cfg$n_values, rightmost.closed = TRUE)
  row_idx <- min(row_idx, length(auto_cfg$n_values))

  # 根据列数选择参数
  col_idx <- findInterval(ncol_data, auto_cfg$n_values, rightmost.closed = TRUE)
  col_idx <- min(col_idx, length(auto_cfg$n_values))

  # 自动设置的参数
  auto_params <- list(
    cellwidth = auto_cfg$cellwidth_values[col_idx],
    cellheight = auto_cfg$cellwidth_values[row_idx],
    fontsize_row = auto_cfg$fontsize_values[row_idx],
    fontsize_col = auto_cfg$fontsize_values[col_idx],
    dpi = max(auto_cfg$dpi_values[row_idx], auto_cfg$dpi_values[col_idx]),
    color = colorRampPalette(c("navy", "white", "firebrick3"))(100)
  )

  # 当 n_values 超过 3000 时不添加 names
  if (nrow_data > .HEATMAP_CONFIG$NAME_THRESHOLD) {
    auto_params$show_rownames <- FALSE
  }
  if (ncol_data > .HEATMAP_CONFIG$NAME_THRESHOLD) {
    auto_params$show_colnames <- FALSE
  }

  # 合并参数：用户参数优先，自动参数作为默认值
  final_params <- user_params
  for (param_name in names(auto_params)) {
    if (!param_name %in% names(final_params) ||
      is.null(final_params[[param_name]]) ||
      (length(final_params[[param_name]]) == 1 && is.na(final_params[[param_name]]))) {
      final_params[[param_name]] <- auto_params[[param_name]]
    }
  }

  return(final_params)
}

#' 校验聚类参数，避免无效聚类
validate_clustering <- function(matrix_data, cluster_rows, cluster_cols) {
  result <- list(cluster_rows = cluster_rows, cluster_cols = cluster_cols)

  if (isTRUE(cluster_rows)) {
    if (nrow(matrix_data) < 2) {
      warning("数据行数不足，无法进行行聚类，将cluster_rows设置为FALSE")
      result$cluster_rows <- FALSE
    } else {
      row_vars <- apply(matrix_data, 1, var, na.rm = TRUE)
      if (sum(row_vars > 0) < 2) {
        warning("数据变异不足，无法进行行聚类，将cluster_rows设置为FALSE")
        result$cluster_rows <- FALSE
      }
    }
  }

  if (isTRUE(cluster_cols)) {
    if (ncol(matrix_data) < 2) {
      warning("数据列数不足，无法进行列聚类，将cluster_cols设置为FALSE")
      result$cluster_cols <- FALSE
    } else {
      col_vars <- apply(matrix_data, 2, var, na.rm = TRUE)
      if (sum(col_vars > 0) < 2) {
        warning("数据变异不足，无法进行列聚类，将cluster_cols设置为FALSE")
        result$cluster_cols <- FALSE
      }
    }
  }

  return(result)
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

  # 检查是否支持该格式
  if (extension %in% names(.HEATMAP_CONFIG$FORMAT_MAPPING)) {
    return(.HEATMAP_CONFIG$FORMAT_MAPPING[[extension]])
  } else {
    # 默认返回png格式
    warning(sprintf("不支持的文件格式 '%s'，将使用PNG格式", extension))
    return("png")
  }
}

#' 数据质量检查和清理函数
#'
#' @param matrix_data 数据矩阵
#' @param remove_zero_variance 是否移除方差为0的行和列（聚类时需要）
#' @return 清理后的数据矩阵
clean_matrix_data <- function(matrix_data, remove_zero_variance = FALSE) {
  # 确保输入是矩阵格式
  if (is.data.frame(matrix_data)) {
    matrix_data <- as.matrix(matrix_data)
  } else if (!is.matrix(matrix_data)) {
    stop("输入数据必须是矩阵或数据框格式")
  }

  # 检查并处理NA、NaN、Inf值
  if (any(is.na(matrix_data)) || any(is.infinite(matrix_data))) {
    warning("发现NA、NaN或Inf值，正在清理数据...")

    # 将NA、NaN、Inf替换为0
    matrix_data[is.na(matrix_data) | is.infinite(matrix_data)] <- 0

    # 检查是否还有问题
    if (any(is.na(matrix_data)) || any(is.infinite(matrix_data))) {
      stop("数据清理后仍存在异常值，无法生成热图")
    }
  }

  # 检查方差为0的行和列（只有聚类时才需要移除）
  if (remove_zero_variance) {
    # 检查方差为0的行（所有值相同）
    row_vars <- apply(matrix_data, 1, var, na.rm = TRUE)
    zero_var_rows <- which(row_vars == 0 | is.na(row_vars))

    if (length(zero_var_rows) > 0) {
      warning(paste("发现", length(zero_var_rows), "行方差为0，将被移除"))
      matrix_data <- matrix_data[-zero_var_rows, , drop = FALSE]
    }

    # 检查方差为0的列（所有值相同）
    col_vars <- apply(matrix_data, 2, var, na.rm = TRUE)
    zero_var_cols <- which(col_vars == 0 | is.na(col_vars))

    if (length(zero_var_cols) > 0) {
      warning(paste("发现", length(zero_var_cols), "列方差为0，将被移除"))
      matrix_data <- matrix_data[, -zero_var_cols, drop = FALSE]
    }
  }

  # 检查矩阵是否为空
  if (nrow(matrix_data) == 0 || ncol(matrix_data) == 0) {
    stop("数据清理后矩阵为空，无法生成热图")
  }

  return(matrix_data)
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
  config <- .HEATMAP_CONFIG

  # 使用自动参数设置函数，先获取参数以确定是否需要聚类
  params <- auto_set_heatmap_parameters(matrix_data, autoset_image_specification, ...)
  # 根据聚类参数决定是否移除方差为0的数据
  needs_clustering <- (isTRUE(params$cluster_rows) || isTRUE(params$cluster_cols))
  matrix_data <- clean_matrix_data(matrix_data, remove_zero_variance = needs_clustering)

  # 获取dpi参数（用于输出设备设置）
  dpi <- if ("dpi" %in% names(params)) params$dpi else config$DEFAULT_DPI

  # 自动检测输出格式（如果未指定）
  if (is.null(output_format) && !is.null(filename)) {
    output_format <- detect_file_extension(filename)
  }

  # 计算所需尺寸
  dimensions <- do.call(calculate_heatmap_dimensions, c(list(matrix_data), params))

  # 针对位图设备限制进行像素尺寸自适应（避免超过设备最大尺寸）
  max_pixels_per_side <- config$MAX_PIXELS_PER_SIDE
  adjust_dpi_if_needed <- function(cur_dpi, width_in, height_in) {
    width_px <- width_in * cur_dpi
    height_px <- height_in * cur_dpi

    # 检查是否超过最大像素限制
    if (width_px > max_pixels_per_side || height_px > max_pixels_per_side) {
      scale_factor <- min(max_pixels_per_side / width_px, max_pixels_per_side / height_px)
      # 至少保持 72 DPI，防止过低导致设备报错或图形过于模糊
      new_dpi <- max(config$MIN_DPI, floor(cur_dpi * scale_factor * 0.95)) # 乘以0.95提供额外安全边际
      message(sprintf("警告：热图尺寸过大，自动调整DPI从 %d 到 %d", cur_dpi, new_dpi))
      return(new_dpi)
    }
    return(cur_dpi)
  }

  # 检查是否需要进行尺寸限制
  max_dimension_inches <- config$MAX_DIMENSION_INCHES
  if (dimensions$width > max_dimension_inches || dimensions$height > max_dimension_inches) {
    scale_factor <- min(max_dimension_inches / dimensions$width, max_dimension_inches / dimensions$height)
    dimensions$width <- dimensions$width * scale_factor
    dimensions$height <- dimensions$height * scale_factor
    message(sprintf(
      "警告：热图尺寸过大，已按比例缩放至 %.2f x %.2f 英寸",
      dimensions$width, dimensions$height
    ))
  }

  # 调试信息：输出 annotation_col 信息
  has_annotation_col <- "annotation_col" %in% names(params) &&
    !is.null(params$annotation_col) &&
    !identical(params$annotation_col, NA) &&
    !(is.data.frame(params$annotation_col) && nrow(params$annotation_col) == 0)

  if (has_annotation_col) {
    cat("=== Column Annotations 调试信息 ===\n")
    if (is.data.frame(params$annotation_col)) {
      cat("Column annotations 维度:", dim(params$annotation_col), "\n")
      cat("Column annotations 列名:", colnames(params$annotation_col), "\n")
    } else {
      cat("Column annotations 类型:", class(params$annotation_col), "\n")
    }
    cat("================================\n")
  } else {
    cat("=== Column Annotations 调试信息 ===\n")
    cat("未提供 column annotations 或为 NA\n")
    cat("================================\n")
  }

  cluster_checked <- validate_clustering(
    matrix_data,
    params$cluster_rows,
    params$cluster_cols
  )
  params$cluster_rows <- cluster_checked$cluster_rows
  params$cluster_cols <- cluster_checked$cluster_cols

  # 如果指定了文件名，设置输出设备
  if (!is.null(filename)) {
    # 确保设备能被安全关闭
    device_opened <- FALSE

    tryCatch(
      {
        if (output_format == "pdf") {
          pdf(
            file = filename,
            width = dimensions$width,
            height = dimensions$height
          )
          device_opened <- TRUE
        } else if (output_format == "png") {
          dpi <- adjust_dpi_if_needed(dpi, dimensions$width, dimensions$height)
          # 确保像素尺寸不超过限制
          width_px <- dimensions$width * dpi
          height_px <- dimensions$height * dpi

          if (width_px > max_pixels_per_side || height_px > max_pixels_per_side) {
            stop(sprintf(
              "热图尺寸过大 (%.0f x %.0f 像素)，无法生成PNG文件。建议：1) 使用PDF格式 2) 减少数据维度 3) 调整cellwidth/cellheight参数",
              width_px, height_px
            ))
          }

          png(
            filename = filename,
            width = width_px,
            height = height_px,
            res = dpi
          )
          device_opened <- TRUE
        } else if (output_format == "tiff") {
          dpi <- adjust_dpi_if_needed(dpi, dimensions$width, dimensions$height)
          tiff(
            filename = filename,
            width = dimensions$width * dpi,
            height = dimensions$height * dpi,
            res = dpi
          )
          device_opened <- TRUE
        } else if (output_format == "svg") {
          svg(
            filename = filename,
            width = dimensions$width,
            height = dimensions$height
          )
          device_opened <- TRUE
        } else if (output_format == "jpg" || output_format == "jpeg") {
          dpi <- adjust_dpi_if_needed(dpi, dimensions$width, dimensions$height)
          jpeg(
            filename = filename,
            width = dimensions$width * dpi,
            height = dimensions$height * dpi,
            res = dpi
          )
          device_opened <- TRUE
        } else {
          stop("不支持的输出格式。请选择: 'pdf', 'png', 'tiff', 'svg', 'jpg', 或 'jpeg'")
        }
      },
      error = function(e) {
        stop(sprintf("无法创建输出设备: %s\n建议使用PDF格式或减少热图尺寸", e$message))
      }
    )

    # 绘制热图并保证关闭设备
    on.exit(
      {
        if (isTRUE(device_opened)) {
          try(dev.off(), silent = TRUE)
        }
      },
      add = TRUE
    )

    # 绘制热图
    heatmap_result <- draw_heatmap_safe(matrix_data, params)

    message(sprintf(
      "热图已保存到: %s (尺寸: %.2f × %.2f 英寸)",
      filename, dimensions$width, dimensions$height
    ))
  } else {
    # 如果没有指定文件名，直接在设备上绘制
    heatmap_result <- draw_heatmap_safe(matrix_data, params)
  }

  return(invisible(heatmap_result))
}

#' 安全绘制热图，必要时禁用聚类
draw_heatmap_safe <- function(matrix_data, params) {
  # 过滤掉 NULL 或空字符串的参数，避免传递给 pheatmap 时出错
  params <- params[!sapply(params, function(x) is.null(x) || (is.character(x) && length(x) == 1 && !nzchar(x)))]

  tryCatch(
    {
      do.call(pheatmap, c(list(matrix_data), params))
    },
    error = function(e) {
      if (grepl("hclust|cluster", e$message, ignore.case = TRUE)) {
        warning("聚类失败，尝试生成无聚类热图...")
        params$cluster_rows <- FALSE
        params$cluster_cols <- FALSE
        do.call(pheatmap, c(list(matrix_data), params))
      } else {
        stop("热图生成失败: ", e$message)
      }
    }
  )
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
    make_option(c("--input"),
      type = "character", default = NULL,
      help = "输入数据文件（txt/csv/tsv/xlsx/xls），首列为行名", metavar = "character"
    ),
    make_option(c("--input_sheet"),
      type = "character", default = NULL,
      help = "当输入为 Excel 时，指定数据矩阵所在的 Sheet 名或序号"
    ),
    make_option(c("--annotation_row"),
      type = "character", default = NULL,
      help = "行注释文件（txt/csv/tsv/xlsx/xls），首列为行名"
    ),
    make_option(c("--annotation_row_sheet"),
      type = "character", default = NULL,
      help = "行注释 Excel Sheet 名或序号"
    ),
    make_option(c("--annotation_col"),
      type = "character", default = NULL,
      help = "列注释文件（txt/csv/tsv/xlsx/xls），首列为列名"
    ),
    make_option(c("--annotation_col_sheet"),
      type = "character", default = NULL,
      help = "列注释 Excel Sheet 名或序号"
    ),
    make_option(c("--output"),
      type = "character", default = NULL,
      help = "输出图文件名，自动根据后缀选择格式（pdf/png/tiff/svg/jpg）"
    ),
    # 常见作图参数
    make_option(c("--main"), type = "character", default = NULL, help = "主标题"),
    make_option(c("--scale"), type = "character", default = "none", help = "按行或列缩放：none,row,column"),
    make_option(c("--cluster_rows"), type = "logical", default = TRUE, help = "是否聚类行"),
    make_option(c("--cluster_cols"), type = "logical", default = TRUE, help = "是否聚类列"),
    make_option(c("--show_rownames"), type = "logical", default = TRUE, help = "显示行名"),
    make_option(c("--show_colnames"), type = "logical", default = TRUE, help = "显示列名"),
    make_option(c("--fontsize_row"), type = "double", default = NA, help = "行名字体大小"),
    make_option(c("--fontsize_col"), type = "double", default = NA, help = "列名字体大小"),
    make_option(c("--cellwidth"), type = "double", default = NA, help = "单元格宽度"),
    make_option(c("--cellheight"), type = "double", default = NA, help = "单元格高度"),
    make_option(c("--treeheight_row"), type = "double", default = 30, help = "行树高度"),
    make_option(c("--treeheight_col"), type = "double", default = 30, help = "列树高度"),
    make_option(c("--legend"), type = "logical", default = TRUE, help = "显示图例"),
    make_option(c("--dpi"), type = "integer", default = NA, help = "输出 DPI，留空自动"),
    make_option(c("--autoset"), type = "logical", default = TRUE, help = "启用自动规格设置")
  )
  opt_parser <- OptionParser(option_list = option_list)
  opt <- parse_args(opt_parser)

  if (is.null(opt$input)) stop("必须提供 --input")

  # 读取数据矩阵
  input_sheet_val <- parse_sheet_param(opt$input_sheet)
  data_df <- read_indexed_table(opt$input, sheet = input_sheet_val, expect_numeric = TRUE)
  data_mat <- as.matrix(data_df)

  # 读取注释
  row_anno <- NULL
  if (!is.null(opt$annotation_row)) {
    row_sheet_val <- parse_sheet_param(opt$annotation_row_sheet)
    row_anno <- read_indexed_table(opt$annotation_row, sheet = row_sheet_val, expect_numeric = FALSE)
  }
  col_anno <- NULL
  if (!is.null(opt$annotation_col)) {
    col_sheet_val <- parse_sheet_param(opt$annotation_col_sheet)
    col_anno <- read_indexed_table(opt$annotation_col, sheet = col_sheet_val, expect_numeric = FALSE)
  }

  # 对齐注释到数据
  if (!is.null(row_anno)) {
    common_rows <- intersect(rownames(data_mat), rownames(row_anno))
    if (length(common_rows) == 0) stop("行注释与数据行名没有交集")
    row_anno <- row_anno[common_rows, , drop = FALSE]
    data_mat <- data_mat[common_rows, , drop = FALSE]
  }
  if (!is.null(col_anno)) {
    common_cols <- intersect(colnames(data_mat), rownames(col_anno))
    if (length(common_cols) == 0) stop("列注释与数据列名没有交集")
    col_anno <- col_anno[common_cols, , drop = FALSE]
    data_mat <- data_mat[, common_cols, drop = FALSE]
  }

  # 组装 pheatmap 参数（用户参数优先，NA 表示自动）
  heatmap_params <- list(
    annotation_row = if (!is.null(row_anno)) row_anno else NA,
    annotation_col = if (!is.null(col_anno)) col_anno else NA,
    scale = opt$scale,
    cluster_rows = isTRUE(opt$cluster_rows),
    cluster_cols = isTRUE(opt$cluster_cols),
    show_rownames = isTRUE(opt$show_rownames),
    show_colnames = isTRUE(opt$show_colnames),
    fontsize_row = opt$fontsize_row,
    fontsize_col = opt$fontsize_col,
    cellwidth = opt$cellwidth,
    cellheight = opt$cellheight,
    treeheight_row = opt$treeheight_row,
    treeheight_col = opt$treeheight_col,
    legend = isTRUE(opt$legend),
    dpi = opt$dpi
  )

  # 生成
  do.call(
    smart_heatmap,
    c(
      list(
        matrix_data = data_mat,
        filename = opt$output,
        autoset_image_specification = isTRUE(opt$autoset)
      ),
      heatmap_params
    )
  )
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

# nrow or ncol 2   4   8   16  64  128 256 512
# cellwidth    120 60  30  20  15  7.5 3.3 1.2
# fontsize     10  10  10  10  10  8   4   1
# dpi          300 300 300 300 300 300 512 1024
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
