pkgs <- c(
  "ggthemes", "ggplot2", "pheatmap", "reshape2", "ggcorrplot", "corrplot",
  "FactoMineR", "plyr", "dplyr", "ropls", "ggrepel", "mixOmics",
  "openxlsx", "optparse", "magrittr", "data.table"
)
suppressPackageStartupMessages(invisible(lapply(pkgs, require, character.only = TRUE)))
Sys.setenv(R_LIBS="/home/data/opt/biosoft/R-422/lib64/R/library")

source('/home/colddata/qinqiang/script/Plot/Heatmap/heatmap_1.r', echo = TRUE, encoding = 'UTF-8')
source('/home/colddata/qinqiang/script/Plot/PCA/pca_1.r', echo = TRUE, encoding = 'UTF-8')
source('/home/colddata/qinqiang/script/Plot/Corrplot/correlation_1.r', echo = TRUE, encoding = 'UTF-8')

option_list <- list(
  make_option(c("-t", "--runtype"), type="character", default="normal",
              help="运行类型, normal 和 zscore", metavar="character"),
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="输入数据文件", metavar="character"),
  make_option(c("-s", "--samples"), type="character", default=NULL,
              help="样本信息文件", metavar="character"), 
  make_option(c("-c", "--compare"), type="character", default=NULL,
              help="比较组信息文件", metavar="character"),
  make_option(c("-d", "--definition"), type="character", default=NULL,
              help="代谢物定义文件", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default="./",
              help="输出目录", metavar="character"),
  make_option("--sort_by_class", type="logical", default=TRUE, metavar="logical",
              help="是否对 class 进行排序"),
  make_option(c("-l", "--log2data"), type="logical", default=FALSE,
              help="是否进行log2转换", metavar="logical"),
  make_option("--each_class_heatmap", type='logical', default=FALSE, metavar="logical",
              help="是否进行每个 Class 画一个 heatmap 图，常用于 Class 很多的情况"),
  make_option("--ncomp", type="integer", default=NULL, metavar="integer",
              help="设置 PLS-DA 的主成分数量 (ncomp)，如果不指定则根据样本数量自动计算")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# 检查必需参数
if (is.null(opt$input) || is.null(opt$samples)) {
  print_help(opt_parser)
  stop("必需的参数 --input 和 --samples 缺失", call.=FALSE)
}

# 定义颜色方案
my_set_colors <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF",
  "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494",
  "#B3B3B3", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D",
  "#666666", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F"
)

#' 绘制类别热图（综合图和分图）
#' @param data_frame 表达矩阵（行名需与compound_def一致）
#' @param compound_def 包含Class列的注释数据框
#' @param output_dir 输出目录路径
#' @param show_all 是否绘制综合热图（默认TRUE）
#' @param show_each 是否绘制每个类别的热图（默认TRUE）
row_class_heatmap <- function(data_frame, samples_df, compound_def, output_dir) {
  
  # 检查输入有效性
  if (!is.data.frame(compound_def) || nrow(compound_def) == 0) {
    warning("compound_def 为空或无效，跳过类别热图绘制")
    return(invisible(TRUE))
  }
  
  if (!all(rownames(data_frame) %in% rownames(compound_def))) {
    stop("表达矩阵的行名必须与注释数据框完全匹配")
  }
  
  # 创建输出目录
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # 准备注释数据
  annotation_row_df <- gen_annotation_row_df(compound_def)
  annotation_col_df <- gen_annotation_col_df(samples_df)
  unique_classes <- unique(compound_def$Class)

  for (class_name in unique_classes) {
    class_samples <- rownames(compound_def)[compound_def$Class == class_name]
    class_data <- data_frame[class_samples, , drop = FALSE]
    
    # 跳过样本不足的类别
    if (nrow(class_data) < 2) {
      warning("跳过类别 [", class_name, "] - 样本数不足 (n=", nrow(class_data), ")")
      next
    }
    
    # 安全文件名处理
    safe_name <- make.names(class_name)
    
    # 输出按 Class 的热图
    smart_heatmap(
      matrix_data = class_data, 
      filename = file.path(output_dir, paste0(safe_name, "_heatmap.jpg")),
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      annotation_row = annotation_row_df,
      annotation_col = annotation_col_df,
      scale = 'row',
      main = paste("Class: ", class_name)
    )
    smart_heatmap(
      matrix_data = class_data, 
      filename = file.path(output_dir, paste0(safe_name, "_heatmap_cluster_row.jpg")),
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      annotation_row = annotation_row_df,
      annotation_col = annotation_col_df,
      scale = 'row',
      main = paste("Class: ", class_name)
    )

    # 如果存在 SubClass，则为当前 Class 下的每个 SubClass 分别绘图，放入 Class 子目录
    if ("SubClass" %in% colnames(compound_def)) {
      class_dir <- file.path(output_dir, safe_name)
      if (!dir.exists(class_dir)) dir.create(class_dir, recursive = TRUE, showWarnings = FALSE)
      class_subclasses <- unique(compound_def[class_samples, "SubClass"]) 
      class_subclasses <- class_subclasses[!is.na(class_subclasses) & class_subclasses != ""]
      for (sub_name in class_subclasses) {
        sub_samples <- rownames(compound_def)[compound_def$Class == class_name & compound_def$SubClass == sub_name]
        sub_data <- data_frame[sub_samples, , drop = FALSE]
        if (nrow(sub_data) < 2) {
          warning("跳过子类别 [", class_name, "/", sub_name, "] - 样本数不足 (n=", nrow(sub_data), ")")
          next
        }
        safe_sub <- make.names(sub_name)
        smart_heatmap(
          matrix_data = sub_data,
          filename = file.path(class_dir, paste0(safe_sub, "_heatmap.jpg")),
          cluster_rows = FALSE,
          cluster_cols = FALSE,
          annotation_row = annotation_row_df,
          annotation_col = annotation_col_df,
          scale = 'row',
          main = paste("Class/SubClass: ", class_name, "/", sub_name)
        )
        smart_heatmap(
          matrix_data = sub_data,
          filename = file.path(class_dir, paste0(safe_sub, "_heatmap_cluster_row.jpg")),
          cluster_rows = TRUE,
          cluster_cols = FALSE,
          annotation_row = annotation_row_df,
          annotation_col = annotation_col_df,
          scale = 'row',
          main = paste("Class/SubClass: ", class_name, "/", sub_name)
        )
      }
    }
  }
  invisible(TRUE)
}


# 检查数据框列类型并转换为数值类型的函数
check_and_convert_numeric <- function(data_frame) {
  print(paste("检查输入数据的列类型..."))
  problematic_columns <- c()

  for (col_name in colnames(data_frame)) {
    col_data <- data_frame[[col_name]]
    
    # 检查是否为数值类型
    if (!is.numeric(col_data)) {
      # 如果是字符类型，先尝试去除空格
      if (is.character(col_data)) {
        col_data_cleaned <- trimws(col_data, which = "both")
        # 将空字符串转换为NA
        col_data_cleaned[col_data_cleaned == ""] <- NA
      } else {
        col_data_cleaned <- col_data
      }
      
      # 尝试转换为数值类型
      col_data_numeric <- suppressWarnings(as.numeric(col_data_cleaned))
      
      # 检查转换是否成功
      if (all(is.na(col_data_numeric)) || length(unique(col_data_numeric[!is.na(col_data_numeric)])) == 0) {
        problematic_columns <- c(problematic_columns, col_name)
        print(paste("警告: 列", col_name, "无法转换为数值类型"))
      } else {
        data_frame[[col_name]] <- col_data_numeric
      }
    }
  }

  # 如果有问题列，输出并退出程序
  if (length(problematic_columns) > 0) {
    print(paste("以下列无法转换为数值类型:"))
    for (col in problematic_columns) {
      print(paste("  -", col))
    }
    print("程序退出")
    quit(status = 1)
  }

  print(paste("所有列类型检查完成，数据已转换为数值类型"))
  return(data_frame)
}

# 将数据框/矩阵中的 0 替换为一个很小的正数（默认 1e-6）
# 仅对数值型列生效，非数值列不做改动
replace_zeros_with_epsilon <- function(x, epsilon = 1e-6, tolerance = 0) {
  if (is.matrix(x)) {
    x[x == 0] <- epsilon
    return(x)
  }
  if (is.data.frame(x)) {
    for (col_name in colnames(x)) {
      col_data <- x[[col_name]]
      if (is.numeric(col_data)) {
        zero_idx <- !is.na(col_data) & if (tolerance > 0) abs(col_data) <= tolerance else col_data == 0
        if (any(zero_idx)) {
          col_data[zero_idx] <- epsilon
          x[[col_name]] <- col_data
        }
      } else if (is.character(col_data)) {
        trimmed <- trimws(col_data)
        suppressWarnings(num_vals <- as.numeric(trimmed))
        zero_idx <- !is.na(num_vals) & if (tolerance > 0) abs(num_vals) <= tolerance else num_vals == 0
        if (any(zero_idx)) {
          if (all(!is.na(num_vals) | is.na(trimmed))) {
            num_vals[zero_idx] <- epsilon
            x[[col_name]] <- num_vals
          } else {
            trimmed[zero_idx] <- as.character(epsilon)
            x[[col_name]] <- trimmed
          }
        }
      } else if (is.factor(col_data)) {
        char_vals <- as.character(col_data)
        trimmed <- trimws(char_vals)
        suppressWarnings(num_vals <- as.numeric(trimmed))
        zero_idx <- !is.na(num_vals) & if (tolerance > 0) abs(num_vals) <= tolerance else num_vals == 0
        if (any(zero_idx)) {
          trimmed[zero_idx] <- as.character(epsilon)
          x[[col_name]] <- factor(trimmed)
        }
      }
    }
    return(x)
  }
  stop("输入必须是 data.frame 或 matrix")
}

# 根据 definition_df 生成注释行数据框
gen_annotation_row_df <- function(definition_df) {
  # 若为 NA 或非数据框或为空，直接返回 NA
  if (!is.data.frame(definition_df) || nrow(definition_df) == 0) {
    return(NA)
  }
  # 若不存在 Class 列，直接返回 NA
  if (!("Class" %in% colnames(definition_df))) {
    return(NA)
  }
  if ("SubClass" %in% colnames(definition_df)) {
    annotation_row_df <- data.frame(
      Class = definition_df$Class,
      SubClass = definition_df$SubClass,
      row.names = rownames(definition_df)
    )
  } else {
    annotation_row_df <- data.frame(
      Class = definition_df$Class,
      row.names = rownames(definition_df)
    )
  }
  return(annotation_row_df)
}

# 根据 samples_df 生成列注释数据框（用于热图列注释）
gen_annotation_col_df <- function(samples_df) {
  # 若为 NA 或非数据框或为空，直接返回 NA
  if (!is.data.frame(samples_df) || nrow(samples_df) == 0) {
    return(NA)
  }
  # 需要存在 group 和 sample 两列
  required_cols <- c("group", "sample")
  if (!all(required_cols %in% colnames(samples_df))) {
    return(NA)
  }
  # 去除缺失或空样本名
  valid_rows <- !is.na(samples_df$sample) & trimws(samples_df$sample) != ""
  if (!any(valid_rows)) {
    return(NA)
  }
  samples_df <- samples_df[valid_rows, required_cols, drop = FALSE]
  # 构造列注释：行名为样本名，列为 group
  annotation_col_df <- data.frame(
    group = samples_df$group,
    row.names = samples_df$sample,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  return(annotation_col_df)
}

# 代谢物分析主函数
metabolite_analysis <- function(samples_file, data_matrix, definition_df = NA, output_dir, log2data = FALSE, each_class_heatmap = FALSE, ncomp = NULL) {
  select_sample_info <- read.table(
    samples_file, 
    sep = "\t", 
    header = T, 
    check.names = F, 
    stringsAsFactors = F
  )

  groups <- select_sample_info$group

  # 按照样本排列顺序
  select_data_frame <- data_matrix[, select_sample_info$sample, drop = FALSE]
  # mutligroup samples 的时候会出现每行都为 0 跑不了的情况
  select_data_frame <- select_data_frame[rowSums(select_data_frame != 1e-6) > 0, ]
  metabolites <- as.matrix(t(select_data_frame))
  metabolites <- metabolites[rowSums(metabolites != 1e-6) > 0, ]
  
  # 创建 annotation_row_df，函数内自处理 NA 情况
  if (is.data.frame(definition_df) && nrow(definition_df) > 0) {
    current_definition_df <- definition_df[rownames(select_data_frame), , drop = FALSE]
  } else {
    current_definition_df <- NA
  }
  annotation_row_df <- gen_annotation_row_df(current_definition_df)
  annotation_col_df <- gen_annotation_col_df(select_sample_info)

  smart_heatmap(
    matrix_data = select_data_frame,
    file.path(output_dir, 'Compounds_heatmap.png'),
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    annotation_row = annotation_row_df,
    annotation_col = annotation_col_df,
    scale = 'row'
  )
  smart_heatmap(
    matrix_data = select_data_frame,
    file.path(output_dir, 'Compounds_heatmap_cluster_row.png'),
    cluster_cols = FALSE,
    cluster_rows = TRUE,
    annotation_row = annotation_row_df,
    annotation_col = annotation_col_df,
    scale = 'row'
  )
  
  # 每个 Class 一张图
  if (each_class_heatmap) {
    each_class_heatmap_dir = file.path(output_dir, 'Each_class_heatmap')
    dir.create(each_class_heatmap_dir)
    row_class_heatmap(
      data_frame = select_data_frame,
      samples_df = select_sample_info,
      compound_def = current_definition_df,
      output_dir = each_class_heatmap_dir
    )
  }

  # 根据样本数量设置 ncomp（如果用户未指定）
  num_samples <- nrow(metabolites)
  num_metabolites <- ncol(metabolites)
  num_groups <- length(unique(groups))
  
  # 检查数据维度是否满足 plsda 的要求（至少需要3行3列）
  if (num_samples < 3 || num_metabolites < 3) {
    warning(paste("数据维度不足（样本数:", num_samples, "代谢物数:", num_metabolites, 
                  "），无法进行 PLS-DA 分析。PLS-DA 需要至少3个样本和3个代谢物。"))
    return(invisible(NULL))
  }
  
  if (is.null(ncomp)) {
    # 一般规则：ncomp 不应超过样本数量的 1/10，且不应超过组数-1，但是 2 个不好看，目前至少画 4 个
    max_ncomp <- min(floor(num_samples / 10), num_groups - 1, num_metabolites - 1)
    
    # 设置合理的 ncomp 值
    if (num_samples < 10) {
      ncomp <- min(4, max_ncomp)
    } else if (num_samples < 20) {
      ncomp <- min(6, max_ncomp)
    } else if (num_samples < 50) {
      ncomp <- min(8, max_ncomp)
    } else {
      ncomp <- min(10, max_ncomp)
    }
    
    # 确保 ncomp 至少为 4（如果样本数允许），但不超过实际限制
    ncomp <- max(4, min(4, ncomp, max_ncomp))
  } else {
    # 用户指定的 ncomp 需要满足限制条件
    max_ncomp <- min(floor(num_samples / 10), num_groups - 1, num_metabolites - 1)
    if (ncomp > max_ncomp) {
      warning(paste("指定的 ncomp (", ncomp, ") 超过最大允许值 (", max_ncomp, 
                    ")，已自动调整为", max_ncomp))
      ncomp <- max_ncomp
    }
    # 确保 ncomp 至少为 2
    if (ncomp < 2) {
      warning(paste("ncomp 必须至少为 2，已自动调整为 2"))
      ncomp <- 2
    }
  }
  
  # 打印调试信息
  print(paste("样本数量:", num_samples, "代谢物数量:", num_metabolites, 
              "组数:", num_groups, "设置的ncomp:", ncomp))
  
  df_plsda <- plsda(metabolites, groups, ncomp = ncomp)

  # scree plot
  scree_df <- as.data.frame(df_plsda$prop_expl_var$X)
  scree_df$comp <- rownames(scree_df)
  colnames(scree_df)[1] <- "value"

  df <- data.frame(x = letters[1:6], y = seq(10, 60, 10))
  # 创建 scree plot 图
  scree_df$comp <- factor(scree_df$comp, levels = scree_df$comp)
  scree_plot <- ggplot(scree_df, aes(x = comp, y = value)) +
    #  geom_bar(stat = "identity", fill = "steelblue") +
    labs(x = "Principal Component", y = "Explained Variance") +
    ggtitle("Scree Plot") +
    theme_minimal() +
    geom_line(aes(group = 1), color = "red") +
    geom_point()

  ggsave(
    file.path(output_dir, "/Metabolite_quantitation_scree_plot.jpeg"), 
    scree_plot, 
    width = ncomp * 0.7 + 1, 
    height = 4
  )

  comp_load_df <- as.data.frame(df_plsda$loadings$X)
  comp_load_df <- cbind(rownames(comp_load_df), comp_load_df)
  colnames(comp_load_df)[1] <- "compound_name"
  
  write.xlsx(
    comp_load_df, 
    file = file.path(output_dir, "pc_loading_value.xlsx"), 
    sheetName = "Sheet1", 
    rowNames = FALSE
  )

  df <- unclass(df_plsda)
  df1 <- as.data.frame(df$variates$X)
  df1$group <- groups
  df1$samples <- rownames(df1)

  explain <- df$prop_expl_var$X
  x_label <- round(explain[1], digits = 3)
  y_label <- round(explain[2], digits = 3)

  num_samples <- nrow(df1)
  plot_width <- 5 + num_samples * 0.2
  plot_height <- 4 + num_samples * 0.2


  p1 <- ggplot(df1, aes(x = comp1, y = comp2, color = group)) +
    theme_bw() +
    geom_point(size = 5) +
    geom_text_repel(
      aes(label = samples),
      size = 5,
      box.padding = 0.35,
      point.padding = 0.5
    ) +
    labs(
      x = paste0("P1 (", x_label * 100, "%)"), 
      y = paste0("P2 (", y_label * 100, "%)")
    ) +
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

  x_min <- x_range[1] - 2
  x_max <- x_range[2] + 2
  y_min <- y_range[1]
  y_max <- y_range[2]
  
  p1 <- p1 + coord_cartesian(
    xlim = c(x_min, x_max), 
    ylim = c(y_min, y_max)
  )

  ggsave(
    file.path(output_dir, "Multigroup_Plsda_Distribution_Graph.jpeg"), 
    p1, 
    width = plot_width, 
    height = plot_height
  )

  # vip 值使用 fpkm, 其他值计算都是用 reads
  # predI 这里不知道为什么有时候需要写 NA 才行
  plsda_model <- opls(
    x = metabolites, 
    y = groups, 
    predI = 1, 
    crossvalI = 5
  )
  
  vip_values <- plsda_model@vipVn
  df.vip <- as.data.frame(vip_values)

  reads_data_with_def <- cbind(select_data_frame, VIP = df.vip$vip_values)
  reads_data_with_def$Metabolite <- rownames(reads_data_with_def)
  reads_data_with_def <- reads_data_with_def[
    , c("Metabolite", setdiff(names(reads_data_with_def), "Metabolite"))
  ]

  # 没有定义跳过
  class_count <- ""
  if (is.data.frame(definition_df)) {
    reads_data_with_def <- merge(
      reads_data_with_def, 
      definition_df, 
      by = "Metabolite", 
      all.x = TRUE
    )
    
    greater_than_one_data_def <- reads_data_with_def[reads_data_with_def$VIP > 1, ]
    # 如果 greater_than_one_data_def 没有数据，就是用 reads_data_with_def
    if (nrow(greater_than_one_data_def) == 0) {
      print("没有 VIP 大于 1，输出全部")
      greater_than_one_data_def <- reads_data_with_def
    }

    # 添加 class count 和对应的名称
    # 创建函数来处理分类统计和导出
    process_class_analysis <- function(data, class_col, output_filename) {
      if (class_col %in% colnames(data)) {
        class_count <- data %>%
          group_by(!!sym(class_col)) %>%
          summarize(
            count = n(),
            compounds = paste(Metabolite, collapse = ", ")
          ) %>%
          filter(!!sym(class_col) != "") %>%
          arrange(desc(count))
        
        write.xlsx(
          class_count,
          file = file.path(output_dir, output_filename),
          sheetName = "Sheet1", 
          rowNames = FALSE, 
          colNames = TRUE
        )
      }
    }
    
    # 处理Class和SubClass分析
    write.xlsx(greater_than_one_data_def, file=file.path(output_dir, 'temp.xlsx'), rowNames=TRUE, colNames=TRUE)    
    process_class_analysis(greater_than_one_data_def, "Class", "Significant_compound_count_by_class.xlsx")
    process_class_analysis(greater_than_one_data_def, "SubClass", "Significant_compound_count_by_subclass.xlsx")
  }

  reads_data_with_def <- reads_data_with_def[
    order(-reads_data_with_def$VIP, na.last = TRUE), 
  ]
  
  write.xlsx(
    reads_data_with_def, 
    file = file.path(output_dir, "Metabolite_quantitation_VIP.xlsx"), 
    sheetName = "Sheet1", 
    rowNames = FALSE, 
    colNames = TRUE
  )
}

# 组间分析函数
comparison_analysis <- function(compare_file, samples_file, reads_df, fpkm_df = NA, definition_df = NA, output_dir, log2data = FALSE, each_class_heatmap = FALSE) {
  sample_info <- read.table(samples_file, sep = "\t", header = T, check.names = F, stringsAsFactors = F)
  comp_info <- read.table(compare_file, sep = "\t", header = T, check.names = F, stringsAsFactors = F)
  
  
  comparisons <- list()
  for (i in seq_along(1:nrow(comp_info))) {
    comparisons <- append(comparisons, list(as.character(comp_info[i, ])))
  }

  plot_types <- c(
    "correlation", "outlier", "overview", "permutation",
    "predict-train", "x-loading", "x-score",
    "x-variance", "xy-score"
  )

  deg_data <- data.frame(
    group = character(0), 
    Total = numeric(0),
    All_Diff_compound = numeric(0), 
    Up = numeric(0), 
    Down = numeric(0),
    Up_ration = character(0),
    Down_ratio = character(0),
    Up_list = character(0), 
    Down_list = character(0)
  )

  # 循环中注意可能需要修改 corssvalI 值
  for (i in seq_along(comparisons)) {
    compare_name <- paste(comparisons[[i]], collapse = "-vs-")
    compare_path <- file.path(output_dir, compare_name)
    dir.create(compare_path, showWarnings = FALSE)
    print(paste0("Processing:", compare_name)) # 打印当前处理的组合名称
    groups <- comparisons[[i]]
    groupA_cols <- sample_info$sample[sample_info$group == groups[1]]
    groupB_cols <- sample_info$sample[sample_info$group == groups[2]]

    # 获取当前两组的样本的表达量数据
    current_samples <- c(groupA_cols, groupB_cols)

    if (is.data.frame(fpkm_df)) {
      current_expression_fpkm_data <- fpkm_df[
        , current_samples, 
        drop = FALSE
      ]
      rows_to_remove <- apply(
        current_expression_fpkm_data, 
        1, 
        function(row) all(row == 1e-6)
      )
      current_expression_fpkm_data <- current_expression_fpkm_data[!rows_to_remove, ]
    }
    
    current_expression_data <- reads_df[, current_samples, drop = FALSE]
    keep_rows <- rowSums(current_expression_data != 1e-6, na.rm = TRUE) > 0
    current_expression_data <- current_expression_data[keep_rows, , drop = FALSE]

    # 创建 annotation_row_df，函数内自处理 NA 情况
    if (is.data.frame(definition_df) && nrow(definition_df) > 0) {
      current_definition_df <- definition_df[rownames(current_expression_data), , drop = FALSE]
    } else {
      current_definition_df <- NA
    }
    current_annotation_row_df <- gen_annotation_row_df(current_definition_df)
    current_annotation_col_df <- gen_annotation_col_df(sample_info[sample_info$sample %in% current_samples, , drop = FALSE])
    
    # p_value_current_expression_data <- reads_data[, current_samples, drop = FALSE]
    # 检查因变量的水平数量
    y_factor <- factor(
      sample_info$group[sample_info$sample %in% current_samples],
      levels = groups
    )
    if (length(unique(y_factor)) < 2) {
      print(paste(
        "y_factor 小于 2 个", 
        compare_name, 
        " 确认 compare_info 和 samples_described.txt 是否完全匹配"
      ))
      next
    }

    # 检查是否有缺失值
    if (any(is.na(current_expression_data))) {
      print(paste("Missing values found in expression data for:", compare_name))
      next # 跳过当前的迭代
    }

    # 进行OPLS分析
    if (is.data.frame(fpkm_df)) {
      transposed_expression_data <- t(current_expression_fpkm_data)
    } else {
      transposed_expression_data <- t(current_expression_data)
    }
    
    crossval_number <- nrow(transposed_expression_data) - 1
    opls_model <- try(
      opls(
        x = transposed_expression_data, 
        y = y_factor, 
        predI = 1, 
        orthoI = 2, 
        crossvalI = crossval_number # crossvalI 默认是 7, crossvalI 需要小于等于两组样本的数量
      ),
      silent = FALSE
    )

    # 检查模型是否成功建立
    if (class(opls_model) == "try-error") {
      print(paste("Failed to build model for:", compare_name))
      next # 跳过当前的迭代
    }

    # 获取VIP值，如果模型失败则赋予 0
    vip_values <- tryCatch(
      {opls_model@vipVn},
      error = function(e) {rep(0, ncol(current_expression_data))}
    )

    # 组间 heatmap 图
    if (is.data.frame(fpkm_df)) {
      smart_heatmap(
        matrix_data = current_expression_fpkm_data,
        file.path(compare_path, paste0(compare_name, "_heatmap.png")),
        cluster_cols = FALSE,
        cluster_rows = FALSE,
        annotation_row = current_annotation_row_df,
        annotation_col = current_annotation_col_df,
        scale = 'row'
      )
      smart_heatmap(
        matrix_data = current_expression_fpkm_data,
        file.path(compare_path, paste0(compare_name, "_heatmap_cluster_row.png")),
        cluster_cols = FALSE,
        cluster_rows = TRUE,
        annotation_row = current_annotation_row_df,
        annotation_col = current_annotation_col_df,
        scale = 'row'
      )
    } else {
      smart_heatmap(
        matrix_data = current_expression_data,
        file.path(compare_path, paste0(compare_name, "_heatmap.png")),
        cluster_cols = FALSE,
        cluster_rows = FALSE,
        annotation_row = current_annotation_row_df,
        annotation_col = current_annotation_col_df,
        scale = 'row'
      )
      smart_heatmap(
        matrix_data = current_expression_data,
        file.path(compare_path, paste0(compare_name, "_heatmap_cluster_row.png")),
        cluster_cols = FALSE,
        cluster_rows = TRUE,
        annotation_row = current_annotation_row_df,
        annotation_col = current_annotation_col_df,
        scale = 'row'
      )
    }

    # 计算FoldChange
    baseMeanA <- rowMeans(current_expression_data[, groupA_cols])
    baseMeanB <- rowMeans(current_expression_data[, groupB_cols])
    FoldChange <- ifelse(baseMeanB > 0, abs(baseMeanA / baseMeanB), 0)

    # 定义一个函数进行成对t检验
    paired_t_test <- function(x, y) {
      # 如果x中所有值都相同
      if (length(unique(x)) == 1) {
        x[1] <- x[1] + 0.000000001
      }

      # 如果y中所有值都相同
      if (length(unique(y)) == 1) {
        y[1] <- y[1] + 0.000000001
      }
      test_result <- t.test(x, y)
      return(test_result$p.value)
    }

    p_values <- apply(
      current_expression_data, 
      1, 
      function(row) paired_t_test(row[groupA_cols], row[groupB_cols])
    )

    current_expression_data_def <- cbind(
      current_expression_data, 
      baseMeanA, 
      baseMeanB, 
      FoldChange, 
      pvalues = p_values
    )
    current_expression_data_def <- as.data.frame(current_expression_data_def)

    # 创建一个长度与current_expression_data行数相同的全零向量
    vip_values_adjusted <- rep(0, nrow(current_expression_data_def))
    # 更新vip_values_adjusted向量中相应的位置
    # 这假设vip_values的名称与current_expression_data的行名对应
    if (!is.null(names(vip_values))) {
      matching_rows <- match(names(vip_values), rownames(current_expression_data_def))
      vip_values_adjusted[matching_rows] <- vip_values
    }
    
    current_expression_data_def <- cbind(
      current_expression_data_def, 
      VIP = vip_values_adjusted
    )

    current_expression_data_def <- current_expression_data_def[
      order(current_expression_data_def$pvalues, na.last = TRUE), 
    ]
    
    current_expression_data_def$Metabolite <- rownames(current_expression_data_def)
    current_expression_data_def <- current_expression_data_def[
      , c("Metabolite", setdiff(names(current_expression_data_def), "Metabolite"))
    ]
    
    current_expression_data_def$padj <- p.adjust(
      current_expression_data_def$pvalues, 
      "BH"
    )

    if (is.data.frame(definition_df) && nrow(definition_df) > 0) {
      current_expression_data_def <- merge(
        current_expression_data_def, 
        definition_df, 
        by = "Metabolite", 
        all.x = TRUE
      )
    }

    current_expression_data_def <- current_expression_data_def[
      order(-current_expression_data_def$VIP, na.last = TRUE), 
    ]
    current_expression_data_def[is.na(current_expression_data_def)] <- 0

    # 将更新后的数据框保存为文本文件
    write.xlsx(
      current_expression_data_def,
      file = file.path(compare_path, paste0(compare_name, "_VIP.xlsx")),
      sheetName = "Sheet1",
      rowNames = FALSE,
      colNames = TRUE
    )

    deg_df <- current_expression_data_def[
      current_expression_data_def$VIP > 1, 
    ]
    deg_up <- nrow(deg_df[deg_df$FoldChange >= 1.2, ])
    deg_up_idlist <- deg_df[deg_df$FoldChange >= 1.2, ]$Metabolite
    deg_down <- nrow(deg_df[deg_df$FoldChange <= 0.8, ])
    deg_down_idlist <- deg_df[deg_df$FoldChange <= 0.8, ]$Metabolite
    
    new_row <- data.frame(
      group = compare_name,
      Total = nrow(current_expression_data_def),
      All_Diff_compound = deg_up + deg_down,
      Up = deg_up,
      Down = deg_down,
      Up_ration = paste0(round(deg_up / nrow(current_expression_data_def) * 100, 2), '%'),
      Down_ratio = paste0(round(deg_down / nrow(current_expression_data_def) * 100, 2), '%'),
      Up_list = paste0(deg_up_idlist, collapse = ","),
      Down_list = paste0(deg_down_idlist, collapse = ",")
    )
    deg_data <- rbind(deg_data, new_row)

    # 每个 Class 一张图
    if (each_class_heatmap && is.data.frame(current_definition_df) && nrow(current_definition_df) > 0) {
      each_class_heatmap_dir <- file.path(compare_path, 'Each_class_heatmap')
      dir.create(each_class_heatmap_dir, showWarnings = FALSE)
      
      # 根据是否有 fpkm 数据选择使用的数据
      if (is.data.frame(fpkm_df)) {
        data_for_class_heatmap <- current_expression_fpkm_data
      } else {
        data_for_class_heatmap <- current_expression_data
      }
      
      # 使注释行名与热图数据完全一致
      if (is.data.frame(definition_df) && nrow(definition_df) > 0) {
        current_definition_df <- definition_df[rownames(data_for_class_heatmap), , drop = FALSE]
      } else {
        current_definition_df <- NA
      }
      
      row_class_heatmap(
        data_frame = data_for_class_heatmap,
        samples_df = sample_info[sample_info$sample %in% current_samples, , drop = FALSE],
        compound_def = current_definition_df,
        output_dir = each_class_heatmap_dir
      )
    }

    for (plot_type in plot_types) {
      file_name <- file.path(
        compare_path, 
        paste0(compare_name, "_OPLS_DA_", plot_type, ".png")
      )

      # 检查模型是否有效
      if (inherits(opls_model, "opls")) {
        png(file_name)
        try({plot(opls_model, typeVc = plot_type)})
        dev.off() # 确保关闭图形设备
      }
    }
  }
  
  write.xlsx(
    deg_data, 
    file = file.path(output_dir, "Significant_differential_metabolite_count_summary.xlsx"), 
    sheetName = "Sheet1", 
    rowNames = FALSE, 
    colNames = TRUE
  )
}

# ================================= 文件预处理 =========================================
zhengtifenxi_dir <- "03_代谢物整体分析"
chayifenxi_dir <- "04_代谢物差异分析"
comparison_analysis_dir <- file.path(chayifenxi_dir, "组间分析")

dir.create("00_Background_materials")
dir.create("01_原始质谱数据")
dir.create("02_代谢物定量图数据")
dir.create(zhengtifenxi_dir)
dir.create(chayifenxi_dir)
dir.create(comparison_analysis_dir)

samples_file <- opt$samples
sample_info <- read.table(samples_file, sep = "\t", header = T, check.names = F, stringsAsFactors = F)

# 读取数据
reads_data <- read.xlsx(opt$input, sheet = 1, rowNames = TRUE)
# reads_data <- fread(opt$input, fileEncoding = "UTF-8", check.names = TRUE)
# rownames(reads_data) <- make.names(reads_data[[1]], unique = TRUE)
# reads_data <- reads_data[, -1, with = FALSE]

# 检查并转换数据框列类型
reads_data <- check_and_convert_numeric(reads_data)
# 替换 0 为 1e-6
reads_data <- replace_zeros_with_epsilon(reads_data)
# 按照 sample 样本列进行排序
reads_data <- reads_data[sample_info$sample]
# 去掉一行全是 1e-6 的
reads_data <- reads_data[rowSums(reads_data != 1e-6) > 0, ]



Compound_def_file <- opt$definition
# 如果需要合并定义则读取单独定义文件，没有则跳过
if (!is.null(opt$definition)) {
  definition_df <- read.xlsx(Compound_def_file, sheet = 1, rowNames = TRUE)
  # definition_df <- fread(Compound_def_file, fileEncoding = "UTF-8", check.names = TRUE)
  # rownames(definition_df) <- make.names(definition_df[[1]], unique = TRUE)
  # definition_df <- definition_df[, -1, with = FALSE]
  definition_df$Metabolite <- rownames(definition_df)

  if (opt$sort_by_class) {
    # 按照 class 排序，reads_data 也按照 class 的顺序排序
    definition_df <- definition_df[order(definition_df$Class), ]
    sorted_ids <- rownames(definition_df)
    reads_data <- reads_data[sorted_ids, ]
  }
  current_definition_df <- definition_df[rownames(reads_data), , drop = FALSE]
  annotation_row_df <- gen_annotation_row_df(current_definition_df)
} else {
  definition_df <- NA
  current_definition_df <- NA
  annotation_row_df <- NA
}

# ====================== 主程序流程 =================================
# 1) 根据 runtype 生成处理后的数据和主图参数
if (opt$runtype == "zscore") {
  data_used <- as.data.frame(t(apply(reads_data, 1, function(x) {
    (x - mean(x)) / (sd(x)**0.5)
  })))
  fpkm_for_pair <- data_used
  fpkm_for_pair <- replace_zeros_with_epsilon(fpkm_for_pair)
} else {
  data_used <- reads_data
  fpkm_for_pair <- FALSE
}

annotation_col_df <- gen_annotation_col_df(sample_info)

# 3) 整体热图
smart_heatmap(
  matrix_data = data_used,
  filename = file.path(zhengtifenxi_dir, "All_metabolites_heatmap.png"),
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  annotation_row = annotation_row_df,
  annotation_col = annotation_col_df,
  scale = 'row',
  main = 'All metabolites heatmap'
)
smart_heatmap(
  matrix_data = data_used,
  filename = file.path(zhengtifenxi_dir, "All_metabolites_heatmap_cluster_row.png"),
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  annotation_row = annotation_row_df,
  annotation_col = annotation_col_df,
  scale = 'row',
  main = 'All metabolites heatmap'
)

# 4) 组平均热图与导出
print(paste0('生成', file.path(zhengtifenxi_dir, "All_metabolites_groupmean_heatmap.png")))
grouped_means <- sapply(unique(sample_info$group), function(g) {
  samples_in_group <- sample_info$sample[sample_info$group == g]
  if (length(samples_in_group) == 1) {
    data_used[[samples_in_group]]
  } else {
    rowMeans(data_used[, samples_in_group, drop = FALSE])
  }
})

grouped_means <- as.data.frame(grouped_means)
colnames(grouped_means) <- unique(sample_info$group)
rownames(grouped_means) <- rownames(reads_data)

smart_heatmap(
  matrix_data = grouped_means,
  filename = file.path(zhengtifenxi_dir, 'All_metabolites_groupmean_heatmap.png'),
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  annotation_row = annotation_row_df,
  annotation_col = NA,
  scale = 'row'
)
smart_heatmap(
  matrix_data = grouped_means,
  filename = file.path(zhengtifenxi_dir, 'All_metabolites_groupmean_heatmap_cluster_row.png'),
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  annotation_row = annotation_row_df,
  annotation_col = NA,
  scale = 'row'
)

write.xlsx(
  grouped_means,
  file = file.path(zhengtifenxi_dir, "All_metabolites_groupmean.xlsx"),
  sheetName = "Sheet1",
  rowNames = TRUE,
  colNames = TRUE
)

# 5) 可选 log2 图（保留原有逻辑差异）
if (opt$log2data) {
  zhengtifenxi_log2heatmap_pic_name <- file.path(
    zhengtifenxi_dir,
    "All_metabolites_log_heatmap.png"
  )
  if (opt$runtype == "normal") {
    log2data_reads_data <- log2(reads_data + 1)
    smart_heatmap(
      matrix_data = log2data_reads_data,
      filename = zhengtifenxi_log2heatmap_pic_name,
      cluster_cols = FALSE,
      cluster_rows = TRUE,
      annotation_row = annotation_row_df,
      annotation_col = annotation_col_df,
      scale = 'row',
      main = 'All metabolites log2 heatmap'
    )
  } else {
    log2data_reads_data <- log2(reads_data + 1)
    heatmap_plot(log2data_reads_data, zhengtifenxi_log2heatmap_pic_name)
  }
}

# 6) 按类别热图（根据数据来源切换）
if (opt$each_class_heatmap) {
  each_class_heatmap_dir <- file.path(zhengtifenxi_dir, 'Each_class_heatmap')
  dir.create(each_class_heatmap_dir)
  row_class_heatmap(
    data_frame = data_used,
    samples_df = sample_info,
    compound_def = current_definition_df,
    output_dir = each_class_heatmap_dir
  )
}

# 7) 相关性图与 PCA（根据数据来源切换）
cor_plot(
  data_used,
  file.path(zhengtifenxi_dir, 'Metabolite_correlation_graph.png')
)

pca_plot(
  data_frame = data_used,
  samples_file = samples_file,
  output_prefix = file.path(zhengtifenxi_dir, 'Metabolite_')
)

# 8) 多组分析（根据数据来源切换）
multigroup_dir <- file.path(chayifenxi_dir, "多组分析")
dir.create(multigroup_dir)

metabolite_analysis(
  samples_file = samples_file,
  data_matrix = data_used,
  definition_df = current_definition_df,
  output_dir = multigroup_dir,
  log2data = opt$log2data,
  each_class_heatmap = opt$each_class_heatmap,
  ncomp = opt$ncomp
)

# 9) 组间分析（zscore 时传递 fpkm；normal 传 FALSE，保持原逻辑）
if (!is.null(opt$compare)) {
  comparison_analysis(
    compare_file = opt$compare,
    samples_file = samples_file,
    reads_df = reads_data,
    fpkm_df = fpkm_for_pair,
    definition_df = current_definition_df,
    output_dir = comparison_analysis_dir,
    log2data = opt$log2data,
    each_class_heatmap = opt$each_class_heatmap
  )
}



