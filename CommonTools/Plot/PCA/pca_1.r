library(FactoMineR)  # PCA分析
library(ggplot2)     # 绘图
library(ggrepel)     # 文本标签防重叠
library(plyr)        # ddply函数
library(openxlsx)    # Excel文件读写

pca_plot <- function(data_frame, samples_file, output_prefix) {
  # 检查输入数据是否有效
  if (is.null(data_frame) || nrow(data_frame) == 0 || ncol(data_frame) == 0) {
    warning("输入数据为空，无法进行PCA分析")
    return(invisible(NULL))
  }
  
  # 确保数据是数值型
  if (!is.matrix(data_frame) && !is.data.frame(data_frame)) {
    warning("输入数据格式不正确，需要矩阵或数据框")
    return(invisible(NULL))
  }
  
  # 转换为矩阵并检查是否为数值型
  data_mat <- as.matrix(data_frame)
  
  # 检查是否包含非数值型数据
  if (!is.numeric(data_mat)) {
    warning("数据包含非数值型数据，无法进行PCA分析")
    return(invisible(NULL))
  }
  
  # 检查并处理NA和Inf值
  if (any(is.na(data_mat)) || any(is.infinite(data_mat))) {
    warning("数据包含NA或Inf值，正在清理...")
    data_mat[is.na(data_mat) | is.infinite(data_mat)] <- 0
  }
  
  # 转置数据
  data_frame <- t(data_mat)
  
  # 检查转置后的数据
  if (nrow(data_frame) == 0 || ncol(data_frame) == 0) {
    warning("转置后数据为空，无法进行PCA分析")
    return(invisible(NULL))
  }
  
  # 检查是否至少有两个样本（行）和两个变量（列）
  if (nrow(data_frame) < 2) {
    warning("样本数不足（至少需要2个），无法进行PCA分析")
    return(invisible(NULL))
  }
  
  if (ncol(data_frame) < 2) {
    warning("变量数不足（至少需要2个），无法进行PCA分析")
    return(invisible(NULL))
  }
  
  # 尝试进行PCA分析
  tryCatch({
    gene.pca <- PCA(data_frame, ncp = 2, scale.unit = TRUE, graph = FALSE)
  }, error = function(e) {
    warning("PCA分析失败: ", e$message)
    return(invisible(NULL))
  })
  
  # 检查PCA结果是否有效
  if (is.null(gene.pca) || is.null(gene.pca$var) || is.null(gene.pca$ind)) {
    warning("PCA分析结果无效")
    return(invisible(NULL))
  }
  
  tryCatch({
    pca_component_compound_table <- as.data.frame(gene.pca[['var']][['cor']])
    pca_component_compound_table_file_name <- paste0(
      output_prefix, 
      "PCA_component_compound.xlsx"
    )
    
    write.xlsx(
      pca_component_compound_table, 
      pca_component_compound_table_file_name, 
      sheetName = "Sheet1", 
      rowNames = TRUE
    )
  }, error = function(e) {
    warning("无法写入PCA结果Excel文件: ", e$message)
  })

  tryCatch({
    pca_sample <- data.frame(gene.pca$ind$coord[, 1:2])
    colnames(pca_sample) <- c("Dim.1", "Dim.2")

    pca_eig1 <- round(gene.pca$eig[1, 2], 2)
    pca_eig2 <- round(gene.pca$eig[2, 2], 2)

    # 检查样本文件是否存在
    if (!file.exists(samples_file)) {
      warning("样本文件不存在: ", samples_file, "，跳过分组信息")
      pca_sample$group <- "All"
    } else {
      tryCatch({
        group <- read.delim(
          samples_file, 
          row.names = 2, 
          sep = "\t", 
          check.names = FALSE, 
          header = TRUE
        )
        group <- group[rownames(pca_sample), , drop = FALSE]

        # Ensure group is a factor and add it to pca_sample
        if (nrow(group) > 0 && ncol(group) > 0) {
          pca_sample$group <- factor(group[[1]])
        } else {
          pca_sample$group <- "All"
        }
      }, error = function(e) {
        warning("无法读取样本分组信息: ", e$message, "，使用默认分组")
        pca_sample$group <<- "All"
      })
    }

    p <- ggplot(data = pca_sample, aes(x = Dim.1, y = Dim.2)) +
      geom_point(aes(color = group), size = 5) +
      theme(
        panel.grid = element_blank(),
        panel.background = element_rect(color = "black", fill = "transparent"),
        legend.key = element_rect(fill = "transparent")
      ) +
      labs(
        x = paste("PC1:", pca_eig1, "%"), 
        y = paste("PC2:", pca_eig2, "%"), 
        color = ""
      ) +
      geom_text_repel(aes(label = rownames(pca_sample), color = group), 
                       show.legend = FALSE)

    # 只有当有多个组时才绘制多边形边界
    if (length(unique(pca_sample$group)) > 1) {
      tryCatch({
        cluster_border <- ddply(
          pca_sample, 
          .(group), 
          function(df) df[chull(df$Dim.1, df$Dim.2), ]
        )
        
        p <- p + geom_polygon(
          data = cluster_border, 
          aes(group = group, fill = group, color = group), 
          alpha = 0.3, 
          show.legend = FALSE
        )
      }, error = function(e) {
        warning("无法绘制分组边界: ", e$message)
      })
    }

    ggsave(
      paste0(output_prefix, "PCA_analysis.jpeg"), 
      p, 
      dpi = 300, 
      width = 10, 
      height = 10
    )
  }, error = function(e) {
    warning("无法生成PCA图: ", e$message)
    return(invisible(NULL))
  })
  
  return(invisible(NULL))
}