library(ggthemes)
library(ggplot2)
library(pheatmap)
library(reshape2)
library(ggcorrplot)
library(corrplot)
library(FactoMineR)
library(plyr)
library(ropls)
library(ggrepel)
library(mixOmics)
library(openxlsx)
library(optparse)
rm(list = ls())


heatmap_plot <- function(data_frame, output_pic, log2data = FALSE, log2pic_fname = NA, ...) {
  data_frame <- data_frame[rowSums(data_frame != 0) > 0, ]
  all.heatmap <- pheatmap(
    data_frame,
    show_rownames = nrow(data_frame) <= 100,
    encoding = "UTF-8",
    cluster_cols = FALSE,
    scale = 'row',
    ...
  )
  ggsave(output_pic, all.heatmap, dpi = 300, width = 20, height = 10, limitsize = FALSE)

  if (log2data) {
    data_frame <- log2(data_frame + 1)
    heatmap_plot(data_frame, log2pic_fname, log2data = FALSE, log2pic_fname = NA, ...)
  }
}

correlation_plot <- function(data_frame, output_dir) {
  data_frame <- na.omit(data_frame)
  data_frame <- data_frame + 0.000000001
  fpkm.m <- as.matrix(data_frame)
  fpkm.cor <- cor(fpkm.m)
  correlation_df <- as.data.frame(fpkm.cor)
  correlation_df$ID <- rownames(correlation_df)
  correlation_df <- correlation_df[, c("ID", setdiff(names(correlation_df), "ID"))]
  write.xlsx(correlation_df, paste(output_dir, "Metabolite_correlation.xlsx", sep = "/"), sheetName = "Sheet1", rowNames = FALSE)
  min(fpkm.cor)
  # 根据样本数量，设置图的大小
  sample_num <- ncol(fpkm.cor)
  correlation_plot_width <- 5 + sample_num * 0.5
  correlation_plot_height <- 4 + sample_num * 0.5
  png(paste(output_dir, "Metabolite_correlation_graph.png", sep = "/"), width = correlation_plot_width, height = correlation_plot_height, units = "in", res = 300)
  corrplot(fpkm.cor, is.corr = F, col = rev(COL2("PiYG")), method = "color", addCoef.col = "black", tl.col = "black", col.lim = c(min(fpkm.cor) - 0.01, max(fpkm.cor)), cl.ratio = 0.1)
  dev.off()
}

pca_plot <- function(reads_data_frame = NA, fpkm_data_frame = NA, samples_file, output_dir) {
  if (is.data.frame(fpkm_data_frame)) {
    data_frame <- log2(fpkm_data_frame + 1)
  } else {
    data_frame <- reads_data_frame
  }

  data_frame <- t(data_frame)
  gene.pca <- PCA(data_frame, ncp = 2, scale.unit = TRUE, graph = FALSE)

  pca_sample <- data.frame(gene.pca$ind$coord[, 1:2])
  colnames(pca_sample) <- c("Dim.1", "Dim.2")

  pca_eig1 <- round(gene.pca$eig[1, 2], 2)
  pca_eig2 <- round(gene.pca$eig[2, 2], 2)

  group <- read.delim(samples_file, row.names = 2, sep = "\t", check.names = FALSE, header = TRUE)
  group <- group[rownames(pca_sample), , drop = FALSE]

  # Ensure group is a factor and add it to pca_sample
  pca_sample$group <- factor(group[[1]])

  p <- ggplot(data = pca_sample, aes(x = Dim.1, y = Dim.2)) +
    geom_point(aes(color = group), size = 5) +
    theme(
      panel.grid = element_blank(),
      panel.background = element_rect(color = "black", fill = "transparent"),
      legend.key = element_rect(fill = "transparent")
    ) +
    labs(x = paste("PCA1:", pca_eig1, "%"), y = paste("PCA2:", pca_eig2, "%"), color = "") +
    geom_text_repel(aes(label = rownames(pca_sample)))

  cluster_border <- ddply(pca_sample, .(group), function(df) df[chull(df$Dim.1, df$Dim.2), ])
  p <- p + geom_polygon(data = cluster_border, aes(group = group, fill = group), color = "black", alpha = 0.3, show.legend = FALSE)

  ggsave(paste(output_dir, "Metabolite_PCA_analysis.jpeg", sep = "/"), p, dpi = 300, width = 10, height = 10)
}

metabolite_analysis <- function(samples_file, reads_data_frame, fpkm_data_frame = NA, definition_df = NA, output_dir, log2data = FALSE) {
  # >>>>>> 多组分析中有多个组的
  # select_sample_info <- read.table("multigroup4_samples_described.txt",sep="\t",header=T,check.names=F,stringsAsFactors = F)
  select_sample_info <- read.table(samples_file, sep = "\t", header = T, check.names = F, stringsAsFactors = F)

  groups <- select_sample_info$group

  if (is.data.frame(fpkm_data_frame)) {
    print("plsda 使用归一化数据")
    select_fpkm <- fpkm[, select_sample_info$sample, drop = FALSE]
    fpkm_t <- t(select_fpkm)
    select_data_frame <- fpkm[, select_sample_info$sample, drop = FALSE]
    metabolites <- as.matrix(fpkm_t)
  } else {
    select_reads <- reads_data_frame[, select_sample_info$sample, drop = FALSE]
    select_reads <- select_reads[rowSums(select_reads != 0) > 0, ]
    reads_data_t <- t(select_reads)
    select_data_frame <- reads_data_frame[, select_sample_info$sample, drop = FALSE]
    metabolites <- as.matrix(reads_data_t)
  }

  # mutligroup samples 的时候会出现每行都为 0 跑不了的情况
  select_data_frame <- select_data_frame[rowSums(select_data_frame != 0) > 0, ]
  metabolites <- metabolites[rowSums(metabolites != 0) > 0, ]
  
  heatmap_plot(
    select_data_frame, paste0(output_dir, '/MultiGroup_heatmap.jpeg'),
    log2data = log2data ,log2pic_fname = paste0(output_dir, '/MultiGroup_log_heatmap.jpeg')
    )

  # 如果代谢物的数量小于 10 用 4，10-20 用 6，20 以上用 10
  if (ncol(metabolites) < 10) {
    ncomp <- 4
  } else if (ncol(metabolites) < 20) {
    ncomp <- 6
  } else {
    ncomp <- 10
  }

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

  ggsave(paste0(output_dir, "/Metabolite_quantitation_scree_plot.jpeg"), scree_plot, width = ncomp * 0.7, height = 4)

  comp_load_df <- as.data.frame(df_plsda$loadings$X)
  comp_load_df <- cbind(rownames(comp_load_df), comp_load_df)
  colnames(comp_load_df)[1] <- "compound_name"
  # write.table(comp_load_df, file="多组分析/pc_loading_value.txt", sep='\t', row.names=FALSE,col.names = TRUE,quote = FALSE)
  # write.xlsx(comp_load_df, file=paste0(multigroup_dir, "/pc_loading_value.xlsx"), sheetName = "Sheet1", rowNames = FALSE)
  write.xlsx(comp_load_df, file = paste0(output_dir, "/pc_loading_value.xlsx"), sheetName = "Sheet1", rowNames = FALSE)

  df <- unclass(df_plsda)

  df1 <- as.data.frame(df$variates$X)
  # df1$group = sample_info$group
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
    geom_point(size = 1.8) +
    geom_text_repel(
      aes(label = samples),
      size = 2,
      box.padding = 0.35,
      point.padding = 0.5
    ) +
    labs(x = paste0("P1 (", x_label * 100, "%)"), y = paste0("P2 (", y_label * 100, "%)")) +
    #stat_ellipse(
    #  data = df1, geom = "polygon", level = 0.95,
    #  linetype = 2, linewidth = 0.5, aes(fill = group),
    #  alpha = 0.2, show.legend = TRUE
    #) +
    # scale_color_discrete() + # 这个有时候画图会出错
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

  # 根据当前值每个增加2
  x_min <- x_range[1] - 2

  x_max <- x_range[2] + 2
  y_min <- y_range[1] # - 1
  y_max <- y_range[2] # + 1
  p1 <- p1 + coord_cartesian(xlim = c(x_min, x_max), ylim = c(y_min, y_max))

  # ggsave(paste0(multigroup_dir, "/Multigroup_Plsda_Distribution_Graph.jpeg"), p1, width = plot_width, height = plot_height)
  ggsave(paste0(output_dir, "/Multigroup_Plsda_Distribution_Graph.jpeg"), p1, width = plot_width, height = plot_height)


  # vip 值使用 fpkm, 其他值计算都是用 reads
  # predI 这里不知道为什么有时候需要写 NA 才行
  plsda_model <- opls(x = metabolites, y = groups, predI = NA, crossvalI = 5)
  vip_values <- plsda_model@vipVn
  df.vip <- as.data.frame(vip_values)

  reads_data_with_def <- cbind(select_data_frame, VIP = df.vip$vip_values)

  reads_data_with_def$Metabolite <- rownames(reads_data_with_def)
  reads_data_with_def <- reads_data_with_def[, c("Metabolite", setdiff(names(reads_data_with_def), "Metabolite"))]

  # 没有定义跳过
  class_count <- ""
  if (is.data.frame(definition_df)) {
    reads_data_with_def <- merge(reads_data_with_def, definition_df, by = "Metabolite", all.x = TRUE)
    greater_than_one_data_def <- reads_data_with_def[reads_data_with_def$VIP > 1, ]
    # 如果 greater_than_one_data_def 没有数据，就是用 reads_data_with_def
    if (nrow(greater_than_one_data_def) == 0) {
      print("没有 VIP 大于 1，输出全部")
      greater_than_one_data_def <- reads_data_with_def
    }

    # 添加 class count 和对应的名称
    # class_count <- aggregate(greater_than_one_data_def$Metabolite, by = list(greater_than_one_data_def$Class), length)
    if ("Class" %in% colnames(greater_than_one_data_def)) {
      class_count <- greater_than_one_data_def %>%
      group_by(Class) %>%
      summarize(
        count = n(),
        compounds = paste(Metabolite, collapse = ", ")
      )
      class_count <- class_count[class_count$Class != "", ]
      class_count <- class_count[order(-class_count$count), ]
      write.xlsx(class_count,
        file = paste0(output_dir, "/Significant_compound_count_by_class.xlsx"),
        sheetName = "Sheet1", rowNames = FALSE, colNames = TRUE
      )
    }

    if ("Subclass" %in% colnames(greater_than_one_data_def)) {
      class_count <- greater_than_one_data_def %>%
      group_by(Subclass) %>%
      summarize(
        count = n(),
        compounds = paste(Metabolite, collapse = ", ")
      )
      class_count <- class_count[class_count$Subclass != "", ]
      class_count <- class_count[order(-class_count$count), ]
      write.xlsx(class_count,
        file = paste0(output_dir, "/Significant_compound_count_by_subclass.xlsx"),
        sheetName = "Sheet1", rowNames = FALSE, colNames = TRUE
      )
    }
    
  }

  reads_data_with_def <- reads_data_with_def[order(-reads_data_with_def$VIP, na.last = TRUE), ]
  write.xlsx(reads_data_with_def, file = paste0(output_dir, "/Metabolite_quantitation_VIP.xlsx"), sheetName = "Sheet1", rowNames = FALSE, colNames = TRUE)
}


zujianfenxi <- function(compare_file, samples_file, reads_data_frame, fpkm_data_frame=NA, definition_df=NA, output_dir, log2data=FALSE) {
  sample_info <- read.table(
    samples_file, sep = "\t", header = T, check.names = F, stringsAsFactors = F)
  comp_info <- read.table(
    compare_file, sep = "\t", header = T, check.names = F, stringsAsFactors = F)
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
    group = character(0), All = numeric(0), Up = numeric(0), Down = numeric(0),
    Up_list = character(0), Down_list = character(0)
  )
  # 循环中注意可能需要修改 corssvalI 值
  for (i in seq_along(comparisons)) {
    compare_name <- paste(comparisons[[i]], collapse = "_vs_")
    compare_path <- paste(output_dir, compare_name, sep = "/")
    dir.create(compare_path, showWarnings = FALSE)
    print(paste("Processing:", compare_name)) # 打印当前处理的组合名称
    groups <- comparisons[[i]]
    groupA_cols <- sample_info$sample[sample_info$group == groups[1]]
    groupB_cols <- sample_info$sample[sample_info$group == groups[2]]

    # 获取当前两组的样本的表达量数据
    current_samples <- sample_info$sample[sample_info$group %in% groups]

    if (is.data.frame(fpkm_data_frame)) {
      
      current_expression_fpkm_data <- fpkm[, current_samples, drop = FALSE]
      rows_to_remove <- apply(current_expression_fpkm_data, 1, function(row) all(row == 0))
      current_expression_fpkm_data <- current_expression_fpkm_data[!rows_to_remove, ]
    }
    current_expression_data <- reads_data_frame[, current_samples, drop = FALSE]
    rows_to_remove <- apply(current_expression_data, 1, function(row) all(row == 0))
    current_expression_data <- current_expression_data[!rows_to_remove, ]
    
    # p_value_current_expression_data <- reads_data[, current_samples, drop = FALSE]
    # 检查因变量的水平数量
    y_factor <- factor(sample_info$group[sample_info$sample %in% current_samples])
    if (length(unique(y_factor)) < 2) {
      print(paste("y_factor 小于 2 个", compare_name))
      next # 跳过当前的迭代
    }

    # 检查是否有缺失值
    if (any(is.na(current_expression_data))) {
      print(paste("Missing values found in expression data for:", compare_name))
      next # 跳过当前的迭代
    }

    # 进行OPLS分析
    if (is.data.frame(fpkm_data_frame)) {
      transposed_expression_data <- t(current_expression_fpkm_data)
    } else {
      transposed_expression_data <- t(current_expression_data)
    }
    crossval_number <- length(row_number(transposed_expression_data)) - 1
    opls_model <- try(
      # crossvalI 默认是 7, crossvalI 需要小于等于两组样本的数量
      opls(x = transposed_expression_data, y = y_factor, predI = 1, orthoI = 2, crossvalI = crossval_number),
      silent = FALSE
    )

    # 检查模型是否成功建立
    if (class(opls_model) == "try-error") {
      print(paste("Failed to build model for:", compare_name))
      next # 跳过当前的迭代
    }

    # 获取VIP值，如果模型失败则赋予 0
    vip_values <- tryCatch(
      {
        opls_model@vipVn
      },
      error = function(e) {
        rep(0, ncol(current_expression_data))
      }
    )

    
    # 组间 heatmap 图
    if (is.data.frame(fpkm_data_frame)) {
      heatmap_plot(current_expression_fpkm_data, paste0(compare_path, "/", compare_name, "_heatmap.jpeg"),
      log2data = log2data, log2pic_fname = paste0(compare_path, "/", compare_name, "_log_heatmap.jpeg"))
    } else {
      heatmap_plot(current_expression_data, paste0(compare_path, "/", compare_name, "_heatmap.jpeg"),
      log2data = log2data, log2pic_fname = paste0(compare_path, "/", compare_name, "_log_heatmap.jpeg"))
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

    # 使用apply函数
    p_values <- apply(current_expression_data, 1, function(row) paired_t_test(row[groupA_cols], row[groupB_cols]))
    # padj <- p.adjust(p_values, "BH")
    # 将 VIP 值和 FoldChange 保存到文件
    current_expression_data_def <- cbind(current_expression_data, baseMeanA, baseMeanB, FoldChange, pvalues = p_values)
    current_expression_data_def <- as.data.frame(current_expression_data_def)

    # 创建一个长度与current_expression_data行数相同的全零向量
    vip_values_adjusted <- rep(0, nrow(current_expression_data_def))
    # 更新vip_values_adjusted向量中相应的位置
    # 这假设vip_values的名称与current_expression_data的行名对应
    if (!is.null(names(vip_values))) {
      matching_rows <- match(names(vip_values), rownames(current_expression_data_def))
      vip_values_adjusted[matching_rows] <- vip_values
    }
    current_expression_data_def <- cbind(current_expression_data_def, VIP = vip_values_adjusted)

    current_expression_data_def <- current_expression_data_def[order(current_expression_data_def$pvalues, na.last = TRUE), ]
    current_expression_data_def$Metabolite <- rownames(current_expression_data_def)
    current_expression_data_def <- current_expression_data_def[, c("Metabolite", setdiff(names(current_expression_data_def), "Metabolite"))]
    current_expression_data_def$padj <- p.adjust(current_expression_data_def$pvalues, "BH")

    if (is.data.frame(definition_df)) {
      current_expression_data_def <- merge(current_expression_data_def, definition_df, by = "Metabolite", all.x = TRUE)
    }

    current_expression_data_def <- current_expression_data_def[order(-current_expression_data_def$VIP, na.last = TRUE), ]
    current_expression_data_def[is.na(current_expression_data_def)] <- 0

    # 将更新后的数据框保存为文本文件
    # write.table(current_expression_data_def, file = paste0('组间分析/', compare_name, '/', compare_name, '_VIP.txt'), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    write.xlsx(current_expression_data_def, file = paste0(compare_path, "/", compare_name, "_VIP.xlsx"), sheetName = "Sheet1", rowNames = FALSE, colNames = TRUE)

    # 计算 deg
    deg_df <- current_expression_data_def[current_expression_data_def$VIP > 1, ]
    deg_up <- nrow(deg_df[deg_df$FoldChange >= 1.2, ])
    deg_up_idlist <- deg_df[deg_df$FoldChange >= 1.2, ]$Metabolite
    deg_down <- nrow(deg_df[deg_df$FoldChange <= 0.8, ])
    deg_down_idlist <- deg_df[deg_df$FoldChange <= 0.8, ]$Metabolite
    new_row <- data.frame(
      group = compare_name, All = deg_up + deg_down, Up = deg_up, Down = deg_down,
      Up_list = paste0(deg_up_idlist, collapse = ","), Down_list = paste0(deg_down_idlist, collapse = ",")
    )
    deg_data <- rbind(deg_data, new_row)

    # 生成图形
    for (plot_type in plot_types) {
      # 创建文件名
      file_name <- paste0(compare_path, "/", compare_name, "_OPLS_DA_", plot_type, ".png")

      # 检查模型是否有效
      if (inherits(opls_model, "opls")) {
        png(file_name)
        try({plot(opls_model, typeVc = plot_type)},silent = TRUE)
        dev.off() # 确保关闭图形设备
      }
    }
  }
  write.xlsx(deg_data, file = paste0(output_dir, "/Differential_metabolite_count_summary.xlsx"), sheetName = "Sheet1", rowNames = FALSE, colNames = TRUE)

  # class 分组计数
  #if (is.data.frame(definition_df)) {
  all_files <- list.files(path = output_dir, recursive = TRUE)
  vip_files <- all_files[grep("VIP", all_files)]

  class_count_list <- list()
  subclass_count_list <- list()

  for (vip_file in vip_files) {
    print(vip_file)
    compare_group_name = gsub("_VIP.xlsx", "", basename(vip_file))
    
    vip_df <- read.xlsx(paste(output_dir, vip_file, sep = "/"), sheet = 1, rowNames = FALSE)

    # volcano$regulation <- as.factor(ifelse(volcano_filter_col < filter_value & abs(volcano$log2FoldChange) >= bs_pos, ifelse(volcano$log2FoldChange >= bs_pos, "Up", "Down"), "NoSignificant"))
    vip_df$regulation <- as.factor(ifelse(vip_df$VIP > 1 & (vip_df$FoldChange > 1.2 | vip_df$FoldChange < 0.8), ifelse(vip_df$FoldChange >= 1.2, "Up", "Down"), "NoSignificant"))
    write.xlsx(vip_df, file=paste0(output_dir,'/',vip_file), sheetName = "Sheet1", rowNames = FALSE, colNames = TRUE)
    
    # ======= enrich start ======== 有 KEGG 才能做 enrich
    if ("KEGG" %in% colnames(vip_df)) {
      up_df <- vip_df$KEGG[vip_df$regulation == "Up"]
      down_df <- vip_df$KEGG[vip_df$regulation == 'Down']
      
      dir.create(paste0('05_Enrich/', compare_group_name))
      enrich_output_prefix <- paste0('05_Enrich/', compare_group_name, '/', compare_group_name)
      up_df_filename <- paste0(enrich_output_prefix, '_Up_Compound_ID.txt')
      down_df_filename <- paste0(enrich_output_prefix, '_Down_Compound_ID.txt')
      write.table(up_df, file = up_df_filename, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
      write.table(down_df, file = down_df_filename, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
      
      script_path = "/home/colddata/qinqiang/script/MetaboliteAnalysis/MetaboliteEnrich/metabolite_enrich.r"
      up_df_cmd = paste0("/opt/biosoft/R-4.2.2/bin/Rscript ", script_path, " --datatable ", up_df_filename, " --outputprefix ", enrich_output_prefix, "_Up")
      down_df_cmd = paste0("/opt/biosoft/R-4.2.2/bin/Rscript ", script_path, " --datatable ", down_df_filename, " --outputprefix ", enrich_output_prefix, "_Down")
      Sys.setenv(R_LIBS="/opt/biosoft/R-4.2.2/lib64/R/library")
      system(up_df_cmd)
      system(down_df_cmd)
    }
    # ======= enrich end =======

    vip_df_vipgt1 <- vip_df[vip_df$VIP > 1 & (vip_df$FoldChange > 1.2 | vip_df$FoldChange < 0.8), ]

    if ("Class" %in% colnames(vip_df_vipgt1)) {
      class_count <- aggregate(Metabolite ~ Class, data = vip_df_vipgt1, FUN = length)
      class_count <- class_count[class_count$Class != "", ]
      names(class_count)[2] <- strsplit(vip_file, "/")[[1]][1] # 修改列名为对应的 count_name
      class_count_list <- append(class_count_list, list(class_count))
    }

    if ("Subclass" %in% colnames(vip_df_vipgt1)) {
      subclass_count <- aggregate(Metabolite ~ Subclass, data = vip_df_vipgt1, FUN = length)
      subclass_count <- subclass_count[subclass_count$Subclass != "", ]
      names(subclass_count)[2] <- strsplit(vip_file, "/")[[1]][1] # 修改列名为对应的 count_name
      subclass_count_list <- append(subclass_count_list, list(subclass_count))
    }
  }

  # 去除重复的行?
  # class_count_list <- lapply(class_count_list, function(df) df[!duplicated(df$class), ])

  # 合并 class_count_list 中所有的 class_count，根据第一列的 class 合并，合并方式为并集
  if (length(class_count_list) != 0) {
    class_count_result <- Reduce(function(x, y) merge(x, y, by = "Class", all = TRUE), class_count_list)
    class_count_result[is.na(class_count_result)] <- 0
    write.xlsx(class_count_result,
      file = paste0(output_dir, "/Significant_compound_count_by_class.xlsx"),
      sheetName = "Sheet1", rowNames = FALSE, colNames = TRUE
    )
  }

  if (length(subclass_count_list) != 0) {
    subclass_count_result <- Reduce(function(x, y) merge(x, y, by = "Subclass", all = TRUE), subclass_count_list)
    subclass_count_result[is.na(subclass_count_result)] <- 0
    write.xlsx(subclass_count_result,
      file = paste0(output_dir, "/Significant_compound_count_by_subclass.xlsx"),
      sheetName = "Sheet1", rowNames = FALSE, colNames = TRUE
    )
  }
}


# =========== READ DATA =============

zhengtifenxi_dir <- "03_代谢物整体分析/"
chayifenxi_dir <- "04_代谢物差异分析/"
zujianfenxi_dir <- paste(chayifenxi_dir, "组间分析", sep = "")

dir.create("00_Background_materials")
dir.create("01_原始质谱数据")
dir.create("02_代谢物定量图数据")
dir.create(zhengtifenxi_dir)
dir.create(chayifenxi_dir)
dir.create(zujianfenxi_dir)

# 读取文件
reads_data <- read.xlsx("All_sample_data.xlsx", sheet = 1, rowNames = TRUE)
sample_info <- read.table("samples_described.txt", sep = "\t", header = T, check.names = F, stringsAsFactors = F)

# 按照 sample 样本列进行排序
reads_data <- reads_data[sample_info$sample]
reads_data <- reads_data[rowSums(reads_data != 0) > 0, ]
# 清理每个元素头尾的空格
# reads_data <- apply(reads_data, 2, function(x) trimws(x, which = c("both")))
if (any(is.na(reads_data))) {
  print("检查数据，有 NA")
  quit()
}

Compound_def_file <- 'Compound_def.xlsx'
# 如果需要合并定义则读取单独定义文件，没有则跳过
if (file.exists(Compound_def_file)) {
  definition_df <- read.xlsx(Compound_def_file, sheet = 1, rowNames = TRUE)
  definition_df$Metabolite <- rownames(definition_df)
}

# ============ RUN ================
# definition_df = NA 或者 definition_df = definition_df
# 整体分析
zhengtifenxi_heatmap_pic_name <- paste(zhengtifenxi_dir, "All_metabolites_heatmap.jpeg", sep = "/")
zhengtifenxi_log2heatmap_pic_name <- paste(zhengtifenxi_dir, "All_metabolites_log_heatmap.jpeg", sep = "/")
heatmap_plot(reads_data, zhengtifenxi_heatmap_pic_name, log2data=FALSE, log2pic_fname=zhengtifenxi_log2heatmap_pic_name)
correlation_plot(reads_data, zhengtifenxi_dir)
pca_plot(reads_data_frame = reads_data, samples_file = "samples_described.txt", output_dir = zhengtifenxi_dir)
# 多组分析
multigroup_dir <- paste(chayifenxi_dir, "多组分析3", sep = "/")
dir.create(multigroup_dir)
metabolite_analysis(samples_file="multigroup_samples_described_3.txt", reads_data, definition_df=definition_df, output_dir=multigroup_dir)
# 组间分析
zujianfenxi(compare_file = "compare_info.txt", samples_file = "samples_described.txt",
            reads_data_frame = reads_data, fpkm_data_frame = FALSE,
            definition_df=definition_df, output_dir = zujianfenxi_dir, log2data=FALSE)


# ========== z-score RUN ============
fpkm <- as.data.frame(t(apply(reads_data, 1, function(x) {
  (x - mean(x)) / (sd(x)**0.5)
})))
heatmap_plot(data_frame = fpkm, output_pic = paste(zhengtifenxi_dir, "All_metabolites_heatmap.jpeg", sep = "/"),log2data = FALSE)
correlation_plot(data_frame = fpkm, output_dir = zhengtifenxi_dir)
pca_plot(reads_data_frame = FALSE, fpkm_data_frame = fpkm, samples_file = "samples_described.txt", output_dir = zhengtifenxi_dir)
# 多组分析
multigroup_dir <- paste(chayifenxi_dir, "多组分析", sep = "/")
dir.create(multigroup_dir)
metabolite_analysis(samples_file="samples_described.txt", reads_data_frame=NA, fpkm_data_frame=fpkm, definition_df=definition_df, output_dir = multigroup_dir)
zujianfenxi(compare_file="compare_info.txt", samples_file="samples_described.txt", reads_data_frame=reads_data, fpkm_data_frame=fpkm, output_dir=zujianfenxi_dir, log2data=FALSE)
