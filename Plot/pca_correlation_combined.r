suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(FactoMineR))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(corrplot))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(ggplotify))

option_list <- list(
  make_option(c("--data_file"),
    type = "character", default = "fpkm_matrix_filtered.txt",
    help = "fpkm_matrix_filtered.txt", metavar = "character"
  ),
  make_option(c("--samples"),
    type = "character", default = "samples_described.txt",
    help = "提供 samples_described.txt 文件", metavar = "character"
  ),
  make_option(c("--output"),
    type = "character", default = "pca",
    help = "输出图片文件名称", metavar = "character"
  ),
  make_option(c("--log2_transform"),
    type = "logical", default = FALSE,
    help = "是否对输入数据进行 log2 转换，默认不转换", metavar = "logical"
  ),
  make_option(c("--connect_lines"),
    type = "logical", default = FALSE,
    help = "同组样本是否连线（默认不连）", metavar = "logical"
  ),
  make_option(c("--add_background"),
    type = "logical", default = FALSE,
    help = "是否绘制组背景（凸包多边形，默认不画）", metavar = "logical"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Assign the first argument to prefix
data_file <- opt$data_file
data_file_extension <- tools::file_ext(data_file)
samples_file <- opt$samples
output_prefix <- opt$output
log2_transform <- opt$log2_transform
connect_lines <- opt$connect_lines
add_background <- opt$add_background

pca_plot <- function(data_frame, samples_file, output_prefix, log2_transform = FALSE, connect_lines = FALSE, add_background = FALSE) {
  if (log2_transform) {
    data_frame <- log2(data_frame + 1)
  }
  data_frame <- t(data_frame)
  pca_model <- PCA(data_frame, ncp = 2, scale.unit = TRUE, graph = FALSE)

  gene_pca_var_contrib_filename <- paste0(output_prefix, '_PCA_contrib.xlsx')
  gene_pca_var_contrib <- as.data.frame(pca_model[["var"]][["contrib"]])
  write.xlsx(x = gene_pca_var_contrib, file = gene_pca_var_contrib_filename, rowNames=TRUE)

  gene_pca_var_cor_filename <- paste0(output_prefix, '_PCA_cor.xlsx')
  gene_pca_var_cor <- as.data.frame(pca_model[["var"]][["cor"]])
  write.xlsx(x = gene_pca_var_cor, file=gene_pca_var_cor_filename, rowNames=TRUE)

  pca_sample <- data.frame(pca_model$ind$coord[, 1:2])
  colnames(pca_sample) <- c("Dim.1", "Dim.2")

  pca_eig1 <- round(pca_model$eig[1, 2], 2)
  pca_eig2 <- round(pca_model$eig[2, 2], 2)

  group <- read.delim(samples_file, row.names = 2, sep = "\t", check.names = FALSE, header = TRUE)
  group <- group[rownames(pca_sample), , drop = FALSE]

  pca_sample$group <- factor(group[[1]])
  pca_sample$sample_order <- seq_len(nrow(pca_sample))

  p <- ggplot(data = pca_sample, aes(x = Dim.1, y = Dim.2)) +
    geom_point(aes(color = group), size = 5) +
    theme(
      panel.grid = element_blank(),
      panel.background = element_rect(color = "black", fill = "transparent"),
      legend.key = element_rect(fill = "transparent")
    ) +
    labs(x = paste("PCA1:", pca_eig1, "%"), y = paste("PCA2:", pca_eig2, "%"), color = "") +
    geom_text_repel(aes(label = rownames(pca_sample)),
                    max.overlaps = Inf,
                    force = 2,
                    min.segment.length = 0,
                    box.padding = 0.5,
                    point.padding = 0.3)

  if (connect_lines) {
    line_df <- pca_sample[order(pca_sample$group, pca_sample$sample_order), , drop = FALSE]
    p <- p + geom_path(data = line_df, aes(group = group, color = group), linewidth = 0.7, alpha = 0.7, show.legend = FALSE)
  }

  if (add_background) {
    cluster_border <- ddply(pca_sample, .(group), function(df) {
      if (nrow(df) < 3) return(df[0, ])
      df[chull(df$Dim.1, df$Dim.2), ]
    })
    if (nrow(cluster_border) > 0) {
      p <- p + geom_polygon(data = cluster_border, aes(group = group, fill = group), color = "black", alpha = 0.3, show.legend = FALSE)
    }
  }
  return(p)
}


correlation_plot <- function(data_frame) {
  data_frame <- data_frame[, sapply(data_frame, is.numeric), drop = FALSE]
  cor_matrix <- cor(data_frame, use = "pairwise.complete.obs")

  correlation_plot_result <- as.ggplot(function() {
    corrplot(
      cor_matrix,
      type = "lower",
      is.corr = TRUE,
      col = rev(COL2("RdBu")),
      method = "circle",
      addCoef.col = "black",
      tl.col = "black",
      cl.ratio = 0.1
    )
  })

  return(correlation_plot_result)
}

if (data_file_extension == "txt" || data_file_extension == "tsv") {
  data_frame <- read.table(data_file, sep = ifelse(data_file_extension == "tsv", "\t", "\t"),
                           header = TRUE, row.names = 1, check.names = FALSE)
} else if (data_file_extension == "csv") {
  data_frame <- read.table(data_file, sep = ",", header = TRUE, row.names = 1, check.names = FALSE)
} else if (data_file_extension == "xlsx" || data_file_extension == "xls") {
  data_frame <- read_excel(data_file, col_names = TRUE, na = "")
  data_frame <- as.data.frame(data_frame)
  rownames(data_frame) <- data_frame[[1]]
  data_frame <- data_frame[-1]
} else {
  stop("Unsupported file format. Please provide a .txt, .csv, .tsv, .xlsx, or .xls file.")
}

# 绘制 PCA 图
pca_plot_result <- pca_plot(
  data_frame = data_frame,
  samples_file = samples_file,
  output_prefix = output_prefix,
  log2_transform = log2_transform,
  connect_lines = connect_lines,
  add_background = add_background
)

correlation_plot_result <- correlation_plot(data_frame)

ggsave(paste0(output_prefix, '_pca.jpg'), pca_plot_result, dpi = 320, width = 10, height = 10)
ggsave(paste0(output_prefix, '_correlation.jpg'), correlation_plot_result, dpi = 320, width = 10, height = 10)

p_final <- pca_plot_result + correlation_plot_result + plot_layout(widths = c(2, 2))
ggsave(paste0(output_prefix, '_pca_correlation_combined.jpg'), p_final, dpi = 320, width = 18, height = 10)
