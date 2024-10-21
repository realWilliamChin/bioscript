library(pheatmap)
library(openxlsx)
library(ggplot2)

# 获取参数
args <- commandArgs(trailingOnly = TRUE)
data_file <- args[1]
output_pic <- args[2]

reads_data <- read.xlsx(data_file, sheet = 1, rowNames = TRUE)
# 数据预处理
reads_data <- reads_data[rowSums(reads_data != 0) > 0, ]
# 清理每个元素头尾的空格
# reads_data <- apply(reads_data, 2, function(x) trimws(x, which = c("both")))
# 检查是否油 NA，油则退出
if (any(is.na(reads_data))) {
  print("检查数据，有 NA")
  quit()
}

heatmap_plot <- function(data_frame, output_pic) {
  all.heatmap <- pheatmap(data_frame, show_rownames = TRUE, cutree_rows = 5, gaps_row = c(5), cluster_cols = FALSE, scale = "row", cluster_rows = TRUE)
  df_row_cluster = data.frame(cluster = cutree(all.heatmap$tree_row, k=5))
  write.table(df_row_cluster, file="df_row_cluster.txt", sep='\t', row.names=TRUE,col.names = TRUE,quote = FALSE)
  ggsave(output_pic, all.heatmap, dpi = 300, width = 40, height = 40, limitsize = FALSE)
}

heatmap_plot(reads_data, output_pic)