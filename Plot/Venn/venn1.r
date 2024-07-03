library(eulerr)
library(purrr)
args=commandArgs(T)
input_file <- args[1]
output_file <- args[2]
venn_dat  <- read.delim(input_file)
head(venn_dat)
colnames(venn_dat)
venn_list <- as.list(venn_dat)
venn_list <- purrr::map(venn_list, na.omit)            # 删除列表中每个向量中的NA
venn_list <- lapply(venn_list, function(x) x[x != ""]) # 删除列表中每个向量中的""空字符串
venn_list <- lapply(venn_list, unique)                 # 移除重复元素
# if (ncol(venn_dat) == 3) {
#   venn.title<-paste0(colnames(venn_dat)[1],"_vs_",colnames(venn_dat)[2],"_vs_",colnames(venn_dat)[3])
#   commonid = intersect(intersect(venn_list[[1]],venn_list[[2]]),venn_list[[3]])
# }else if (ncol(venn_dat) == 2) {
#   venn.title<-paste0(colnames(venn_dat)[1],"_vs_",colnames(venn_dat)[2])
#   commonid = intersect(venn_list[[1]],venn_list[[2]])
# }

commonid <- reduce(venn_list, intersect)
common_id_outfile <- gsub(".jpeg", "_common_id.txt", output_file)
write.table(commonid, common_id_outfile, quote = FALSE, row.names = FALSE, col.names = FALSE)
venn.title <- paste(colnames(venn_dat),collapse = " vs ")

png(output_file, width = 1000, height = 1000, res = 320)
plot(euler(
  venn_list,
  shape = "circle"),                    # 图案的形状
  quantities = list(type = c("counts"), cex = 1.2),  # 显示类型和数字大小
  # 默认情况下不显示组名标签，如果需要可以取消注释并调整大小
  # labels = list(cex = 1.0),            # 组名标签的大小
  edges = list(col = "black", lwd = 2),  # 图形边缘的颜色和线宽
  fills = list(fill = c("#1E90FF", "#FF8C00", "#80FF00"), alpha = 0.7),  # 填充的颜色和透明度
  legend = list(side = "right", cex = 1.0),  # 图例的位置和字体大小
  main = venn.title,  # 图形标题
  # factor_names = FALSE  # 如果需要隐藏组名，可以取消注释  
)

# 关闭图形设备
dev.off()








