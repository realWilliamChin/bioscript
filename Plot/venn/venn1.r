library(eulerr)
library(purrr)  
venn_dat  <- read.delim("/home/colddata/qinqiang/script/Rscript/data2.txt")
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
write.table(commonid, "commonid.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
venn.title <- paste(colnames(venn_dat),collapse = "_vs_")

png("venn.png",width = 1000, height = 1000,res=320)
plot(euler(
  venn_list,
  shape = "circle"),                    # 图案的形状，椭圆ellipse 或圆circle
  # quantities = list(type = c("percent","counts"),cex=1),          # 显示类型，百分比和数字，数字大小
  quantities = list(type = c("counts"),cex=1),          # 显示类型，百分比和数字，数字大小
  # labels=list(cex=0),                   # 组名标签的大小
  edges = list(col = "black", lex = 1), # 图形边缘的颜色和大小
  fills = list(fill = c("#1E90FF", "#FF8C00", "#80FF00"),alpha=0.7), # 填充的颜色和透明度
  legend = list(side = "right",cex=1.0),
  main=venn.title,
  # 图例的位置
  # factor_names = FALSE
)
dev.off()








