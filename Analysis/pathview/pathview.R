# setwd("/home/colddata/qinqiang/Project/2024_03_06_Loropetalum_chinensis/2024_03_13_kegg_pathway")
library(tidyverse)
library(scales)
library(pathview)
library(optparse)
rm(list=ls())
# setwd("/home/colddata/qinqiang/Project/2024_03_06_Loropetalum_chinensis/2024_03_13_kegg_pathway")

option_list = list(
  make_option(c("-r", "--regulation"), type="character", default=NULL, 
              help="kegg_regulation file", metavar="character"),
  make_option(c("-p", "--passedpath"), type="character", default=NULL, 
              help="passed_path 文件", metavar="character")
) 
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Check that all required arguments are provided
if (is.null(opt$regulation)){
  print_help(opt_parser)
  stop("请输入 regulation 文件", call.=FALSE)
}else if (is.null(opt$passedpath)){
  print_help(opt_parser)
  stop("请输入 passed path 文件", call.=FALSE)
}

regulation_file=opt$regulation
passed_path_file=opt$passedpath

kid.dt<-read.table(regulation_file,header=T)
ko.data<-kid.dt[,2]
names(ko.data)<-kid.dt[,1]
head(kid.dt)
path.dt<-read.table(passed_path_file,header=F,sep="\t")
#rownames(path.dt)
head(path.dt)
pathway.names<-path.dt[,1]
#head(pathway.names)
pathway.names<-gsub("ko",'',pathway.names)

# pathview(
#   gene.data = ko.data,
#   pathway.id = pathway.names,
#   species = "ko",
#   out.suffix = "ko.data",
#   kegg.native = T,
#   kegg.dir='/home/colddata/qinqiang/script/Rscript/pathview/ceshi'
# )

node.font <- list(fontname = "Calibri")
edge.font <- list(fontname = "Calibri")
for (pw_id in pathway.names) {
  # 使用 tryCatch 来捕获错误
  tryCatch({
    # 设置输出文件的路径
    # out_file <- file.path(output_dir, paste0(pw_id, ".ko.data.png"))
    
    # 调用 pathview 函数
    pathview(
      gene.data = ko.data,
      pathway.id = pw_id,
      species = "ko",
      out.suffix = "ko.data",
      kegg.native = TRUE,
      node.attrs = node.font,
      edge.attrs = edge.font,
      kegg.dir = '/home/colddata/qinqiang/script/Analysis/pathview/kegg_files'
    )
    
    # 如果成功，可以打印一条消息（可选）
    cat("Successfully plotted pathway:", pw_id, "\n")
    
  }, error = function(e) {
    # 如果出现错误，打印错误消息并跳过该通路
    cat("Error plotting pathway:", pw_id, ":", e$message, "\n")
  })
}



head(path.dt)
colnames(path.dt)<-c('ID','Name')
library(fs)
dir()
for (i in path.dt[,1]) {
  old_name <- paste0(i,'.ko.data.png')
  new_name <- paste0(i,'_',path.dt$Name[path.dt$ID == i],'.ko.data.png')
  
  # 使用tryCatch捕获错误
  result <- tryCatch({
    file.rename(old_name,new_name)
  }, error = function(e) {
    # 打印错误信息，用于调试
    print(paste("重命名文件", old_name, "时出错:", e))
    # 返回一个值，表示操作失败
    return(FALSE)
  })
  
  # 如果重命名成功，result为NULL，否则为FALSE
  if (is.null(result)) {
    print(paste("文件", old_name, "重命名为", new_name, "成功"))
  }
}
# i<-''
#library(fs)
# 移动文件
# dir()
# for (i in grep('.ko.data.png',list.files(),value = T)) {
#   file_move(i,paste0('./CHX_vs_HBX/',i))
# }
#file_move()

#file.rename('ko04928.ko.data.png.png','ko04928.ko.data.png')
#pathview(gene.data = ko.data, pathway.id = "00901", species = "ko", out.suffix = "ko.data", kegg.native = T)
#a<-'temp'
#dir.create(a)
# library(ggplot2)
# head(mpg)
# p<-ggplot(mpg,aes(x=displ,y=hwy))+geom_point()
# p+facet_grid(drv~.)
# p+facet_grid(.~drv)


# df<-data.frame(x=1:10,y=1:10,group=rep(c('a','b'),each=5))
# df
