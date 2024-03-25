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
pathview(gene.data = ko.data, pathway.id = pathway.names, species = "ko", out.suffix = "ko.data", kegg.native = T,
          kegg.dir='/home/colddata/qinqiang/script/Rscript/pathview/kegg_files')

head(path.dt)
colnames(path.dt)<-c('ID','Name')
library(fs)
dir()
for (i in path.dt[,1]) {
  old_name <- paste0(i,'.ko.data.png')
  new_name <- paste0(i,'_',path.dt$Name[path.dt$ID == i],'.ko.data.png')
  file.rename(old_name,new_name)
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
