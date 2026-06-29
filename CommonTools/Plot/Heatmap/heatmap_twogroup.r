# setwd("D:/pll/R_work/heatmap/")
library(openxlsx)
library(pheatmap)
library(ggplot2)
library(optparse)
rm(list=ls())

option_list = list(
  make_option(c("-f", "--datafile"), type="character", default=NULL, 
              help="两组之间比对的 fpkm 表达量 excel file", metavar="character"),
  make_option(c("-o", "--outputfile"), type="character", default=NULL, 
              help="输出的图片名称", metavar="character")
) 
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Check that all required arguments are provided
if (is.null(opt$datafile)){
  print_help(opt_parser)
  stop("请输入 datafile 文件", call.=FALSE)
}

datafile=opt$datafile
output=opt$outputfile
# data_file="Ba_vs_BHI_heatmap.xlsx"
data <- read.xlsx(datafile, rowNames = T, sheet = 1)
data_new<-as.data.frame(t(apply(data,1,as.numeric)))
colnames(data_new)<-colnames(data)
data<-data_new
data_samples <- read.xlsx(datafile, rowNames = T, sheet = 2)
#unique(data_samples$group)
data_samples
ann_column<-data_samples[1]
class(ann_column)
#ann_column = data.frame(group=data_samples$group)
#ann_column
#rownames(ann_column)<-rownames(data_samples)
#ann_column
uniq_colors<-unique(data_samples$colors)
uniq_colors
#ann_colors = list(group=c(Ba="blue",BHI="cyan"))
ann_colors = list(group=uniq_colors)
ann_colors
data_samples$group[data_samples$colors==uniq_colors[1]][1]
data_samples$group[data_samples$colors==uniq_colors[2]][1]
color_names<-c(data_samples$group[data_samples$colors==uniq_colors[1]][1],data_samples$group[data_samples$colors==uniq_colors[2]][1])
#ann_colors = list(group=c(data_samples[1,1]="blue",BHI="cyan"))
names(ann_colors[['group']])<-color_names
mode(ann_colors)
#as.data.frame(ann_colors)
#show_rownames <- if (nrow(data) > 100) FALSE else TRUE
#ann_column<-as.data.frame(ann_column)
tmp.plot <- pheatmap(data, cluster_cols = FALSE,
                     clustering_distance_rows = "correlation",
                     scale = "row",
                     colorRampPalette(c("green", "black", "red"))(50),
                     fontsize_row = 6,
                     annotation_col = ann_column,
                     annotation_colors = ann_colors,
                     show_rownames = TRUE,
                     annotation_names_col = FALSE
)
p1_height <- 3 + (nrow(data) / 8)
p1_width <- 3 + ncol(data) * 0.8
if ((p1_width - p1_height) > 5 * p1_height) {
  p1_width <- p1_height * 5
} else if ((p1_height - p1_width) > 3 * p1_width) {
  p1_width <- p1_height / 3
}
ggsave(output,tmp.plot,dpi=320,width=p1_width,height=p1_height,limitsize=FALSE)
#?scale_color_manual
