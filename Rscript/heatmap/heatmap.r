# setwd("D:/pll/R_work/heatmap/")
library(openxlsx)
library(pheatmap)
library(ggplot2)
library(optparse)
rm(list=ls())

option_list = list(
  make_option(c("-f", "--datafile"), type="character", default=NULL, 
              help="两组之间比对的 fpkm 表达量 excel file", metavar="character")
) 
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Check that all required arguments are provided
if (is.null(opt$datafile)){
  print_help(opt_parser)
  stop("请输入 datafile 文件", call.=FALSE)
}

datafile=opt$datafile
# data_file="Ba_vs_BHI_heatmap.xlsx"
data <- read.xlsx(datafile, rowNames = T, sheet = 1)
data_samples <- read.xlsx(datafile, rowNames = T, sheet = 2)
#unique(data_samples$group)
data_samples
ann_column<-data_samples[1]
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
ann_colors

if (nrow(data) > 100) {
  tmp.plot<-pheatmap(data, cluster_col = FALSE,
                     clustering_distance_rows = "correlation",
                     scale="row",
                     colorRampPalette(c("green", "black", "red"))(50),
                     fontsize_row=6,
                     annotation_col = ann_column,
                     annotation_colors = ann_colors,
                     show_rownames = FALSE
  )
} else {
  tmp.plot<-pheatmap(data, cluster_col = FALSE,
                     clustering_distance_rows = "correlation",
                     scale="row",
                     colorRampPalette(c("green", "black", "red"))(50),
                     fontsize_row=6,
                     annotation_col = ann_column,
                     annotation_colors = ann_colors,
                     show_rownames = TRUE
  )
}
ggsave(paste0(datafile,".jpeg"),tmp.plot,dpi=320)
#?scale_color_manual
