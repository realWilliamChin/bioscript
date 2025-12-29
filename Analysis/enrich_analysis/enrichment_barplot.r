# setwd("d:/pll/R_work/clusterProfile/")
library(ggplot2)
library(openxlsx)
library(ggthemes)
library(optparse)
rm(list=ls())

option_list = list(
  make_option(c("-f", "--datafile"), type="character", default=NULL, 
              help="输入文件，输入列包括 ID，Ontology，SubOntology，RichFactor", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="输出文件", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Check that all required arguments are provided
if (is.null(opt$datafile)){
  print_help(opt_parser)
  stop("请输入 datafile 文件", call.=FALSE)
}else if (is.null(opt$output)){
  print_help(opt_parser)
  stop("请输入 output 文件名", call.=FALSE)
}

datafile=opt$datafile
output_file=opt$output

plot_title <- gsub(".xlsx", "", datafile)
plot_title <- basename(plot_title)
go_data<-read.xlsx(datafile)
ontology_value <- go_data[1,4]
go_data<-as.data.frame(go_data)
go_data<-go_data[order(go_data$SubOntology),]
go_data$ID<-factor(go_data$ID,levels=go_data$ID)

# 确保数值列是正确的类型
go_data$RichFactor <- as.numeric(go_data$RichFactor)
go_data$pvalue <- as.numeric(go_data$pvalue)

# 检查pvalue是否为数值类型
if (!is.numeric(go_data$pvalue)) {
  stop("pvalue列必须是数值类型，请检查输入数据")
}

mode(go_data$RichFactor)
go_data$star <- cut(go_data$pvalue, breaks = c(0, 0.001, 0.01, 0.05, Inf),
    labels = c("***        ", "**     ", "*   ", ""))
go_p <- ggplot(go_data, aes(x=factor(ID, levels=ID), y=RichFactor, fill=SubOntology)) + 
  geom_bar(position="identity",stat="identity")+
  geom_text(aes(label = star), position = position_dodge(width = 1), size = 7) +  # 添加星号
  # geom_bar(position="identity",stat="identity",aes(fill=SubOntology))+
  # facet_grid(SubOntology~.,scale="free")+
  theme_bw()+
  # theme(strip.text=element_text(face="bold",size=11),
  #       axis.text=element_text(size=9),
  #       strip.background = element_rect(fill="lightgrey",linewidth=1))+
  xlab(ontology_value)+
  ylab("RichFactor")+
  ggtitle(plot_title) +
  coord_flip() +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 18, face="bold"))
  # labs(fill=ontology_value)

# height 根据 datafile 的行数调整
go_p_height <- 2 + 0.4 * nrow(go_data)
# Add a minimum height for plots with few bars
if (nrow(go_data) < 5) {
  go_p_height <- 6 # Set a minimum height
}
ggsave(output_file,go_p,height = go_p_height,width=14,dpi=320)
