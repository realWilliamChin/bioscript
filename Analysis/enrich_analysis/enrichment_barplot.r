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
mode(go_data$RichFactor)
go_data$star <- cut(go_data$pvalue, breaks = c(0, 0.001, 0.01, 0.05, Inf), 
    labels = c("***        ", "**     ", "*   ", ""))
go_p<-ggplot(go_data,aes(x=ID, y=RichFactor)) + 
  geom_bar(position="identity",stat="identity",aes(fill=SubOntology), color="black", width = 0.5, height = 0.5)+
  geom_text(aes(label = star), size = 7) +  # 添加星号
  # geom_bar(position="identity",stat="identity",aes(fill=SubOntology))+
  # facet_grid(SubOntology~.,scale="free")+
  theme_minimal()+
  # theme(strip.text=element_text(face="bold",size=11),
  #       axis.text=element_text(size=9),
  #       strip.background = element_rect(fill="lightgrey",linewidth=1))+
  xlab(ontology_value)+
  coord_flip() +
  ggtitle(plot_title) +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 14, face = "bold"),
        legend.title = element_text(face = "bold")) +
  labs(fill=ontology_value)

go_p
# height 根据 datafile 的行数调整
go_p_height <- 2 + 0.4 * nrow(go_data)
# Add a minimum height for plots with few bars
if (nrow(go_data) < 5) {
  go_p_height <- 6 # Set a minimum height
}
ggsave(output_file,go_p,height = go_p_height,width=14,dpi=320)
