library(ggplot2)
library(ggsignif)
library(ggpubr)
library(RColorBrewer)
library(openxlsx)

setwd("C:/Users/Analysis/OneDrive/Work/zhiwusuo/00_Script/02_Rscript/practice")
file='alpha_diversity.txt'
plot_data <- read.csv(file, sep='\t')
plot_data
y_var <- 'Insimpson'
p <- ggplot(data=plot_data)+
  geom_boxplot(mapping=aes(x=group,y=Insimpson,colour=group),
               alpha = 0.5,
               size = 1.5,
               width = 0.6)+
  geom_jitter(mapping=aes(x=group,y=Insimpson,colour=group),
              alpha = 0.3, size = 3)+
  #scale_color_manual(limits=c("A", "F"),
   #                  values=c("#85B22E", "#5F80B4"))+
  geom_signif(mapping=aes(x=group,y=Insimpson),
              comparisons = list(c("MS", "GZ")),
              map_signif_level=T,
              tip_length=c(0,0,0,0,0,0,0,0,0,0,0,0),
              #y_position = c(40,41,42,39,38,40),
              size=1,
              textsize=4,
              test = "t.test")+
  theme_classic(
    base_line_size = 1
  )+
  labs(title=y_var, x="group", y=y_var)+
  theme(plot.title = element_text(size = 15,
                                  colour = "black",
                                  hjust = 0.5),
        axis.title.y = element_text(size = 15, 
                                    # family = "myFont", 
                                    color = "black",
                                    face = "bold", 
                                    vjust = 1.9, 
                                    hjust = 0.5, 
                                    angle = 90),
        legend.title = element_text(color="black", # 修改图例的标题
                                    size=15, 
                                    face="bold"),
        legend.text = element_text(color="black", # 设置图例标签文字
                                   size = 10, 
                                   face = "bold"),
        axis.text.x = element_text(size = 13,  # 修改X轴上字体大小，
                                   color = "black", # 颜色
                                   face = "bold", #  face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
                                   vjust = 0.5, # 位置
                                   hjust = 0.5, 
                                   angle = 0), #角度
        axis.text.y = element_text(size = 13,  
                                   color = "black",
                                   face = "bold", 
                                   vjust = 0.5, 
                                   hjust = 0.5, 
                                   angle = 0) 
  )
p
outfile=paste0(y_var, '.svg')
ggsave(outfile, p, dpi=300)
