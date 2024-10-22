library('microeco')
library(openxlsx)
library(ggplot2)
library(optparse)
rm(list=ls())

#TODO: 加入参数解析功能

my_set_colors <- c(
  "#f44336", "#e91e63", "#9c27b0", "#3f51b5", "#2196f3",
  "#00bcd4", "#009688", "#4caf50", "#8bc34a", "#cddc39",
  "#ffeb3b", "#ffc107", "#ff5722", "#795548", "#9e9e9e",
  "#607d8b", "#c8bfe7", "#b97a57", "#ffaec9", "#1ee1c4"
)

feature_table<-read.table("Species.txt",row.names=1,header = T,check.names=F)
sample_table<-read.table("samples_described.txt",header = T,check.names=F)
rownames(sample_table)<-sample_table$sample
tax_table<-read.table("Species_taxon.txt",row.names=1,header = T,check.names=F)
dataset <- microtable$new(  
  sample_table = sample_table,  
  otu_table = feature_table,  
  tax_table = tax_table
)

lefse_res<-trans_diff$new(dataset=dataset,
                          method = "lefse",
                          group = "group",
                          alpha=0.05,
                          lefse_subgroup=NULL,
                          p_adjust_method = "none"
  
)
res.dt<-lefse_res$res_diff
write.xlsx(res.dt, file = "LDA_result.xlsx", sheetName = "Sheet1", rowNames = FALSE)

p1<-lefse_res$plot_diff_bar(add_sig = T,
                        use_number = 1:20,
                        width=0.8,
                        group_order = unique(sample_table$group),
                        color_values = my_set_colors,
                        )
ggsave("LDA_bar.png", p1, dpi=320,width=8,height=10)

p2<-lefse_res$plot_diff_abund(add_sig = T,
                          color_values = my_set_colors,
                          use_number = 1:20,
                        #width=0.8,
                          group_order = unique(sample_table$group),
)
ggsave("Abundance_bar.png",p2, dpi=320,width=8,height=13)

p3<-lefse_res$plot_diff_cladogram(color = my_set_colors,
                              use_taxa_num = 500,
                              use_feature_num = 20,
                              clade_label_level = 5,
                              #clade_label_size = 100,
                              #select_show_labels = "c__Bacilli",
                              node_size_offset = 1,
                              #annotation_shape = 21,
                              #annotation_shape_size = 4.5,
                              group_order = unique(sample_table$group)
                        #clade_label_size =1
                      )

ggsave("Clade_tree.png",p3,dpi=320,width=12,height=12)
