library('microeco')
library(openxlsx)
library(ggplot2)
setwd("/home/colddata/qinqiang/Project/2024_08_12_hongjiyinzu_sdzlyy/microeco2")
rm(list=ls())
feature_table<-read.table("Species.txt",row.names=1,header = T,check.names=F)
sample_table<-read.table("samples_described.txt",header = T,check.names=F)
rownames(sample_table)<-sample_table$SampleID
tax_table<-read.table("Species_taxon.txt",row.names=1,header = T,check.names=F)
dataset <- microtable$new(  
  sample_table = sample_table,  
  otu_table = feature_table,  
  tax_table = tax_table
)
#?trans_diff
#?p.adjust

#trans_diff$new()
lefse_res<-trans_diff$new(dataset=dataset,
                          method = "lefse",
                          group = "Group",
                          alpha=0.05,
                          lefse_subgroup=NULL,
                          p_adjust_method = "none"
  
)
#lefse_res$plot_diff_abund()
#colnames(lefse_res$res_abund)
res.dt<-lefse_res$res_diff
write.xlsx(res.dt, file = "LDA_result.xlsx", sheetName = "Sheet1", rowNames = FALSE)
#colnames(res.dt)
#unique(res.dt$Group)
#colnames(res.dt)
#range(res.dt$LDA)
#head(lefse_res$res_diff)
#colnames(lefse_res$res_diff)
#dim(lefse_res$res_diff)
#class(lefse_res$res_diff)
#summary()
#lefse_res$plot_diff_bar(
#    
#)
#?png()
# png("LDA_bar.png",width=800,height=800)
p1<-lefse_res$plot_diff_bar(add_sig = T,
                        use_number = 1:20,
                        width=0.8,
                        group_order = c("R","NR"),
                        color_values = c("red","green")
                        )
ggsave("LDA_bar.png", p1, dpi=320)
# dev.off()
# png("Abundance_bar.png",width=800,height=800)
p2<-lefse_res$plot_diff_abund(add_sig = T,
                          color_values = c("red","green"),
                          use_number = 1:20,
                        #width=0.8,
                          group_order = c("R","NR"),
)
ggsave("Abundance_bar.png",p2, dpi=320)
#dev.off()


#png("Clade_tree.png",width=800,height=800)

p3<-lefse_res$plot_diff_cladogram(color = c("red","green"),
                              use_taxa_num = 500,
                              use_feature_num = 20,
                              clade_label_level = 5,
                              #clade_label_size = 100,
                              #select_show_labels = "c__Bacilli",
                              node_size_offset = 1,
                              #annotation_shape = 21,
                              #annotation_shape_size = 4.5,
                              group_order = c("R","NR")
                        #clade_label_size =1
                      )
ggsave("Clade_tree.png",p3,dpi=320)
#dev.off()
#ggsave('Clade_tree.pdf',p)


#library(ggtree)

#lefse_res$plot_diff_cladogram(
#)


