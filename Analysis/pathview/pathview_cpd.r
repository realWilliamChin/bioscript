setwd("d:/pll/R_work/pathview/")
library(tidyverse)
library("pathview")
rm(list=ls())

# sim.cpd.data<-sim.mol.data(mol.type = "cpd",nmol = 3000)
# head(sim.cpd.data)
test_cpd <- read.table("cpd_data.txt",header=F)
head(test_cpd)
cpd_data <- test_cpd[,2]
cpd_data
names(cpd_data) <- test_cpd[,1]
head(cpd_data)
cpd_path <- read.table("cpd_path.txt",header=F)
cpd_path
cpd.path.name <- gsub('hsa','',cpd_path[,1])
cpd.path.name
pathview(cpd.data = cpd_data, pathway.id = cpd.path.name, species = "ko", out.suffix = "cpd.data", kegg.native = T,cpd.idtype = "kegg",kegg.dir='/home/colddata/qinqiang/script/Rscript/pathview/kegg_files')
