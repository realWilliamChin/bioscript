setwd("d:/pll/R_work/pathview/")
library(tidyverse)
library(scales)
library(pathview)
rm(list=ls())
set.seed(2)
test.data <- data.frame(values = replicate(1, sample(0:10,1000, rep=TRUE)))
head(test.data)
ggplot(test.data, aes(x=values)) +
  geom_histogram(fill="grey70") + 
  stat_ecdf(aes(y=..y..*100),geom="smooth") + 
  scale_y_continuous(sec.axis=sec_axis(trans = ~./100 , name="percentage")) +
  theme_bw()
percent


head(gse16873.d)
eco.dat.kegg <- sim.mol.data(mol.type="gene",id.type="kegg",species="eco",nmol=3000)
head(eco.dat.kegg)
class(eco.dat.kegg)
names(eco.dat.kegg)
mode(eco.dat.kegg)
is.data.frame(eco.dat.kegg)
names(ko.data)
ko.data=sim.mol.data(mol.type="gene.ko", nmol=5000)
?sim.mol.data
length(ko.data)
head(ko.data)
class(ko.data)
ko.names<-names(ko.data)
ko.data<-rep(c(-1,1),each=2500)
names(ko.data)<-ko.names
head(ko.data)
pathview(gene.data = ko.data, pathway.id = "04112", species = "ko", out.suffix = "ko.data", kegg.native = T)
p <- pathview(gene.data = ko.data, pathway.id = "04112", species = "ko", out.suffix = "ko.data", kegg.native = T)
?pathview
p
getwd()
dir()
head(gse16873.d)
gsub('ko','',c("04112","04113"))
pathview(gene.data = ko.data, pathway.id = c("04112","04113"), species = "ko", out.suffix = "ko.data", kegg.native = T)
ko.names<-c('K04503','K10151','K10152','K02089','K02091')
ko.data<-rep(1,5)
names(ko.data)<-ko.names
ko.data
getwd()
tmp.dt<-read.table("cell_cycle2.txt",header=F)
tmp.dt[,1]
ko.data<-tmp.dt[,2]
names(ko.data)<-tmp.dt[,1]
pathview(gene.data = ko.data, pathway.id = c("04110","04210"), species = "ko", out.suffix = "ko.data", kegg.native = T)
###################################程序正式开始，之前是测试
setwd("d:/pll/R_work/pathview/")
rm(list=ls())
kid.dt<-read.table("CHX_vs_HBX_DEG_data.result",header=T)
ko.data<-kid.dt[,2]
names(ko.data)<-kid.dt[,1]
head(kid.dt)
path.dt<-read.table("CHX_vs_HBX_KEGG_path.txt",header=F,sep="\t")
#rownames(path.dt)
head(path.dt)
pathway.names<-path.dt[,1]
#head(pathway.names)
pathway.names<-gsub("ko",'',pathway.names)
?pathview
#head(pathway.names)
#pathview(gene.data = ko.data, pathway.id = pathway.names, species = "ko", out.suffix = "ko.data", kegg.native = T)

pathview(gene.data = ko.data, pathway.id = pathway.names, species = "ko", out.suffix = "ko.data", kegg.native = T, kegg.dir='./kegg_files/')

#pathview(gene.data = ko.data, pathway.id = '00901', species = "ko", out.suffix = "ko.data", kegg.native = T, kegg.dir='./kegg_files/')


head(path.dt)
colnames(path.dt)<-c('ID','Name')
#path.dt$Name[path.dt$ID=='ko00010']
#grep('.ko.data.png',list.files(),value = T)
library(fs)
dir()
for (i in path.dt[,1]) {
  old_name <- paste0(i,'.ko.data.png')
  new_name <- paste0(i,'_',path.dt$Name[path.dt$ID == i],'.ko.data.png')
  file.rename(old_name,new_name)
  
}
i<-''
#library(fs)
dir()
for (i in grep('.ko.data.png',list.files(),value = T)) {
  file_move(i,paste0('./CHX_vs_HBX/',i))
}
#file_move()

#file.rename('ko04928.ko.data.png.png','ko04928.ko.data.png')
#pathview(gene.data = ko.data, pathway.id = "00901", species = "ko", out.suffix = "ko.data", kegg.native = T)
#a<-'temp'
#dir.create(a)
library(ggplot2)
head(mpg)
p<-ggplot(mpg,aes(x=displ,y=hwy))+geom_point()
p+facet_grid(drv~.)
p+facet_grid(.~drv)


df<-data.frame(x=1:10,y=1:10,group=rep(c('a','b'),each=5))
df
