setwd("/home/colddata/qinqiang/work/2023_07_03_potato_vcf_r")
library(CMplot)
data("pig60K")
head(pig60K)
library(openxlsx)
library(qqman)
potato<-read.xlsx("potato_snp_MI_specific_SNP_basicinfo.xlsx")
#potato<-read.table("potato_all_SNP_basicinfo_process.txt",sep='\t',header=T)
colnames(potato)
colnames(potato)[1]<-"SNP"
potato<-potato[grep('chr',potato$Chrom),]
CMplot(potato,plot.type="d",bin.size=1e6,col=c("darkgreen", "yellow", "red"),file="jpg",file.name="pecific",dpi=300,file.output=TRUE, verbose=TRUE)