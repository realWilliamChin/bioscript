setwd("C:/Users/Analysis/OneDrive/Work/zhiwusuo/00_Script/02_Rscript/GWAS/谷子")
getwd()
library("GWASpoly")
# Read in the arguments
library("optparse")
option_list = list(
  make_option(c("-t", "--traits"), type="integer", default=NULL, 
              help="traits 数量，从phenofile中查看", metavar="integer"),
  make_option(c("-n", "--samplesNumber"), type="integer", default=NULL, 
              help="样本的数量", metavar="integer"),lt=NULL, 
  make_option(c("-o", "--outPrefix"), type="character", default="default", 
              help="文件前缀", metavar="character")
) 
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

prefix<-opt$outPrefix
traits<-opt$
#######设置样本数
N<-opt$samplesNumber
rm(list=ls())

prefix<-"mature"
traits<-13
N<-203
#VCF2dosage(VCF.file="patato_snp_sp5.vcf",dosage.file = "patato_snp_dp1_sp5.csv",geno.code="GT",ploidy=4,samples=NULL,min.DP=1, max.missing=0.5, min.minor=5)
#设置miss率=0.5（即丢掉在一半以上的样品中缺失的snp位点），该值越小则保存的snp位点越少
#VCF2dosage(VCF.file="heading_stage.recode.vcf",dosage.file = "heading_dp1.csv",geno.code="GT",ploidy=2,samples=NULL,min.DP=1, max.missing=0.5, min.minor=5)

#设置miss率为0.05
VCF2dosage(
  VCF.file = paste0(prefix, "_stage.recode.vcf"),
  dosage.file = paste0(prefix, "_dp1_miss_0.05.csv"),
  geno.code="GT",ploidy=2,samples=NULL,min.DP=1, max.missing=0.05, min.minor=5)

#VCF2dosage(VCF.file="heading_stage.recode.vcf",dosage.file = "heading_dp5_miss_0.05.csv",geno.code="GT",ploidy=2,samples=NULL,min.DP=0.05, max.missing=0.05, min.minor=5)

#VCF2dosage(VCF.file="heading_stage.recode.vcf",dosage.file = "heading_dp2.csv",geno.code="GT",ploidy=2,samples=NULL,min.DP=2, max.missing=0.5, min.minor=5)


#genofile <- system.file("extdata", "new_potato_geno.csv", package = "GWASpoly")
#phenofile <- system.file("extdata", "new_potato_pheno.csv", package = "GWASpoly")

genofile<-paste0(prefix, "_dp1_miss_0.05.csv")
phenofile<-paste0(prefix, "_stage_phe_all.csv")
#######注意修改n.traits这个参数，确定表型参数个数
data <- read.GWASpoly(ploidy=2, pheno.file=phenofile, geno.file=genofile,
                      format="numeric", n.traits=traits, delim=",")
class(data)
#data.loco <- set.K(data,LOCO=TRUE,n.core=2)
###为了做后面的qtl，所以这里选择参数LOCO=F
data.original <- set.K(data,LOCO=FALSE,n.core=2)


#params <- set.params(geno.freq = 1 - 2/N, fixed = "env", fixed.type = "factor",MAF = 0.01)

params <- set.params(geno.freq = 1 - 2/N, fixed = "env", fixed.type = "factor",MAF = 0.01,n.PC=5)

#GWASpoly的traits选项，默认为phenofile中的所有traits，也可以单独指定一个trait
#data.loco.scan <- GWASpoly(data=data.loco,models=c("additive","1-dom"),
#                   params=params,n.core=2)
#summary(data.loco.scan)

#file.remove("tmp.txt")
#write.GWASpoly(data.loco.scan, trait="disease", filename="tmp.txt", what = "scores", delim = "\t")


#data.loco.scan
data.original.scan <- GWASpoly(data.original,models=c("additive","1-dom"),
                               params=params,n.core=2)
#?GWASpoly
library(ggplot2)
for(i in paste0("trait",1:traits)){
  qq.p<-qq.plot(data.original.scan,trait= i) + ggtitle(label="Original")
  outfile<-paste0(prefix,"_",i,"_qq.jpeg")
  ggsave(outfile,qq.p,dpi=350)
}
rm(i)
#qq.plot(data.original.scan,trait="trait1") + ggtitle(label="Original")
#qq.plot(data.original.scan,trait="trait6") + ggtitle(label="Original")


#head(data.loco.scan)
#qq.plot(data.original.scan,trait="disease") + ggtitle(label="Original")
#qq.plot(data.loco.scan,trait="trait1") + ggtitle(label="LOCO")
#qq.plot(data.loco.scan,trait="trait6") + ggtitle(label="LOCO")
#data2 <- set.threshold(data.loco.scan,method="Bonferroni",level=0.999)

#data2 <- set.threshold(data.loco.scan,method="Bonferroni",level=0.05)

#data2 <- set.threshold(data.loco.scan,method="M.eff",level=0.05)

#data2 <- set.threshold(data.original.scan,method="Bonferroni",level=0.05)
data2 <- set.threshold(data.original.scan,method="M.eff",level=0.05)
data2
#########上面这一步产生-log10(p)的阈值，会显示在后面的曼哈顿图中，要记录下来，这个样品中是5.56
class(data2)
summary(data2)
#paste0("trait",1:10)
#p <- manhattan.plot(data2,traits=paste0("trait",1:10))
#p
p <- manhattan.plot(data2)+ theme(axis.text.x = element_text(angle=90,vjust=0.5))
ggsave(paste0(prefix, "_traits_all_man.jpeg"),p,dpi=350, width = 10, height = 15)
rm(p)

for(i in paste0("trait",1:traits)){
  man.p<-manhattan.plot(data2,traits=i)+ theme(axis.text.x = element_text(angle=90,vjust=0.5))
  outfile<-paste0(prefix,"_",i,"_man.jpeg")
  ggsave(outfile,man.p,dpi=350)
}
rm(i)
#p + theme(axis.text.x = element_text(angle=90,vjust=0.5))

#chr_name<-c(paste0("chr0",1:9),"chr10","chr11","chr12")
# 1 3 4 6 7 8 9 
#for (i in 1:12) {
# outfile<-paste0(chr_name[i],".jpeg")
#chr_num=paste0("chr",i)
#  p.chr<-manhattan.plot(data2,traits="disease",chrom=chr_name[i])
#ggsave(outfile,p.chr)
#}
#p.chr<-manhattan.plot(data2,traits="disease",chrom="chr12")

# 1 3 5 6 8 9 11
#manhattan.plot(data2,traits="disease",chrom="chr11")

LD.p <- LD.plot(data2, max.loci=1000)
ggsave(paste0(prefix, "_LD.jpeg"),LD.p,dpi=350)

#p + xlim(0,30) 

#qtl <- get.QTL(data=data2,traits="disease",models=c("additive","1-dom"),bp.window=5e6)
qtl <- get.QTL(data=data2,models=c("additive","1-dom"),bp.window=5e6)

class(qtl)
knitr::kable(qtl)
dim(qtl)
write.table(qtl,paste0(prefix, "_qtl_result.txt"),quote = F,sep="\t",row.names=F)
#?fit.QTL
#fit.ans <- fit.QTL(data=data2,trait="trait1",
#                   qtl=qtl[,c("Marker","Model")],
#                   fixed=data.frame(Effect="env",Type="factor"))
#class(fit.ans)
#write.table(fit.ans,"fit.ans_heading_result_v3.txt",quote = F,sep="\t",row.names=F)
#knitr::kable(fit.ans,digits=3)

#write.GWASpoly(data2, trait="trait1", filename="heading_trait1_gwas.txt", what = "scores", delim = "\t")
for(i in paste0("trait",1:traits)){
  outfile<-paste0(prefix,"_",i,"_gwas.txt")
  write.GWASpoly(data2, trait=i, filename=outfile, what = "scores", delim = "\t")
}


