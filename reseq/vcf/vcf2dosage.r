setwd("C:/Users/Analysis/OneDrive/Work/zhiwusuo/01_Work/2023_03_30_马铃薯/target_gene")
getwd()
library("GWASpoly")
rm(list=ls())
dir()
#VCF2dosage(VCF.file="patato_snp_sp5.vcf",dosage.file = "patato_snp_dp1_sp5.csv",geno.code="GT",ploidy=4,samples=NULL,min.DP=1, max.missing=0.5, min.minor=5)
samples = ["sample112","sample142","sample162","sample138","sample47","sample32","sample123","sample129", "sample188"]
VCF2dosage(VCF.file="potato_snp_MI_specific.vcf",dosage.file = "potato_snp_MI_specific.csv",geno.code="GT",ploidy=4,samples=NULL,min.DP=1, max.missing=0, min.minor=5)
