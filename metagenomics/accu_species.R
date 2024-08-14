library(ggplot2)  
# setwd("d:/pll/R_work/tmp/")
rm(list=ls())
library(vegan)  
library(ggplot2)  
otu <- read.delim('Species.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE,header = T)
otu=t(otu)
nrow(otu)
sp <- specaccum(otu, method = 'random')
sp.min<-1000*trunc(min((sp$richness)/1000))
sp.max<-1000*ceiling(max((sp$richness)/1000))
seq(sp.min,sp.max,by=1000)
summary(sp)
png("species_accumulation.png",width = 2000, height = 1600,res=320)
plot(sp, ci.type = 'poly', col = 'blue', lwd = 2, ci.lty = 0, ci.col = 'white', bty="L",ylim=c(sp.min,sp.max),xaxt = "n",yaxt="n",xlab="Sample Number",ylab="Species Count")
axis(1, at = 1:nrow(otu), labels = 1:nrow(otu))
axis(2, at = seq(sp.min,sp.max,by=500), labels = seq(sp.min,sp.max,by=500))
boxplot(sp, col = 'yellow', add = TRUE, pch = '+')
dev.off()
tiff("species_accumulation.tiff",width = 2000, height = 1600,res=320)
plot(sp, ci.type = 'poly', col = 'blue', lwd = 2, ci.lty = 0, ci.col = 'white', bty="L",ylim=c(sp.min,sp.max),xaxt = "n",yaxt="n",xlab="Sample Number",ylab="Species Count")
axis(1, at = 1:nrow(otu), labels = 1:nrow(otu))
axis(2, at = seq(sp.min,sp.max,by=500), labels = seq(sp.min,sp.max,by=500))
boxplot(sp, col = 'yellow', add = TRUE, pch = '+')
dev.off()








