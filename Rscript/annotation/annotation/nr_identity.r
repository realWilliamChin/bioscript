library(ggplot2)
library(ggthemes)
setwd("d:/pll/R_work/annotation")
data.nr<-read.table("Trifolium_repens_Baisanye_unigene_CDS_nr_diamond_uniq.blast",sep="\t",header=F)
dim(data.nr)
x<-data.nr[,3]
dim(data.nr[,3])
tmp<-hist(x)
bin.start<-tmp$breaks[-length(tmp$breaks)]
bin.end<-tmp$breaks[-1]
identity<-paste(bin.start,bin.end,sep="-")
count<-tmp$counts
id.count<-data.frame(Identity=identity,Count=count)
newlegend <- paste(id.count$Identity, " (", round(id.count$Count/sum(id.count$Count)* 100, 2), "%)", sep = "")
id.p<-ggplot(id.count, aes(x = "", y = Count, fill = Identity)) + geom_bar(stat = "identity") + coord_polar(theta = "y") +
       theme_base()+theme(panel.grid = element_blank(), axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank())+
       scale_fill_discrete(breaks = id.count$Identity, labels = newlegend)+ggtitle("Identity distribution")
ggsave("nr_id_distribution.jpeg",id.p,dpi=320,width=10)

evalue<-data.nr[,11]
evalue.hist<-hist(evalue,breaks=10)
e.start<-evalue.hist$breaks[-length(evalue.hist$breaks)]
e.start<-gsub(".00000000000001",'',as.character(e.start))

e.end<-evalue.hist$breaks[-1]
e.end<-gsub(".00000000000001",'',as.character(e.end))
e_value<-paste(e.start,e.end,sep="~~")
e_count<-evalue.hist$counts
e.count<-data.frame(e_value=e_value,Count=e_count)
e.p<-ggplot(e.count, aes(x = "", y = Count, fill = e_value)) + geom_bar(stat = "identity") + coord_polar(theta = "y") +
  theme_base()+theme(panel.grid = element_blank(), axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank())+
  scale_fill_discrete(breaks = e.count$e_value)+ggtitle("E value distribution")
ggsave("nr_evalue_distribution.jpeg",e.p,dpi=320,width=10)


sp<-as.vector(data.nr[,13])
sp<-data.nr[,13]
gsub('.+\\[',"",sp,perl=TRUE)->sp
gsub('\\]',"",sp,perl=TRUE)->sp
gsub("\n","",sp,perl=TRUE)->sp

head(sp)
table(as.character(sp))
sp.dt<-as.data.frame(table(as.character(sp)))
head(sp.dt)

######################################nr annotation count picture species#####
species<-read.table("species_count.txt",sep="\t",header=F,check.names = F)
#species
colnames(species)<-c("Species","Count")
sp.p<-ggplot(species, aes(x = "", y = Count, fill = Species)) + geom_bar(stat = "identity") + coord_polar(theta = "y") +
  theme_base()+theme(panel.grid = element_blank(), axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank())+
  scale_fill_discrete(breaks = species$Species)+ggtitle("Species distribution")
ggsave("nr_sp_distribution.jpeg",sp.p,dpi=320,width=10)
####################################nr annotation identity picture###########################################

id.nr<-read.table("identity.txt",sep="\t",header=F)
id<-id.nr[,1]
id.count<-as.data.frame(summary(cut(as.vector(id),10)))
dim(id.count)
id.count$Identity<-rownames(id.count)
id.count
colnames(id.count)[1]<-"Count"
newlegend <- paste(id.count$Identity, " (", round(id.count$Count/sum(id.count$Count)* 100, 2), "%)", sep = "")
id.p<-ggplot(id.count, aes(x = "", y = Count, fill = Identity)) + geom_bar(stat = "identity") + coord_polar(theta = "y") +
  theme_base()+theme(panel.grid = element_blank(), axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank())+
  scale_fill_discrete(breaks = id.count$Identity,labels=newlegend)+ggtitle("Identity distribution")
id.p
ggsave("nr_id_distribution.jpeg",id.p,dpi=320,width=10)
#################################nr evalue count picture#############################################
e.nr<-read.table("evalue.txt",sep="\t",header=F)
evalue<-e.nr[,1]
e.count<-as.data.frame(summary(cut(as.vector(evalue),10)))
e.count
rownames(e.count)[1]="(0,1e-06]"
e.count$Evalue<-rownames(e.count)
colnames(e.count)[1]<-"Count"
#colnames(e.count)[2]<-"Evalue"
e.p<-ggplot(e.count, aes(x = "", y = Count, fill = Evalue)) + geom_bar(stat = "identity") + coord_polar(theta = "y") +
  theme_base()+theme(panel.grid = element_blank(), axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank())+
  scale_fill_discrete(breaks = e.count$Evalue)+ggtitle("E value distribution")

ggsave("nr_evalue_distribution.jpeg",e.p,dpi=320,width=10)





