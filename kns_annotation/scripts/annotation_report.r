##### input files
#NR_ID.list,GO_ID.list,KEGG_ID.list,Swisss_ID.list: one column: gene ID
#all_gene_id.txt :one column with all gene ID of the Unigenes
#KEGG_clean.txt : KEGG annotation file
#GO_BP.txt, GO_CC.txt,GO_ MF.txt : two columns with geneID,GO ID
#COG count picture COG_count.txt
#nr species count: species_count.txt
#nr identity count: identity.txt
#nr evalue count: evalue.txt
################################ venn ########################################
NR<-read.table("NR_ID.list",sep="\t",header=F,stringsAsFactors = F)
GO<-read.table("GO_ID.list",sep="\t",header=F,stringsAsFactors = F)
KEGG<-read.table("KEGG_ID.list",sep="\t",header=F,stringsAsFactors = F)
Swiss<-read.table("Swiss_ID.list",sep="\t",header=F,stringsAsFactors = F)
#install.packages("VennDiagram")
library(VennDiagram)
library(RColorBrewer)
library(ggplot2)
library(ggthemes)
head(KEGG)
#########annotation venn picture
venn.diagram(x=list(NR=NR[,1],SwissProtein =Swiss[,1],KEGG=KEGG[,1],GO=GO[,1]),
             filename="annotation_venn.jpeg",resolution = 320,
             fill = c("cornflowerblue", "green", "yellow", "darkorchid1"),
             col = "transparent",
             fontfamily = "serif", 
             cex = 1.5,  # ????????????
             cat.cex = 2,
             #cat.pos = c(-27, 27, 135,315),
             width=4000,
             #rotation.degree = 270,
             fontface = "bold",  # ??????
             )
all<-read.table("all_gene_id.txt",header=F)
#Database	Total unigenes	Annotated  uingene	Percentage output the summary of annotated results
perc<-100*round(c(nrow(KEGG),nrow(GO),nrow(NR),nrow(Swiss))/nrow(all),digits = 3)
perc<-paste(perc,"%",sep='')
anno.summary<-data.frame(Database=c("KEGG","GO","NR","Swiss Protein"),Total=rep(nrow(all),4),Annotated=c(nrow(KEGG),nrow(GO),nrow(NR),nrow(Swiss)),Percentage=perc)
anno.summary
write.table(anno.summary,"annotation_summary.txt",sep="\t",quote=F,row.names = F)

##############kegg terms levelB count picture
data.kegg<-read.table("KEGG_clean.txt",sep="\t",stringsAsFactors = F,header=F,quote = '')
count.kegg<-data.frame(table(data.kegg[,3:4]))
count.kegg<-count.kegg[count.kegg[,3]>0,]
head(count.kegg)
count.kegg$V4<-gsub("A\\d+:","",count.kegg$V4,perl=TRUE)
count.kegg$V3<-gsub("\\d+:","",count.kegg$V3,perl=TRUE)

colnames(count.kegg)<-c("LevelB","LevelA","Count")
colnames(count.kegg)
n_levels <- length(unique(count.kegg$LevelB))

# 设置图片的高度，可以根据需要调整乘数以获得合适的尺寸
height <- n_levels * 0.2 # 例如，每个级别0.5英寸的高度

# 绘制图表
kegg.p <- ggplot(count.kegg, aes(y = LevelB, x = Count, fill= LevelA)) +
          geom_bar(stat="identity") +
          theme_base() +
          scale_y_discrete(limits=count.kegg$LevelB)

# 打印图表
print(kegg.p)

# 保存图表，使用动态计算的高度
ggsave("kegg_count.jpeg", kegg.p, dpi=300, width=15, height=height)
##########GO stat and drawing pictures ############
#data.go<-read.table("go_stat.txt",sep="\t",header=F,stringsAsFactors = F)
#colnames(data.go)<-c("Term_name","Count","GO")
#head(data.go)
#data.go$Term_name<-gsub("GO:\\d+_","",data.go$Term_name,perl=TRUE)
#go.p<-ggplot(data.go,aes(x=Term_name,y=Count,fill=GO))+geom_bar(stat="identity")+
 # scale_x_discrete(limits=data.go$Term_name)+theme_base()+theme(axis.text.x=element_text(angle=60,hjust=1,vjust=1))
#ggsave("GO_count.jpeg",go.p,dpi=300,width=15)
################################################## GO stat and drawing pictures Top20 Terms for each GO categories(BP,CC,MF)
BP<-read.table("swiss_GO_BP_ID.txt",sep="\t",header=F,stringsAsFactors = F,quote = '')
#head(BP)
#table(head(BP[,2]))
BP.count<-data.frame(table(BP[,2]))
BP.20 <- head(BP.count[order(-BP.count[,2]),],n=20)

CC<-read.table("swiss_GO_CC_ID.txt",sep="\t",header=F,stringsAsFactors = F,quote='')
CC.count<-data.frame(table(CC[,2]))
CC.20 <- head(CC.count[order(-CC.count[,2]),],n=20)
#CC.20
#?read.table
MF<-read.table("swiss_GO_MF_ID.txt",sep="\t",header=F,stringsAsFactors = F,quote ="")
MF.count<-data.frame(table(MF[,2]))
MF.20 <- head(MF.count[order(-MF.count[,2]),],n=20)

GO.stat<-rbind(BP.20,CC.20,MF.20)

#rep(c("BP","CC","MF"),each=20)
GO.stat$GO<-rep(c("BP","CC","MF"),each=20)
colnames(GO.stat)[1:2]<-c("Term","Count")
GO.stat$Term<-gsub("GO:\\d+_","",GO.stat$Term,perl=TRUE)
#GO.stat$Count<-log2(GO.stat$Count+1)
#GO.stat
go.p<-ggplot(GO.stat,aes(x=Term,y=Count,fill=GO))+geom_bar(stat="identity")+
  scale_x_discrete(limits=GO.stat$Term)+theme_base()+theme(axis.text.x=element_text(angle=75,hjust=1,vjust=1))
ggsave("GO_count.jpeg",go.p,dpi=300,width=15,height=10)
######################COG distribution picture
cog<-read.table("COG_count.txt",sep="\t",header=F,stringsAsFactors = F,check.names = F)
colnames(cog)<-c("Category","Group","Function","Count")
cog.p<-ggplot(cog,aes(x=Category,y=Count,fill=Category))+geom_bar(stat="identity")+scale_fill_discrete(labels=cog$Function)+theme_base()
#cog.p
ggsave("cog_function_count.jpg",cog.p,dpi=320,width=18,height=12)

######################################################
######################################nr annotation count picture species#####
species<-read.table("species_count_Top10.txt",sep="\t",header=F,check.names = F)
#species
colnames(species)<-c("Species","Count")
sp.p<-ggplot(species, aes(x = "", y = Count, fill = Species)) + geom_bar(stat = "identity") + coord_polar(theta = "y") +
  theme_base()+theme(panel.grid = element_blank(), axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank())+
  scale_fill_discrete(breaks = species$Species)+ggtitle("Species distribution")
ggsave("nr_sp_distribution_Top10.jpeg",sp.p,dpi=320,width=10)
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
# e.nr<-read.table("evalue.txt",sep="\t",header=F)
# evalue<-e.nr[,1]
# e.count<-as.data.frame(summary(cut(as.vector(evalue),10)))
# e.count
# rownames(e.count)[1]="(0,1e-06]"
# e.count$Evalue<-rownames(e.count)
# colnames(e.count)[1]<-"Count"
# #colnames(e.count)[2]<-"Evalue"
# e.p<-ggplot(e.count, aes(x = "", y = Count, fill = Evalue)) + geom_bar(stat = "identity") + coord_polar(theta = "y") +
#   theme_base()+theme(panel.grid = element_blank(), axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank())+
#   scale_fill_discrete(breaks = e.count$Evalue)+ggtitle("E value distribution")
# 
# ggsave("nr_evalue_distribution.jpeg",e.p,dpi=320,width=10)

