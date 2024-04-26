library(VennDiagram)
library(RColorBrewer)
library(ggplot2)
library(ggthemes)
H<-read.table("Species_H.txt",sep="\t",header=F,stringsAsFactors = F)
LC<-read.table("Species_LC.txt",sep="\t",header=F,stringsAsFactors = F)
venn.diagram(x=list(H=H[,1],LC=LC[,1]),
	direct.area=FALSE,
	filename="annotation_venn.jpeg",resolution = 320,
	fill = c("blue", "red"),
	col = "transparent",
	fontfamily = "serif",
	cex = 1.5,
	cat.cex = 2,
	#cat.pos = c(-27, 27, 135,315),
	width=4000,
	#rotation.degree = 270,
	fontface = "bold",
)
