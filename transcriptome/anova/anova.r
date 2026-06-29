library(optparse)
library(readxl)
rm(list=ls())

option_list = list(
  make_option(c("-f", "--fpkmfile"), type="character", default=NULL,
              help="输入文件", metavar="character"),
  make_option(c("-s", "--samplesfile"), type="character", default=NULL,
              help="输入文件", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="输出文件", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Check that all required arguments are provided
if (is.null(opt$fpkmfile)){
  print_help(opt_parser)
  stop("请输入 fpkmfile 文件", call.=FALSE)
}else if (is.null(opt$samplesfile)){
  print_help(opt_parser)
  stop("请输入 samplesfile 文件名", call.=FALSE)
}else if (is.null(opt$output)){
  print_help(opt_parser)
  stop("请输入 output 文件名", call.=FALSE)
}

fpkmfile=opt$fpkmfile
samplesfile=opt$samplesfile
output_file=opt$output

if (endsWith(fpkmfile, ".xlsx")) {
  data <- read_xlsx(fpkmfile,1)
} else if (endsWith(fpkmfile, ".txt")) {
  data <- read.table(fpkmfile, sep="\t", header=TRUE,stringsAsFactors = F,quote="",check.names=F)
} else {
  stop("Unsupported file format. Supported formats are: .xlsx and .txt")
}

data<-as.data.frame(data)
head(data)
rownames(data)<-data$GeneID
data<-data[,-1]
data<-subset(data,rowSums(data)>0)
group<-read.table(samplesfile,sep="\t",header=T,check.names = F)
colnames(data)
tmp.group<-group$group[match(colnames(data),group$sample)]
p<-NULL
for (i in rownames(data)) {
  tmp.count<-as.numeric(data[i,])
  tmp.dt<-data.frame(tmp.count,group=tmp.group)
  colnames(tmp.dt)[1]<-i
  tmp.dt$group<-as.factor(tmp.dt$group)
  mod <- aov( tmp.dt[,1]~group , data = tmp.dt)
  p<-c(p,summary(mod)[[1]][5]$Pr[1])
}
data$p_value<-p
data<-data[order(data$p_value),]
p.adjust(data$p_value,method = "BH")->p.aj
data$BH_p_value<-p.aj
df.rname<-data.frame(GeneID=rownames(data))
data<-cbind(df.rname,data)
write.table(data,output_file,sep="\t",quote=F,row.names = F)