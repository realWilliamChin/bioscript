suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(openxlsx)) # 读取.xlsx文件
suppressPackageStartupMessages(library(ggplot2)) # 柱状图和点状图
suppressPackageStartupMessages(library(stringr)) # 基因ID转换
suppressPackageStartupMessages(library(enrichplot)) # GO,KEGG,GSEA
suppressPackageStartupMessages(library(clusterProfiler)) # GO,KEGG,GSEA
suppressPackageStartupMessages(library(GOplot)) # 弦图，弦表图，系统聚类图
suppressPackageStartupMessages(library(DOSE))
suppressPackageStartupMessages(library(ggnewscale))
suppressPackageStartupMessages(library(topGO)) # 绘制通路网络图
suppressPackageStartupMessages(library(circlize)) # 绘制富集分析圈图
suppressPackageStartupMessages(library(ComplexHeatmap)) # 绘制图例
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyr))

option_list <- list(
  make_option(c("--inputidfile"),
    type = "character", default = NULL,
    help = "输入 GeneID list 文件，文件名后面需要是 _ID.txt, header 为 GeneID", metavar = "character"
  ),
  make_option(c("--genego"),
    type = "character", default = NULL,
    help = "提供 swiss 注释出来的 gene_go.txt 文件", metavar = "character"
  ),
  make_option(c("--keggclean"),
    type = "character", default = NULL,
    help = "提供 kegg 注释出来的 KEGG_clean.txt 文件", metavar = "character"
  ),
  make_option(c("--outputdir"),
    type = "character", default = NULL,
    help = "输出结果文件的目录，输出 EnrichmentGO.xlsx 和 EnrichmentKEGG.xlsx", metavar = "character"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check that all required arguments are provided
if (!file.exists(opt$genego)) {
  print_help(opt_parser)
  stop("请提供 gene go 文件", call. = FALSE)
} else if (!file.exists(opt$keggclean)) {
  print_help(opt_parser)
  stop("请提供 kegg 注释出来的 KEGG_clean.txt 文件", call. = FALSE)
}

if (!dir.exists(opt$outputdir)) {
  dir.create(opt$outputdir)
  cat("文件夹已创建：", opt$outputdir, "\n")
} else {
  cat("文件夹已存在，自动覆盖已有的文件：", opt$outputdir, "\n")
}

if (substr(opt$outputdir, nchar(opt$outputdir), nchar(opt$outputdir)) != "/") {
  opt$outputdir <- paste0(opt$outputdir, "/")
}

# Assign the first argument to prefix
id_list_file <- opt$inputidfile
gene_go <- opt$genego
kegg_clean <- opt$keggclean
output_dir <- opt$outputdir

# 计算比值的函数
calculate_ratio <- function(ratio_string) {
  parts <- strsplit(ratio_string, "/")[[1]] # 分割字符串
  numerator <- as.numeric(parts[1]) # 转换分子为数值
  denominator <- as.numeric(parts[2]) # 转换分母为数值
  return(numerator / denominator) # 返回比值
}

go_enrich <- function(deg_id_file, output_dir) {
  gene_list <- read.delim(deg_id_file, stringsAsFactors = FALSE, header = T)
  names(gene_list)[1] <- c("gene_id")
  gene_select <- gene_list$gene_id
  go_rich_bp <- enricher(
    gene = gene_select,
    TERM2GENE = go_anno[go_anno$Ontology == "biological_process", ][c("ID", "gene_id")],
    TERM2NAME = go_anno[go_anno$Ontology == "biological_process", ][c("ID", "Description")],
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    qvalueCutoff = 1,
    minGSSize = 1,
    maxGSSize = 1000
  )

  go_rich_mf <- enricher(
    gene = gene_select,
    TERM2GENE = go_anno[go_anno$Ontology == "molecular_function", ][c("ID", "gene_id")],
    TERM2NAME = go_anno[go_anno$Ontology == "molecular_function", ][c("ID", "Description")],
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    qvalueCutoff = 1,
    minGSSize = 1,
    maxGSSize = 1000
  )

  go_rich_cc <- enricher(
    gene = gene_select,
    TERM2GENE = go_anno[go_anno$Ontology == "cellular_component", ][c("ID", "gene_id")],
    TERM2NAME = go_anno[go_anno$Ontology == "cellular_component", ][c("ID", "Description")],
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    qvalueCutoff = 1,
    minGSSize = 1,
    maxGSSize = 1000
  )

  value_counts = 0
  if (length(go_rich_bp@result) != 0) {
    bp_result <- go_rich_bp@result
    bp_result$GeneRatio <- sapply(bp_result$GeneRatio, calculate_ratio)
    bp_result$BgRatio <- sapply(bp_result$BgRatio, calculate_ratio)
    bp_result$RichFactor <- bp_result$GeneRatio / bp_result$BgRatio
    bp_result$Ontology <- rep("BP", nrow(bp_result))
    # go_rich.bp <- as.data.frame(bp_result)
    go_rich.total <- as.data.frame(bp_result)
    value_counts = 1
  }

  if (length(go_rich_cc@result) != 0) {
    cc_result <- go_rich_cc@result
    cc_result$GeneRatio <- sapply(cc_result$GeneRatio, calculate_ratio)
    cc_result$BgRatio <- sapply(cc_result$BgRatio, calculate_ratio)
    cc_result$RichFactor <- cc_result$GeneRatio / cc_result$BgRatio
    cc_result$Ontology <- rep("CC", nrow(cc_result))
    # go_rich.cc <- as.data.frame(go_rich_cc)
    if (is.data.frame(go_rich.total)) {
      go_rich.total <- rbind(go_rich.total, cc_result)
    } else {
      go_rich.total <- as.data.frame(cc_result)
    }
    value_counts = 1
  }

  if (length(go_rich_mf@result) != 0) {
    mf_result <- go_rich_mf@result
    mf_result$GeneRatio <- sapply(mf_result$GeneRatio, calculate_ratio)
    mf_result$BgRatio <- sapply(mf_result$BgRatio, calculate_ratio)
    mf_result$RichFactor <- mf_result$GeneRatio / mf_result$BgRatio
    mf_result$Ontology <- rep("MF", nrow(mf_result))
    # go_rich.mf <- as.data.frame(go_rich_mf)
    if (is.data.frame(go_rich.total)) {
      go_rich.total <- rbind(go_rich.total, mf_result)
    } else {
      go_rich.total <- as.data.frame(mf_result)
    }
    value_counts = 1
  }
  
  # go_rich.total <- rbind(bp_result, cc_result, mf_result)

  # 检测是否存在 go_rich.total 的一个 dataframe 且不为空
  # TODO: 检测不到，还是会出错
  if (value_counts == 0) {
    print("没有富集到任何的GO信息\n")
    return()
  }

  go_rich.total <- go_rich.total[, c("ID", "Description", "GeneRatio", "BgRatio", "RichFactor", "pvalue", "p.adjust", "qvalue", "geneID", "Count", "Ontology")]

  go_enrich_file_name <- paste0(output_dir, sub("_ID.txt", "_EnrichmentGO.xlsx", basename(deg_id_file)))
  cat("正在输出文件",go_enrich_file_name, "\n")
  write.xlsx(go_rich.total, go_enrich_file_name)
}

kegg_enrich <- function(deg_id_file, kegg2gene, kegg2name, output_dir) {
  ###### 导入目标基因列表################################
  gene <- read.table(deg_id_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  gene <- as.factor(gene$V1)
  # 富集分析
  enrichment <- enricher(gene,
    TERM2GENE = kegg2gene, TERM2NAME = kegg2name, pvalueCutoff = 1,
    pAdjustMethod = "BH", qvalueCutoff = 1, minGSSize = 1, maxGSSize = 1000
  )
  if (length(enrichment@result) == 0) {
    cat("没有富集到任何的KEGG信息\n")
    return()
  }
  # # 绘制条形图, 绘制气泡图
  # barplot(enrichment)
  # dotplot(enrichment)
  enrich_result <- enrichment@result
  enrich_result$GeneRatio <- sapply(enrich_result$GeneRatio, calculate_ratio)
  enrich_result$BgRatio <- sapply(enrich_result$BgRatio, calculate_ratio)
  enrich_result$RichFactor <- enrich_result$GeneRatio / enrich_result$BgRatio
  enrich_result$Ontology <- rep("KEGG", nrow(enrich_result))
  enrich_result <- na.omit(enrich_result)  # 在 KO_id_pathway_name.txt 中没包含所有的 ko，有的对上了没有描述，所以要去掉这些
  go_enrich_file_name <- paste0(
    output_dir,
    sub("_ID.txt", "_EnrichmentKEGG.xlsx", basename(deg_id_file))
  )
  enrich_result <- enrich_result[, c("ID", "Description", "GeneRatio", "BgRatio", "RichFactor", "pvalue", "p.adjust", "qvalue", "geneID", "Count", "Ontology")]
  write.xlsx(enrich_result, go_enrich_file_name)
  # enrichplot::cnetplot(go_rich_bp,circular=FALSE,colorEdge = T,edge = F,color_category = "red", color_gene = "green",cex_category = 1,cex_label_category = 1)
}

# 加载GO注释描述文件,GOid\tdescriptiom\tGOclass(BP CC MF)
go_class <- read.delim(
  "/home/colddata/qinqiang/script/lib/go_term.list",
  header = FALSE,
  stringsAsFactors = FALSE
)
names(go_class) <- c("ID", "Description", "Ontology")

# 加载背景库文件，gene\tGOid
go_anno <- read.delim(gene_go, header = FALSE, stringsAsFactors = FALSE)
names(go_anno) <- c("gene_id", "ID")
# 合并背景与GO文件
go_anno <- merge(go_anno, go_class, by = "ID", all.x = TRUE)

kegg_clean <- read.table(kegg_clean,
  header = F, sep = "\t", quote = "",
  col.names = c("GeneID", "KEGG_Pathway", "Metabolism_Category", "General_Metabolism", "K_Number", "Protein", "EC_Number")
)
kegg2gene <- kegg_clean[c("KEGG_Pathway", "GeneID")]
kegg2gene$KEGG_Pathway <- sub(":.*", "", kegg2gene$KEGG_Pathway)
kegg2name <- read.table("/home/colddata/qinqiang/script/lib/KO_id_pathway_name.txt", header = F, sep = "\t")

go_enrich(id_list_file, output_dir)
kegg_enrich(id_list_file, kegg2gene, kegg2name, output_dir)
# suppressPackageStartupMessages(library(clusterProfiler))
# suppressPackageStartupMessages(library(openxlsx)) # 读取.xlsx文件
# suppressPackageStartupMessages(library(ggplot2)) # 柱状图和点状图
# suppressPackageStartupMessages(library(stringr)) # 基因ID转换
# suppressPackageStartupMessages(library(enrichplot)) # GO,KEGG,GSEA
# suppressPackageStartupMessages(library(clusterProfiler)) # GO,KEGG,GSEA
# suppressPackageStartupMessages(library(GOplot)) # 弦图，弦表图，系统聚类图
# suppressPackageStartupMessages(library(DOSE))
# suppressPackageStartupMessages(library(ggnewscale))
# suppressPackageStartupMessages(library(topGO)) # 绘制通路网络图
# suppressPackageStartupMessages(library(circlize)) # 绘制富集分析圈图
# suppressPackageStartupMessages(library(ComplexHeatmap)) # 绘制图例
# suppressPackageStartupMessages(library(optparse))
# suppressPackageStartupMessages(library(tidyr))

# option_list <- list(
#   make_option(c("--inputidfile"),
#     type = "character", default = NULL,
#     help = "输入 GeneID list 文件", metavar = "character"
#   ),
#   make_option(c("--genego"),
#     type = "character", default = NULL,
#     help = "提供 swiss 注释出来的 gene_go.txt 文件", metavar = "character"
#   ),
#   make_option(c("--keggclean"),
#     type = "character", default = NULL,
#     help = "提供 kegg 注释出来的 KEGG_clean.txt 文件", metavar = "integer"
#   ),
#   make_option(c("--outputdir"),
#     type = "character", default = NULL,
#     help = "输出结果文件的目录", metavar = "character"
#   )
# )
# opt_parser <- OptionParser(option_list = option_list)
# opt <- parse_args(opt_parser)

# # Check that all required arguments are provided
# if (!file.exists(opt$genego)) {
#   print_help(opt_parser)
#   stop("请提供 gene go 文件", call. = FALSE)
# } else if (!file.exists(opt$keggclean)) {
#   print_help(opt_parser)
#   stop("请提供 kegg 注释出来的 KEGG_clean.txt 文件", call. = FALSE)
# }

# if (!dir.exists(opt$outputdir)) {
#   dir.create(opt$outputdir)
#   cat("文件夹已创建：", opt$outputdir, "\n")
# } else {
#   cat("文件夹已存在，自动覆盖已有的文件：", opt$outputdir, "\n")
# }

# if (substr(opt$outputdir, nchar(opt$outputdir), nchar(opt$outputdir)) != "/") {
#   opt$outputdir <- paste0(opt$outputdir, "/")
# }

# # Assign the first argument to prefix
# id_list_file <- opt$inputidfile
# gene_go <- opt$genego
# kegg_clean <- opt$keggclean
# output_dir <- opt$outputdir

# # 计算比值的函数
# calculate_ratio <- function(ratio_string) {
#   parts <- strsplit(ratio_string, "/")[[1]] # 分割字符串
#   numerator <- as.numeric(parts[1]) # 转换分子为数值
#   denominator <- as.numeric(parts[2]) # 转换分母为数值
#   return(numerator / denominator) # 返回比值
# }

# go_enrich <- function(deg_id_file, output_dir) {
#   gene_list <- read.delim(deg_id_file, stringsAsFactors = FALSE, header = T)
#   names(gene_list)[1] <- c("gene_id")
#   gene_select <- gene_list$gene_id
#   go_rich_bp <- enricher(
#     gene = gene_select,
#     TERM2GENE = go_anno[go_anno$Ontology == "biological_process", ][c("ID", "gene_id")],
#     TERM2NAME = go_anno[go_anno$Ontology == "biological_process", ][c("ID", "Description")],
#     pvalueCutoff = 1,
#     pAdjustMethod = "BH",
#     qvalueCutoff = 1,
#     minGSSize = 1,
#     maxGSSize = 1000
#   )

#   go_rich_mf <- enricher(
#     gene = gene_select,
#     TERM2GENE = go_anno[go_anno$Ontology == "molecular_function", ][c("ID", "gene_id")],
#     TERM2NAME = go_anno[go_anno$Ontology == "molecular_function", ][c("ID", "Description")],
#     pvalueCutoff = 1,
#     pAdjustMethod = "BH",
#     qvalueCutoff = 1,
#     minGSSize = 1,
#     maxGSSize = 1000
#   )

#   go_rich_cc <- enricher(
#     gene = gene_select,
#     TERM2GENE = go_anno[go_anno$Ontology == "cellular_component", ][c("ID", "gene_id")],
#     TERM2NAME = go_anno[go_anno$Ontology == "cellular_component", ][c("ID", "Description")],
#     pvalueCutoff = 1,
#     pAdjustMethod = "BH",
#     qvalueCutoff = 1,
#     minGSSize = 1,
#     maxGSSize = 1000
#   )

#   value_counts = 0
#   if (length(go_rich_bp@result) != 0) {
#     bp_result <- go_rich_bp@result
#     bp_result$GeneRatio <- sapply(bp_result$GeneRatio, calculate_ratio)
#     bp_result$BgRatio <- sapply(bp_result$BgRatio, calculate_ratio)
#     bp_result$RichFactor <- bp_result$GeneRatio / bp_result$BgRatio
#     bp_result$Ontology <- rep("BP", nrow(bp_result))
#     # go_rich.bp <- as.data.frame(bp_result)
#     go_rich.total <- as.data.frame(bp_result)
#     value_counts = 1
#   }

#   if (length(go_rich_cc@result) != 0) {
#     cc_result <- go_rich_cc@result
#     cc_result$GeneRatio <- sapply(cc_result$GeneRatio, calculate_ratio)
#     cc_result$BgRatio <- sapply(cc_result$BgRatio, calculate_ratio)
#     cc_result$RichFactor <- cc_result$GeneRatio / cc_result$BgRatio
#     cc_result$Ontology <- rep("CC", nrow(cc_result))
#     # go_rich.cc <- as.data.frame(go_rich_cc)
#     if (is.data.frame(go_rich.total)) {
#       go_rich.total <- rbind(go_rich.total, cc_result)
#     } else {
#       go_rich.total <- as.data.frame(cc_result)
#     }
#     value_counts = 1
#   }

#   if (length(go_rich_mf@result) != 0) {
#     mf_result <- go_rich_mf@result
#     mf_result$GeneRatio <- sapply(mf_result$GeneRatio, calculate_ratio)
#     mf_result$BgRatio <- sapply(mf_result$BgRatio, calculate_ratio)
#     mf_result$RichFactor <- mf_result$GeneRatio / mf_result$BgRatio
#     mf_result$Ontology <- rep("MF", nrow(mf_result))
#     # go_rich.mf <- as.data.frame(go_rich_mf)
#     if (is.data.frame(go_rich.total)) {
#       go_rich.total <- rbind(go_rich.total, mf_result)
#     } else {
#       go_rich.total <- as.data.frame(mf_result)
#     }
#     value_counts = 1
#   }
  
#   # go_rich.total <- rbind(bp_result, cc_result, mf_result)

#   # 检测是否存在 go_rich.total 的一个 dataframe 且不为空
#   # TODO: 检测不到，还是会出错
#   if (value_counts == 0) {
#     print("没有富集到任何的GO信息\n")
#     return()
#   }

#   go_rich.total <- go_rich.total[, c("ID", "Description", "GeneRatio", "BgRatio", "RichFactor", "pvalue", "p.adjust", "qvalue", "geneID", "Count", "Ontology")]

#   go_enrich_file_name <- paste0(output_dir, sub("_ID.txt", "_EnrichmentGO.xlsx", basename(deg_id_file)))
#   cat("正在输出文件",go_enrich_file_name, "\n")
#   write.xlsx(go_rich.total, go_enrich_file_name)
# }

# kegg_enrich <- function(deg_id_file, kegg2gene, kegg2name, output_dir) {
#   ###### 导入目标基因列表################################
#   gene <- read.table(deg_id_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
#   gene <- as.factor(gene$V1)
#   # 富集分析
#   enrichment <- enricher(gene,
#     TERM2GENE = kegg2gene, TERM2NAME = kegg2name, pvalueCutoff = 1,
#     pAdjustMethod = "BH", qvalueCutoff = 1, minGSSize = 1, maxGSSize = 1000
#   )
#   if (length(enrichment@result) == 0) {
#     cat("没有富集到任何的KEGG信息\n")
#     return()
#   }
#   # # 绘制条形图, 绘制气泡图
#   # barplot(enrichment)
#   # dotplot(enrichment)
#   enrich_result <- enrichment@result
#   enrich_result$GeneRatio <- sapply(enrich_result$GeneRatio, calculate_ratio)
#   enrich_result$BgRatio <- sapply(enrich_result$BgRatio, calculate_ratio)
#   enrich_result$RichFactor <- enrich_result$GeneRatio / enrich_result$BgRatio
#   enrich_result$Ontology <- rep("KEGG", nrow(enrich_result))
#   enrich_result <- na.omit(enrich_result)  # 在 KO_id_pathway_name.txt 中没包含所有的 ko，有的对上了没有描述，所以要去掉这些
#   go_enrich_file_name <- paste0(
#     output_dir,
#     sub("_ID.txt", "_EnrichmentKEGG.xlsx", basename(deg_id_file))
#   )
#   enrich_result <- enrich_result[, c("ID", "Description", "GeneRatio", "BgRatio", "RichFactor", "pvalue", "p.adjust", "qvalue", "geneID", "Count", "Ontology")]
#   write.xlsx(enrich_result, go_enrich_file_name)
#   # enrichplot::cnetplot(go_rich_bp,circular=FALSE,colorEdge = T,edge = F,color_category = "red", color_gene = "green",cex_category = 1,cex_label_category = 1)
# }

# # 加载GO注释描述文件,GOid\tdescriptiom\tGOclass(BP CC MF)
# go_class <- read.delim(
#   "/home/colddata/qinqiang/script/lib/go_term.list",
#   header = FALSE,
#   stringsAsFactors = FALSE
# )
# names(go_class) <- c("ID", "Description", "Ontology")

# # 加载背景库文件，gene\tGOid
# go_anno <- read.delim(gene_go, header = FALSE, stringsAsFactors = FALSE)
# names(go_anno) <- c("gene_id", "ID")
# # 合并背景与GO文件
# go_anno <- merge(go_anno, go_class, by = "ID", all.x = TRUE)

# kegg_clean <- read.table(kegg_clean,
#   header = F, sep = "\t", quote = "",
#   col.names = c("GeneID", "KEGG_Pathway", "Metabolism_Category", "General_Metabolism", "K_Number", "Protein", "EC_Number")
# )
# kegg2gene <- kegg_clean[c("KEGG_Pathway", "GeneID")]
# kegg2gene$KEGG_Pathway <- sub(":.*", "", kegg2gene$KEGG_Pathway)
# kegg2name <- read.table("/home/colddata/qinqiang/script/lib/KO_id_pathway_name.txt", header = F, sep = "\t")

# go_enrich(id_list_file, output_dir)
# kegg_enrich(id_list_file, kegg2gene, kegg2name, output_dir)
