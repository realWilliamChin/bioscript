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
suppressPackageStartupMessages(library(ggthemes))
suppressPackageStartupMessages(library(dplyr))


option_list <- list(
  make_option(c("--datatable"),
    type = "character", help = "输入文件（ sample ID 是列名，C number 是行名）", metavar = "character"
  ),
  make_option(c("--outputprefix"),
    type = "character", help = "输出文件 .xlsx 结尾的结果", metavar = "character"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Assign the first argument to prefix
data_table <- opt$datatable
output_prefix <- opt$outputprefix
group_output_dir <- dirname(output_prefix)
group_name <- basename(output_prefix)
enrichment_graphs_dir <- file.path(group_output_dir, 'Enrichment_Graphs')
m_graphs_dir <- file.path(group_output_dir, paste0(group_name, '_M_Graphs'))

dir.create(m_graphs_dir)
dir.create(enrichment_graphs_dir)


metabolite_enrich <- function(data_table, output_prefix) {
  compound_matrix_df <- read.delim('/home/colddata/qinqiang/script/MetaboliteAnalysis/MetaboliteEnrich/Compound_matrix.txt', sep="\t", header=T, check.names = F, stringsAsFactors = F)

  # 计算比值的函数
  calculate_ratio <- function(ratio_string) {
    parts <- strsplit(ratio_string, "/")[[1]] # 分割字符串
    numerator <- as.numeric(parts[1]) # 转换分子为数值
    denominator <- as.numeric(parts[2]) # 转换分母为数值
    return(numerator / denominator) # 返回比值
  }

  gene_list <- read.delim(data_table, stringsAsFactors = FALSE, header = T)
  names(gene_list)[1] <- c("c_number")
  gene_select <- gene_list$c_number
  
  # 标准化基因ID格式（统一转换为大写，去除空格）
  gene_select <- toupper(trimws(gene_select))
  gene_select <- gene_select[!is.na(gene_select) & gene_select != ""]
  
  # 标准化compound_matrix_df中的c_number格式
  compound_matrix_df$c_number <- toupper(trimws(compound_matrix_df$c_number))
  
  # 检查有多少基因ID可以匹配
  matched_genes <- intersect(gene_select, compound_matrix_df$c_number)
  unmatched_genes <- setdiff(gene_select, compound_matrix_df$c_number)
  
  if (length(matched_genes) == 0) {
    warning(paste0("错误：输入的基因ID无法在Compound_matrix.txt中找到匹配。\n",
                   "输入的基因ID: ", paste(gene_select, collapse = ", "), "\n",
                   "请检查基因ID格式是否正确（应为C开头的化合物编号，如C00279）。"))
    message(paste0(output_prefix, "：没有可匹配的基因ID，跳过富集分析。"))
    return(invisible(NULL))
  }
  
  if (length(unmatched_genes) > 0) {
    warning(paste0("警告：以下基因ID无法匹配: ", paste(unmatched_genes, collapse = ", "), "\n",
                   "已匹配的基因ID数量: ", length(matched_genes), "/", length(gene_select)))
  }
  
  # 只使用可以匹配的基因ID进行富集分析
  gene_select_matched <- matched_genes

  go_rich <- enricher(
    gene = gene_select_matched,
    TERM2GENE = compound_matrix_df[c("entry", "c_number")],
    TERM2NAME = compound_matrix_df[c("entry", "name")],
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    qvalueCutoff = 1,
    minGSSize = 1,
    maxGSSize = 1000
  )

  # 检查enricher是否返回NULL
  if (is.null(go_rich)) {
    warning(paste0("错误：enricher()返回NULL，无法进行富集分析。\n",
                   "输入的基因ID: ", paste(gene_select_matched, collapse = ", "), "\n",
                   "请检查Compound_matrix.txt文件格式是否正确。"))
    message(paste0(output_prefix, "：富集分析失败，跳过后续步骤。"))
    return(invisible(NULL))
  }
  
  # 检查结果是否为空
  if (nrow(go_rich@result) == 0) {
    warning(paste0("警告：富集分析完成，但没有找到任何富集结果。\n",
                   "输入的基因ID: ", paste(gene_select_matched, collapse = ", ")))
    message(paste0(output_prefix, "：没有富集结果，跳过后续步骤。"))
    return(invisible(NULL))
  }

  go_rich_result <- go_rich@result
  go_rich_result$GeneRatio <- sapply(go_rich_result$GeneRatio, calculate_ratio)
  go_rich_result$BgRatio <- sapply(go_rich_result$BgRatio, calculate_ratio)
  go_rich_result$RichFactor <- go_rich_result$GeneRatio / go_rich_result$BgRatio

  go_rich_result$entry <- go_rich_result$ID

  # 去掉 entry c_number name 列
  go_rich_result <- left_join(go_rich_result, compound_matrix_df, by = "entry")
  go_rich_result <- go_rich_result[ , !(names(go_rich_result) %in% c("entry", "c_number", "name"))]
  # 去重
  go_rich_result <- distinct(go_rich_result)
  go_rich_result <- rename(go_rich_result, Compound_ID = geneID, KEGG_Ortholog = kog)

  # 创建文件夹
  pic_path = '/home/colddata/qinqiang/script/MetaboliteAnalysis/MetaboliteEnrich/M_entry_pictures/'

  write.xlsx(go_rich_result, file.path(m_graphs_dir, paste0(group_name, "_enrich_result.xlsx")))

  output_m_id_df <- go_rich_result[go_rich_result$pvalue < 0.05, ]
  # 循环 go_rich_result$ID
  for (i in output_m_id_df$ID) {
    cmd = paste0("cp -r /home/colddata/qinqiang/script/MetaboliteAnalysis/MetaboliteEnrich/M_entry_pictures/", i, " ", m_graphs_dir, sep='')
    system(cmd)
  }

  #dt <- read.xlsx(output_file, sheet = 1)
  dt <- go_rich_result
  dt$MID <- paste(dt$ID, dt$Description, sep='_')

  dt <- dt[order(dt$pvalue), ]
  dt <- dt[dt$Count >= 1 & dt$pvalue < 0.05, ]
  dt <- head(dt[order(dt$pvalue), ], n = 15)
  dt <- dt[order(-dt$pvalue), ]
  dt[, c('ID', 'Count', 'pvalue')]
  if (nrow(dt) > 0) {
    p <- ggplot(dt, aes(y = factor(MID, levels = MID), x = RichFactor, size = Count, colour = pvalue)) +
      geom_point() +
      scale_size_continuous(range = c(2, 8)) +
      scale_y_discrete(limits = dt$MID) +
      scale_colour_gradient(low = "red", high = "green") +
      ggtitle(paste0(group_name, "_Enrich_Bubble")) +
      ylab("Metabolite Entry") +
      theme_base()
    ggsave(file.path(enrichment_graphs_dir, paste0(group_name, "_enrich_Bubble.png")), p, dpi = 320, width = 20, height = 10)
  } else {
    message(paste0(output_prefix, "没有显著富集的。"))
  }

  dt <- head(dt[order(dt$pvalue), ], n = 10)
  dt <- dt[order(dt$ID), ]

  dt$star <- cut(dt$pvalue, breaks = c(0, 0.001, 0.01, 0.05, Inf), 
    labels = c("***        ", "**     ", "*   ", ""))
  bar.p <- ggplot(dt, aes(x= factor(MID, levels = MID), y=RichFactor)) + 
    geom_bar(stat="identity", position = "identity", fill='skyblue')+
    geom_text(aes(label = star), position = position_dodge(width = 1), size = 7) +  # 添加星号
    # geom_bar(position="identity",stat="identity",aes(fill=SubOntology))+
    # facet_grid(SubOntology~.,scale="free")+
    # theme_bw()+
    xlab("MID Term") +  # 设置 x 轴标签
    ylab("RichFactor") +    # 设置 y 轴标签
    coord_flip() +
    ggtitle(paste0(group_name, "_Enrich_bargraph")) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 14, face = "bold"),
      legend.title = element_text(face = "bold"))
  # height 根据 datafile 的行数调整
  p_height <- 2 + 0.4 * nrow(dt)
  ggsave(file.path(enrichment_graphs_dir, paste0(group_name, "_enrich_barplot.png")), bar.p, dpi = 320, width = 20, height = p_height)
}

metabolite_enrich(data_table, output_prefix)
