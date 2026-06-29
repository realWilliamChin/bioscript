20260622
r443
R
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)
library(pheatmap)
library(ggpubr)
library(tibble)
library(stringr)
library(harmony)
library(scGate)
library(reticulate)
library(celldex)
library(SummarizedExperiment)
library(SingleR)

library(clusterProfiler)
library(patchwork)
library(org.Hs.eg.db)  # 人类数据库，如果是小鼠改为 org.Mm.eg.db
library(AnnotationDbi)
library(enrichplot)
library(RColorBrewer)
library(ggrepel)
library(tidyverse)
library(readr)
library(gridExtra)  # 用于组合多个图形
library(grid)       # 用于自定义图形布局
library(ggthemes)

#一、h5ad文件转换为rds文件及其保存
setwd("/home/colddata/qinqiang/Project/Test/singleCellPipeline_test/pipelineTest1/")
use_python("~/miniconda3/envs/test_scanpy_python311/bin/python")  # 指定 Python 路径（如 conda 环境）
ad <- import("anndata")

adataRDS_dir <- file.path(getwd(), 'adataRDS')
adata<- ad$read_h5ad(adataRDS_dir, 'adataFiltrd_celltypist.h5ad'))
counts <- t(adata$raw$X) # 转置矩阵（Seurat 需要 基因×细胞）
#dim(counts)
#[1]  33282 138672

# 提取细胞和基因的元数据
cell_meta <- as.data.frame(adata$obs)
gene_meta <- as.data.frame(adata$var)
rownames(counts) <- rownames(gene_meta) # 基因名作为行名
colnames(counts) <- rownames(cell_meta)  # 细胞ID作为列名

seurat_obj <- CreateSeuratObject(counts = counts, meta.data = as.data.frame(adata$obs))
dim(seurat_obj@meta.data)
[1] 138672     43

saveRDS(seurat_obj, file = file.path(adataRDS_dir, 'adataFiltrd_celltypist.rds'))
#seurat_obj <- readRDS(file.path(adataRDS_dir, 'adataFiltrd_celltypist.rds'))

#二、‘04_样本或组汇总’部分的基因富集
Idents(seurat_obj) <- seurat_obj$group
# 1. 标准化（必须！生成 data 层）
seurat_obj <- NormalizeData(seurat_obj)
# 2. 找高可变基因（必须）
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

timestamp <- format(Sys.time(), "%Y%m%d")
unique_groups <- unique(unique(seurat_obj@meta.data$group))
group_pairs <- combn(unique_groups, 2)
pair_list <- lapply(1:ncol(group_pairs), function(i) group_pairs[, i])
cleaned_kegg <- readRDS("/home/colddata/qinqiang/Project/2025_Project/2025_12_26_昆明第一医院_人_皮肤样本_单细胞/zhgq/kegg_hierarchy_20260131.rds")

enrichGO_plotting_20260205 <- function(enrichGO_topN, comb, direction, workDir) {
	new_comb <- str_replace(comb, ": ", "_")
	new_comb <- str_replace_all(new_comb, " ", "-")
	
	combName=paste(new_comb, direction, sep='_')
	color_map <- c("BP" = "#33a02c", "CC" = "#e31a1c", "MF" = "#1f78b4")
	
	# 绘制横向条形图（分面展示BP/CC/MF，y轴包含截断Term+基因数量）
	enrichGO_plot <- ggplot(enrichGO_topN, aes(x = log_pval, y = reorder(Description, log_pval), fill = ONTOLOGY)) +
	geom_col(width = 0.7) +  # 条形宽度
	facet_wrap(~ONTOLOGY, ncol = 1, scales = "free_y") +  # 单列分面，y轴自由缩放
	scale_fill_manual(values = color_map) +
	labs(x = "-log10(p-value)", y = "GO Term",
		title = paste0(comb, " (", direction, ")\n Top 10 GO Term Enrichment")
	) +
	theme_bw() +
	theme(
	panel.grid = element_blank(),        # 去除网格线
	axis.text.y = element_text(size = 10),# 调整y轴字体大小（适配截断后标签）
	axis.text.x = element_text(size = 10),# x轴标签大小
	strip.text = element_text(size = 12, face = "bold"),# 分面标题
	legend.position = "none",           # 隐藏图例（分面已区分）
	plot.title = element_text(hjust = 0.5, size = 14, face = "bold")# 标题居中
	)
		
	ggsave(file.path(workDir, paste0(combName, '_GO.pdf')), enrichGO_plot, width = 14, height = 10, dpi = 300)
	ggsave(file.path(workDir, paste0(combName, '_GO.png')), enrichGO_plot, width = 14, height = 10, dpi = 300)
}

extract_pathway_name <- function(pathway_string) {
  # 使用正则表达式移除常见的后缀模式（如 [PATH:ko...] 或 (PATH:ko...)）
  # 同时也会处理掉字符串开头可能存在的数字和空格
  cleaned_name <- gsub("^\\d+\\s+|\\s*\\[[^]]*\\]$|\\s*\\([^)]*\\)$", "", pathway_string)
  return(trimws(cleaned_name))
}

enrichGO_topN_picking_20260205 <- function(goEnrich, top_n) {
	go_bp_top10 <- goEnrich@result %>%
	filter(ONTOLOGY == "BP") %>%
	arrange(p.adjust) %>%
	slice_head(n = top_n)
	
	go_cc_top10 <- goEnrich@result %>%
	filter(ONTOLOGY == "CC") %>%
	arrange(p.adjust) %>%
	slice_head(n = top_n)
	
	go_mf_top10 <- goEnrich@result %>%
	filter(ONTOLOGY == "MF") %>%
	arrange(p.adjust) %>%
	slice_head(n = top_n)
	
	# 合并三分类结果并计算-log10(p.adjust)
	go_all_top10 <- bind_rows(go_bp_top10, go_cc_top10, go_mf_top10) %>%
	mutate(log_pval = -log10(p.adjust))
	
	return(go_all_top10)
}

markerFinding_combSeuratObjExtracting <- function(celltypeSeuratObj, cond1, cond2) {
	combSeuratInfosLst <- list()
	cellType_degs <- FindMarkers(object = celltypeSeuratObj,
		ident.1 = cond1,
		ident.2 = cond2,
		min.pct = 0.1,
		logfc.threshold = 0.25,
		test.use = "wilcox",
		only.pos = FALSE,
		verbose = FALSE
	)
	combSeuratInfosLst$degs <- cellType_degs

	cond1_celltype_cells <- rownames(celltypeSeuratObj@meta.data)[celltypeSeuratObj$group == cond1]
	cond2_celltype_cells <- rownames(celltypeSeuratObj@meta.data)[celltypeSeuratObj$group == cond2]
	comb_cells <- c(cond1_celltype_cells, cond2_celltype_cells)
	comb_seuratObj <- subset(celltypeSeuratObj, cells = comb_cells)
	combSeuratInfosLst$seuratObj <- comb_seuratObj
	
	return(combSeuratInfosLst)
}

kegg_enrichPlotting_20260205 <- function(genes, kegg_hierarchy, comb, direc, workDir,
	universe_genes = NULL, pvalue_cutoff = 0.05, qvalue_cutoff = 0.2) {
	gene2ko <- kegg_hierarchy %>%
		filter(!is.na(Gene_Symbol) & !is.na(KO_ID)) %>%
		dplyr::select(Gene_Symbol, KO_ID) %>%
		distinct()
	
	ko2pathway <- kegg_hierarchy %>%
		dplyr::select(KO_ID, Level_C) %>%
		distinct() %>%
		mutate(Pathway_ID = paste0("ko", str_extract(Level_C, "\\d+")), Pathway_Description = Level_C
		) %>%
		filter(!is.na(Pathway_ID))
	
	# 从KEGG层级数据中提取通路信息
	pathway_info <- kegg_hierarchy %>%
		mutate(Pathway_Num = str_extract(Level_C, "\\d+")) %>%
		filter(!is.na(Pathway_Num)) %>%
		dplyr::select(Pathway_Num, Level_A, Level_B, Level_C) %>%
		distinct()
	
	# 执行富集分析
	kegg_enrich <- perform_kegg_enrichment(genes = rownames(genes),
		gene2ko = gene2ko, ko2pathway = ko2pathway,
		#universe_genes = gene_list$gene_id,
		pvalue_cutoff = 0.05, qvalue_cutoff = 0.2)
	
	if (!is.null(kegg_enrich)) {
	   kegg_enrich_df <- as.data.frame(kegg_enrich)
	   kegg_enrich_df <- kegg_enrich_df %>%
	   	mutate(Pathway_Num = str_extract(ID, "\\d+"),Pathway_ID = ID)
	   
	   timestamp <- format(Sys.time(), "%Y%m%d")
	   new_string <- str_replace(comb, ": ", "_")
	   new_string <- str_replace_all(new_string, " ", "-")
	   combName=paste(new_string, direc, sep='_')
	   
	   if (dim(kegg_enrich_df)[1] ==0) {
		  write.csv(kegg_enrich_df, file.path(workDir, paste0(paste(combName, "KEGG", timestamp, sep='_'), ".csv")), row.names=FALSE)
	      } else {
	      # 合并富集结果与层级信息
	      kegg_enrich_df <- kegg_enrich_df %>%
		left_join(pathway_info, by = "Pathway_Num") %>%
		# 计算富集分数
		mutate(GeneRatio_num = as.numeric(str_split(GeneRatio, "/", simplify = TRUE)[,1]) /
		  as.numeric(str_split(GeneRatio, "/", simplify = TRUE)[,2]),
		BgRatio_num = as.numeric(str_split(BgRatio, "/", simplify = TRUE)[,1]) /
		  as.numeric(str_split(BgRatio, "/", simplify = TRUE)[,2]),
		Fold_Enrichment = GeneRatio_num / BgRatio_num,
		Rich_Factor = Count / as.numeric(str_split(BgRatio, "/", simplify = TRUE)[,2])
		  ) %>%
		# 排序
		arrange(p.adjust) %>% mutate(Rank = row_number()) # 添加排名
		
	      kegg_enrich_df <- kegg_enrich_df %>%
		mutate(Pathway_Name = sapply(Description, extract_pathway_name)) %>%
		mutate(LevelA = sapply(Level_A, extract_pathway_name)) %>%
		arrange(pvalue)
	
	     write.csv(kegg_enrich_df[colnames(kegg_enrich_df)[-seq(13,22)]], file.path(workDir, paste0(paste(combName, "KEGG", timestamp, sep='_'), ".csv")), row.names=FALSE)
	
	     top20_enrich_df <- kegg_enrich_df %>%
			arrange(p.adjust) %>% head(20)
	     	
	     p_complete <- ggplot(top20_enrich_df, 
				aes(x = GeneRatio_num, y = reorder(Pathway_Name, GeneRatio_num))) +
				geom_point(aes(size = Count, color = pvalue), alpha = 0.8) +
				scale_color_gradientn(colors = rev(brewer.pal(11, "Spectral")), name = "P-value", guide = guide_colorbar(reverse = TRUE, order = 1)) +
				scale_size_continuous(range = c(3, 8), name = "Gene Count",
		                breaks = scales::pretty_breaks(n = 4)(range(top20_enrich_df$Count)),
			        guide = guide_legend(order = 2)
		             ) +
		           theme_bw(base_size = 12) +
		           theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
			      axis.text.y = element_text(size = 10),
			      axis.title = element_text(size = 12, face = "bold"),
			      strip.text = element_text(size = 11, face = "bold", color = "white"),
			      strip.background = element_rect(fill = "#2E86AB", color = "#2E86AB"),
			      legend.position = "right",
			      legend.box = "vertical",
			      legend.title = element_text(face = "bold"),
			      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
			      plot.subtitle = element_text(hjust = 0.5, size = 12),
			      panel.grid.major = element_line(color = "gray90", size = 0.3),
			      panel.grid.minor = element_blank()
			   ) +
		           labs(title = paste("Top", 20, "KEGG Pathways by Significance"),
			   #subtitle = paste("Organism:", 'hsa'),
			    x = "Gene Ratio", y = "Pathway Description"
			   ) + ggforce::facet_col(~ LevelA, scales = "free_y", space = "free")
	
	      fileName <- paste0("top", 20, "_KEGG_", timestamp)
	      ggsave(filename = file.path(workDir, paste0(paste(combName, fileName, sep='_'), ".pdf")), width = 14, height = 10, dpi = 300)
	      ggsave(filename = file.path(workDir, paste0(paste(combName, fileName, sep='_'), ".png")), width = 14, height = 10, dpi = 300)
	   }
	}
}

perform_kegg_enrichment <- function(genes, gene2ko, ko2pathway, 
      universe_genes = NULL, 
      pvalue_cutoff = 0.05,
      qvalue_cutoff = 0.2) {
  # 将基因转换为KO
  gene_ko_mapping <- gene2ko %>%
    filter(Gene_Symbol %in% genes)
  
  ko_list <- unique(gene_ko_mapping$KO_ID)
  cat("输入基因对应的KO数量:", length(ko_list), "\n")
  
  # 准备背景基因集（如果提供）
  if (!is.null(universe_genes)) {
    universe_ko <- gene2ko %>%
      filter(Gene_Symbol %in% universe_genes) %>%
      pull(KO_ID) %>%
      unique()
    cat("背景KO数量:", length(universe_ko), "\n")
  } else {
    universe_ko <- NULL
  }
  
  # 执行富集分析
  enrich_result <- enricher(
    gene = ko_list,
    pvalueCutoff = pvalue_cutoff,
    pAdjustMethod = "BH",
    qvalueCutoff = qvalue_cutoff,
    universe = universe_ko,
    minGSSize = 10,
    maxGSSize = 2000,
    TERM2GENE = ko2pathway %>% dplyr::select(Pathway_ID, KO_ID),
    TERM2NAME = ko2pathway %>% dplyr::select(Pathway_ID, Pathway_Description)
  )
  return(enrich_result)
}

#此次新加部分！！！
seuratObjExtracting_20260618 <- function(seuratObj, targetLst) {
	grpsObj <-seuratObj %>%
		subset(
		subset = group %in% pair
	)		
	Idents(grpsObj) <- "group"
	
	if (is.null(grpsObj[["RNA"]]$data) || nrow(grpsObj[["RNA"]]$data) == 0) {
		cat("data层为空，正在进行数据归一化...\n")
		grpsObj <- NormalizeData(grpsObj, 
                              normalization.method = "LogNormalize", 
                              scale.factor = 10000,
                              verbose = FALSE)
		cat("数据归一化完成\n")
	}
	return(grpsObj)
}

heatmap_plotting_20260622 <- function(seurat_obj, comb, pair, topN_genes) {
	timestamp <- format(Sys.time(), "%Y%m%d")
	#new_string <- str_replace(comb, ": ", "_")
	#new_string <- str_replace_all(new_string, " ", "-")
	#combName=paste(new_string, sep='_')
	#combName <- str_replace(combName, "_", ": ")

	#seurat_obj <- ScaleData(seurat_obj, features =  rownames(seurat_obj))
	seurat_obj@meta.data$group <-  factor(seurat_obj@meta.data$group, levels = unlist(pair))  #重新分组去掉原先的levels的痕迹，很重要！！！
	cell_keep <- c(
		WhichCells(comb_seuratObj, expression = group == pair[1]),
		WhichCells(comb_seuratObj, expression = group == pair[2])
	)
	
	p <- DoHeatmap(seurat_obj, group.by = "group", slot = "scale.data",  cells = cell_keep, size = 3, angle = 0, features = topN_genes) + scale_fill_gradientn(colors = c("lightblue", "firebrick3")) +theme(plot.title = element_text(size = 14)) + labs(title = comb)
	return(p)
}

volcano_plotting_20260623 <- function(degs_data, comb, workDir) {	
	degs_data$regulation <- ifelse(degs_data$p_val_adj < 0.05 & abs(degs_data$avg_log2FC) > 0.5,
           ifelse(degs_data$avg_log2FC > 0.5, "Up", "Down"),
                  "NoSignificant")
	all_up <- degs_data[degs_data$regulation == "Up", ]
	all_down <- degs_data[degs_data$regulation == "Down",]
	degs_data[degs_data$regulation == "Down",]
	up_total <- nrow(all_up)
	down_total <- nrow(all_down)
	# 提取Top10最显著上调/下调（按padj从小到大，显著性最高）
	top10_up <- all_up[order(all_up$p_val_adj), ][1:10, ]
	top10_down <- all_down[order(all_down$p_val_adj), ][1:10, ]
	label_gene <- rbind(top10_up, top10_down)
	
	# 修正极小padj，避免无穷大绘图报错
	degs_data$p_val_adj <- ifelse(degs_data$p_val_adj < 1e-300, 1e-300, degs_data$p_val_adj)
	max(degs_data$avg_log2FC)
	FC_max<-ceiling(max(abs(degs_data$avg_log2FC)))
	
	# 绘图
	p <- ggplot(data = degs_data, aes(x = avg_log2FC, y = -log10(p_val_adj), colour = regulation)) +
		geom_point(alpha = 0.7, size = 1.5) +
		geom_vline(xintercept = c(-0.5, 0.5), lty = 4, col = "black", lwd = 0.8) +
	geom_hline(yintercept = -log10(0.05), lty = 4, col = "black", lwd = 0.8) +
	# 现在up_total/down_total已提前定义，不会报错，图例正常渲染
	scale_color_manual(values = c("green", "gray60", "red"),
              labels = c(paste0("Up(", up_total, ")"),
              "Not significant",
		paste0("Down(", down_total, ")"))) +
		# 标注基因（仅隐藏标签图例，不影响主颜色图例）
		geom_label_repel(data = label_gene,
                   aes(label = Gene_Symbol),
                   size = 2.5,
                   max.overlaps = 100,
                   box.padding = 0.3,
                   label.padding = 0.1,
                   segment.color = "black",
                   show.legend = FALSE) +
		labs(x = "log2(Fold Change)",
		y = "-log10(Adjusted P-value)",
		title = comb,
		color = "Gene expression status") +scale_x_continuous(limits = c(-FC_max, FC_max))+
		theme_bw() + # theme_base报错替换为theme_bw()
		# 去掉全部网格线
	theme(
	panel.grid.major = element_blank(),  # 主网格线
	panel.grid.minor = element_blank(),   # 次网格线
	plot.title = element_text(hjust = 0.5),
	legend.position = "right",   
	legend.key.size = unit(1, "cm"), # 图例点大小微调
	legend.key.height = unit(0.5, "cm"),
	legend.spacing.y = unit(0, "cm")  
	)
	
	timestamp <- format(Sys.time(), "%Y%m%d")
	new_string <- str_replace(comb, ": ", "_")
	new_string <- str_replace_all(new_string, " ", "-")
	combName=paste(new_string, sep='_')
	
	ggsave(filename = file.path(workDir, paste0(paste(combName, 'volcano', timestamp, sep='_'), ".png")), plot = p, width = 8, height = 6, units = "in")
}

setwd('01_Primary_clustering')
directions <- c('Up', 'Down', 'Total')
diffDir=file.path(getwd(), '04_样本或组汇总/diffExpre')

for (pair in pair_list) {
	pairDir <- file.path(diffDir,  paste(str_replace(pair[1], ' ', ''), str_replace(pair[2], ' ', ''), sep='_vs_'))
	if (!dir.exists(pairDir)) {
		dir.create(pairDir, recursive = TRUE)
	} else {
		print("目录已存在。")
	}
	print(pairDir)

	targetObj <- seuratObjExtracting_20260618(seurat_obj, pair)
	combSeuratObj <- markerFinding_combSeuratObjExtracting(targetObj, pair[1], pair[2])
	
	write.csv(rownames_to_column(combSeuratObj$degs, var = "Gene_Symbol"), file.path(pairDir, paste0(paste('diffExpre', timestamp, sep='_'), '.csv')), row.names=FALSE)
	#saveRDS(combSeuratObj$seuratObj, file = file.path(pairDir, paste0(paste('seuratObj', timestamp, sep='_'), ".rds")))

	for (e in directions) {
	   degs <- combSeuratObj$degs %>%
	          mutate(
			direction = case_when(
			p_val_adj < 0.05 & avg_log2FC > 0.5 ~ "Up",
			p_val_adj < 0.05 & avg_log2FC < -0.5 ~ "Down",
			TRUE ~ "Not significant"
			),
			significance = ifelse(p_val_adj < 0.05, "Significant", "Not significant")
		)
		
	if (e != "Total") {
		genes <- degs %>%
			filter(direction == e)
		} else {
			genes <- degs
		}

	top20_down <- degs %>%
		filter(direction == "Down") %>%
		arrange(p_val_adj) %>% head(20)
	
	top20_up <- degs %>%
		filter(direction == "Up") %>%
		arrange(p_val_adj) %>% head(20)
	
	comb <- paste(str_replace(pair[1], ' ', ''), str_replace(pair[2], ' ', ''), sep=' vs ')
	
	comb_seuratObj <- ScaleData(combSeuratObj$seuratObj, features =  rownames(combSeuratObj$seuratObj))
	#避免运行heatmap绘图函数时报如下错误：Error in DoHeatmap(seurat_obj, features = topN_genes, group.by = "group",  :No requested features found in the scale.data layer for the RNA assay.
	
	timestamp <- format(Sys.time(), "%Y%m%d")
	pt <- heatmap_plotting_20260622(comb_seuratObj, comb, pair, c(rownames(top20_up), rownames(top20_down)))
	ggsave(filename = file.path(pairDir, paste0(str_replace_all(comb, ' ', '-'), '_heatmap_' , timestamp, '.png')), plot = pt, width = 8, height = 6, units = "in", dpi=300)
	degs <- rownames_to_column(degs, var = "Gene_Symbol")
	volcano_plotting_20260623(degs, comb, pairDir)
	
	for (e in directions) {
		print(c(e, dim(genes)))
		if (dim(genes)[1] !=0) {
			gene_ids <- tryCatch({
			gene_ids <- mapIds(org.Hs.eg.db, keys = rownames(genes), column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
			}, 
			error = function(e) {
			message("捕获到一个错误: ", e$message)
				return(NULL)
				}
			)
			    
			if (!is.null(gene_ids)) {
				goEnrich <- enrichGO(gene = gene_ids,
				OrgDb = org.Hs.eg.db,
				keyType = "ENTREZID",
				ont = "ALL",          # 富集类型：CC=细胞组分
				pAdjustMethod = "BH",
				qvalueCutoff = 0.05,
				readable = TRUE
				)
				
				GO_dir <- file.path(pairDir, 'GO_enrichment1')
				dir.create(GO_dir)
				print(c('goEnrich', dim(genes)))
				timestamp <- format(Sys.time(), "%Y%m%d")
				#new_string <- str_replace(comb, ": ", "_")
				#new_string <- str_replace_all(new_string, " ", "-")
				combName=paste(comb, e, sep='_')
				write.csv(as.data.frame(goEnrich), file.path(GO_dir, paste0(paste(combName, 'GO', timestamp, sep='_'), ".csv")), row.names=FALSE)
				
				if (dim(goEnrich)[1] >0) {
					enrichGO_top10 = enrichGO_topN_picking_20260205(goEnrich, 10)
						if (dim(enrichGO_top10)[1] >0) {
						enrichGO_plotting_20260205(enrichGO_top10, comb, e, GO_dir)
						}
					}
				KEGG_dir <- file.path(pairDir, 'KEGG_enrichment1')
				dir.create(KEGG_dir)
				kegg_enrichPlotting_20260205(genes, cleaned_kegg, comb, e, KEGG_dir)
			}
		}
	}
}


#三、细胞通讯
library(CellChat)
library(Seurat)
library(patchwork)
library(future)
library(NMF)
library(ggalluvial)
library(ComplexHeatmap)
library(stringr)
library(magick)

Groups <- c('HD', 'NoLN', 'Class_III', 'Class_III-IV', 'Class_IV-V')
targetClusters <- c('0', '1', '2', '3', '4', '6')

cluster2annot_mapping <- function(df) {
  df <- as.data.frame(df)
  annot_map <- df[[2]]          # 第2列 = 注释名
  names(annot_map) <- df[[1]]  # 第1列 = 聚类编号
  return(annot_map)
}

# 你的注释表（2列：leiden_r025 + 注释）
annot_df <- data.frame(
  cluster = c("0", "2", "1", "3", "4", "6"),
  label   = c("Tem/Tcm", "Tem/Tcm", "Naive T", "Tem", "Naive T/Tcm", "Effector_CD4/Temra/Tem")
)

# 自动生成 annot_map
cluster2annot_map <- cluster2annot_mapping(annot_df)

# 核心：提取聚类并转换
clusters <- as.character(grpObj$leiden_r025)  # 转为字符
annot_result <- cluster2annot_map[clusters]   # 匹配注释

# 关键：保留细胞名称（解决Seurat报错）
names(annot_result) <- colnames(grpObj)

# 用 Seurat 官方函数添加
grpObj <- AddMetaData(grpObj, metadata = annot_result, col.name = "r025_annot")


RDS_dir <- '/home/colddata/qinqiang/Project/Test/singleCellPipeline_test/pipelineTest1/adataRDS'

object.list <- list()
cellchatLst <- list()
for (grp in Groups) {
	if (!dir.exists(RDS_dir)) {
			dir.create(RDS_dir, recursive = TRUE)
		} else {
			print("目录已存在。")
		}
	
	grpObj <-seurat_obj %>%
		subset(
		subset = leiden_r025 %in% targetClusters & group == grp
	)
	
	clusters <- as.character(grpObj$leiden_r025)  # 转为字符
	annot_result <- cluster2annot_map[clusters]   # 匹配注释
	names(annot_result) <- colnames(grpObj)
	grpObj <- AddMetaData(grpObj, metadata = annot_result, col.name = "r025_annot")
			
	grpObj$r025_annot <- factor(grpObj$r025_annot)
	meta <- grpObj@meta.data
	data.input <- GetAssayData(grpObj, assay = "RNA", layer = "counts")
	common_cells <- intersect(colnames(data.input), rownames(meta))
	data.input <- data.input[, common_cells]
	meta <- meta[common_cells, ]
	labels <- grpObj$r025_annot
	meta <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe of the cell labels
	cellchat <- createCellChat(object = data.input, meta = meta)
	cellchat@idents <- factor(grpObj$r025_annot)

	CellChatDB <- CellChatDB.human
	CellChatDB.use <- CellChatDB
	cellchat@DB <- CellChatDB.use
	cellchat <- subsetData(cellchat)
	cellchat <- identifyOverExpressedGenes(cellchat)
	
	cellchat <- identifyOverExpressedInteractions(cellchat)
	cellchat <- computeCommunProb(cellchat)
	
	cellchat <- filterCommunication(cellchat, min.cells = 10)
	cellchat <- computeCommunProbPathway(cellchat)
	cellchat <- aggregateNet(cellchat)
	
	saveRDS(cellchat, file = file.path(RDS_dir, paste(grp, sep='_', 'cellchat.rds')))
	pathways.show.all <- cellchat@netP$pathways
	length(pathways.show.all)
	
	#write.csv(pathways.show.all, file.path(cellchat_dir, 'pathways.csv'), row.names = FALSE, fileEncoding = "UTF-8")
	
	object.list[[grp]] <- grpObj
	cellchatLst[[grp]] <- cellchat
}

cellchat_dir <- '/home/colddata/qinqiang/Project/Test/singleCellPipeline_test/pipelineTest1/05_Cell_communication_analysis'
dir.create(file.path(cellchat_dir, 'Graph_data'), recursive = TRUE)
dir.create(file.path(cellchat_dir, 'Table_data'), recursive = TRUE)

cellchatLst <- list()
for (grp in names(object.list)) {
	cellchatLst[[grp]] <- readRDS(file.path(RDS_dir, paste(grp, sep='_', 'cellchat.rds')))
}
#上面4行cellchatLst相关命令可以省略！！！

cellchatMergdLst <- list()
for (pair in pair_list) {
	print(paste0(pair[1], ' ', pair[2], '\n'))
	grps2 <- paste0(pair[1], '_', pair[2])
	grp_comp=c(pair[1], pair[2])
	grpComp=paste0(str_replace(pair[1], ' ', ''), '_vs_', str_replace(pair[2], ' ', ''))
	
	print(grps2)
	print(grpComp)
	objectLst <- list()
	for (grp in pair) {
		objectLst[[grp]] <- cellchatLst[[grp]]
	}
	cellchatMergd <- mergeCellChat(objectLst, add.names = names(objectLst))
	cellchatMergd@meta$datasets = factor(cellchatMergd@meta$datasets, levels = pair) # set factor level
	
	cells1 <- rownames(cellchatMergd@net[[1]]$weight)
	cells2 <- rownames(cellchatMergd@net[[2]]$weight)
		# 判断是否完全一致
	is_aligned <- identical(sort(cells1), sort(cells2))
	if (!is_aligned) {
		cat("⚠️  细胞群未对齐，正在自动安全对齐（不影响原有分析）...\n")
		cellchatMergd <- liftCellChat(cellchatMergd)
	}
		
	cellchatMergd
	cellchatMergdLst[[grps2]] <- cellchatMergd
	
	grpCompDir <- file.path(cellchat_dir, 'Graph_data', grpComp)
	if (!dir.exists(grpCompDir)) {
		dir.create(grpCompDir, recursive = TRUE)
	} else {
		print("目录已存在。")
	}
	tblCompDir <- file.path(cellchat_dir, 'Table_data')
	if (!dir.exists(tblCompDir)) {
		dir.create(tblCompDir, recursive = TRUE)
	} else {
		print("目录已存在。")
	}
	
	iteractionsName=paste0('The_number_of_', grpComp, '_interactions')
	png(file.path(grpCompDir, paste0('01_', iteractionsName, '_network.png')), width  = 8, height = 6, units= "in", res=300, pointsize=12)
	netVisual_diffInteraction(cellchatMergd, weight.scale = T)
	dev.off()
    
	strengthName=paste0('The_number_of_', grpComp, '_interactions_strength')
	png(file.path(grpCompDir, paste0('01_', strengthName, '_network.png')), width  = 8, height = 6, units= "in", res=300, pointsize=12)
	 netVisual_diffInteraction(cellchatMergd, weight.scale = T, measure = "weight")
	dev.off()
	#The width of edges represent the relative number of interactions or interaction strength. Red (or blue) colored edges represent increased (or decreased) signaling in the second dataset compared to the first one.
	
	gg1 <- rankNet(cellchatMergd, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
	gg2 <- rankNet(cellchatMergd, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)
	ggsave(file.path(grpCompDir, paste0('02_', grpComp, '_all_clusters_infromationFlow.png')), gg1 + gg2, width = 14, height = 10, dpi = 300)
	
	grp1 <- names(objectLst)[1]
	grp2 <- names(objectLst)[2]
	sourceClusters <- unique(unique(cellchatLst[[grp1]]@idents), unique(cellchatLst[[grp2]]@idents))
	targetClusters <- unique(unique(cellchatLst[[grp1]]@idents), unique(cellchatLst[[grp2]]@idents))
	print(targetClusters)
	lrPair_plotting_20260603(sourceClusters, targetClusters, grp_comp, cellchatMergd, file.path(grpCompDir), file.path(tblCompDir))
    }

lrPair_plotting_20260603 <- function(sourceClusters, targetClusters, groupsPair, cellChatData, plotDir, tableDir) {
    for (e in sourceClusters) {
        target_list <- as.character(targetClusters[targetClusters != e])
        print(target_list)
        
        tryCatch({
	    gg1 <- netVisual_bubble(cellChatData, sources.use = e, targets.use = target_list,  comparison = c(1, 2), max.dataset = 1, title.name = paste0('Increased signaling in ', groupsPair[1]), angle.x = 45, remove.isolate = T, font.size.title = 14)
	    # 2. 【关键】获取真正的 Y 轴标签数（去重后）
	    n_labels <- length(unique(gg1$data$interaction_name_2))
	    # 3. 自动计算最佳高度（每个标签占 0.45~0.55 英寸）
	    height_auto <- n_labels * 0.3
	    height_auto <- max(height_auto, 6)  # 最小高度 6 英寸
	    png(file.path(plotDir, paste0('03_source', str_replace_all(e, '/', ''), '_Incr', str_replace_all(str_replace(grpComp, '_vs_', '_'), ' ', ''), '_communicationProbsChange_dotplot.png')), width= 8, height =height_auto, units  = "in", res= 300, pointsize = 12)
	    print(gg1) #直接用gg1时，在for循环中不能生成图片保存！！！
	    dev.off()
            file.path(plotDir, paste0('03_source', str_replace_all(e, '/', ''), '_Incr', str_replace_all(str_replace(grpComp, '_vs_', '_'), ' ', ''), '_communicationProbsChange_dotplot.png'))
	    }, error=function(e) message( "跳过气泡图1: ", e$message))
    	
        tryCatch({
	    gg2 <- netVisual_bubble(cellChatData, sources.use = e, targets.use = target_list,  comparison = c(1, 2), max.dataset = 2, title.name = paste0('Decreased signaling in ',  groupsPair[1]), angle.x = 45, remove.isolate = T, font.size.title = 14)
	    n_labels <- length(unique(gg2$data$interaction_name_2))
	    height_auto <- n_labels * 0.3
	    height_auto <- max(height_auto, 6)  # 最小高度 6 英寸
	    png(
		file.path(plotDir, paste0('03_source', str_replace_all(e, '/', ''), '_Decr', str_replace_all(str_replace(grpComp, '_vs_', '_'), ' ', ''), '_communicationProbsChange_dotplot.png')), width= 8, height =height_auto, units  = "in", res= 300, pointsize = 12)
	    print(gg2)
	    dev.off()
	    }, error=function(e) message("跳过气泡图1"))
        
        tryCatch({
	    gg1 <- netVisual_bubble(cellChatData, targets.use = e, sources.use = target_list,  comparison = c(1, 2), max.dataset = 1, title.name = paste0('Increased signaling in ', groupsPair[1]), angle.x = 45, remove.isolate = T, font.size.title = 14)
            n_labels <- length(unique(gg1$data$interaction_name_2))
	    height_auto <- n_labels * 0.3
	    height_auto <- max(height_auto, 6)  # 最小高度 6 英寸
            png(file.path(plotDir, paste0('03_target', str_replace_all(e, '/', ''), '_Incr', str_replace_all(str_replace(grpComp, '_vs_', '_'), ' ', ''), '_communicationProbsChange_dotplot.png')), width= 8, height =height_auto, units  = "in", res= 300, pointsize = 12)
        print(gg1)
        dev.off()
        }, error=function(e) message("跳过气泡图1"))
        
        tryCatch({
	    gg2 <- netVisual_bubble(cellChatData, targets.use = e, sources.use = target_list,  comparison = c(1, 2), max.dataset = 2, title.name = paste0('Decreased signaling in ',  groupsPair[1]), angle.x = 45, remove.isolate = T, font.size.title = 14)
        n_labels <- length(unique(gg2$data$interaction_name_2))
	    height_auto <- n_labels * 0.3
	    height_auto <- max(height_auto, 6)  # 最小高度 6 英寸
        png(file.path(plotDir, paste0('03_target', str_replace_all(e, '/', ''), '_Decr', str_replace_all(str_replace(grpComp, '_vs_', '_'), ' ', ''), '_communicationProbsChange_dotplot.png')), width= 8, height =height_auto, units  = "in", res= 300, pointsize = 12)
        print(gg2)
        dev.off()
        }, error=function(e) message("跳过气泡图1"))
    	
        # extract the ligand-receptor pairs with upregulated ligands and upregulated receptors in NL, i.e.,downregulated in LS
        pos.dataset = groupsPair[1]
        features.name = paste0(pos.dataset, ".merged")
        tryCatch({
                cellChatData <- identifyOverExpressedGenes(cellChatData, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.05)
                net <- netMappingDEG(cellChatData, features.name = features.name)
                write.csv(net, file.path(tableDir, paste0(str_replace_all(grpComp, ' ', ''), '_netMappingDEG.csv')), row.names = FALSE, fileEncoding = "UTF-8")
    	        
                net_up <- subsetCommunication(cellChatData, net = net, datasets = groupsPair[1], ligand.logFC = 0.05, receptor.logFC = NULL)
                net_down <- subsetCommunication(cellChatData, net = net, datasets = groupsPair[1], ligand.logFC = -0.05, receptor.logFC = NULL)
    		
                pairLR.use.up = net_up[, "interaction_name", drop = F]
                tryCatch({
	            gg1 <- netVisual_bubble(cellChatData, pairLR.use = pairLR.use.up, sources.use = e, targets.use = target_list, comparison = c(1, 2),  angle.x = 45, remove.isolate = T, font.size.title = 14, title.name = paste0('Up-regulated signaling in ', groupsPair[1]))
                n_labels <- length(unique(gg1$data$interaction_name_2))
                height_auto <- n_labels * 0.3
                height_auto <- max(height_auto, 6)  # 最小高度 6 英寸
                png(file.path(plotDir, paste0('04_source', str_replace_all(e, '/', ''), '_Incr', str_replace_all(str_replace(grpComp, '_vs_', '_'), ' ', ''), '_genesExpr_dotplot.png')), width= 8, height =height_auto, units  = "in", res= 300, pointsize = 12)
                print(gg1)
                dev.off()
                }, error=function(e) message("跳过气泡图1"))
    	        
                pairLR.use.down = net_down[, "interaction_name", drop = F]
                tryCatch({
	            gg2 <- netVisual_bubble(cellChatData, pairLR.use = pairLR.use.down, sources.use = e, targets.use = target_list, comparison = c(1, 2),  angle.x = 45, remove.isolate = T, font.size.title = 14, title.name = paste0('Down-regulated signaling in ',  groupsPair[1]))
                n_labels <- length(unique(gg2$data$interaction_name_2))
                height_auto <- n_labels * 0.3
                height_auto <- max(height_auto, 6)  # 最小高度 6 英寸
                png(file.path(plotDir, paste0('04_source', str_replace_all(e, '/', ''), '_Decr', str_replace_all(str_replace(grpComp, '_vs_', '_'), ' ', ''), '_genesExpr_dotplot.png')), width= 8, height =height_auto, units  = "in", res= 300, pointsize = 12)
                print(gg2)
                dev.off()
                 }, error=function(e) message("跳过气泡图1"))
    	        
                tryCatch({
	            gg1 <- netVisual_bubble(cellChatData, pairLR.use = pairLR.use.up, targets.use = e, sources.use = target_list, comparison = c(1, 2),  angle.x = 45, remove.isolate = T, font.size.title = 14, title.name = paste0('Up-regulated signaling in ', groupsPair[1]))
                n_labels <- length(unique(gg1$data$interaction_name_2))
                height_auto <- n_labels * 0.3
                height_auto <- max(height_auto, 6)  # 最小高度 6 英寸
                png(file.path(plotDir, paste0('04_targets', str_replace_all(e, '/', ''), '_Incr', str_replace_all(str_replace(grpComp, '_vs_', '_'), ' ', ''), '_genesExpr_dotplot.png')), width= 8, height =height_auto, units  = "in", res= 300, pointsize = 12)
                print(gg1)
                dev.off()
                }, error=function(e) message("跳过气泡图1"))
    	        
                #pairLR.use.down = net_down[, "interaction_name", drop = F]
                tryCatch({
	            gg2 <- netVisual_bubble(cellChatData, pairLR.use = pairLR.use.down, targets.use = e, sources.use = target_list, comparison = c(1, 2),  angle.x = 45, remove.isolate = T, font.size.title = 14, title.name = paste0('Down-regulated signaling in ',  groupsPair[1]))
                n_labels <- length(unique(gg2$data$interaction_name_2))
                height_auto <- n_labels * 0.3
                height_auto <- max(height_auto, 6)  # 最小高度 6 英寸
                png(file.path(plotDir, paste0('04_targets', str_replace_all(e, '/', ''), '_Decr', str_replace_all(str_replace(grpComp, '_vs_', '_'), ' ', ''), '_genesExpr_dotplot.png')), width= 8, height =height_auto, units  = "in", res= 300, pointsize = 12)
                print(gg2)
                dev.off()
                }, error=function(e) message("跳过气泡图1"))
	}, error = function(e) message("⚠️ identifyOverExpressedGenes 跳过：", e$message))
    }
}


cellchat_dir <- '/home/colddata/qinqiang/Project/Test/singleCellPipeline_test/pipelineTest1/05_Cell_communication_analysis_re1'
dir.create(file.path(cellchat_dir, 'Graph_data'), recursive = TRUE)
dir.create(file.path(cellchat_dir, 'Table_data'), recursive = TRUE)

all_pairs <- expand.grid(grp1 = names(cellchatLst), grp2 = names(cellchatLst))
all_pairs <- all_pairs[all_pairs$grp1 != all_pairs$grp2, ]
allPairs_list <- split(all_pairs, seq(nrow(all_pairs)))

grpcomp2cellchatLst <- list()
grpcomp2infoFlowLst <- list()
for (pair in allPairs_list) {
	grps2 <- paste0(pair$grp1, '_', pair$grp2)
	grpComp=paste0(str_replace(pair$grp1, ' ', ''), '_vs_', str_replace(pair$grp2, ' ', ''))
	
	print(grps2)
	print(grpComp)
	objectLst <- list()
	objectLst <- cellchatLst[unlist(pair)]
	cellchatMergd <- mergeCellChat(objectLst, add.names = names(objectLst))
	cellchatMergd@meta$datasets = factor(cellchatMergd@meta$datasets, levels = pair) # set factor level
	
	cells1 <- rownames(cellchatMergd@net[[1]]$weight)
	cells2 <- rownames(cellchatMergd@net[[2]]$weight)
		# 判断是否完全一致
	is_aligned <- identical(sort(cells1), sort(cells2))
	if (!is_aligned) {
		cat("⚠️  细胞群未对齐，正在自动安全对齐（不影响原有分析）...\n")
		cellchatMergd <- liftCellChat(cellchatMergd)
	}
	grpcomp2cellchatLst[[grpComp]] <- cellchatMergd
			
	gg1 <- rankNet(cellchatMergd, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
	grpcomp2infoFlowLst[[grpComp]] <- gg1
}

combGraphs_dir <- file.path(cellchat_dir, "Graph_data/combGraphs/")
if (!dir.exists(combGraphs_dir)) dir.create(combGraphs_dir, recursive = TRUE)

# 全部对比名称
all_comp <- names(grpcomp2infoFlowLst)
# 提取前缀（_vs_前面部分）
prefix_all <- sapply(strsplit(all_comp, "_vs_"), `[`, 1)

# 按前缀分类
group_split <- split(all_comp, prefix_all)

# 遍历5大前缀分组
for (prefix in names(group_split)) {
  comp_vec <- group_split[[prefix]]
  plot_list <- list()
  
  # 读取每组4张气泡图（假设你已提前导出对象，这里以绘图对象举例）
  for (comp in comp_vec) {
    # 这里替换成你实际绘图代码，得到单张气泡图p
    p <- grpcomp2infoFlowLst[[comp]]
    plot_list[[comp]] <- p
  }
  
  # 4图横向拼接
  combine_p <- wrap_plots(plot_list, nrow = 1, ncol = 4) #+
    #plot_annotation(title = paste0("Reference group: ", prefix))
  
  # 保存大图
  ggsave(
    filename = file.path(combGraphs_dir, paste0(prefix,"_combine.png")),
    plot = combine_p,
    width = 32,
    height = 8,
    dpi = 300
  )
}

tableDir <- file.path(cellchat_dir, 'Table_data')
tableGrpDir <- file.path(tableDir, 'Table_by_group')
dir.create(tableGrpDir, recursive = TRUE)

cellchatMergdAll <- mergeCellChat(cellchatLst, add.names = names(cellchatLst))
cellchatMergdAll@meta$datasets = factor(cellchatMergdAll@meta$datasets, levels = names(cellchatLst)) # set factor level

#下面为替代20260625脚本中生成‘write.csv(net_grp_sel, file.path(tableGrpDir, paste0(grp, '_netMapping.csv')), row.names = FALSE, fileEncoding = "UTF-8")’部分的代码
for (grp in names(cellchatLst)) {
	pos.dataset = grp
	features.name = paste0(pos.dataset, ".merged")
	cellchatMergdAll_NoLN <- identifyOverExpressedGenes(cellchatMergdAll, group.dataset = "datasets", pos.dataset = grp, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.05)
	net_grp <- netMappingDEG(cellchatMergdAll_NoLN, features.name = features.name)
	print(c(grp, dim(net_grp)))
	net_grp_sel <- subset(net_grp, net_grp$datasets==grp)
	
	write.csv(net_grp_sel, file.path(tableGrpDir, paste0(grp, '_netMapping.csv')), row.names = FALSE, fileEncoding = "UTF-8")
}

gg1 <- compareInteractions(cellchatMergdAll, show.legend = F, group = c(1,2,3,4,5), x.lab.rot=45)
gg1 + theme(axis.text.x=element_text(angle=45, hjust=1))
gg2 <- compareInteractions(cellchatMergdAll, show.legend = F, measure = "weight", group = c(1,2,3,4,5), x.lab.rot=45)
gg2 + theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave(file.path(cellchat_dir, 'Graph_data', '00_compareInteractions.png'), gg1 + gg2, dpi = 300)
#theme中的angle参数不起作用！！！

num.link <- sapply(cellchatLst, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})

weight.MinMax <- c(min(unlist(num.link)), max(unlist(num.link)))
gg <- list()
for (i in 1:length(cellchatLst)) {
  cellchatLst[[i]] <- netAnalysis_computeCentrality(cellchatLst[[i]], slot.name = "netP")
  gg[[i]]  <- netAnalysis_signalingRole_scatter(cellchatLst[[i]], title = names(cellchatLst)[i], weight.MinMax = weight.MinMax)  + theme(
  panel.border = element_rect(
    colour = "black",   # 边框颜色
    fill = NA,          # 不填充背景
    linewidth = 1       # 边框粗细
  ))
}

gg <- patchwork::wrap_plots(plots = gg) &
  theme(plot.margin = margin(10, 15, 10, 15))
ggsave(file.path(cellchat_dir, 'Graph_data', '00_clustersSignal_dotplot.png'), gg, width = 14, height = 10, dpi = 300)

for (ref_grp in names(group_split)) {
  comp_vec <- group_split[[ref_grp]]
  img_list <- list()
  print(comp_vec)
  
  for (i in seq_along(comp_vec)) {
    comp_name <- comp_vec[i]
    cc_obj <- grpcomp2cellchatLst[[comp_name]]
    
    # 保存单张 PNG（注意：这里不需要再赋值给 p）
    png_file <- file.path(combGraphs_dir, paste0('01_', comp_name, '_network.png'))
    png(png_file, width = 1000, height = 800, res = 150)
    par(mar = c(5,4,4,2) + 0.1, oma = c(0,0,2,0))
    netVisual_diffInteraction(
      cc_obj,
      weight.scale = T, measure = "weight",
      #title.name = comp_name
    )
    title(main = comp_name, line = 1, cex.main = 1.5, font.main = 1)
    dev.off()
    
    # 读取图片到 magick 对象
    img_list[[i]] <- image_read(png_file)
  }
  
  # 横向拼接所有图片（一行）
  combined_img <- image_append(do.call(c, img_list), stack = FALSE)
  
  # 添加标题（可选）：用 magick 的 image_annotate 或 ggplot 叠加
  # 简单做法：用 image_annotate 直接在图片上方加文字
  combined_img <- image_annotate(
    combined_img,
    text = '', #paste0("Reference group: ", ref_grp),
    gravity = "north",
    size = 40,
    color = "black",
    location = "+0+10"
  )
  
  # 保存最终拼接图（PNG）
  image_write(combined_img, 
              path = file.path(combGraphs_dir, paste0(ref_grp, "_diffNet.png")),
              format = "png", density = 300)

  if (file.exists(file.path(combGraphs_dir, paste0(ref_grp, "_diffNet.png")))) {
	filesDel <- list.files(
		path = combGraphs_dir,
		pattern = "01_",  # 例如 "\\.csv$" 匹配所有 .csv 文件
		full.names = TRUE                # 返回完整路径，否则只返回文件名
		)
	file.remove(filesDel)
	}
}


#四、细胞轨迹分析（工作站）
#服务器
r443

R
library(monocle3)
library(Seurat)
library(ggplot2)
library(dplyr)

setwd('../')
getwd()
[1] "/home/colddata/qinqiang/Project/Test/singleCellPipeline_test/pipelineTest1"

psudotimeDir <- file.path(getwd(),  '03_Celltrack_analysis')
#psudotimeDir
[1] "/home/colddata/qinqiang/Project/Test/singleCellPipeline_test/pipelineTest1/03_Celltrack_analysis"

dir.create(psudotimeDir, recursive = TRUE)

seurat_obj <- readRDS("adataRDS/adataFiltrd_celltypist.rds")
targetClusters <- c('0', '1', '2', '3', '4', '6')
targetObj <-seurat_obj %>%
	subset(
	subset = leiden_r025 %in% targetClusters
	)

#dim(targetObj@meta.data)
#[1] 123845     43

# 构建Monocle3所需的三个核心数据[6,8](@ref)
expression_matrix <- GetAssayData(targetObj, assay = 'RNA', layer = 'counts')  # 表达矩阵
cell_metadata <-targetObj@meta.data  # 细胞元数据
gene_annotation <- data.frame(gene_short_name = rownames(expression_matrix))  # 基因注释
rownames(gene_annotation) <- rownames(expression_matrix)

cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
# 数据预处理
cds <- preprocess_cds(cds, num_dim = 50)
cds <- align_cds(cds)
cds <- reduce_dimension(cds)

#细胞聚类
cds <- cluster_cells(cds)
# 学习轨迹图
cds <- learn_graph(cds)

rootCluster='1'
root_cells <- colnames(cds)[colData(cds)$leiden_r025 == rootCluster]
cds <- order_cells(cds, root_cell = root_cells)
Pseudotime <- pseudotime(cds, reduction_method='UMAP')
head(Pseudotime)
GSM8708828_cell_0 GSM8708828_cell_1 GSM8708828_cell_2 GSM8708828_cell_3
      0.006322076       0.031199069       0.003320957       0.116897061
GSM8708828_cell_4 GSM8708828_cell_5
      0.006652090       0.039881129

p <- plot_cells(cds,
           color_cells_by = "pseudotime",
	   group_cells_by = "leiden_r025", 
           label_cell_groups=TRUE,
           label_leaves=FALSE,
	   label_groups_by_cluster = TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5)

ggsave(file.path(psudotimeDir, paste0('umap_trajectoryPseudotime.png')), p, width = 14, height = 10, dpi = 300)

p1 <- plot_cells(cds,
           color_cells_by = "leiden_r025",
	   group_cells_by = "cluster", 
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=3)

ggsave(file.path(psudotimeDir, paste0('umap_trajectoryCluster.png')), p1, width = 14, height = 10, dpi = 300)

cluster2pseudotime <- data.frame(Cell_Barcode = targetObj$barcode, sample_ID=targetObj$sampleID, seurat_ClusterID = targetObj$leiden_r025,  monocle3_ClusterID = clusters(cds, reduction_method = 'UMAP'), Pseudotime=pseudotime(cds))
#dim(cluster2pseudotime)
#[1] 123845      5

write.csv(cluster2pseudotime, file.path(psudotimeDir,  paste0('cluster2pseudotime.csv')), row.names = FALSE, fileEncoding = "UTF-8")

grp2monocleCDS_list <- list()
for (grp in unique(targetObj$group)) {
	print(grp)
	targetCells <- rownames(targetObj@meta.data)[as.character(targetObj$group) == grp]
	targetCells_seuratObj <- subset(targetObj, cells = targetCells)
	
	# 构建Monocle3所需的三个核心数据[6,8](@ref)
	expression_matrix <- GetAssayData(targetCells_seuratObj, assay = 'RNA', layer = 'counts')  # 表达矩阵
	cell_metadata <- targetCells_seuratObj@meta.data  # 细胞元数据
	gene_annotation <- data.frame(gene_short_name = rownames(expression_matrix))  # 基因注释
	rownames(gene_annotation) <- rownames(expression_matrix)
	
	cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
	# 数据预处理
	cds <- preprocess_cds(cds, num_dim = 50)
	cds <- align_cds(cds)
	cds <- reduce_dimension(cds)

	#细胞聚类
	cds <- cluster_cells(cds)
	
	# 学习轨迹图
	cds <- learn_graph(cds)
	
	root_cells <- colnames(cds)[colData(cds)$leiden_r025 == 1]
	cds <- order_cells(cds, root_cell = root_cells)

	Pseudotime <- pseudotime(cds, reduction_method='UMAP')
	psudotimeDir <- file.path(getwd(),  '03_Celltrack_analysis')
	p <- plot_cells(cds,
           color_cells_by = "pseudotime",
	   group_cells_by = "leiden_r025", 
           label_cell_groups=TRUE,
           label_leaves=FALSE,
	   label_groups_by_cluster = TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5)
	
	ggsave(file.path(psudotimeDir, paste0(str_replace(grp, ' ', ''), '_umap_trajectoryPseudotime.png')), p, width = 14, height = 10, dpi = 300)
	
	p1 <- plot_cells(cds,
           color_cells_by = "leiden_r025",
	   group_cells_by = "cluster", 
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=3)
	
	ggsave(file.path(psudotimeDir, paste0(str_replace(grp, ' ', ''), '_umap_trajectoryCluster.png')), p1, width = 14, height = 10, dpi = 300)
	
	cluster2pseudotime <- data.frame(Cell_Barcode = targetCells_seuratObj$barcode, sample_ID=targetCells_seuratObj$sampleID, seurat_ClusterID = targetCells_seuratObj$leiden_r025,  monocle3_ClusterID = clusters(cds, reduction_method = 'UMAP'), Pseudotime=pseudotime(cds))
	write.csv(cluster2pseudotime, file.path(psudotimeDir,  paste0(str_replace(grp, ' ', ''), '_cluster2pseudotime.csv')), row.names = FALSE, fileEncoding = "UTF-8")

	grp2monocleCDS_list[[str_replace(grp, ' ', '')]] <- cds
}

grp2monocleCDS_lst <- list()
for (grp in names(grp2monocleCDS_list)) {
	grp2monocleCDS_lst[[str_replace(grp, ' ', '')]] <- grp2monocleCDS_list[grp]
}

save(grp2monocleCDS_lst,file = file.path(getwd(), 'adata_files', "grp2monocleCDS_lst.Rdata"))


#五、singleR细胞注释
getwd()
[1] "/home/colddata/qinqiang/Project/Test/singleCellPipeline_test/pipelineTest1/01_Primary_clustering"

otherAnnotDir <- file.path(getwd(),  '03_细胞注释/otherAnnot')
dir.create(otherAnnotDir, recursive = TRUE)

monaco <- MonacoImmuneData(legacy=TRUE)
seurat_obj <- NormalizeData(seurat_obj, 
                           normalization.method = "LogNormalize",
                           scale.factor = 1e4,  # CPM
                           verbose = FALSE)

sce <- as.SingleCellExperiment(seurat_obj, assay = "RNA")
pred <- SingleR(
  test = sce,
  ref = monaco,
  labels = monaco$label.fine,
  genes = "de",  # 使用差异表达基因，加快计算
  de.method = "wilcox",  # Wilcoxon检验，速度快
  de.n = 50,  # 每个细胞类型使用前50个差异基因
  aggr.ref = FALSE,  # 不聚合参考数据
  #BPPARAM = BiocParallel::SerialParam()  # 单核运行
  #clusters = pbmc$seurat_clusters  # 按聚类注释
)

seurat_obj$SingleR_Label <- pred$labels
seurat_obj$SingleR_Score <- pred$scores

saveRDS(seurat_obj, file = file.path(adataRDS_dir, 'adataFiltrd_singleR.rds'))

#monacoAnnon_df=data.frame(table(seurat_obj$SingleR_Label))
#colnames(monacoAnnon_df)=c('group', 'count')
r025_monacoAnnon_df=as.data.frame.matrix(table(seurat_obj$leiden_r025, seurat_obj$SingleR_Label))
r025_monacoAnnon_df=rownames_to_column(r025_monacoAnnon_df, var = "clusterID")
#write.csv(monacoAnnon_df, file.path(otherAnnotDir, 'monacoAnnon.csv'), row.names = FALSE)
write.csv(r025_monacoAnnon_df, file.path(otherAnnotDir, 'resolution025_clusters_monacoAnnon.csv'), row.names = FALSE)

sample_monacoAnnon_df=as.data.frame.matrix(table(seurat_obj$sampleID, seurat_obj$SingleR_Label))
sample_monacoAnnon_df=rownames_to_column(sample_monacoAnnon_df, var = "sampleID")
write.csv(sample_monacoAnnon_df, file.path(otherAnnotDir, 'sample_monacoAnnon.csv'), row.names = FALSE)

sample_cellCount_table <- table(seurat_obj$sampleID, seurat_obj$SingleR_Label)
cellProp_by_sample <- prop.table(t(sample_cellCount_table), margin = 2) * 100
cellProp_by_sample_df <- as.data.frame.matrix(round(cellProp_by_sample, 2))
#cellProp_by_sample_df$HD3
 [1]  0.71  0.01  0.05  0.01  9.74  0.00  0.00  1.35  0.05  0.02 51.11  2.35
[13]  0.05  0.00  0.00  0.61  0.00  0.00  0.00  0.00  8.59  2.88  0.23  5.47
[25]  6.46  5.13  3.82  1.38

cellProp_by_sample_df <- rownames_to_column(cellProp_by_sample_df, var = "cellType")
write.csv(cellProp_by_sample_df, file.path(otherAnnotDir, 'cellProp_by_sample.csv'), row.names = FALSE)

