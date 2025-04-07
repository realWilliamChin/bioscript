library(ggplot2)
library(ggthemes)
library(gdata)
library(svglite)
library(openxlsx)

infiles <- grep("_EnrichmentGO.xlsx", dir(), value = T)
infiles
crd <- getwd()

barplot_dir <- file.path(crd, "Pathway_enrichment_bar_plot")
bubble_dir <- file.path(crd, "Pathway_enrichment_bubble_plot")
raw_data <- file.path(crd, "Pathway_enrichment_raw_data")

if (!dir.exists(barplot_dir)) {
  dir.create(barplot_dir, recursive = TRUE)
}
if (!dir.exists(bubble_dir)) {
  dir.create(bubble_dir, recursive = TRUE)
}
if (!dir.exists(raw_data)) {
  dir.create(raw_data, recursive = TRUE)
}

for (i in infiles) {
  out.name <- gsub("_EnrichmentGO.xlsx", "", i)
  out.name
  dt <- read.xlsx(i, sheet = 1)
  dt$GOID <- paste(dt$ID, dt$Description, sep='_')

  dt.bp <- dt[dt$Ontology == "BP", ]
  dt.bp <- dt.bp[order(dt.bp$pvalue), ]
  dt.bp <- dt.bp[dt.bp$Count >= 3 & dt.bp$pvalue < 0.05, ]
  dt.bp <- head(dt.bp[order(dt.bp$pvalue), ], n = 15)
  dt.bp <- dt.bp[order(-dt.bp$pvalue), ]
  dt.bp[, c('ID', 'Count', 'pvalue')]
  if (nrow(dt.bp) > 0) {
    p <- ggplot(dt.bp, aes(y = factor(GOID, levels = GOID), x = RichFactor, size = Count, colour = pvalue)) +
      geom_point() +
      scale_size_continuous(range = c(2, 8))+
      scale_y_discrete(limits = dt.bp$GOID) +
      scale_colour_gradient(low = "red", high = "green") +
      ggtitle(paste0(out.name, "_GOBP", sep = "")) +
      ylab(colnames(data)[1]) +
      theme_base()
    ggsave(paste0(bubble_dir, "/", out.name, "_enrich_GOBP_Bubble.png", sep = ""), p, dpi = 320, width = 20, height = 10)
  } else {
    message(paste0(out.name, "没有显著富集的GOBP通路。"))
  }

  dt.cc <- dt[dt$Ontology == "CC", ]
  dt.cc <- dt.cc[order(dt.cc$pvalue), ]
  dt.cc <- dt.cc[dt.cc$Count >= 3 & dt.cc$pvalue < 0.05, ]
  dt.cc <- head(dt.cc[order(dt.cc$pvalue), ], n = 15)
  dt.cc <- dt.cc[order(-dt.cc$pvalue), ]
  dt.cc[, c('ID', 'Count', 'pvalue')]
  if (nrow(dt.cc) > 0) {
    p <- ggplot(dt.cc, aes(y = factor(GOID, levels = GOID), x = RichFactor, size = Count, colour = pvalue)) +
      geom_point() +
      scale_size_continuous(range = c(2, 8))+
      scale_y_discrete(limits = dt.cc$GOID) +
      scale_colour_gradient(low = "red", high = "green") +
      ggtitle(paste0(out.name, "_GOCC", sep = "")) +
      ylab(colnames(data)[1]) +
      theme_base()
    ggsave(paste0(bubble_dir, "/", out.name, "_enrich_GOCC_Bubble.png", sep = ""), p, dpi = 320, width = 20, height = 10)
  } else {
    message(paste0(out.name, "没有显著富集的GOCC通路。"))
  }

  dt.mf <- dt[dt$Ontology == "MF", ]
  dt.mf <- dt.mf[order(dt.mf$pvalue), ]
  dt.mf <- dt.mf[dt.mf$Count >= 3 & dt.mf$pvalue < 0.05, ]
  dt.mf <- head(dt.mf[order(dt.mf$pvalue), ], n = 15)
  dt.mf <- dt.mf[order(-dt.mf$pvalue), ]
  dt.mf[, c('ID', 'Count', 'pvalue')]
  if (nrow(dt.mf) > 0) {
    p <- ggplot(dt.mf, aes(y = factor(GOID, levels = GOID), x = RichFactor, size = Count, colour = pvalue)) +
      geom_point() +
      scale_size_continuous(range = c(2, 8))+
      scale_y_discrete(limits = dt.mf$GOID) +
      scale_colour_gradient(low = "red", high = "green") +
      ggtitle(paste0(out.name, "_GOMF", sep = "")) +
      ylab(colnames(data)[1]) +
      theme_base()
    ggsave(paste0(bubble_dir, "/", out.name, "_enrich_GOMF_Bubble.png", sep = ""), p, dpi = 320, width = 20, height = 10)
  } else {
    message(paste0(out.name, "没有显著富集的GOMF通路"))
  }

  kegg_file <- gsub("_EnrichmentGO.xlsx", "_EnrichmentKEGG.xlsx", i)
  dt.kegg <- read.xlsx(kegg_file, check.names = F, sheet = 1)
  dt.kegg$KEGGID <- paste(dt.kegg$ID, dt.kegg$Description, sep='_')
  dt.kegg <- dt.kegg[order(dt.kegg$pvalue), ]
  dt.kegg <- dt.kegg[dt.kegg$Count >= 3 & dt.kegg$pvalue < 0.05, ]
  dt.kegg <- head(dt.kegg[order(dt.kegg$pvalue), ], n = 15)
  dt.kegg <- dt.kegg[order(-dt.kegg$pvalue), ]
  dt.kegg[, c('ID', 'Count', 'pvalue')]
  if (nrow(dt.kegg) > 0) {
    p <- ggplot(dt.kegg, aes(y = factor(KEGGID, levels = KEGGID), x = RichFactor, size = Count, colour = pvalue)) +
      geom_point() + 
      scale_size_continuous(range = c(2, 8))+
      scale_y_discrete(limits = dt.kegg$KEGGID) +
      scale_colour_gradient(low = "red", high = "green") +
      ggtitle(paste0(out.name, "_KEGG", sep = "")) +
      ylab(colnames(data)[1]) +
      theme_base()
    ggsave(paste0(bubble_dir, "/", out.name, "_enrich_KEGG_Bubble.png", sep = ""), p, dpi = 320, width = 20, height = 10)
  } else {
    message(paste0(out.name, "没有显著富集的KEGG通路。"))
  }

  dt.bp <- head(dt.bp[order(dt.bp$pvalue), ], n = 10)
  dt.cc <- head(dt.cc[order(dt.cc$pvalue), ], n = 10)
  dt.mf <- head(dt.mf[order(dt.mf$pvalue), ], n = 10)
  dt.bp <- dt.bp[order(dt.bp$ID), ]
  dt.cc <- dt.cc[order(dt.cc$ID), ]
  dt.mf <- dt.mf[order(dt.mf$ID), ]
  
  dt.total <- rbind(dt.bp, dt.cc, dt.mf)
  
  dt.total$star <- cut(dt.total$pvalue, breaks = c(0, 0.001, 0.01, 0.05, Inf), 
    labels = c("***        ", "**     ", "*   ", ""))
  #dt.total$GOID <- factor(dt.total$GOID, levels = dt.total$GOID)
  bar.p <- ggplot(dt.total, aes(x= factor(GOID, levels = GOID), y=RichFactor, fill=Ontology)) + 
    geom_bar(stat="identity", position = "identity")+
    geom_text(aes(label = star), position = position_dodge(width = 1), size = 7) +  # 添加星号
    # geom_bar(position="identity",stat="identity",aes(fill=SubOntology))+
    # facet_grid(SubOntology~.,scale="free")+
    theme_bw()+
    xlab("GO Term") +  # 设置 x 轴标签
    ylab("RichFactor") +    # 设置 y 轴标签
    coord_flip() +
    ggtitle(paste0(out.name, "_Enrich_Graph", sep = "")) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 14, face = "bold"),
      legend.title = element_text(face = "bold"))
  # height 根据 datafile 的行数调整
  go_p_height <- 2 + 0.4 * nrow(dt.total)
  ggsave(paste0(barplot_dir, "/", out.name, "_enrich_GObarplot.png", sep = ""), bar.p, dpi = 320, width = 20, height = go_p_height)
}

xlsx_files <- list.files(path = crd, pattern = "\\.xlsx$", full.names = TRUE)
destination_folder <- file.path(crd, "Pathway_enrichment_raw_data")
for (file in xlsx_files) {
  destination_file <- file.path(destination_folder, basename(file))
  success <- file.rename(file, destination_file)
  
  if (success) {
    a <- 1
  } else {
    message(paste("文件", basename(file), "移动失败。"))
  }
}