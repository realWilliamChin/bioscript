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

plot_go_bubble <- function(dt, ontology, out.name, bubble_dir) {
  suffix <- switch(ontology, BP = "GOBP", CC = "GOCC", MF = "GOMF")
  dt.sub <- dt[dt$Ontology == ontology, ]
  dt.sub <- dt.sub[order(dt.sub$pvalue), ]
  dt.sub <- dt.sub[dt.sub$Count >= 3 & dt.sub$pvalue < 0.05, ]
  dt.sub <- head(dt.sub[order(dt.sub$pvalue), ], n = 15)
  dt.sub <- dt.sub[order(-dt.sub$pvalue), ]
  if (nrow(dt.sub) > 0) {
    p <- ggplot(dt.sub, aes(y = factor(GOID, levels = GOID), x = RichFactor, size = Count, colour = pvalue)) +
      geom_point() +
      scale_size_continuous(range = c(2, 8))+
      scale_y_discrete(limits = dt.sub$GOID) +
      scale_colour_gradient(low = "red", high = "green") +
      ggtitle(paste0(out.name, "_", suffix, sep = "")) +
      ylab("GO Term") +
      theme_base()
    ggsave(paste0(bubble_dir, "/", out.name, "_enrich_", suffix, "_Bubble.png", sep = ""), p, dpi = 320, width = 20, height = 10)
  } else {
    message(paste0(out.name, "没有显著富集的", suffix, "通路。"))
  }
  return(dt.sub)
}

for (i in infiles) {
  out.name <- gsub("_EnrichmentGO.xlsx", "", i)
  out.name
  dt <- read.xlsx(i, sheet = 1)
  dt$GOID <- paste(dt$ID, dt$Description, sep='_')
  dt.bp <- plot_go_bubble(dt, "BP", out.name, bubble_dir)
  dt.cc <- plot_go_bubble(dt, "CC", out.name, bubble_dir)
  dt.mf <- plot_go_bubble(dt, "MF", out.name, bubble_dir)

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
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),  # 设置标题字体大小
      axis.title = element_text(size = 20, face = "bold"),  # 设置 x 轴标签字体大小
      axis.text.y = element_text(size = 16, face = "bold"),  # 设置 y 轴标签字体大小
      legend.title = element_text(size = 20, face = "bold"),  # Ontology 图例标题字体大小
      legend.text = element_text(size = 18, face = "bold")  # GO_BP/CC/MF 的图例字体大小
    )
  # height 根据 datafile 的行数调整
  go_p_height <- 2 + 0.4 * nrow(dt.total)
  ggsave(paste0(barplot_dir, "/", out.name, "_enrich_GObarplot.png", sep = ""), bar.p, dpi = 320, width = 24, height = go_p_height)
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