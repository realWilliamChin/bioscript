library(ggplot2)

# 输入 multideseq.py --degvalue 的参数值
degvalue <- 2
bs_pos <- log2(degvalue)
bs_neg <- -bs_pos

prefix <- 'L_vs_Z_volcano'
input_file <- paste0(prefix,'.txt')
outfile_jpeg <- paste0(prefix,".jpeg")
outfile_tiff <- paste0(prefix,'.tiff')
outfile_pdf <- paste0(prefix, '.pdf')

volcano = read.table(input_file, sep = "\t", header = T, check.names = F)

# 根据条件筛选添加
# volcano$label = ifelse(volcano$pvalue < 0.05 & abs(log2(volcano$FC)) >= 1,volcano$GeneID,"")

p.volcano <- ggplot(
  data = volcano, 
  aes(
    x = log2FoldChange,
    y = -log10(padj), 
    colour = regulation
    )
  ) +
  geom_point() +
  scale_color_manual(values = c("green", "grey", "red")) +
  geom_vline(xintercept = c(bs_neg, bs_pos), lty = 4, col = "black", lwd = 0.8) +
  geom_hline(yintercept = -log10(0.05), lty = 4, col = "black", lwd = 0.8) +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank()) +
  geom_text_repel(
    data = volcano,
    aes(
      x = log2FoldChange,
      y = -log10(padj),
      label = Genesymbol),
    size = 4,
    color = 'black',
    box.padding = unit(0.5, "lines"), # 增加标签之间的间距
    point.padding = unit(0.5, "lines"), # 调整点与标签之间的距离
    segment.color = "black",
    show.legend = FALSE,
    max.overlaps = Inf,
    #force = 5, 
    #force_pull = 1
  )
p.volcano

ggsave(outfile_jpeg, p.volcano, dpi = 300, width = 12, height = 12)
ggsave(outfile_tiff, p.volcano, dpi = 300, width = 12, height = 12)
ggsave(outfile_pdf, p.volcano, dpi = 300, width = 12, height = 12)


