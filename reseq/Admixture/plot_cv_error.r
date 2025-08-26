library(ggplot2)
file <- 'Setaria_italica_snp_01_05.cv.error'
out_file <- 'Setaria_italica_snp_01_05_plot.jpeg'
df <- read.table(file, sep = " ", head=T)

df_plot <- ggplot(df, aes(x=K,y=CV_Error))+geom_point()+geom_line()
ggsave(out_file, df_plot, dpi=300)