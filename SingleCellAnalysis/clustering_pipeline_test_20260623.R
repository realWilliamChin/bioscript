20260611
#r443
cd /home/colddata/qinqiang/Project/Test/clustering_test
R

library(ggplot2)
library(ggdendro)
library(readxl)
library(showtext)
library(factoextra)
library(cluster)
library(ComplexHeatmap)
library(pheatmap)
library(dplyr)
library(FactoMineR)
library(ggplotify)
library(grid)
library(NbClust)
library(writexl)
library(openxlsx)
library(tibble)

set.seed(123)

testDir <- file.path(getwd(), 'clustering_pipeline_test')
df=read_excel('test-data.xlsx')
df=as.data.frame(df)
rownames(df) <- df$实验用编号
df <- df[, -which(names(df) == "实验用编号")]

#挑选用于聚类的指标变量
clustering_df <- df %>% select('白蛋白*g/L', 'TC *mmol/L', 'HDL*mmol/L', 'LDL*mmol/L', 'TG*mmol/L', '尿酸*umol/L', 'Cr*mmol/L', '糖化血红蛋白', '右房大小', '左房前后径', '左室前后径', '射血分数', 'E', 'A', 'e')

for (e in colnames(clustering_df)) {
	clustering_df[[e]]=as.numeric(clustering_df[[e]])
}

# 数据+标准化（PCA/聚类用同一套数据）
clustering_df_scale <- scale(clustering_df) #（行=样本，列=代谢物）

# 层次聚类(Ward)版本，做对照
nb_ward <- NbClust(clustering_df_scale, min.nc=2, max.nc=10, method="ward.D2", index="all")

index_names <- colnames(nb_ward$Best.nc)
opt_k_vals <- as.integer(nb_ward$Best.nc[1, ])
# 筛选核心指标结果
index_k_table <- data.frame(
  Index = index_names,
  Opt_K = opt_k_vals,
  row.names = NULL
)

# 定义核心指标清单（K-Means+PCA 场景优先）
core_index <- c("Hartigan", "CH", "Silhouette", "DB", "TraceW", "Ball")
core_k_table <- subset(index_k_table, Index %in% core_index)

# 统计每个K获得的票数（投票分布）
vote_count <- table(index_k_table$Opt_K)

bestK=3
hc <- hclust(dist(clustering_df_scale), method = "ward.D2")
ward_group <- cutree(hc, k = bestK)
ward_group <- as.factor(ward_group)
clustering_df.pca <- PCA(clustering_df_scale, graph = F)

# PCA绘图，habillage填入聚类分组，带椭圆
grp_pca <- fviz_pca_ind(
  clustering_df.pca,
  mean.point = F,
  label = "none",
  habillage = ward_group,
  palette = "lancet",
  addEllipses = TRUE,
  ellipse.type = "t",  # 凸包，几乎无样本量限制
  ellipse.level = 0.95,      # 降低置信水平
  ggtheme = theme_bw(base_size = 14)
)

ggsave(file.path(testDir, "01_Cluster_PCA_ward.png"), grp_pca, width=9, height=7, dpi=300)

ward_pca <- fviz_pca_ind(
  clustering_df.pca,
  mean.point = F,
  label = "none",
  habillage = ward_group,
  palette = c("#00AFBB", "#E7B800", "#FC4E07"),
  addEllipses = TRUE,
  ellipse.level = 0.95
)

# 1. 先确认 ward_group 分组数和颜色数一致
# 你当前是 3 个 cluster，所以 palette 用 3 个颜色
pal_clust <- c("#00AFBB", "#E7B800", "#FC4E07")

# 2. 绘制 PCA 图
grp_pca <- fviz_pca_ind(
  clustering_df.pca,
  mean.point = F,
  label = "none",
  habillage = df$年龄阶段,
  palette = "lancet",
  addEllipses = TRUE,
  ellipse.type = "t",  # 凸包，几乎无样本量限制
  ellipse.level = 0.95,      # 降低置信水平
  ggtheme = theme_bw(base_size = 14)
)

ggsave(file.path(testDir, "02_Cluster_PCA_ward_by_age.png"), grp_pca, width=9, height=7, dpi=300)

grp_pca <- fviz_pca_ind(
  clustering_df.pca,
  mean.point = F,
  label = "none",
  habillage = df$分组,
  palette = "lancet",
  addEllipses = TRUE,
  ellipse.type = "t",  # 凸包，几乎无样本量限制
  ellipse.level = 0.95,      # 降低置信水平
  ggtheme = theme_bw(base_size = 14)
)

ggsave(file.path(testDir, "02_Cluster_PCA_ward_by_disease.png"), grp_pca, width=9, height=7, dpi=300)

df <- df %>%
	mutate(
		Sex = case_when(
		性别 == '女' ~ "Female",
		性别 == '男' ~ "Male",
		)
	)

grp_pca <- fviz_pca_ind(
  clustering_df.pca,
  mean.point = F,
  label = "none",
  habillage = df$Sex,
  palette = "lancet",
  addEllipses = TRUE,
  ellipse.type = "t",  # 凸包，几乎无样本量限制
  ellipse.level = 0.95,      # 降低置信水平
  ggtheme = theme_bw(base_size = 14)
)

ggsave(file.path(testDir, "02_Cluster_PCA_ward_by_sex.png"), grp_pca, width=9, height=7, dpi=300)

colnames(clustering_df_scale) <- c('Albumin', 'TC', 'HDL', 'LDL', 'Triglycerides', 'Uric Acid', 'Creatinine', 'Hemoglobin A1c', 'RAS', 'LADD', 'LVAD', 'EjectionFraction', 'E', 'A', 'e')

dist_mat <- dist(clustering_df_scale, method = "euclidean")
hc <- hclust(dist_mat, method = "ward.D2")
cluster_result <- cutree(hc, k = bestK)

#4. 构建样本注释（热图标聚类颜色） =====================
annotation_row <- data.frame(
  Cluster = factor(paste0("C", cluster_result))
)
rownames(annotation_row) <- rownames(clustering_df_scale)

#5. 画热图 + 聚类树（论文标准图） =====================
p1 <- pheatmap(t(clustering_df_scale),  # 转置：列=样本，行=代谢物
         scale = "none", # 已标准化
         clustering_method = "ward.D2",
         clustering_distance_cols = "euclidean",
         clustering_distance_rows = "euclidean",
         annotation_col = annotation_row,  # 聚类颜色条
         treeheight_col = 40,    # 样本聚类树高度
         treeheight_row = 30,    # 代谢物聚类树高度
         show_rownames = T,
         show_colnames = T,
         silent = TRUE,       # 关键：不直接弹出图
         fontsize = 12,
         color = colorRampPalette(c("#4575B4","white","#D73027"))(100),
         main = paste0("heatmap | best_k = ", bestK))

# 3. 保存图片（PNG 高清 300dpi）
png(file.path(testDir, "03_Heatmap_by_cluster.png"), width = 3200, height = 2000, res = 300)
# 先画背景（最底层）
grid.rect(gp = gpar(fill = "gray95", col = NA))
# 再把 pheatmap 画在背景上（字体不缩放！）
print(p1, newpage = FALSE)
dev.off()  # 保存完成

vote_count_df <- as.data.frame(vote_count)
colnames(vote_count_df) <- c('Cluster_k', 'Vote_Num')
write_xlsx(vote_count_df, file.path(testDir, "04_Clustering_full_index.xlsx"))

colnames(core_k_table) <- c('Index_Name', 'Recommend_K')
write_xlsx(core_k_table, file.path(testDir, "04_Clustering_core_index.xlsx"))

df$Cluster <- ward_group
df_new <- df[-21]
df_new <- df_new[order(df_new$分组, df_new$Cluster), ]
df_new <- rownames_to_column(df_new, var = "实验用编号")

df_srt <- df_new[, c(1, 17, 22, 2:16, 18:21)]
write_xlsx(df_srt, file.path(testDir, "04_SampleID_to_clusterID.xlsx"))
