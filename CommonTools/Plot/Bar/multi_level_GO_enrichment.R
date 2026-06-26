suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(argparse)
})

# 创建参数解析器
parser <- ArgumentParser(description = "多层级GO富集分析可视化")

# 添加参数
parser$add_argument("--input_file", "-i", 
                   help = "输入文件路径", 
                   required = TRUE)
parser$add_argument("--plot_title", "-t", 
                   help = "图表标题", 
                   required = TRUE)
parser$add_argument("--out_name", "-o", 
                   help = "输出文件名", 
                   required = TRUE)

# 解析参数
args <- parser$parse_args()

# 获取参数
input_file <- args$input_file
plot_title <- args$plot_title
out_name <- args$out_name

Fon <- 'sans'
input_df <- read.table(input_file,header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
# min_bar_width <- (max(input_df$Count)) * 0.005
# input_df$Count[input_df$Count == 0] <- min_bar_width

# 1. 数据准备 ------------------------------------------------------

# 计算 Description 最长字符数
max_desc_len <- max(nchar(as.character(input_df$Description)), na.rm = TRUE)

# 基于最长 Description 文本宽度动态计算 y 轴扩展值
# 使用 grid::textGrob 与 grid::convertWidth 估算以英寸为单位的文本宽度
max_desc_str <- input_df$Description[which.max(nchar(as.character(input_df$Description)))]
if (length(max_desc_str) == 0 || is.na(max_desc_str)) {
  max_desc_str <- ""
}
if (!requireNamespace("grid", quietly = TRUE)) {
  stop("Package 'grid' is required for width calculation. Please install it.")
}
# ggplot2 中 size=4 等价于 fontsize=4 * .pt（points），.pt ≈ 72.27/25.4
fontsize_pts <- 4 * ggplot2::.pt
desc_width_in <- grid::convertWidth(
  grid::grobWidth(grid::textGrob(label = max_desc_str, gp = grid::gpar(fontsize = fontsize_pts))),
  unitTo = "in",
  valueOnly = TRUE
)
# 估算面板宽度（英寸）：按最终组合布局的宽度占比计算（1.3 / (1.3 + 1.3 + 0.04 + 0.04)）
figure_width_in <- 26 # 与 ggsave 一致
panel_prop <- 1.3 / (1.3 + 1.3 + 0.04 + 0.04)
panel_width_in <- figure_width_in * panel_prop
# y 轴扩展量换算为数据单位（与 xmax 成比例）
y_extension_factor <- ifelse(panel_width_in > 0, desc_width_in / panel_width_in, 0)

## 按 Regulation 分组，原表已按 Comparison 分组，直接使用
plot_data_up <- input_df %>%
  filter(Regulation == 'Up') %>%
  mutate(Count = -Count,                     # 左侧镜像显示
         Group = "plot_data_up",
         Comparison_Description = paste0(Comparison, ": ", Description)) %>%
  distinct(Comparison_Description, .keep_all = TRUE) %>%  # 去重，确保每个 Comparison_Description 只有一行
  mutate(Comparison_Description = factor(Comparison_Description, levels = unique(Comparison_Description)))

plot_data_down <- input_df %>%
  filter(Regulation == 'Down') %>%
  mutate(Group = "plot_data_down",
         Comparison_Description = paste0(Comparison, ": ", Description)) %>%
  distinct(Comparison_Description, .keep_all = TRUE) %>%  # 去重，确保每个 Comparison_Description 只有一行
  mutate(Comparison_Description = factor(Comparison_Description, levels = unique(Comparison_Description)))


## 各侧分别识别 Comparison 并生成独立配色
comparisons_up <- sort(unique(plot_data_up$Comparison))
fill_colors_up <- scales::hue_pal()(length(comparisons_up))
names(fill_colors_up) <- comparisons_up

comparisons_down <- sort(unique(plot_data_down$Comparison))
fill_colors_down <- scales::hue_pal()(length(comparisons_down))
names(fill_colors_down) <- comparisons_down

# 计算各侧 Comparison 的分组边界索引，用于绘制分隔线
levels_up <- levels(plot_data_up$Comparison_Description)
comp_seq_up <- sub(":.*$", "", levels_up)
rle_up <- rle(comp_seq_up)
boundary_up <- head(cumsum(rle_up$lengths), -1)

levels_down <- levels(plot_data_down$Comparison_Description)
comp_seq_down <- sub(":.*$", "", levels_down)
rle_down <- rle(comp_seq_down)
boundary_down <- head(cumsum(rle_down$lengths), -1)

# x 轴最大值
up_x_max <- max(abs(plot_data_up$Count), na.rm = TRUE)
down_x_max <- max(abs(plot_data_down$Count), na.rm = TRUE)
xmax <- max(up_x_max, down_x_max)
message(sprintf("DEBUG: 最大 Count 值: %.3f", xmax))
ext_val <- xmax * y_extension_factor

## 创建左图（Up）：按Comparison分组的柱状图
plot_data_up$star <- cut(plot_data_up$p.adjust, breaks = c(0, 0.001, 0.01, 0.05, Inf), labels = c("***", "**", "*", ""))
# 计算相对偏移量，基于最大 Count 值的百分比
left_plot <- ggplot(plot_data_up, aes(x = Comparison_Description, y = Count, fill = Comparison)) +
  geom_bar(stat = "identity", width = 0.9, na.rm=TRUE) +
  geom_text(aes(label = Description), hjust = 1.05, vjust = 0.5, size = 4, color = "black") +  # Description
  geom_label(aes(label = star, y = Count), nudge_y = xmax * 0.02, color="black", alpha=0, label.size = 0, fill='white') +  # ***
  geom_segment(aes(y = Inf, yend = Inf, x = -Inf, xend = Inf), color = "black", linewidth = 0.5) +
  scale_y_continuous(limits = c(-xmax - ext_val, 0), labels = abs, expand = c(0, 0)) +
  scale_fill_manual(values = fill_colors_up, drop = FALSE) +
  coord_flip(clip = "off") +
  labs(title = "Up regulated", x = "", y = "") +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 14, vjust = -0.5, face='bold', family=Fon),
        plot.title = element_text(size = 18, hjust = 0.5, face='bold', family=Fon),
        plot.margin = margin(2, 2, 2, 2),
        legend.position = "none",
        panel.grid.major = element_blank(), # Remove major gridlines
        panel.grid.minor = element_blank()
        )
left_plot

## 创建右图（Down）：按Comparison分组的柱状图
plot_data_down$star <- cut(plot_data_down$p.adjust, breaks = c(0, 0.001, 0.01, 0.05, Inf), labels = c("***", "**", "*", ""))
# 计算相对偏移量，基于最大 Count 值的百分比
right_plot <- ggplot(plot_data_down, aes(x = Comparison_Description, y = Count, fill = Comparison)) +
  geom_bar(stat = "identity", width = 0.9) +
  geom_text(aes(label = Description), hjust = -0.05, vjust = 0.5, size = 4.2, color = "black") +
  geom_label(aes(label = star, y = Count), nudge_y = -xmax * 0.02, color="black", alpha=0, label.size = 0, fill='white') +
  scale_y_continuous(limits = c(0, xmax + ext_val), expand = c(0, 0)) +
  scale_fill_manual(values = fill_colors_down, drop = FALSE) +
  coord_flip(clip = "off") +
  labs(title = "Down regulated", x = "", y = "") +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_line(linewidth = 0.5, color='black'),
        axis.text.x = element_text(size = 14, vjust = -0.5, face='bold', family=Fon),
        plot.title = element_text(size = 18, hjust = 0.5, face='bold', family=Fon),
        plot.margin = margin(2, 2, 2, 2),
        legend.position = "none",
        panel.grid.major = element_blank(), # Remove major gridlines
        panel.grid.minor = element_blank()
        )

## 为左右两侧添加颜色注释列（annotation_row 风格，按 Comparison 聚合整块）
starts_up <- cumsum(c(1, head(rle_up$lengths, -1)))
ends_up <- cumsum(rle_up$lengths)
blocks_up <- data.frame(
  Comparison = rle_up$values,
  ymin = starts_up - 0.5,
  ymax = ends_up + 0.5,
  y_mid = (starts_up + ends_up) / 2,
  stringsAsFactors = FALSE
)
anno_up <- ggplot(blocks_up) +
  geom_rect(aes(xmin = 0.5, xmax = 1.5, ymin = ymin, ymax = ymax, fill = Comparison)) +
  geom_text(aes(x = 1, y = y_mid, label = Comparison), size = 5, color = "black", angle = 90) +
  scale_y_continuous(limits = c(0.5, length(levels_up) + 0.5), expand = c(0, 0)) +
  scale_fill_manual(values = fill_colors_up, drop = FALSE) +
  scale_x_continuous(limits = c(0.5, 1.5), expand = c(0, 0)) +
  theme_void() +
  theme(
    plot.margin = margin(2, 0, 2, 0),
    legend.position = "none",
    plot.title = element_text(face='bold', family=Fon)
    )

levels_down <- levels(plot_data_down$Comparison_Description)
comp_seq_down <- sub(":.*$", "", levels_down)
rle_down <- rle(comp_seq_down)
starts_down <- cumsum(c(1, head(rle_down$lengths, -1)))
ends_down <- cumsum(rle_down$lengths)
blocks_down <- data.frame(
  Comparison = rle_down$values,
  ymin = starts_down - 0.5,
  ymax = ends_down + 0.5,
  y_mid = (starts_down + ends_down) / 2,
  stringsAsFactors = FALSE
)
anno_down <- ggplot(blocks_down) +
  geom_rect(aes(xmin = 0.5, xmax = 1.5, ymin = ymin, ymax = ymax, fill = Comparison)) +
  geom_text(aes(x = 1, y = y_mid, label = Comparison), size = 5, color = "black", angle = 90) +
  scale_y_continuous(limits = c(0.5, length(levels_down) + 0.5), expand = c(0, 0)) +
  scale_fill_manual(values = fill_colors_down, drop = FALSE) +
  scale_x_continuous(limits = c(0.5, 1.5), expand = c(0, 0)) +
  theme_void() +
  theme(
    plot.margin = margin(2, 0, 2, 0),
    legend.position = "none",
    plot.title = element_text(face='bold', family=Fon)
    )

## 合并左右图：各自保留独立图例
final_plot <- anno_up + left_plot + right_plot + anno_down + 
  plot_layout(nrow = 1, widths = c(0.04, 1.3, 1.3, 0.04)) +
  plot_annotation(title = plot_title,theme = theme(plot.title = element_text(hjust = 0.5, size = 20, family=Fon, face='bold')))


# 保存图表（高度根据两侧条目数自适应）
num_terms_left <- nlevels(plot_data_up$Comparison_Description)
num_terms_right <- nlevels(plot_data_down$Comparison_Description)
plot_height <- max(max(num_terms_left, num_terms_right) * 0.22 + 3, 10)
ggsave(out_name, width = 26, height = plot_height, dpi = 300, limitsize = FALSE)

