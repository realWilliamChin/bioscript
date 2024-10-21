# 读取数据文件
data <- read.table("21mer_out.histo", header = FALSE, sep = " ")

# 获取第一列和第二列数据
x_values <- data$V1  # 第一列作为 x 轴的位置
y_values <- data$V2  # 第二列作为 y 轴的值（直方图高度）


# 保存为 JPEG 文件
jpeg("test.jpeg", quality = 100)  # 打开 JPEG 输出，quality 参数设置图像质量（0-100，100 为最高质量）

# 绘制直方图
barplot(height = y_values, names.arg = x_values,
        xlab = "X Values", ylab = "Y Values", main = "Custom Histogram",
        col = "skyblue", border = "black",
        ylim = c(0, max(y_values) * 1.2))  # 设置 y 轴的范围略大于最大值


dev.off()  # 关闭绘图设备，保存图像为 test.jpeg