# 计算在肿瘤条件下，显著性基因之间的相关性，并用环形和弦图和相关性矩阵热图可视化。识别基因间的协同表达关系

rm(list = ls())  # 清除工作空间

if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}
# corrplot（用于绘制相关性矩阵图）、circlize（用于绘制环形图）
pacman::p_load(corrplot, circlize, here)


# 读取数据文件
data = read.table(here("3_outputs/diffGeneExp.txt"), header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)  # 读取带有标题的表达数据文件
rownames(data)  # 查看行名（基因名）)
# 选择"Treat"组的数据
group = gsub("(.*)\\_(.*)", "\\2", colnames(data))
group
data = data[, group == "Tumor", drop = FALSE]
rt = t(data)  # 转置数据矩阵，使基因为列，样本为行

# 计算相关性矩阵
cor1 = cor(rt)  # 计算基因之间的相关性矩阵

# 设置相关性矩阵的颜色
col = c(rgb(1, 0, 0, seq(1, 0, length = 32)), rgb(0, 1, 0, seq(0, 1, length = 32)))  # 颜色渐变从红色到绿色
cor1[cor1 == 1] = 0  # 将完全相关（值为1）的部分设为0，以便在绘图时不显示
c1 = ifelse(c(cor1) >= 0, rgb(1, 0, 0, abs(cor1)), rgb(0, 1, 0, abs(cor1)))  # 根据相关性值设置颜色
col1 = matrix(c1, nc = ncol(rt))  # 将颜色矩阵按数据矩阵的列数重新排列

# 绘制环形相关性图
pdf(file = here("3_outputs/circos.pdf"), width = 7, height = 7)  # 创建PDF文件，指定宽度和高度
par(mar = c(2, 2, 2, 4))  # 设置边距
circos.par(gap.degree = c(3, rep(2, nrow(cor1) - 1)), start.degree = 180)  # 设置环形图的参数
chordDiagram(cor1, grid.col = rainbow(ncol(rt)), col = col1, transparency = 0.5, symmetric = TRUE)  # 绘制环形图
par(xpd = TRUE)  # 允许图形超出绘图区域
colorlegend(col, vertical = TRUE, labels = c(1, 0, -1), xlim = c(1.1, 1.3), ylim = c(-0.4, 0.4))  # 添加颜色图例
dev.off()  # 关闭PDF设备
circos.clear()  # 清除circos设置

# 绘制相关性矩阵的热图
pdf(file = here("3_outputs/corrplot.pdf"), width = 7, height = 7)  # 创建PDF文件，指定宽度和高度
corrplot(cor1,  # 绘制相关性图
         method = "pie",  # 使用饼图表示相关性
         order = "hclust",  # 按照层次聚类的顺序排列
         type = "upper",  # 只显示上三角矩阵
         col = colorRampPalette(c("green", "white", "red"))(50)  # 颜色渐变从绿色到白色再到红色
)
dev.off()  # 关闭PDF设备
