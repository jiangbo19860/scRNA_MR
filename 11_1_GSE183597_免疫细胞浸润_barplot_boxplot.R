# 堆积条形图：可视化免疫细胞组成比例。1个输入：免疫细胞组成结果文件inputFile = here("3_outputs/CIBERSORT-Results.txt")，输出2个图：堆积条形图（免疫细胞浸润组成图）和箱线图（免疫细胞比例组间差异图）

rm(list = ls())  # 清除工作空间
pacman::p_load(reshape2, ggpubr, corrplot, here)  # 加载必要的R包
# 添加：定义当天日期（格式为YYYYMMDD，如20250807）
today_date <- format(Sys.Date(), "%Y%m%d")

# 设置输入文件路径
inputFile = here("3_outputs/20250810_CIBERSORT-Results.txt")  # 输入的免疫细胞组成结果文件


# 读取免疫细胞组成结果文件
rt = read.table(inputFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

# 将样本分为对照组和治疗组（基于行名进行分类）
con = grepl("_Control", rownames(rt), ignore.case = TRUE)  # 匹配行名中包含"Control"的样本
Tumor = grepl("_Tumor", rownames(rt), ignore.case = TRUE)  # 匹配行名中包含"Tumor"的样本
conData = rt[con, ]  # 提取对照组数据
TumorData = rt[Tumor, ]  # 提取治疗组数据
conNum = nrow(conData)  # 计算对照组样本数量
TumorNum = nrow(TumorData)  # 计算治疗组样本数量
data = t(rbind(conData, TumorData))  # 合并并转置数据，使其适合绘制堆积条形图

# 绘制堆积条形图
# 原代码：
# pdf(file = here("3_outputs/免疫细胞浸润barplot.pdf"), width = 13, height = 7.5)

# 修改后：
pdf(file = here("3_outputs", paste0(today_date, "_ImmuneCellInfiltration__Barplot.pdf")), width = 13, height = 7.5)

col = rainbow(nrow(data), s = 0.7, v = 0.7)  # 为每个免疫细胞类型生成颜色
par(las = 1, mar = c(8, 5, 4, 16), mgp = c(3, 0.6, 0), cex.axis = 1)  # 设置绘图参数。mgp = c(3, 0.1, 0))：3表示轴标题与轴线的距离，0.1表示刻度标签与轴线的距离，0表示刻度线与轴线的距离
a1 = barplot(data, col = col, xaxt = "n", yaxt = "n", ylab = "Relative Percent", cex.lab = 1.2)  # 绘制条形图
a2 = axis(2, tick = FALSE, labels = FALSE)  # 添加y轴刻度
axis(2, a2, paste0(a2 * 100, "%"))  # 将y轴刻度转换为百分比形式
par(srt = 0, xpd = TRUE)  # 设置文本方向和扩展绘图区域
rect(xleft = a1[1] - 0.5, ybottom = -0.01, xright = a1[conNum] + 0.5, ytop = -0.06, col = "#0088FF")  # 绘制对照组矩形
text(a1[conNum] / 2, -0.035, "Control", cex = 1.2)  # 标注对照组
rect(xleft = a1[conNum] + 0.5, ybottom = -0.01, xright = a1[length(a1)] + 0.5, ytop = -0.06, col = "#FF5555")  # 绘制治疗组矩形
text((a1[length(a1)] + a1[conNum]) / 2, -0.035, "Tumor", cex = 1.2)  # 标注治疗组
ytick2 = cumsum(data[, ncol(data)])  # 计算堆积条形图的顶部位置
ytick1 = c(0, ytick2[-length(ytick2)])  # 计算底部位置
legend(par('usr')[2] * 0.98, par('usr')[4], legend = rownames(data), col = col, pch = 15, bty = "n", cex = 1)  # 添加图例
dev.off()  # 关闭PDF输出

################## 箱线图 ##################
# 转换数据以适合ggplot2绘图
Type = gsub("(.*)\\_(.*)", "\\2", rownames(rt))  # 提取对照组和治疗组的标签
data = cbind(as.data.frame(t(data)), Type)  # 将组别信息加入数据
data = melt(data, id.vars = c("Type"))  # 将数据重构为长格式，以便使用ggplot2
colnames(data) = c("Type", "Immune", "Expression")  # 重命名列

# 绘制箱线图
group = levels(factor(data$Type))  # 获取组别信息
bioCol = c("#0066FF", "#FF0000", "#6E568C", "#7CC767", "#223D6C", "#D20A13", "#FFD121", "#088247", "#11AA4D")  # 设置颜色
bioCol = bioCol[1:length(group)]  # 为组别分配颜色
boxplot = ggboxplot(data, x = "Immune", y = "Expression", color = "Type",
                    xlab = "",  # 设置x轴标签为空
                    ylab = "Fraction",  # 设置y轴标签
                    legend.title = "Type",  # 图例标题
                    add = "point",  # 在箱线图上添加数据点
                    width = 0.8,  # 设置箱线宽度
                    palette = bioCol) +
  rotate_x_text(50) +  # 旋转x轴标签
  stat_compare_means(aes(group = Type),  # 进行组间的统计学比较，并添加显著性符号
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "")),
                     label = "p.signif")

# 原代码：
# pdf(file = here("3_outputs/免疫细胞浸润immune.diff2.pdf"), width = 8, height = 6)

# 修改后：
pdf(file = here("3_outputs", paste0(today_date, "_ImmuneCellInfiltration_Boxplot.pdf")), width = 8, height = 6)
print(boxplot)  # 打印并保存图形
dev.off()  # 关闭PDF输出
