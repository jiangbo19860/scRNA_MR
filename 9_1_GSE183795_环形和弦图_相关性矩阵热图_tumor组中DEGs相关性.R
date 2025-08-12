# 基因表达量的相关性分析：Chord Diagram（环形和弦图）和Correlation Heatmap（相关性矩阵热图）。
# 1个输入：显著差异基因（如 p<0.05）在不同样本中的表达量数据：3_outputs/20250809_original_diffGeneExp_p05_24genes.txt，行名为基因名，列名为样本名。
# Pearson相关性分析计算两个基因的表达量之间的线性相关程度，取值范围为[-1, 1]：cor_matrix <- cor(rt)。通过基因间的相关性，探索功能关联（如协同表达的基因可能参与同一生物学通路）。
# 层次聚类（Hierarchical Clustering）计算基因间的距离（基于相关性矩阵），将表达模式相似的基因聚为一类，最终在热图中按聚类结果排序，使相似基因相邻。order = "hclust"（在corrplot函数中，表示按照层次聚类（Hierarchical Clustering） 的结果对基因进行排序，使表达模式相似（即相关性较高）的基因在热图中相邻排列。）

rm(list = ls())  # 清除工作空间

if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}
# 加载所需包：corrplot（相关性热图）、circlize（环形和弦图）、here（路径管理）
pacman::p_load(corrplot, circlize, here)

# 获取当天日期（格式：年-月-日，如20250806）
today_date <- format(Sys.Date(), "%Y%m%d")

# 读取数据文件
data <- read.table(here("3_outputs/20250810_original_diffGeneExp_p05_26genes.txt"), header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)  # header = TRUE识别表头，将第一行识别为列名（例如样本名）。"\t" 表示制表符（Tab 键），即文件是 Tab 分隔的文本文件（.txt），.csv则sep = ","。R 默认会将列名中的特殊字符（如空格、减号、括号等）修改为点号（.），避免命名错误。设置为check.names = FALSE时，保留文件中原始的列名（即使包含特殊字符，如 Sample-1、Tumor (1) 等）。row.names = 1使用文件的第 1 列作为行名（索引）。如果文件第 1 列是基因名（如 TP53、BRCA1），则设置 row.names = 1 后，这些基因名将成为数据框的行名，方便后续按基因筛选或分析。
colnames(data)  # 查看列名（样本名），确认数据格式正确)

# 筛选肿瘤组（Tumor）样本
group <- gsub("(.*)\\_(.*)", "\\2", colnames(data))  # 从样本名提取分组（如"Tumor"）
group
data_tumor <- data[, group == "Tumor", drop = FALSE]  # 仅保留肿瘤组数据

# 数据转置（样本为行，基因为列，用于计算基因间相关性）
rt <- t(data_tumor)

# 计算基因间的Pearson相关性矩阵
cor_matrix <- cor(rt)

# 相关性可视化参数设置
# 1. cor_colors <- c(正相关颜色, 负相关颜色)，最终 cor_colors 是一个包含 64 种颜色 的向量（32 种红色渐变 + 32 种绿色渐变），用于后续相关性可视化（如环形和弦图的连接线颜色）。
cor_colors <- c(
  rgb(1, 0, 0, seq(1, 0, length = 32)),  # rgb()：用于生成 RGB 颜色值，参数格式为 rgb(red, green, blue, alpha)，其中：red/green/blue：取值范围 [0, 1]，分别表示红、绿、蓝三原色的强度（1 为最强，0 为无）。 正相关：红色渐变（强度从高到低）。alpha：取值范围 [0, 1]，表示透明度（1 为完全不透明，0 为完全透明）。seq(1, 0, length = 32)：生成从 1 到 0 的 32 个等间隔数值（如 1.00、0.97、…、0.03、0.00），用于控制颜色强度或透明度的渐变。red = 1，green = 0，blue = 0：基础色为纯红色（#FF0000）。alpha = seq(1, 0, length = 32)：透明度从1逐渐变为0（即颜色从完全不透明的纯红色逐渐过渡到完全透明的红色）。
  rgb(0, 1, 0, seq(0, 1, length = 32))   # 负相关：绿色渐变（强度从低到高）
)
# 2. 去除自相关（对角线值为1的部分，避免环形图中基因自身连接）
cor_matrix[cor_matrix == 1] <- 0
# 3. 为环形图生成颜色矩阵
flat_cor <- c(cor_matrix)  # 相关性矩阵展平为向量，cor_matrix 是 n×n 的方阵（如 3 个基因则为 3×3 矩阵），存储基因间的相关系数。c(cor_matrix) 会按列优先的顺序将矩阵元素 “展开” 为一个向量。
# ifelse(条件, 满足条件时的值, 不满足条件时的值)
color_vector <- ifelse(flat_cor >= 0,
                       rgb(1, 0, 0, abs(flat_cor)),  # 正相关：红色（透明度反映强度）。
                       rgb(0, 1, 0, abs(flat_cor)))  # 负相关：绿色（透明度反映强度）
color_matrix <- matrix(color_vector, ncol = ncol(cor_matrix))  # 重塑为矩阵


# 绘制环形和弦图chordDiagram（含日期的文件名，调整基因名字体大小）
circos_filename <- paste0(today_date, "_circos_correlation.pdf")
pdf(file = here("3_outputs", circos_filename), width = 7, height = 7)
par(mar = c(2, 2, 2, 4))

# 初始化环形图参数
circos.par(
  gap.degree = c(3, rep(2, nrow(cor_matrix) - 1)),
  start.degree = 180
)

# 绘制和弦图，关闭默认标签轨道
chordDiagram(
  cor_matrix,
  grid.col = rainbow(ncol(rt)),
  col = color_matrix,
  transparency = 0.5,
  symmetric = TRUE,
  annotationTrack = "grid",  # 仅绘制网格，不绘制默认标签
  preAllocateTracks = list(track.height = uh(5, "mm"))  # 预留标签轨道高度
)

# 自定义基因名标签，控制字体大小
circos.trackPlotRegion(
  track.index = 1,  # 使用预留的轨道
  bg.border = NA,  # 去除轨道边框
  panel.fun = function(x, y) {
    sector.index <- get.cell.meta.data("sector.index")  # 获取当前基因名
    # 在轨道中心添加基因名，cex控制字体大小（0.6为示例，可更小）
    circos.text(
      x = CELL_META$xcenter,  # 标签x位置（扇区中心）
      y = 0.4,  # 标签y位置（距离内环的距离）
      labels = sector.index,
      facing = "inside",  # 文字方向（顺时针）
      cex = 0.6,  # 字体大小（核心参数，值越小字体越小）
      col = "black"  # 字体颜色
    )
  }
)

# 添加颜色图例
par(xpd = TRUE)
colorlegend(
  cor_colors,
  vertical = TRUE,
  labels = c("1", "0", "-1"),
  xlim = c(1.1, 1.2),  # 进一步缩小x轴范围（比原来更窄）
  ylim = c(-0.25, 0.25), # 缩小y轴范围（比原来更矮）
  cex = 0.7,  # 标签字体缩小
  lwd = 2     # 颜色条边缘线条调细
)

dev.off()
circos.clear()

# 绘制相关性矩阵热图（含日期的文件名）
heatmap_filename <- paste0(today_date, "_corrplot_heatmap.pdf")  # 新文件名：日期+功能
pdf(file = here("3_outputs", heatmap_filename), width = 7, height = 7)
corrplot(
  cor_matrix,
  method = "pie",  # 用扇形面积表示相关性强度
  order = "hclust",  # 按层次聚类排序基因（相似表达模式的基因相邻）
  type = "upper",  # 仅显示上三角矩阵（避免重复）
  col = colorRampPalette(c("green", "white", "red"))(50),  # 颜色：绿色（负）→白色（0）→红色（正）
  title = "Expression correlation matrix of significant genes in tumor samples",  # 热图标题
  cex.main = 0.9,  # 缩小标题字体
  mar = c(0, 0, 2, 0),  # 调整边距（增加顶部边距避免标题被截）
  tl.col = "black",  # 核心修改1：基因标签文本颜色为黑色
)
dev.off()

# 输出结果提示
cat("环形和弦图已保存至：", here("3_outputs", circos_filename), "\n")
cat("相关性热图已保存至：", here("3_outputs", heatmap_filename), "\n")
