
# 分析非癌细胞亚群在不同分组（group）中的比例差异，并通过森林图可视化结果。步骤可概括为：筛选非癌细胞和目标区域（core/edge）的元信息 → 2. 构建亚群标签和区域 - 患者组合标签 → 3. 计算各亚群在不同组合中的比例 → 4. 通过线性回归lm分析亚群比例与 group 的关联（提取 t 值和 p 值） → 5. 绘制森林图展示关联强度和显著性，并通过灰色背景标记主要细胞类型。

rm(list = ls())

# 加载所需R包：
# here：路径管理；Seurat：单细胞数据分析；dplyr：数据处理；
# ggplot2/cowplot/patchwork：可视化；ggrepel：标签避免重叠；hdf5r：处理HDF5文件；reshape2用于数据重塑（Data Reshaping），其中的melt ()：将宽格式数据转换为长格式，dcast ()：将长格式数据转换为宽格式（数据框输出）tidyr中的pivot_longer ()：将宽格式转换为长格式（替代 reshape2::melt ()）
pacman::p_load(here, Seurat, dplyr, ggplot2, cowplot, patchwork, ggrepel, hdf5r, reshape2, tidyr)

# 显示当前工作目录的路径（基于here包的项目路径管理）
here()

# 加载细胞亚群数据（可能包含细胞聚类、注释等信息）
load("All_cell_subcluster.Rds")
# 加载之前保存的区域相关元信息（meta2和meta3，包含细胞的区域、患者等标签）
load("region.Rds")

# 从meta2中筛选细胞：排除"middle"区域，且排除"Cancer"细胞类型
# 目的：聚焦于core/edge区域的非癌细胞（如免疫细胞、基质细胞等）
meta <- subset(meta2, region != "middle" & cell != "Cancer")

# 统计筛选后的数据中，每个患者（PatientNumber）在core/edge区域的细胞数量分布
table(meta$PatientNumber, meta$region)

# 给meta添加新列"Cluster"：从All.meta中提取对应细胞的dbCluster（聚类ID）
meta$Cluster <- All.meta[rownames(meta), "dbCluster"]

# 给meta添加新列"cell"：从All.meta中提取对应细胞的细胞类型（可能是更详细的注释）
meta$cell <- All.meta[rownames(meta), "cell"]

# 给meta添加新列"sub"：将细胞类型（cell）和聚类ID（Cluster）组合为亚群标签（如"B_cell_1"）
meta$sub <- paste0(meta$cell, "_", meta$Cluster)

# 给meta添加新列"region2"：将区域（region）和患者编号（PatientNumber）组合（如"core_1"）
# 目的：区分不同患者的同一区域（避免混淆不同患者的core/edge）
meta$region2 <- paste0(meta$region, "_", meta$PatientNumber)

# 统计每个亚群（sub）在每个"区域-患者"组合（region2）中的细胞数量，转换为数据框
data_sub <- as.data.frame(table(meta$sub, meta$region2))
# 结果列：Var1=亚群（sub），Var2=region2，Freq=细胞数量

# 统计每个患者在core/edge区域的总细胞数（用于后续计算比例）
data <- as.data.frame(table(meta$PatientNumber, meta$region))

# 将data_sub从长格式转换为宽格式：行=region2，列=亚群（sub），值=细胞数量
data_plot <- spread(data_sub, Var1, Freq)
# 将行名设置为region2（便于后续匹配）
rownames(data_plot) <- data_plot$Var2
# 移除原Var2列（已作为行名）
data_plot <- data_plot[, -1]
# 添加"total"列：每个region2的总细胞数（从data中提取，与region2匹配）
data_plot$total <- data$Freq

# 计算每个亚群在对应region2中的比例（亚群细胞数 / 该region2总细胞数）
for (i in 1:ncol(data_plot)) {
  data_plot[, i] <- as.numeric(data_plot[, i]) / (data_plot$total)
}

# 移除"total"列（比例计算完成，不再需要）
data_plot <- data_plot[, -ncol(data_plot)]
# 添加"group"列：将样本分为两组（0和1，各5个样本，可能是对照/处理或不同表型分组）
data_plot$group <- as.numeric(c(rep(0, 5), rep(1, 5)))

# 输出meta中所有独特的亚群标签（sub），用于确认细胞亚群名称（后续循环需用到）
dput(unique(meta$sub))

# 定义需要分析的细胞亚群列表（与上述输出的sub匹配，确保覆盖所有目标亚群）
cell_name <- c("B_cell_1", "EC_2", "B_cell_4", "B_cell_5", "B_cell_3", "T_cell_1", 
               "Fibroblast_5", "B_cell_2", "Myeloid_3", "T_cell_3", "Fibroblast_3", 
               "T_cell_9", "T_cell_5", "T_cell_2", "B_cell_7", "B_cell_6", "T_cell_6", 
               "EC_4", "Myeloid_2", "Myeloid_5", "Myeloid_11", "Myeloid_12", 
               "T_cell_8", "Alveolar_1", "T_cell_7", "Fibroblast_2", "B_cell_8", 
               "EC_6", "Alveolar_2", "Myeloid_9", "Myeloid_7", "Myeloid_1", 
               "EC_3", "Epithelial_2", "Alveolar_4", "Fibroblast_6", "Alveolar_7", 
               "T_cell_4", "Alveolar_3", "Fibroblast_4", "Myeloid_10", "Myeloid_4", 
               "Fibroblast_1", "Myeloid_6", "Epithelial_1", "Myeloid_8", "EC_5", 
               "Alveolar_5", "EC_1", "Fibroblast_7")

# 创建空数据框，用于存储每个亚群的统计结果（t值和p值）
data_plot2 <- data.frame(cell = NULL, t = NULL, p = NULL)

# 循环分析每个细胞亚群与group的关联（通过线性回归）
for (i in cell_name) {
  # 提取当前亚群的比例数据和group分组
  data_tmp <- data_plot[, c(i, "group")]
  # 构建线性回归模型：group ~ 亚群比例（分析亚群比例是否随group变化）
  x <- lm(group ~ ., data = data_tmp)
  # 提取回归系数的t值（衡量效应大小的统计量）
  t <- summary(x)$coefficients[, 3]
  # 提取回归系数的p值（显著性）
  p <- summary(x)$coefficients[, 4]
  # 将当前亚群的名称、t值、p值存入临时数据框
  data_tmp2 <- data.frame(cell = i, t = t[2], p = p[2])  # t[2]/p[2]对应亚群比例的系数
  # 合并到结果数据框
  data_plot2 <- rbind(data_plot2, data_tmp2)
}

# 加载ggpubr包，用于增强绘图功能
library(ggpubr)

# 对p值进行分组标记：
# a：p < 0.05（显著）；b：0.05 < p < 0.1（接近显著）；c：p > 0.1（不显著）
data_plot2$pvalue <- ifelse(data_plot2$p > 0.1, "c", ifelse(data_plot2$p > 0.05, "b", "a"))

# 从细胞亚群名称中提取主要细胞类型（如从"B_cell_1"中提取"B_cell"）
data_plot2$major <- data.frame(
  sapply(data_plot2$cell, function(x) unlist(strsplit(x, '_'))[1]), 
  stringsAsFactors = F
)[, 1]

# 加载cowplot包，用于组合图形
library(cowplot)

# 绘制森林图（展示各亚群的t值及显著性）：
# x轴为细胞亚群，y轴为t值（反映亚群比例与group的关联强度）
ggplot(data_plot2, aes(x = cell, y = t, fill = pvalue, shape = pvalue)) + 
  # 添加从y=0到t值的线段（展示t值大小）
  geom_segment(
    aes(y = 0, x = cell, yend = t, xend = cell), 
    color = "#36648B"  # 线段颜色
  ) +
  # 添加点（大小3，圆形带边框，边框颜色固定，填充色按pvalue分组）
  geom_point(
    stat = 'identity', 
    size = 3, shape = 21,  # 21为带边框的圆形
    stroke = 1, color = "#36648B"  # 边框粗细和颜色
  ) +
  # 自定义填充色：显著（a）深蓝色，接近显著（b）中蓝色，不显著（c）白色
  scale_fill_manual(
    values = c("#000080", "#4F94CD", "#FFFFFF"),
    name = "P value",  # 图例名称
    breaks = c("a", "b", "c"),  # 分组标签
    labels = c("p < 0.05", "p < 0.1", "P > 0.1")  # 图例显示文本
  ) +
  theme_classic() +  # 使用经典主题（无网格线）
  theme(
    axis.text.x = element_text(angle = 90),  # x轴标签旋转90度（避免重叠）
    panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid")  # 添加上下左右边框
  ) +
  scale_x_discrete(position = "top") +  # x轴标签显示在顶部
  scale_y_reverse() +  # y轴方向反转（从下到上值减小）
  geom_hline(yintercept = 0, lty = 5) +  # 添加y=0的虚线（参考线）
  geom_hline(yintercept = 1.7, lty = 5, color = "gray") +  # 添加t值=1.7的灰色虚线（可能是临界值）
  geom_hline(yintercept = 2.2, lty = 1, color = "gray")  # 添加t值=2.2的灰色实线（可能是显著临界值）

# 绘制增强版森林图（在上述基础上添加细胞类型分组的灰色背景）：
ggplot(data_plot2, aes(x = cell, y = t, fill = pvalue, shape = pvalue)) + 
  geom_segment(aes(y = 0, x = cell, yend = t, xend = cell), color = "#36648B") +
  geom_point(stat = 'identity', size = 3, shape = 21, stroke = 1, color = "#36648B") +
  scale_fill_manual(values = c("#000080", "#4F94CD", "#FFFFFF"),
                    name = "P value",
                    breaks = c("a", "b", "c"),
                    labels = c("p <0.05", "p<0.1", "P>0.1")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90),
        panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid")) +
  scale_x_discrete(position = "top") +
  scale_y_reverse() +
  geom_hline(yintercept = 0, lty = 5) +
  geom_hline(yintercept = 1.7, lty = 5, color = "gray") +
  geom_hline(yintercept = 2.2, lty = 1, color = "gray") +
  # 添加灰色背景矩形，标记B细胞亚群（B_cell_1到B_cell_8）
  annotate("rect",
           xmin = "B_cell_1", xmax = "B_cell_8",
           ymin = -Inf, ymax = Inf,
           alpha = 0.5, fill = "gray") +
  # 标记上皮细胞亚群（Epithelial_1到Epithelial_2）
  annotate("rect",
           xmin = "Epithelial_1", xmax = "Epithelial_2",
           ymin = -Inf, ymax = Inf,
           alpha = 0.5, fill = "gray") +
  # 标记髓系细胞亚群（Myeloid_1到Myeloid_2）
  annotate("rect",
           xmin = "Myeloid_1", xmax = "Myeloid_2",
           ymin = -Inf, ymax = Inf,
           alpha = 0.5, fill = "gray") +
  # 补充标记更多髓系细胞亚群（Myeloid_1到Myeloid_9）
  annotate("rect",
           xmin = "Myeloid_1", xmax = "Myeloid_9",
           ymin = -Inf, ymax = Inf,
           alpha = 0.5, fill = "gray")

# 将图形保存为PDF文件，尺寸为高6英寸、宽10英寸，保存路径为"Fig/Fig6b.pdf"
ggsave(filename = "Fig/Fig6b.pdf", height = 6, width = 10)

# 该森林图的核心特征：
# 横轴（x 轴）：表示细胞亚群（如B_cell_1、T_cell_3等），按主要细胞类型（通过灰色背景矩形标注，如 B 细胞、上皮细胞、髓系细胞）分组。
# 纵轴（y 轴）：表示线性回归的t值（反映细胞亚群比例与group分组的关联强度，绝对值越大，关联越强）。
# 线段：从纵轴 0 点延伸到对应t值的竖线，直观展示每个亚群的效应大小。
# 点：位于线段末端的点，颜色（填充色）根据p值分组（深蓝色：p<0.05；中蓝色：p<0.1；白色：p>0.1），表示关联的显著性。
# 参考线：
# 虚线（y=0）：表示无关联（t 值为 0）。
# 灰色虚线（y=1.7）和实线（y=2.2）：可能是t值的临界值（用于判断统计显著性）。