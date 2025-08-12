# heatmap: 差异基因表达量与免疫细胞浸润比例之间的相关性，并通过热图可视化这种关联。
# 2个输入：差异基因表达数据expFile = here("3_outputs/20250809_original_diffGeneExp_p05_24genes.txt")，免疫细胞浸润结果文件immFile = here("3_outputs/CIBERSORT-Results.txt")。
# 1个输出：3_outputs/immune_gene_correlation_heatmap.pdf

rm(list = ls())  # 清除工作空间
pacman::p_load(limma, reshape2, tidyverse, ggplot2, here)  # 加载必要的R包

# 设置输入文件路径
expFile = here("3_outputs/20250810_original_diffGeneExp_p05_26genes.txt")  # 差异基因表达数据文件
immFile = here("3_outputs/20250810_CIBERSORT-Results.txt")  # 免疫细胞浸润结果文件

# 读取基因表达数据文件
rt = read.table(expFile, header = TRUE, sep = "\t", check.names = FALSE)
rt = as.matrix(rt)  # 将数据转换为矩阵
rownames(rt) = rt[, 1]  # 将第一列作为行名
exp = rt[, 2:ncol(rt)]  # 提取表达数据（去除第一列基因名称）
dimnames = list(rownames(exp), colnames(exp))  # 创建行名和列名的维度名称
data = matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)  # 将数据转换为数值矩阵
data = avereps(data)  # 对重复数据进行平均

# 筛选处理组样本数据
group = gsub("(.*)\\_(.*)", "\\2", colnames(data))  # 提取样本组别信息
data = data[, group == "Tumor", drop = FALSE]  # 仅保留处理组的样本
data = t(data)  # 转置矩阵，使样本作为行，基因作为列

# 读取免疫细胞浸润数据文件，并对样本进行匹配
immune = read.table(immFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
sameSample = intersect(row.names(data), row.names(immune))  # 获取两个数据集中共有的样本
data = data[sameSample, , drop = FALSE]  # 保留相同样本的基因表达数据
immune = immune[sameSample, , drop = FALSE]  # 保留相同样本的免疫细胞数据

# 计算每个基因与每个免疫细胞类型的相关性
outTab = data.frame()  # 初始化空数据框
for (cell in colnames(immune)) {  # 遍历免疫细胞类型
  if (sd(immune[, cell]) == 0) { next }  # 如果免疫细胞数据的标准差为0，跳过该细胞
  for (gene in colnames(data)) {  # 遍历每个基因
    x = as.numeric(immune[, cell])  # 提取免疫细胞数据
    y = as.numeric(data[, gene])  # 提取基因表达数据
    corT = cor.test(x, y, method = "spearman", exact = FALSE)  # 进行Spearman相关性检验
    cor = corT$estimate  # 提取相关系数
    pvalue = corT$p.value  # 提取p值
    text = ifelse(pvalue < 0.001, "***", ifelse(pvalue < 0.01, "**", ifelse(pvalue < 0.05, "*", "")))  # 根据p值标注显著性符号
    outTab = rbind(outTab, cbind(Gene = gene, Immune = cell, cor, text, pvalue))  # 将结果保存到outTab中
  }
}

# 生成相关性热图（修改文件名）
outTab$cor = as.numeric(outTab$cor)  # 将相关系数转换为数值型

# 定义当天日期（格式：YYYYMMDD）
today_date <- format(Sys.Date(), "%Y%m%d")

# 专业英文文件名：日期 + 免疫细胞与差异基因相关性热图
pdf_file <- paste0(today_date, "_immune_gene_correlation_heatmap.pdf")

# 保存为PDF（含日期和专业英文名称）
pdf(file = here("3_outputs", pdf_file), width = 7, height = 5)

ggplot(outTab, aes(Immune, Gene)) +
  geom_tile(aes(fill = cor), colour = "grey", linewidth = 1) +
  scale_fill_gradient2(low = "#5C5DAF", mid = "white", high = "#EA2E2D") +
  geom_text(aes(label = text), col = "black", size = 3) +
  theme_minimal() +
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8, face = "bold"),
        axis.text.y = element_text(size = 8, face = "bold")) +
  labs(fill = paste0("***  p<0.001", "\n", "**  p<0.01", "\n", " *  p<0.05", "\n", "\n", "Correlation")) +
  scale_x_discrete(position = "bottom")

dev.off()
