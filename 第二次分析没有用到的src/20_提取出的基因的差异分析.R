#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("pheatmap")
#install.packages("reshape2")
#install.packages("ggpubr")
# 清除工作空间
rm(list = ls())  # 清除工作空间
here()
# 加载所需的R包
library(limma)       # 用于基因表达数据的差异分析
library(pheatmap)    # 用于绘制热图
library(reshape2)    # 用于数据重塑
library(ggpubr)      # 用于绘制盒图

# 读取表达数据文件并进行处理
rt = read.table(here("3_outputs/snpExp.txt"), header = TRUE, sep = "\t", check.names = FALSE)
colnames(rt)
rownames(rt)
# 提取所有样本列名（排除第一列"ID"，因为第一列是基因名）
sample_cols <- colnames(rt)[-1]

# 统计以"_Control"结尾的列数
control_count <- sum(grepl("_Control$", sample_cols)) # 16列
# 统计以"_Tumor"结尾的列数
tumor_count <- sum(grepl("_Tumor$", sample_cols)) # 16列

rt = as.matrix(rt)  # 将数据转换为矩阵
rownames(rt) = rt[,1]  # 将第一列作为行名
exp = rt[, 2:ncol(rt)]  # 去除第一列，只保留表达数据
dimnames = list(rownames(exp), colnames(exp))  # 设置行列名
data = matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)  # 转换为数值矩阵
# 检查是否存在重复行名
has_duplicate <- any(duplicated(rownames(data)))

data = avereps(data)  # 合并重复的基因
exp = data  # 保存处理后的数据

# 提取样本的分组信息（假设样本名中包含分组信息）
Type = gsub("(.*)\\_(.*)", "\\2", colnames(data))  # (.*) 是 “贪婪模式”—尽可能多地匹配字符，直到遇到最后一个符合条件的下划线为止。
print(Type)

# 进行差异分析
sigVec = c()  # 存储显著性基因的标签
sigGeneVec = c()  # 存储显著性基因的名称
for (i in row.names(data)) {
  test = wilcox.test(data[i, ] ~ Type)  # 对每个基因进行Wilcoxon秩和检验
  pvalue = test$p.value  # 提取p值
  # 根据p值判断显著性水平，并使用符号表示
  Sig = ifelse(pvalue < 0.001, "***", ifelse(pvalue < 0.01, "**", ifelse(pvalue < 0.05, "*", "")))
  if (pvalue < 0.05) {
    sigVec = c(sigVec, paste0(i, Sig))  # 将显著性基因名称和符号添加到sigVec
    sigGeneVec = c(sigGeneVec, i)  # 添加显著性基因到sigGeneVec
  }
}

# 过滤出显著性基因的数据并保存到文件
data = data[sigGeneVec, ]  # 提取显著性基因的表达数据
outTab = rbind(ID = colnames(data), data)  # 添加列名作为第一行
write.table(outTab, file = here("3_outputs/diffGeneExp.txt"), sep = "\t", quote = FALSE, col.names = FALSE)  # 保存为文件
row.names(data) = sigVec  # 用显著性基因的标签替换行名

# 绘制热图展示差异基因
names(Type) = colnames(data)  # 将样本类型信息添加为列名
Type = as.data.frame(Type)  # 转换为数据框格式
pdf(file = here("3_outputs/heatmap.pdf"), width = 7, height = 4.5)  # 输出到PDF文件
pheatmap(data,
         annotation = Type,  # 添加分组注释
         color = colorRampPalette(c(rep("blue", 2), "white", rep("red", 2)))(100),  # 设置热图的颜色渐变
         cluster_cols = FALSE,  # 不对列进行聚类
         cluster_rows = TRUE,  # 对行进行聚类
         scale = "row",  # 按行标准化
         show_colnames = FALSE,  # 不显示列名
         show_rownames = TRUE,  # 显示行名
         fontsize = 7,  # 设置字体大小
         fontsize_row = 7,
         fontsize_col = 7)
dev.off()  # 关闭PDF设备

# 将表达数据转换为适合ggplot2的格式
exp = as.data.frame(t(exp))  # 转置表达数据
exp = cbind(exp, Type = Type)  # 添加分组信息
data = melt(exp, id.vars = c("Type"))  # 将数据重塑为长格式
colnames(data) = c("Type", "Gene", "Expression")  # 重命名列名

# 绘制基因表达的箱线图
p = ggboxplot(data, x = "Gene", y = "Expression", color = "Type",  # 按Type分组显示不同颜色
              xlab = "",  # 不显示x轴标签
              ylab = "Gene expression",  # 设置y轴标签
              legend.title = "Type",  # 设置图例标题
              palette = c("blue", "red"),  # 设置颜色
              add = "point",  # 添加散点
              width = 0.8)  # 设置箱线图宽度
p = p + rotate_x_text(60)  # 旋转x轴标签，便于显示
p1 = p + stat_compare_means(aes(group = Type),  # 添加p值比较
                            method = "wilcox.test",  # 使用Wilcoxon检验
                            symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),  # 设置显著性符号
                            label = "p.signif")  # 显示显著性标记

dev.off()  # 关闭所有图形设备
# 输出箱线图到PDF文件
pdf(file = here("3_outputs/boxplot.pdf"), width = 9, height = 5)
print(p1)
dev.off()  # 关闭PDF设备

# 查看筛选出的显著性基因数量
cat("显著性基因数量：", length(sigGeneVec), "\n")
# 若输出为 0，说明无显著基因，需进一步排查
print(p1)
