## input是归一化为CPM并合并和去除批次效应的GSE表达csv，和从MR结果中得到的SNP对应的gene symbol的txt文档。
## 输出snpExp.txt，在3_outputs文件夹下。

rm(list = ls())  # 清空工作空间
library(limma)  # 加载limma包，用于处理重复基因

# 1. 设置文件路径（确保路径正确）
expFile <- "/Users/lijiangbo/scRNA_MR/1_data/pancreatic_GEO/merged_GSE119794_GSE171485.csv"  # CSV格式的表达矩阵
geneFile <- "/Users/lijiangbo/scRNA_MR/3_outputs/汇总有显著意义的MR结果/sig_genes.txt"

# 2. 读取表达数据（CSV文件用sep=","，而非"\t"）
rt <- read.csv(expFile, header = TRUE, check.names = FALSE)  # 用read.csv读取CSV文件，自动识别逗号分隔
# 检查数据是否读取正确（查看前2行和前3列）
cat("表达数据预览：\n")
print(head(rt[, 1:3], 2))

# 3. 数据预处理（处理重复基因、过滤低表达基因）
# 转换为矩阵，第一列是基因名
rt_matrix <- as.matrix(rt[, -1, drop = FALSE])  # 去除基因名列，保留表达值
rownames(rt_matrix) <- rt$gene_name  # 用gene_name列设置行名（基因名）

# 处理重复基因（同一基因取平均表达量）
unique_matrix <- avereps(rt_matrix)

# 过滤低表达基因（平均表达量>0）
filtered_matrix <- unique_matrix[rowMeans(unique_matrix) > 0, , drop = FALSE]

# 4. 读取目标基因列表并筛选交集
gene_list <- read.table(geneFile, header = TRUE, sep = "\t", check.names = FALSE)  # gene.txt通常是制表符分隔
target_genes <- as.vector(gene_list[, 1])  # 提取基因列表
head(gene_list)

# 找到表达数据与目标基因的交集
common_genes <- intersect(target_genes, rownames(filtered_matrix))
if (length(common_genes) == 0) {
  stop("目标基因列表与表达数据中没有共同基因，请检查基因名是否一致！")
}

# 5. 提取目标基因的表达数据
gene_exp <- filtered_matrix[common_genes, , drop = FALSE]

# 6. 生成输出文件（格式与原始需求一致）
# 第一行添加样本名（ID行），后续行是基因表达值
out_tab <- rbind(ID = colnames(gene_exp), gene_exp)
colnames(out_tab)
head(out_tab)
str(out_tab)

unmatched_genes <- setdiff(target_genes, rownames(filtered_matrix))
if (length(unmatched_genes) > 0) {
  cat("未匹配到的基因名称：\n")
  print(unmatched_genes)
}

output_path <- here("3_outputs/snpExp.txt")
write.table(out_tab, file = output_path, sep = "\t", quote = FALSE, col.names = FALSE)


# 7. 输出结果信息
cat("提取完成！结果保存为：snpExp.txt\n")
cat("目标基因数量：", length(target_genes), "\n")
cat("成功匹配的基因数量：", length(common_genes), "\n")
cat("输出文件维度：", nrow(out_tab)-1, "个基因，", ncol(out_tab), "个样本\n")

