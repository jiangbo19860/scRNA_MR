# 从GSE的基因表达数据中提取与SNP相关的目标基因的表达量，输出日期_snpExp_genes_samples.txt，

## input是归一化为CPM并合并和去除批次效应的GSE表达csv，和从MR结果中得到的SNP对应的gene symbol的txt文档。
## 输出snpExp.txt，在3_outputs文件夹下。

rm(list = ls())  # 清空工作空间
pacman::p_load(here, data.table, limma, tibble)  # 加载必要的包

# 1. 设置文件路径（确保路径正确）
expFile <- here("1_data/pancreatic_GEO/CPM_converted/merged_CPM.csv")
geneFile <- here("3_outputs/20250804_sig_genes_44.txt")

# 2. 读取表达数据和目标基因列表
rt <- read.csv(expFile, header = TRUE, check.names = FALSE)
print(head(rt[, 1:3], 2)) # 查看前2行和前3列
gene_list <- read.table(geneFile, header = FALSE, sep = "\t", check.names = FALSE)  # gene.txt通常是制表符分隔
head(gene_list)

target_genes <- as.vector(gene_list[, 1])  # 从数据框gene_list中提取所有行和第 1 列的数据

# 3. 数据预处理（处理重复基因、过滤低表达基因）
# 转换为矩阵，第一列是基因名
rt_matrix <- as.matrix(rt[, -1, drop = FALSE])  # 逗号前的空值表示 “选择所有行”。逗号后的-1表示 “排除第 1 列”（保留从第 2 列开始的所有列）。整体含义：选择 rt 中所有行和除第 1 列之外的所有列。drop = FALSE：确保结果始终保持数据框（或矩阵）格式，即使筛选后只剩下一列也不会自动转换为向量。
rownames(rt_matrix) <- rt$gene_name  # 用gene_name列设置行名（基因名）

# 处理重复基因（同一基因取平均表达量）
unique_matrix <- avereps(rt_matrix)

# 过滤低表达基因（平均表达量>0）
filtered_matrix <- unique_matrix[rowMeans(unique_matrix) > 0, , drop = FALSE]

# 找到表达数据与目标基因的交集
common_genes <- intersect(target_genes, rownames(filtered_matrix))
if (length(common_genes) == 0) {
  stop("目标基因列表与表达数据中没有共同基因，请检查基因名是否一致！")
}

# 5. 提取目标基因的表达数据
gene_exp <- filtered_matrix[common_genes, , drop = FALSE] # 第一个逗号前表示选择filtered_matrix中行名与common_genes 向量中的基因名相匹配的行（即筛选出交集基因的表达数据）。第二个位置（两个逗号之间）：为空表示选择所有列（因为表达矩阵的列是样本，这里需要保留所有样本的表达量）。行：基因名，列：样本名
colnames(gene_exp)
class(gene_exp)  # 确认类型为矩阵

# 6. 将矩阵转换为数据框并处理基因名列
# 6.1 矩阵转数据框（保留行名和列名）
gene_df <- as.data.frame(gene_exp, stringsAsFactors = FALSE)

gene_df <- gene_df %>%
  rownames_to_column(var = "ID")  # 将行名（基因名）转换为数据框的一列
dim(gene_df)  # 输出数据框的维度，确认转换成功

# 8. 处理未匹配的基因（可选，用于日志输出）
unmatched_genes <- setdiff(target_genes, rownames(filtered_matrix))
if (length(unmatched_genes) > 0) {
  cat("\n未匹配到的基因名称：\n")
  print(unmatched_genes)
}

# 9. 生成带日期和维度的输出文件名
n_genes <- nrow(gene_df)  # 基因数量
n_samples <- ncol(gene_df) - 1  # 样本数量（减去ID列）
today_date <- format(Sys.Date(), "%Y%m%d")
output_filename <- paste0(today_date, "_snpExp_", n_genes, "genes_", n_samples, "samples.txt")
output_path <- here("3_outputs", output_filename)

# 10. 保存结果文件
write.table(gene_df, file = output_path, sep = "\t", quote = FALSE, row.names = FALSE)

# 11. 输出结果信息
cat("\n文件保存完成！\n")
cat("输出路径：", output_path, "\n")
cat("包含基因数：", n_genes, "\n")
cat("包含样本数：", n_samples, "\n")
