# 从GSE的基因表达数据中提取与SNP相关的目标基因的表达量，输出日期_snpExp_genes_samples.txt，

## input是归一化为CPM并合并和去除批次效应的GSE表达csv，和从MR结果中得到的SNP对应的gene symbol的txt文档。
## 输出snpExp.txt，在3_outputs文件夹下。

rm(list = ls())  # 清空工作空间
pacman::p_load(here, data.table, limma, tibble, dplyr)  # 加载必要的包

expFile <- here("1_data/GEO/GSE183795/GSE183795_normalized_matrix_final.txt")
geneFile <- here("3_outputs/20250804_sig_genes_44.txt")

# 2. 读取表达数据和目标基因列表（修正版）
# 读取表达数据txt文件（制表符分隔）
# 重新读取表达数据，明确指定非数值字符为NA
rt <- read.delim(
  expFile,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("NA", "N/A", "?", "", " ")  # 将常见非数值字符识别为NA
)
# 查看前2行和前3列（如果文件列数较少，会自动显示实际存在的列）
print(head(rt[, 1:min(3, ncol(rt))], 2)) # min(3, ncol(rt))：取 “3” 和 “总列数” 中的较小值。rt[, 1:min(3, ncol(rt))]：提取rt中所有行、且列数范围为 “第 1 列到上一步得到的列数” 的数据（即最多前 3 列）。head(..., 2)：从上述提取的子数据框中，再提取前 2 行（只看前 2 行数据）。head()函数的作用是提取数据的前 n 行，其语法为head(x, n)，其中：x是要提取数据的对象（如数据框、矩阵、向量等）；n是一个整数，表示要提取的行数（默认值为 6，即默认提取前 6 行）。
colnames(rt)
# 将行名转换为第一列，列名为"gene_name"
rt <- rt %>%
  rownames_to_column(var = "gene_name")

# 验证结果：查看转换后的前2行和前4列（含新的gene_name列）
print(head(rt[, 1:min(4, ncol(rt))], 2))

# 读取基因列表txt文件
gene_list <- read.delim(geneFile, header = FALSE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
head(gene_list)

target_genes <- as.vector(gene_list[, 1])  # 从数据框gene_list中提取所有行和第 1 列的数据

# 3. 数据预处理
head(rt$gene_name)
# 3.2 提取表达量列（第4列及以后）并转换为矩阵
exp_columns <- 2:ncol(rt)
rt_matrix <- as.matrix(rt[, exp_columns, drop = FALSE])

# 3.3 设置矩阵行名
rownames(rt_matrix) <- rt$gene_name
head(rownames(rt_matrix))  # 再次确认行名正确

# 3.4 处理重复基因（指定ID参数）
unique_matrix <- avereps(rt_matrix, ID = rownames(rt_matrix))

# 确保矩阵为数值型
if (!is.numeric(unique_matrix)) {
  unique_matrix <- matrix(as.numeric(unique_matrix),
                          nrow = nrow(unique_matrix),
                          ncol = ncol(unique_matrix),
                          dimnames = dimnames(unique_matrix))
}

# 处理NA值
unique_matrix[is.na(unique_matrix)] <- 0

# 3.5 过滤低表达基因
filtered_matrix <- unique_matrix[rowMeans(unique_matrix) > 0, , drop = FALSE]

# 验证结果
class(filtered_matrix)  # 应为"matrix"
is.numeric(filtered_matrix)  # 应为TRUE
dim(filtered_matrix)  # 输出过滤后的基因数和样本数

# 4. 计算目标基因与过滤后表达矩阵的交集（新增步骤）
common_genes <- intersect(target_genes, rownames(filtered_matrix))

# 检查交集是否为空
if (length(common_genes) == 0) {
  stop("目标基因列表与表达数据中没有共同基因，请检查基因名是否一致！")
}

# 5. 提取目标基因的表达数据
gene_exp <- filtered_matrix[common_genes, , drop = FALSE] # 第一个逗号前表示选择filtered_matrix中行名与common_genes 向量中的基因名相匹配的行（即筛选出交集基因的表达数据）。第二个位置（两个逗号之间）：为空表示选择所有列（因为表达矩阵的列是样本，这里需要保留所有样本的表达量）。行：基因名，列：样本名
head(colnames(gene_exp))
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

# 检查gene_df中列名是否以_Control或_Tumor结尾，并统计数量
# 排除第一列"ID"（基因名列），只检查样本列
sample_cols <- colnames(gene_df)[-1]

# 统计以_Control结尾的列数
control_count <- sum(grepl("_Control$", sample_cols))

# 统计以_Tumor结尾的列数
tumor_count <- sum(grepl("_Tumor$", sample_cols))

# 输出结果
cat("\n样本列统计结果：\n")
cat("以_Control结尾的列数：", control_count, "\n")
cat("以_Tumor结尾的列数：", tumor_count, "\n")
cat("总样本列数：", length(sample_cols), "\n")

# 10. 保存结果文件
write.table(gene_df, file = output_path, sep = "\t", quote = FALSE, row.names = FALSE)

# 11. 输出结果信息
cat("\n文件保存完成！\n")
cat("输出路径：", output_path, "\n")
cat("包含基因数：", n_genes, "\n")
cat("包含样本数：", n_samples, "\n")
