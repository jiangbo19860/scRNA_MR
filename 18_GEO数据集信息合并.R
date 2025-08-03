rm(list = ls())
# 加载所需包（用于批次效应校正，可选但推荐）
# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
# }
# BiocManager::install(c("sva", "limma"), update = FALSE)
library(sva)    # 批次效应校正
library(limma)  # 重复基因处理

# 1. 读取数据（确保之前已加载GSE119794和GSE171485）
GSE119794 <- read.csv("/Users/lijiangbo/scRNA_MR/1_data/pancreatic_GEO/CPM_converted/GSE119794_CPM.csv", header=TRUE)
GSE171485 <- read.csv("/Users/lijiangbo/scRNA_MR/1_data/pancreatic_GEO/CPM_converted/GSE171485_CPM.csv", header=TRUE)

# 2. 统一基因名列名（确保两个数据集的基因名列名一致）
colnames(GSE119794)[1] <- "gene_name"  # 将GSE119794的"gene_id"改为"gene_name"，与GSE171485统一

# 3. 处理重复基因（同一基因名可能有多个条目，取平均值）
# 处理GSE119794
gse119794_matrix <- as.matrix(GSE119794[, -1, drop=FALSE])  # 提取表达值矩阵（去除基因名列）
rownames(gse119794_matrix) <- GSE119794$gene_name          # 基因名设为行名
gse119794_unique <- avereps(gse119794_matrix)              # 重复基因取平均

# 处理GSE171485
gse171485_matrix <- as.matrix(GSE171485[, -1, drop=FALSE])
rownames(gse171485_matrix) <- GSE171485$gene_name
gse171485_unique <- avereps(gse171485_matrix)

# 4. 筛选共同基因（取两个数据集的基因交集）
common_genes <- intersect(rownames(gse119794_unique), rownames(gse171485_unique))
cat("共同基因数量：", length(common_genes), "\n")  # 输出交集基因数，确认筛选有效, 17501

# 5. 按共同基因筛选并添加批次前缀（避免样本名冲突）
# GSE119794样本名添加"GSE119794_"前缀
colnames(gse119794_unique) <- paste0("GSE119794_", colnames(gse119794_unique))
gse119794_filtered <- gse119794_unique[common_genes, , drop=FALSE]  # 仅保留共同基因, drop=FALSE是防止筛选结果被自动降维。当您从矩阵或数据框中筛选子集时，R 默认会自动降维：如果筛选后结果只有1 行或 1 列，会被转换为向量（一维结构）。如果筛选后结果只有1 个元素，会被转换为标量（单个值）。drop = FALSE的作用是禁用这种自动降维，强制结果保持原始的矩阵 / 数据框结构（二维），即使筛选后只有 1 行或 1 列。

# GSE171485样本名添加"GSE171485_"前缀
colnames(gse171485_unique) <- paste0("GSE171485_", colnames(gse171485_unique))
gse171485_filtered <- gse171485_unique[common_genes, , drop=FALSE]  # 仅保留共同基因

# 6. 列绑定合并数据集
merged_matrix <- cbind(gse119794_filtered, gse171485_filtered)
dim(merged_matrix)  # 输出合并后的矩阵维度，确认合并成功, 17501 32.

# 7. 批次效应校正（关键步骤，消除两个数据集间的技术差异）
# 定义批次信息（1=GSE119794，2=GSE171485）
batch <- c(
  rep(1, ncol(gse119794_filtered)),  # GSE119794的样本批次为1
  rep(2, ncol(gse171485_filtered))   # GSE171485的样本批次为2
)

# 校正批次效应（若样本量小，使用par.prior=FALSE）
if (length(common_genes) >= 500) {  # 足够多的基因时使用参数先验
  merged_corrected <- ComBat(merged_matrix, batch = batch, par.prior = TRUE)
} else {
  merged_corrected <- ComBat(merged_matrix, batch = batch, par.prior = FALSE)
}

# 8. 转换为数据框并添加基因名列
merged_final <- data.frame(
  gene_name = rownames(merged_corrected),
  merged_corrected,
  check.names = FALSE  # 保留原始列名（含特殊字符）
)

colnames(merged_final)
head(merged_final)
str(merged_final)

# 9. 保存结果到指定文件夹
# 定义目标路径（确保文件夹存在）
target_path <- "/Users/lijiangbo/scRNA_MR/1_data/pancreatic_GEO/merged_GSE119794_GSE171485.csv"

# 保存文件
write.csv(merged_final, file = target_path, row.names = FALSE)

# 提示保存成功
cat("整合后的数据已保存到：", target_path, "\n")
cat("整合后的数据维度：", nrow(merged_final), "个基因，", ncol(merged_final)-1, "个样本\n") # 16244个基因， 32 个样本。
