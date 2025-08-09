# 加载所需包
library(limma)       # 处理重复基因
library(biomaRt)     # 获取基因长度（FPKM转CPM需要）

# 1. 设置文件路径
raw_counts_file <- "/Users/lijiangbo/scRNA_MR/1_data/pancreatic_GEO/GSE119794_10T10C/GSE119794_mRNA_samples_updated.csv"
fpkm_file <- "/Users/lijiangbo/scRNA_MR/1_data/pancreatic_GEO/GSE171485_6T6C/GSE171485_expression_with_gene_names.csv"
output_dir <- "/Users/lijiangbo/scRNA_MR/1_data/pancreatic_GEO/CPM_converted"  # 输出文件夹

# 2. 创建输出文件夹（若不存在）
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 3. 读取原始数据
## 读取raw counts数据（GSE119794）
raw_counts <- read.csv(raw_counts_file, header = TRUE, check.names = FALSE)
colnames(raw_counts)[1] <- "gene_name"  # 确保基因名列名为gene_name
cat("GSE119794原始数据维度：", nrow(raw_counts), "个基因，", ncol(raw_counts)-1, "个样本\n")

## 读取FPKM数据（GSE171485）
fpkm_data <- read.csv(fpkm_file, header = TRUE, check.names = FALSE)
colnames(fpkm_data)[1] <- "gene_name"   # 确保基因名列名为gene_name
cat("GSE171485原始数据维度：", nrow(fpkm_data), "个基因，", ncol(fpkm_data)-1, "个样本\n")

# 4. 处理重复基因（同一基因取平均）
## 处理raw counts
raw_matrix <- as.matrix(raw_counts[, -1, drop = FALSE])
rownames(raw_matrix) <- raw_counts$gene_name
raw_unique <- avereps(raw_matrix)  # 重复基因取平均
cat("GSE119794去重后基因数：", nrow(raw_unique), "\n")

## 处理FPKM
fpkm_matrix <- as.matrix(fpkm_data[, -1, drop = FALSE])
rownames(fpkm_matrix) <- fpkm_data$gene_name
fpkm_unique <- avereps(fpkm_matrix)  # 重复基因取平均
cat("GSE171485去重后基因数：", nrow(fpkm_unique), "\n")

# 5. raw counts转换为CPM
## 计算每个样本的总reads数
sample_totals <- colSums(raw_unique)
## 转换公式：CPM = (counts / 样本总reads数) * 1e6
raw_to_cpm <- function(counts_mat, totals) {
  t(t(counts_mat) / totals) * 1e6
}
gse119794_cpm <- raw_to_cpm(raw_unique, sample_totals)

# 6. FPKM转换为CPM（需要基因长度，通过biomaRt获取）
## 关键公式：CPM = FPKM × 基因长度（kb）
## 6.1 获取人类基因长度（使用Ensembl数据库）
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_symbols <- rownames(fpkm_unique)  # FPKM数据中的基因名

## 获取基因长度（转录本长度的平均值，单位：kb）
gene_lengths <- getBM(
  attributes = c("hgnc_symbol", "transcript_length"),
  filters = "hgnc_symbol",
  values = gene_symbols,
  mart = ensembl
)

## 计算每个基因的平均转录本长度（kb）
gene_lengths <- aggregate(
  transcript_length ~ hgnc_symbol,
  data = gene_lengths,
  FUN = mean
)
gene_lengths$transcript_length_kb <- gene_lengths$transcript_length / 1000  # 转换为kb
colnames(gene_lengths) <- c("gene_name", "transcript_length", "length_kb")

## 6.2 匹配基因长度与FPKM数据
### 筛选在FPKM数据中存在的基因
matched_lengths <- gene_lengths[gene_lengths$gene_name %in% rownames(fpkm_unique), ]
### 按FPKM数据的基因顺序排序
matched_lengths <- matched_lengths[match(rownames(fpkm_unique), matched_lengths$gene_name), ]
### 检查匹配情况
unmatched <- sum(is.na(matched_lengths$gene_name))
cat("GSE171485中无法匹配长度的基因数：", unmatched, "\n")

## 6.3 移除无法匹配长度的基因
fpkm_filtered <- fpkm_unique[!is.na(matched_lengths$gene_name), ]
length_filtered <- matched_lengths$length_kb[!is.na(matched_lengths$gene_name)]
cat("GSE171485匹配长度后基因数：", nrow(fpkm_filtered), "\n")

## 6.4 转换FPKM为CPM
gse171485_cpm <- fpkm_filtered * length_filtered  # 应用公式

# 7. 筛选两个数据集的共同基因
common_genes <- intersect(rownames(gse119794_cpm), rownames(gse171485_cpm))
cat("两个数据集的共同基因数：", length(common_genes), "\n")

# 8. 合并CPM数据并添加前缀
gse119794_final <- gse119794_cpm[common_genes, ]
colnames(gse119794_final) <- paste0("GSE119794_", colnames(gse119794_final))

gse171485_final <- gse171485_cpm[common_genes, ]
colnames(gse171485_final) <- paste0("GSE171485_", colnames(gse171485_final))

merged_cpm <- cbind(gse119794_final, gse171485_final)

# 9. 保存结果
## 保存合并后的CPM数据
write.csv(
  data.frame(gene_name = rownames(merged_cpm), merged_cpm, check.names = FALSE),
  file = file.path(output_dir, "merged_CPM.csv"),
  row.names = FALSE
)

## 分别保存两个数据集的CPM数据
write.csv(
  data.frame(gene_name = rownames(gse119794_cpm), gse119794_cpm, check.names = FALSE),
  file = file.path(output_dir, "GSE119794_CPM.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(gene_name = rownames(gse171485_cpm), gse171485_cpm, check.names = FALSE),
  file = file.path(output_dir, "GSE171485_CPM.csv"),
  row.names = FALSE
)

cat("转换完成！结果保存至：", output_dir, "\n")
cat("合并后CPM数据维度：", nrow(merged_cpm), "个基因，", ncol(merged_cpm), "个样本\n")
