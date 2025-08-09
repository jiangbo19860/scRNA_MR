rm(list = ls())  # 清空工作空间
here()
# 1. 设置文件路径
file1 <- "/Users/lijiangbo/scRNA_MR/1_data/pancreatic_GEO/CPM_converted/GSE119794_CPM.csv"
file2 <- "/Users/lijiangbo/scRNA_MR/1_data/pancreatic_GEO/CPM_converted/GSE171485_CPM.csv"
file3 <- "/Users/lijiangbo/scRNA_MR/1_data/pancreatic_GEO/CPM_converted/GSE196009_CPM.csv"

# 2. 读取数据并提取基因名
read_genes <- function(file_path) {
  # 读取CSV文件，假设基因名列名为"gene_name"
  data <- read.csv(file_path, header = TRUE, check.names = FALSE)
  # 提取基因名并去重
  genes <- unique(as.character(data$gene_name))
  return(genes)
}

# 提取每个文件的基因列表
genes1 <- read_genes(file1)
genes2 <- read_genes(file2)
genes3 <- read_genes(file3)

# 3. 计算基因数量
cat("各文件基因总数：\n")
cat("GSE119794: ", length(genes1), "\n")
cat("GSE171485: ", length(genes2), "\n")
cat("GSE196009: ", length(genes3), "\n\n")

# 4. 两两比较基因重叠数量
## GSE119794与GSE171485
overlap_12 <- intersect(genes1, genes2)
cat("GSE119794与GSE171485的共同基因数：", length(overlap_12), "\n")  # 16244

## GSE119794与GSE196009
overlap_13 <- intersect(genes1, genes3)
cat("GSE119794与GSE196009的共同基因数：", length(overlap_13), "\n")  # 11086

## GSE171485与GSE196009
overlap_23 <- intersect(genes2, genes3)
cat("GSE171485与GSE196009的共同基因数：", length(overlap_23), "\n\n")  # 10066

# 5. 三个文件共同的基因数量
overlap_all <- intersect(intersect(genes1, genes2), genes3)
cat("三个文件共有的基因数：", length(overlap_all), "\n")  # 10032

# 6. 输出共同基因列表
## 保存两两共同基因
write.table(data.frame(gene_name = overlap_12),
            "overlap_GSE119794_GSE171485.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

write.table(data.frame(gene_name = overlap_13),
            "overlap_GSE119794_GSE196009.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

write.table(data.frame(gene_name = overlap_23),
            "overlap_GSE171485_GSE196009.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

## 保存三个文件共有的基因
output_file <- here("1_data/pancreatic_GEO/CPM_converted/overlap_all_three_datasets.txt")
write.table(data.frame(gene_name = overlap_all),
            output_file,
            sep = "\t", row.names = FALSE, quote = FALSE)

# 打印实际保存路径
current_dir <- getwd()
full_path <- file.path(current_dir, output_file)
cat("三个文件共有的基因列表已保存至：", full_path, "\n")

# 合并GSE119794和GSE171485 ---------------------------------------------------
# 7. 新增功能：按共有基因筛选并合并前两个数据集
## 7.1 读取完整的表达数据（包含表达值）
read_expression_data <- function(file_path) {
  data <- read.csv(file_path, header = TRUE, check.names = FALSE)
  # 确保基因名为字符型，避免匹配问题
  data$gene_name <- as.character(data$gene_name)
  return(data)
}

# 读取前两个数据集的完整表达数据
data1 <- read_expression_data(file1)
data2 <- read_expression_data(file2)

## 7.2 筛选出仅包含三个数据集共有的基因的行
data1_filtered <- data1[data1$gene_name %in% overlap_all, ]
data2_filtered <- data2[data2$gene_name %in% overlap_all, ]

# 检查筛选后的基因数量是否一致
cat("筛选后GSE119794的基因数：", nrow(data1_filtered), "\n")
cat("筛选后GSE171485的基因数：", nrow(data2_filtered), "\n\n")

## 7.3 按基因名排序，确保合并时顺序一致
data1_filtered <- data1_filtered[order(data1_filtered$gene_name), ]
data2_filtered <- data2_filtered[order(data2_filtered$gene_name), ]

## 7.4 合并两个数据集（按列合并，保留所有样本）
# 检查基因名是否完全匹配
if (all(data1_filtered$gene_name == data2_filtered$gene_name)) {
  merged_data <- cbind(data1_filtered, data2_filtered[, -1, drop = FALSE])
  cat("数据集合并成功！合并后的数据维度：", nrow(merged_data), "行，", ncol(merged_data), "列\n")
} else {
  # 若基因名不完全匹配，使用merge按基因名合并（更稳健）
  merged_data <- merge(data1_filtered, data2_filtered, by = "gene_name", all = FALSE)
  cat("通过基因名匹配合并数据，最终保留基因数：", nrow(merged_data), "\n")
}

## 7.5 保存合并后的数据集
merged_output <- here("1_data/pancreatic_GEO/CPM_converted/merged_GSE119794_GSE171485_common_genes.csv")
write.csv(merged_data, file = merged_output, row.names = FALSE, quote = FALSE)

# 打印合并文件的保存路径
full_merged_path <- file.path(current_dir, merged_output)
cat("合并后的数据集已保存至：", full_merged_path, "\n")

