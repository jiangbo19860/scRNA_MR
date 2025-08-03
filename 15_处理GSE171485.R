rm(list = ls())  # 清空工作空间
# 安装并加载必要的包
if (!require("pacman")) install.packages("pacman")
pacman::p_load(GEOquery, readr, dplyr, stringr)

# 1. 设置保存路径（根据实际需求修改）
save_dir <- "/Users/lijiangbo/scRNA_MR/1_data/pancreatic_GEO"
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

# 2. 下载补充文件中的表达矩阵（CSV.gz格式）
# 补充文件URL（从GEO页面获取）
expr_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE171nnn/GSE171485/suppl/GSE171485_PDAC-tissue-ExpressionMatrix.csv.gz"
expr_file <- file.path(save_dir, "GSE171485_PDAC-tissue-ExpressionMatrix.csv.gz")

# 下载文件
if (!file.exists(expr_file)) {
  download.file(
    url = expr_url,
    destfile = expr_file,
    mode = "wb"
  )
}

# 3. 读取表达矩阵
expr_matrix <- read_csv(expr_file, show_col_types = FALSE)
# 查看矩阵结构（确认第一列为基因/探针ID）
head(expr_matrix)
# 重命名第一列为"ID"（用于后续匹配基因名）
colnames(expr_matrix)[1] <- "ID"

colnames(expr_matrix)

# 整理表达矩阵：第一列为基因名，保留样本表达值
final_expr <- expr_matrix %>%
  select(
    gene_name = gene_short_name,  # 提取基因名作为第一列
    everything(),
    -ID, -tss_id, -locus  # 移除不需要的列（ID、tss_id、locus）
  )

# 处理重复基因（若存在，取平均值）
final_expr <- final_expr %>%
  group_by(gene_name) %>%
  summarise(across(everything(), mean), .groups = "drop")

# 定义样本对应关系数据框（原始列名、GSM编号、样本类型）
sample_mapping <- data.frame(
  original_name = c(paste0("PDAC-", 1:6), paste0("CT-", 1:6)),
  gsm_id = c(
    "GSM5226232", "GSM5226233", "GSM5226234",
    "GSM5530619", "GSM5530620", "GSM5530621",  # PDAC对应GSM
    "GSM5226229", "GSM5226230", "GSM5226231",
    "GSM5530616", "GSM5530617", "GSM5530618"   # CT对应GSM
  ),
  sample_type = c(rep("Tumor", 6), rep("Control", 6))  # 样本类型
)

# 生成新列名（GSM编号_原始名_类型）
sample_mapping$new_name <- with(
  sample_mapping,
  paste0(gsm_id, "_", original_name, "_", sample_type)
)

# 保留非样本列（ID、gene_short_name等），替换样本列名
non_sample_cols <- colnames(expr_matrix)[1:4]  # 前4列是非样本列
new_colnames <- c(non_sample_cols, sample_mapping$new_name)

# 应用新列名到表达矩阵
colnames(expr_matrix) <- new_colnames

# 查看替换后的列名
colnames(expr_matrix)

# 保存结果
output_file <- file.path(save_dir, "GSE171485_expression_with_gene_names.csv")
write_csv(final_expr, output_file)

message("处理完成！结果保存至：", output_file)
