# 安装必要的包
if (!require("pacman")) install.packages("pacman")
pacman::p_load(GEOquery, readxl, dplyr, stringr)

# 定义保存路径（根据实际需求修改）
save_dir <- "/Users/lijiangbo/scRNA_MR/1_data/pancreatic_GEO"
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

# 补充文件URL（从GEO页面获取）
supp_file_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE196nnn/GSE196009/suppl/GSE196009_PK_tissue.xlsx"

# 下载补充文件
download.file(
  url = supp_file_url,
  destfile = file.path(save_dir, "GSE196009_PK_tissue.xlsx"),
  mode = "wb"  # 二进制模式下载Excel文件
)
# 读取Excel文件（根据实际sheet名调整，默认第一个sheet）
expr_matrix <- read_excel(
  path = file.path(save_dir, "GSE196009_PK_tissue.xlsx"),
  sheet = 1,  # 若sheet名已知，可替换为sheet名称（如"Expression"）
  skip = 0    # 若有表头，无需跳过行；若有注释行，需调整skip值
)

# 查看矩阵结构（确认第一列是否为基因ID/探针ID）
head(expr_matrix)
colnames(expr_matrix)[1] <- "gene"

colnames(expr_matrix)

# 定义原始列名与GSM编号、样本类型的对应关系
sample_info <- data.frame(
  original_name = c("Bi-1T", "Bi-2T", "Bi-3N", "Bi-3T", "Bi-4T", "Bi-5T", "Bi-6T", "Bi-8N", "Bi-8T", "Bi-13N", "Bi-13T", "P-1N", "P-1T", "P-8T", "P-12N", "P-12T", "P-17T", "P-22N", "P-22T"),
  gsm_id = c("GSM5857882", "GSM5857883", "GSM5857884", "GSM5857885", "GSM5857886", "GSM5857887", "GSM5857888", "GSM5857889", "GSM5857890", "GSM5857891", "GSM5857892", "GSM5857893", "GSM5857894", "GSM5857895", "GSM5857896", "GSM5857897", "GSM5857898", "GSM5857899", "GSM5857900"),
  sample_type = ifelse(grepl("N$", c("Bi-1T", "Bi-2T", "Bi-3N", "Bi-3T", "Bi-4T", "Bi-5T", "Bi-6T", "Bi-8N", "Bi-8T", "Bi-13N", "Bi-13T", "P-1N", "P-1T", "P-8T", "P-12N", "P-12T", "P-17T", "P-22N", "P-22T")), "Control", "Cancer")
)

# 生成新列名（GSM编号_类型）
sample_info$new_name <- paste0(sample_info$gsm_id, "_", sample_info$sample_type)

# 替换表达矩阵的列名（保留第一列"gene"不变）
new_colnames <- c("gene", sample_info$new_name)

# 应用新列名到表达矩阵
colnames(expr_matrix) <- new_colnames

# 查看替换后的列名
colnames(expr_matrix)

# 提取除第一列（基因名）外的所有样本列名
sample_cols <- colnames(expr_matrix)[-1]

# 统计Cancer组数量（包含"_Cancer"的列）
cancer_count <- sum(grepl("_Cancer", sample_cols))

# 统计Control组数量（包含"_Control"的列）
control_count <- sum(grepl("_Control", sample_cols))

# 输出结果
cat("Cancer组样本数量：", cancer_count, "\n")
cat("Control组样本数量：", control_count, "\n")

colnames(expr_matrix) <- gsub("Cancer", "Tumor", colnames(expr_matrix))
colnames(expr_matrix)
# 保存处理后的表达矩阵
write.csv(
  expr_matrix,
  file = file.path(save_dir, "GSE196009_PK_tissue_processed.csv"),
  row.names = FALSE
)
