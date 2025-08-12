# 免疫细胞组成分析，使用CIBERSORT算法，数据读取、质量控制、预处理，核心算法运行和结果保存。注意：输入基因表达数据的质量（无重复/缺失基因名），以及正确适配 CIBERSORT 对输入格式的要求（文件路径而非内存数据）。
rm(list = ls())  # 清除工作空间
pacman::p_load(
  here,       # 用于处理文件路径
  data.table,  # 用于数据处理
  preprocessCore,  # 用于数据标准化
  e1071,  # 用于支持向量机等算法
  readr  # 用于读取CSV文件
)

# 确认工作路径
here()

# 1. 读取输入文件
inputFile_path = here("1_data/GEO/GSE183795/GSE183795_normalized_matrix_final.txt")
inputFile <- read.delim(
  inputFile_path,
  header = TRUE,
  sep = "\t",
  row.names = NULL,  # 暂不设置行名，避免重复报错
  check.names = FALSE,
  stringsAsFactors = FALSE
)
head(inputFile)

# 2. 数据有效性检查
# 检查是否有重复基因名（第一列）
gene_names <- inputFile[, 1]
cat("重复基因名数量：", sum(duplicated(gene_names)), "\n")

# 检查是否有缺失或空的基因名
cat("缺失基因名数量：", sum(is.na(gene_names)), "\n")
cat("空基因名数量：", sum(gene_names == ""), "\n")

# 检查是否有空行（所有表达值均为NA的行）
if (nrow(inputFile) > 0) {
  cat("全NA空行数量：", sum(apply(inputFile[, -1, drop = FALSE], 1, function(x) all(is.na(x)))), "\n")
} else {
  stop("输入文件读取后为空，请检查文件路径或文件内容！")
}

# 3. 处理可能的问题（去重和清除空行）
# 处理重复基因（取平均表达量）
pacman::p_load(limma)  # 加载limma包用于基因表达量平均
exp_matrix <- as.matrix(inputFile[, -1, drop = FALSE])  # 提取表达矩阵。[, -1]是排除第一列（基因名）。drop = FALSE是确保结果仍是矩阵/数据框。
unique_matrix <- avereps(exp_matrix, ID = gene_names)  # 对重复基因取平均

# 重新构建数据框并设置基因名为行名
inputFile <- as.data.frame(unique_matrix)
rownames(inputFile) <- rownames(unique_matrix)  # 用去重后的基因名作为行名

# 去除全NA的行
if (nrow(inputFile) > 0) {
  inputFile <- inputFile[!apply(inputFile, 1, function(x) all(is.na(x))), , drop = FALSE]
}

# 再次确认数据非空
if (nrow(inputFile) == 0) {
  stop("数据处理后为空，请检查原始数据质量！")
}

# 4. 保存处理后的文件为临时txt（关键修正：为CIBERSORT准备文件路径）
temp_mixture <- here("1_data/temp_mixture.txt")  # 临时文件路径
write.table(
  inputFile,
  file = temp_mixture,
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,  # 保留基因名作为行名
  col.names = TRUE   # 保留样本名作为列名
)

# 5. 加载CIBERSORT算法脚本
source("2_src/00_GEOimmune.CIBERSORT.R")

# 6. 运行CIBERSORT分析（注意：传递文件路径而非数据框）
outTab <- CIBERSORT(
  sig_matrix = here("1_data/ref.txt"),  # 参考基因表达矩阵路径
  mixture_file = temp_mixture,          # 传入临时文件路径（字符串格式）
  perm = 1000,
  QN = TRUE
)

# 删除临时文件（清理空间）
file.remove(temp_mixture)

# 7. 结果处理
# 筛选P-value小于1的有效样本
if (!is.null(outTab) && nrow(outTab) > 0) {
  outTab <- outTab[outTab[,"P-value"] < 1, , drop = FALSE]

  # 去除最后三列统计信息
  if (ncol(outTab) > 3) {
    outTab <- as.matrix(outTab[, 1:(ncol(outTab) - 3), drop = FALSE])
  }

  # 添加列名并保存结果（核心修改：文件名中加入当天日期）
  outTab <- rbind(id = colnames(outTab), outTab)

  # 获取当天日期（格式为"YYYY-MM-DD"），并拼接文件名
  today_date <- format(Sys.Date(), format = "%Y%m%d")
  output_filename <- paste0(today_date, "_CIBERSORT-Results", ".txt")  # 生成带日期的文件名
  output_path <- here("3_outputs", output_filename)  # 完整路径

  write.table(
    outTab,
    file = output_path,
    sep = "\t",
    quote = FALSE,
    col.names = FALSE
  )

  cat("分析完成！结果已保存至：", output_path, "\n")
} else {
  warning("CIBERSORT未返回有效结果，无法生成输出文件！")
}
