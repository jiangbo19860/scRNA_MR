rm(list = ls())  # 清空工作环境
# 加载必要的包
library(data.table)
library(stringr)

GSE183795_series_matrix1 <- read.delim("~/scRNA_MR/1_data/GEO/GSE183795/GSE183795_series_matrix1.txt", header=FALSE)

# 提取series_matrix1中第一行和第十行的数据（注意R中索引从1开始）
row1 <- GSE183795_series_matrix1[1, ]  # 第一行
row10 <- GSE183795_series_matrix1[10, ]  # 第十行

# 转换为数据框并确保长度一致
# 提取有效列（排除可能存在的NA列）
valid_cols <- !is.na(row1) & !is.na(row10)
df <- data.frame(
  Row1_Info = as.character(row1[valid_cols]),
  Row10_Content = as.character(row10[valid_cols]),
  stringsAsFactors = FALSE
)

# 去掉df的第一行
df <- df[-1, , drop = FALSE]

# 去掉第二列中的"tissue: "部分
df$Row10_Content <- sub("tissue: ", "", df$Row10_Content)

unique(df$Row10_Content)  # 查看Row10_Content中的唯一值
counts <- table(df$Row10_Content)
counts

df$Group <- ifelse(df$Row10_Content %in% c("adjacent non-tumor", "Normal pancreas"),
                   "Control",
                   df$Row10_Content)

# 2. 统计新分组的数量
new_counts <- table(df$Group)

# 3. 转换为数据框并查看结果
new_counts_df <- as.data.frame(new_counts)
colnames(new_counts_df) <- c("Group", "Count")
print(new_counts_df)


# rename expression columns -----------------------------------------------
# 读取表达矩阵（环境中已有df）
GSE183795_normalized_matrix <- read.delim(
  "~/scRNA_MR/1_data/GEO/GSE183795/GSE183795_normalized_matrix.txt",
  row.names = 1,
  stringsAsFactors = FALSE
)

# 从环境中已有的df提取信息（第一列：样本ID；第三列：分组）
hussain_df <- df[grepl("^Hussain_", df[, 1]), ]  # 筛选Hussain样本
hussain_ids <- as.character(hussain_df[, 1])
hussain_groups <- as.character(hussain_df[, 3])

# ----------------------
# 1. 处理Hussain样本：按中间数字部分匹配
# ----------------------
# 提取df中Hussain样本的中间关键部分（xx_xxxx）
extract_middle_part <- function(id) {
  sub("^Hussain_(\\d+_\\d+)_.*$", "\\1", id)
}
df_middle_parts <- sapply(hussain_ids, extract_middle_part, USE.NAMES = FALSE)
middle_to_group <- data.frame(
  middle_part = df_middle_parts,
  group = hussain_groups,
  stringsAsFactors = FALSE
)

# ----------------------
# 2. 处理非Hussain样本：沿用原始逻辑
# ----------------------
process_non_hussain <- function(col) {
  col_processed <- gsub(pattern = "^X", replacement = "", col)
  col_processed <- gsub(pattern = "\\.CEL$", replacement = "", col_processed)
  col_processed <- gsub(pattern = "\\.", replacement = "-", col_processed)
  return(col_processed)
}

# 从df提取非Hussain样本的映射（第一列→第三列）
non_hussain_df <- df[!grepl("^Hussain_", df[, 1]), ]
non_hussain_map <- data.frame(
  original_id = as.character(non_hussain_df[, 1]),
  group = as.character(non_hussain_df[, 3]),
  stringsAsFactors = FALSE
)

# ----------------------
# 3. 批量处理所有列名并匹配分组
# ----------------------
expr_cols <- colnames(GSE183795_normalized_matrix)
new_colnames <- character(length(expr_cols))

for (i in seq_along(expr_cols)) {
  col <- expr_cols[i]

  if (grepl("^Hussain_", col)) {
    # 处理Hussain样本：匹配中间数字部分
    matrix_middle <- sub("^Hussain_(\\d+_\\d+)_.*$", "\\1", col)
    match_pos <- match(matrix_middle, middle_to_group$middle_part)
    new_colnames[i] <- ifelse(
      !is.na(match_pos),
      paste0(col, "_", middle_to_group$group[match_pos]),
      col
    )
  } else {
    # 处理非Hussain样本：沿用原始逻辑
    processed_col <- process_non_hussain(col)
    match_pos <- match(processed_col, non_hussain_map$original_id)
    new_colnames[i] <- ifelse(
      !is.na(match_pos),
      paste0(col, "_", non_hussain_map$group[match_pos]),
      col
    )
  }
}

# ----------------------
# 4. 查看匹配结果
# ----------------------
# 查看Hussain样本匹配示例
hussain_indices <- grepl("^Hussain_", expr_cols)
matched_hussain <- which(hussain_indices & new_colnames != expr_cols)
if (length(matched_hussain) > 0) {
  cat("Hussain样本匹配示例：\n")
  sample_matched <- sample(matched_hussain, min(3, length(matched_hussain)))
  for (i in sample_matched) {
    cat(sprintf("  %s -> %s\n", expr_cols[i], new_colnames[i]))
  }
}

# 查看非Hussain样本匹配示例
non_hussain_indices <- !grepl("^Hussain_", expr_cols)
matched_non_hussain <- which(non_hussain_indices & new_colnames != expr_cols)
if (length(matched_non_hussain) > 0) {
  cat("\n非Hussain样本匹配示例：\n")
  sample_matched <- sample(matched_non_hussain, min(3, length(matched_non_hussain)))
  for (i in sample_matched) {
    cat(sprintf("  %s -> %s\n", expr_cols[i], new_colnames[i]))
  }
}

# 查看未匹配列名
unmatched <- which(new_colnames == expr_cols)
if (length(unmatched) > 0) {
  cat("\n未匹配的列名共", length(unmatched), "个：\n")
  print(expr_cols[unmatched][1:min(5, length(unmatched))])
}

# ----------------------
# 5. 检查除前4列外的列名后缀及计数（新增功能）
# ----------------------
# 替换列名
colnames(GSE183795_normalized_matrix) <- new_colnames
final_cols <- colnames(GSE183795_normalized_matrix)

# 提取前4列之后的列
if (length(final_cols) > 3) {
  target_cols <- final_cols[4:length(final_cols)]

  # 检查是否以_Tumor或_Control结尾
  valid_suffix <- grepl("_(Tumor|Control)$", target_cols)
  invalid_cols <- target_cols[!valid_suffix]

  # 统计两组数量
  tumor_count <- sum(grepl("_Tumor$", target_cols))
  control_count <- sum(grepl("_Control$", target_cols))

  # 输出检查结果
  cat("\n===== 列名后缀检查结果 =====", "\n")
  cat("前3列之后的总列数：", length(target_cols), "\n")
  cat("以_Tumor结尾的列数：", tumor_count, "\n")
  cat("以_Control结尾的列数：", control_count, "\n")

  if (length(invalid_cols) > 0) {
    cat("后缀不符合要求的列名（共", length(invalid_cols), "个）：\n")
    print(head(invalid_cols, 10))  # 显示前10个
  } else {
    cat("所有列名后缀均符合要求（仅_Tumor或_Control）！\n")
  }
} else {
  cat("\n注意：总列数不足4列，无法进行后缀检查。\n")
}

# ----------------------
# 6. 保存结果
# ----------------------
output_file <- "~/scRNA_MR/1_data/GEO/GSE183795/GSE183795_normalized_matrix_final.txt"
write.table(GSE183795_normalized_matrix, output_file, sep = "\t", quote = FALSE)
cat("\n结果已保存至:", output_file, "\n")
