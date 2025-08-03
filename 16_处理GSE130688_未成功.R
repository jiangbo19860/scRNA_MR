# 清空工作空间
rm(list = ls())

# 1. 安装并加载必要的包
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  data.table,   # 高效数据处理
  stringr,      # 字符串处理
  biomaRt,      # 基因ID转换
  utils         # 解压文件
)

# 2. 设置路径与解压文件
tar_path <- "/Users/lijiangbo/scRNA_MR/1_data/pancreatic_GEO/GSE130688/GSE130688_RAW.tar"
extract_dir <- file.path(dirname(tar_path), "extracted_data")
dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)

# 强制解压文件（覆盖旧文件）
unlink(extract_dir, recursive = TRUE, force = TRUE)
untar(tar_path, exdir = extract_dir)

# 3. 读取样本文件列表并去重（核心：确保GSM编号唯一）
sample_files <- list.files(
  extract_dir,
  pattern = "*.txt",
  full.names = TRUE,
  recursive = TRUE
)

# 提取每个文件的GSM编号
sample_gsm <- stringr::str_extract(basename(sample_files), "GSM\\d+")

# 检查并去除重复的GSM样本（保留第一个出现的文件）
if (any(duplicated(sample_gsm))) {
  duplicated_gsm <- unique(sample_gsm[duplicated(sample_gsm)])
  warning("发现重复的GSM编号：", paste(duplicated_gsm, collapse = ", "))
  unique_indices <- !duplicated(sample_gsm)
  sample_files <- sample_files[unique_indices]
  sample_gsm <- sample_gsm[unique_indices]  # 更新GSM列表
  message("去重后保留样本数量：", length(sample_files))
} else {
  message("未发现重复样本，共", length(sample_files), "个样本")
}

# 4. 读取样本数据（过滤低表达+去重gene_id）
# 查看第一个样本结构
sample1 <- fread(sample_files[1], nrows = 5)
message("样本文件列名：", paste(colnames(sample1), collapse = ", "))

# 批量读取并处理样本
expr_list <- lapply(seq_along(sample_files), function(i) {
  file <- sample_files[i]
  gsm <- sample_gsm[i]  # 使用去重后的GSM编号

  # 读取单个样本（只保留需要的列）
  dt <- fread(file, select = c("gene_id", "expected_count"))

  # 过滤低表达基因（expected_count > 0）
  dt <- dt[expected_count > 0]

  # 去除重复gene_id（取表达量平均值）
  dt <- dt[, .(expected_count = mean(expected_count)), by = gene_id]

  # 重命名表达量列为GSM编号（确保唯一）
  setnames(dt, "expected_count", gsm)

  return(dt)
})

# 5. 合并样本为表达矩阵（内连接，无重复列名）
expr_matrix <- Reduce(function(x, y) merge(x, y, by = "gene_id", all = FALSE), expr_list)
setDT(expr_matrix)
message("合并后表达矩阵维度：", nrow(expr_matrix), "行（基因）×", ncol(expr_matrix), "列（样本+gene_id）")

# 6. 关联基因名称（XLOC→基因名）
get_gene_names <- function(xloc_list) {
  tryCatch({
    # 连接Ensembl数据库（人类）
    mart <- useMart(
      biomart = "ensembl",
      dataset = "hsapiens_gene_ensembl",
      host = "uswest.ensembl.org"
    )

    # 取所有XLOC进行匹配（当前基因数量较少，可全量匹配）
    xloc_subset <- xloc_list

    # 匹配基因名
    gene_map <- getBM(
      attributes = c("ensembl_transcript_id", "external_gene_name"),
      filters = "ensembl_transcript_id",
      values = xloc_subset,
      mart = mart
    )

    # 补充模糊匹配
    if (nrow(gene_map) < length(xloc_subset) * 0.3) {
      message("模糊匹配补充基因名...")
      xloc_fragments <- stringr::str_remove(xloc_subset, "XLOC_\\d+_?")
      gene_map_supp <- getBM(
        attributes = c("ensembl_transcript_id", "external_gene_name"),
        filters = "external_gene_name",
        values = xloc_fragments,
        mart = mart
      )
      gene_map <- rbind(gene_map, gene_map_supp) %>% unique(by = "ensembl_transcript_id")
    }

    # 整理结果
    setDT(gene_map)
    setnames(gene_map, c("ensembl_transcript_id", "external_gene_name"), c("gene_id", "gene_name"))
    return(gene_map[!is.na(gene_name) & gene_name != ""])

  }, error = function(e) {
    warning("数据库连接失败，用XLOC作为基因名。错误：", e$message)
    return(data.table(gene_id = xloc_list, gene_name = xloc_list))
  })
}

# 提取XLOC并获取基因名
xloc_list <- expr_matrix$gene_id
gene_annot <- get_gene_names(xloc_list)

# 7. 合并基因名到表达矩阵（第一列为gene_name）
final_matrix <- merge(expr_matrix, gene_annot, by = "gene_id", all.x = TRUE)
# 填充未匹配的基因名（用XLOC本身）
final_matrix[is.na(gene_name), gene_name := gene_id]
# 调整列顺序（基因名在前，删除原始gene_id）
setcolorder(final_matrix, c("gene_name", setdiff(colnames(final_matrix), c("gene_name", "gene_id"))))
final_matrix[, gene_id := NULL]

# 8. 保存结果
output_file <- file.path(dirname(tar_path), "GSE130688_gene_expression_matrix.csv")
fwrite(final_matrix, output_file, na = "")
message("最终矩阵已保存至：", output_file)

# 9. 查看示例（显示前5行和前6列）
message("\n最终矩阵前5行示例：")
print(head(final_matrix[, 1:6]))
