# 加载必要的包
library(data.table)
library(stringr)

# 设置工作目录
setwd("/Users/lijiangbo/scRNA_MR/1_data/pancreatic_GEO/GSE119794")
raw_dir <- "GSE119794_RAW"

# 定义改进的读取函数（处理非数字表达量）
read_expression_data <- function(file_list) {
  data_list <- list()
  problem_files <- c()
  non_numeric_records <- list()  # 记录非数字的条目，便于排查

  for (file in file_list) {
    # 读取文件
    df <- tryCatch({
      fread(file, header = FALSE, sep = "\t",
            col.names = c("gene_id", "expression_str"),  # 先按字符串读取表达量
            stringsAsFactors = FALSE)
    }, error = function(e) {
      return(NULL)
    })

    if (is.null(df) || !is.data.table(df) || ncol(df) != 2) {
      problem_files <- c(problem_files, basename(file))
      next
    }

    # 提取样本名
    sample_name <- str_extract(basename(file), "WGC\\d+[RM]")

    # 关键改进：安全转换表达量为数值（允许小数，非数字值标记为NA）
    df[, expression := as.numeric(expression_str)]  # 用as.numeric替代as.integer，支持小数
    # 记录非数字的条目（便于后续检查）
    non_num <- df[is.na(expression) & expression_str != "", .(gene_id, invalid_value = expression_str)]
    if (nrow(non_num) > 0) {
      non_numeric_records[[sample_name]] <- non_num
    }

    # 保留有效列
    df[, expression_str := NULL]  # 删除原始字符串列
    setnames(df, "expression", sample_name)

    data_list[[sample_name]] <- df
  }

  # 输出非数字值的警告（可选：仅显示每个文件的问题数量）
  if (length(non_numeric_records) > 0) {
    cat("警告：部分表达量无法转换为数字（已标记为NA）：\n")
    for (sample in names(non_numeric_records)) {
      cat(sprintf("- %s：%d条无效记录\n", sample, nrow(non_numeric_records[[sample]])))
    }
  }

  if (length(problem_files) > 0) {
    cat("格式异常文件（已跳过）：\n", paste(problem_files, collapse = "\n"), "\n")
  }

  return(data_list)
}

# 分类文件（mRNA和miRNA）
txt_files <- list.files(path = raw_dir, pattern = "\\.txt$", full.names = TRUE)
mRNA_files <- txt_files[str_detect(basename(txt_files), "WGC\\d+R\\.txt$")]
miRNA_files <- txt_files[str_detect(basename(txt_files), "WGC\\d+M\\.txt$")]

# 处理mRNA数据
cat("===== 处理mRNA数据 =====", "\n")
mRNA_list <- read_expression_data(mRNA_files)
mRNA_merged <- Reduce(function(x, y) merge(x, y, by = "gene_id", all = TRUE), mRNA_list)
mRNA_merged[is.na(mRNA_merged)] <- 0  # NA替换为0

# 处理miRNA数据
cat("\n===== 处理miRNA数据 =====", "\n")
miRNA_list <- read_expression_data(miRNA_files)
miRNA_merged <- Reduce(function(x, y) merge(x, y, by = "gene_id", all = TRUE), miRNA_list)
miRNA_merged[is.na(miRNA_merged)] <- 0  # NA替换为0

# 保存结果
fwrite(mRNA_merged, "GSE119794_mRNA_samples.csv", row.names = FALSE)
fwrite(miRNA_merged, "GSE119794_miRNA_samples.csv", row.names = FALSE)
cat("\n结果已保存，NA值已替换为0\n")

colnames(mRNA_merged)

# 1. 定义正常组和肿瘤组的基础WGC编号（去重，避免重复匹配）
wgc_normal <- unique(c("WGC022997", "WGC024529", "WGC024539", "WGC026140", "WGC026144",
                       "WGC026146", "WGC027958", "WGC027960", "WGC027966", "WGC027968"))

wgc_tumor <- unique(c("WGC022998", "WGC024528", "WGC024538", "WGC026141", "WGC026145",
                      "WGC026147", "WGC027959", "WGC027961", "WGC027967", "WGC027969"))

# 2. 提取mRNA_merged的样本列名（排除gene_id）
sample_cols <- setdiff(colnames(mRNA_merged), "gene_id")

# 3. 匹配编号并添加分组后缀
new_colnames <- sapply(sample_cols, function(col) {
  # 提取基础编号（去除末尾的"R"）
  base_wgc <- str_remove(col, "R$")

  # 判断分组并添加后缀
  if (base_wgc %in% wgc_normal) {
    return(paste0(col, "_Control"))  # 正常组添加_Control
  } else if (base_wgc %in% wgc_tumor) {
    return(paste0(col, "_Tumor"))    # 肿瘤组添加_Tumor
  } else {
    return(col)  # 未匹配的列名保持不变（理论上不存在）
  }
})

# 4. 更新mRNA_merged的列名
colnames(mRNA_merged) <- c("gene_id", new_colnames)

# 5. 查看结果
cat("更新后的列名：\n")
print(colnames(mRNA_merged))

write.csv(mRNA_merged, "GSE119794_mRNA_samples_updated.csv", row.names = FALSE)
