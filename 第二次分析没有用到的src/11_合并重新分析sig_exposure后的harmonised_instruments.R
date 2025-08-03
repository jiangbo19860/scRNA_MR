rm(list = ls())

# 加载所需包
pacman::p_load(tidyverse, here)

# 定义文件夹路径列表
folder_paths <- c(
  "/Users/lijiangbo/scRNA_MR/3_outputs/ebi-a-GCST90001526_ebi-a-GCST90018893",
  "/Users/lijiangbo/scRNA_MR/3_outputs/ebi-a-GCST90001528_ebi-a-GCST90018893",
  "/Users/lijiangbo/scRNA_MR/3_outputs/ebi-a-GCST90001532_ebi-a-GCST90018893",
  "/Users/lijiangbo/scRNA_MR/3_outputs/ebi-a-GCST90001571_ebi-a-GCST90018893",
  "/Users/lijiangbo/scRNA_MR/3_outputs/ebi-a-GCST90001670_ebi-a-GCST90018893",
  "/Users/lijiangbo/scRNA_MR/3_outputs/ebi-a-GCST90001671_ebi-a-GCST90018893",
  "/Users/lijiangbo/scRNA_MR/3_outputs/ebi-a-GCST90001672_ebi-a-GCST90018893",
  "/Users/lijiangbo/scRNA_MR/3_outputs/ebi-a-GCST90001704_ebi-a-GCST90018893",
  "/Users/lijiangbo/scRNA_MR/3_outputs/ebi-a-GCST90001716_ebi-a-GCST90018893",
  "/Users/lijiangbo/scRNA_MR/3_outputs/ebi-a-GCST90001719_ebi-a-GCST90018893",
  "/Users/lijiangbo/scRNA_MR/3_outputs/ebi-a-GCST90001723_ebi-a-GCST90018893",
  "/Users/lijiangbo/scRNA_MR/3_outputs/ebi-a-GCST90001736_ebi-a-GCST90018893",
  "/Users/lijiangbo/scRNA_MR/3_outputs/ebi-a-GCST90001752_ebi-a-GCST90018893",
  "/Users/lijiangbo/scRNA_MR/3_outputs/ebi-a-GCST90001759_ebi-a-GCST90018893",
  "/Users/lijiangbo/scRNA_MR/3_outputs/ebi-a-GCST90001762_ebi-a-GCST90018893",
  "/Users/lijiangbo/scRNA_MR/3_outputs/ebi-a-GCST90001768_ebi-a-GCST90018893",
  "/Users/lijiangbo/scRNA_MR/3_outputs/ebi-a-GCST90001808_ebi-a-GCST90018893",
  "/Users/lijiangbo/scRNA_MR/3_outputs/ebi-a-GCST90001823_ebi-a-GCST90018893",
  "/Users/lijiangbo/scRNA_MR/3_outputs/ebi-a-GCST90001827_ebi-a-GCST90018893",
  "/Users/lijiangbo/scRNA_MR/3_outputs/ebi-a-GCST90001829_ebi-a-GCST90018893",
  "/Users/lijiangbo/scRNA_MR/3_outputs/ebi-a-GCST90001834_ebi-a-GCST90018893",
  "/Users/lijiangbo/scRNA_MR/3_outputs/ebi-a-GCST90001894_ebi-a-GCST90018893",
  "/Users/lijiangbo/scRNA_MR/3_outputs/ebi-a-GCST90001898_ebi-a-GCST90018893",
  "/Users/lijiangbo/scRNA_MR/3_outputs/ebi-a-GCST90001987_ebi-a-GCST90018893",
  "/Users/lijiangbo/scRNA_MR/3_outputs/ebi-a-GCST90002116_ebi-a-GCST90018893"
)

# 统计文件信息和碱基列状态
file_stats <- tibble(
  folder_path = character(),
  file_exists = logical(),
  row_count = integer(),
  allele_cols_valid = logical()  # 标记碱基列是否存在且格式正常
)

# 存储所有清洗后的数据
all_data <- list()

# 定义需要处理的碱基列
target_cols <- c(
  "effect_allele.exposure",
  "other_allele.exposure",
  "effect_allele.outcome",
  "other_allele.outcome"
)

# 读取并清洗每个文件
for (folder in folder_paths) {
  file_path <- file.path(folder, "harmonised_data.csv")
  file_exists <- file.exists(file_path)
  row_count <- 0
  allele_cols_valid <- FALSE  # 标记当前文件的碱基列是否有效

  if (file_exists) {
    # 读取文件时强制将所有列转为字符型（避免因子类型干扰）
    df <- read_csv(
      file_path,
      show_col_types = FALSE,
      col_types = cols(.default = col_character())  # 关键：所有列先读为字符型
    )
    row_count <- nrow(df)

    # 检查目标列是否存在
    missing_cols <- setdiff(target_cols, colnames(df))
    if (length(missing_cols) == 0) {
      allele_cols_valid <- TRUE

      # 清洗碱基列：统一为大写字母，去除空格和特殊字符
      df <- df %>%
        mutate(
          # 处理暴露相关等位基因
          effect_allele.exposure = str_trim(effect_allele.exposure) %>% str_to_upper(),
          other_allele.exposure = str_trim(other_allele.exposure) %>% str_to_upper(),
          # 处理结局相关等位基因
          effect_allele.outcome = str_trim(effect_allele.outcome) %>% str_to_upper(),
          other_allele.outcome = str_trim(other_allele.outcome) %>% str_to_upper(),
          # 将空字符串转为NA（便于识别缺失）
          across(all_of(target_cols), ~ifelse(.x == "", NA, .x))
        )

      # 验证碱基列是否包含预期值（A/T/C/G/NA）
      valid_bases <- c("A", "T", "C", "G", NA)
      for (col in target_cols) {
        invalid_values <- setdiff(unique(df[[col]]), valid_bases)
        if (length(invalid_values) > 0) {
          cat("警告：", file_path, "的", col, "存在非碱基值：", paste(invalid_values, collapse = ","), "\n")
        }
      }

      cat("已读取并清洗：", file_path, "（行数：", row_count, "）\n")
    } else {
      cat("警告：", file_path, "缺失以下碱基列：", paste(missing_cols, collapse = ","), "\n")
    }

    # 存储清洗后的数据
    all_data[[folder]] <- df
  } else {
    cat("警告：文件不存在，已跳过：", file_path, "\n")
  }

  # 更新统计信息
  file_stats <- file_stats %>%
    add_row(
      folder_path = folder,
      file_exists = file_exists,
      row_count = row_count,
      allele_cols_valid = allele_cols_valid
    )
}

# 汇总统计信息
total_files <- nrow(file_stats)
existing_files <- sum(file_stats$file_exists)
valid_allele_files <- sum(file_stats$allele_cols_valid)
total_rows <- sum(file_stats$row_count)

cat("\n=== 汇总统计 ===\n")
cat("总文件夹数量：", total_files, "\n")
cat("存在harmonised_data.csv的文件夹数量：", existing_files, "\n")
cat("碱基列（4列）完整且格式正常的文件数量：", valid_allele_files, "\n")
cat("所有文件的总行数：", total_rows, "\n")

# 合并数据
combined_df <- bind_rows(all_data)
cat("数据合并完成，合并后实际行数：", nrow(combined_df), "\n")

# 拆分exposure列（如果存在）
if ("exposure" %in% colnames(combined_df)) {
  combined_df <- combined_df %>%
    separate(
      col = exposure,
      into = c("exposure_name", "exposure_id"),
      sep = " \\|\\| ",
      remove = FALSE,
      fill = "right"
    )
  cat("已将exposure列拆分为exposure_name和exposure_id\n")
} else {
  cat("警告：合并后的数据中未找到'exposure'列，跳过拆分步骤\n")
}

# 验证合并后的碱基列
if (all(target_cols %in% colnames(combined_df))) {
  cat("\n=== 合并后碱基列验证 ===\n")
  for (col in target_cols) {
    unique_vals <- unique(combined_df[[col]])
    cat(col, "的唯一值：", paste(unique_vals, collapse = ","), "\n")
  }
} else {
  cat("\n警告：合并后的数据缺失部分碱基列，无法验证\n")
}

# 保存结果
output_file <- here("3_outputs", "combined_harmonised_data.csv")
write_csv(combined_df, output_file)
cat("合并结果已保存至：", output_file, "\n")

# 保存统计信息
stats_file <- here("3_outputs", "harmonised_data_stats.csv")
write_csv(file_stats, stats_file)
cat("统计信息已保存至：", stats_file, "\n")
