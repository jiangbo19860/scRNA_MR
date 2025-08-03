rm(list = ls())

# 加载所需包
pacman::p_load(tidyverse, here)
here()

# 定义文件夹路径列表（替换为你的实际路径）
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

# 存储每个文件的effect_allele.exposure列类型
type_records <- tibble(
  file_path = character(),
  col_type = character()
)

# 读取文件并记录列类型
for (folder in folder_paths) {
  file_path <- file.path(folder, "exposure_instruments.csv")

  if (file.exists(file_path)) {
    df <- read_csv(file_path, show_col_types = FALSE)

    # 记录列类型（若列存在）
    if ("effect_allele.exposure" %in% colnames(df)) {
      col_type <- class(df$effect_allele.exposure)[1]  # 获取数据类型
    } else {
      col_type <- "列不存在"
    }

    # 添加到记录
    type_records <- type_records %>%
      add_row(file_path = file_path, col_type = col_type)

  } else {
    type_records <- type_records %>%
      add_row(file_path = file_path, col_type = "文件不存在")
  }
}

# 打印所有文件的列类型
cat("=== effect_allele.exposure列类型汇总 ===\n")
print(type_records, n = nrow(type_records))

# 检测并打印类型不一致的文件
unique_types <- unique(type_records$col_type)
if (length(unique_types) > 1) {
  cat("\n=== 类型不一致的文件 ===\n")
  print(type_records, n = nrow(type_records))
} else {
  cat("\n所有文件的effect_allele.exposure列类型一致（均为", unique_types, "）\n")
}

# 统一类型后合并数据（可选）
all_data <- list()
for (folder in folder_paths) {
  file_path <- file.path(folder, "exposure_instruments.csv")
  if (file.exists(file_path)) {
    df <- read_csv(file_path, show_col_types = FALSE)
    # 统一转为字符型
    if ("effect_allele.exposure" %in% colnames(df)) {
      df <- df %>% mutate(effect_allele.exposure = as.character(effect_allele.exposure))
    }
    all_data[[folder]] <- df
  }
}
combined_df <- bind_rows(all_data)
cat("\n数据合并完成，总行数：", nrow(combined_df), "\n")

# 将exposure列按" || "拆分（注意前后空格）
if ("exposure" %in% colnames(combined_df)) {
  combined_df <- combined_df %>%
    separate(
      col = exposure,
      into = c("exposure_name", "exposure_id"),  # 拆分后的两列名称
      sep = " \\|\\| ",  # 匹配" || "（包含前后空格）
      remove = FALSE,    # 保留原始exposure列
      fill = "right"     # 若拆分失败，右侧列填充NA
    )
  cat("已将exposure列拆分为exposure_name和exposure_id_1\n")
} else {
  stop("合并后的数据中未找到'exposure'列，请检查文件格式")
}

# 保存结果
output_file <- here("3_outputs", "combined_exposure_instruments.csv")
write_csv(combined_df, output_file)
cat("结果已保存至：", output_file, "\n")

# 显示前几行验证
cat("\n前5行数据预览：\n")
print(head(combined_df %>% select(exposure, exposure_name, exposure_id, SNP)))
