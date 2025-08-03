# 提取、复制并合并显著关联的SNP文件（含替换TRUE为T及文件名含唯一值数量）
rm(list = ls())  # 清空工作空间

# 加载所需包
if (!require("pacman")) {
  install.packages("pacman")
}
pacman::p_load(here, fs, readr, dplyr, stringr, purrr)
here()
# 1. 设置路径和创建带序号的输出文件夹（若已存在则递增序号）
today_date <- format(Sys.Date(), "%Y%m%d")
base_output_name <- paste0(today_date, "_sig_SNPs")
output_dir <- here("3_outputs", base_output_name)

# 若文件夹已存在，添加序号（如20250802_sig_SNPs_1）
counter <- 1
while (dir_exists(output_dir)) {
  output_dir <- here("3_outputs", paste0(base_output_name, "_", counter))
  counter <- counter + 1
}
dir_create(output_dir)
message("输出文件夹: ", output_dir)

# 定义源文件夹路径
source_dir <- here("3_outputs/20250801_MR")
message("源文件夹路径: ", source_dir)

# 2. 读取显著结果文件并提取id.exposure
sig_results <- read_csv(
  here("3_outputs/20250802_MR_forest/significant_results_p0.05.csv"),
  show_col_types = FALSE
)

# 提取唯一的id.exposure值
target_ids <- sig_results %>%
  pull(id.exposure) %>%
  str_trim() %>%
  unique() %>%
  na.omit()

if (length(target_ids) == 0) {
  stop("在significant_results_p0.05.csv中未找到有效的id.exposure值")
}
message("提取到的id.exposure数量: ", length(target_ids))
message("部分id示例: ", paste(head(target_ids), collapse = ", "))

# 3. 获取源文件夹中所有CSV文件并匹配目标ID
all_csv_files <- dir_ls(
  path = source_dir,
  type = "file",
  regexp = "\\.csv$"
)

if (length(all_csv_files) == 0) {
  stop("在源文件夹中未找到任何CSV文件")
}
message("源文件夹中的CSV文件数量: ", length(all_csv_files))

# 提取文件名用于匹配
file_names <- basename(all_csv_files)

# 构建匹配模式
pattern <- paste0("^table\\.SNP_.*(", paste(target_ids, collapse = "|"), ")")

# 筛选包含目标ID的文件
matched_idx <- str_detect(file_names, pattern)
matched_files <- all_csv_files[matched_idx]

if (length(matched_files) == 0) {
  stop("未找到与id.exposure匹配的文件，请检查ID是否正确")
}

# 4. 复制匹配的文件到输出文件夹
file_copy(
  path = matched_files,
  new_path = output_dir,
  overwrite = TRUE
)

message("成功复制 ", length(matched_files), " 个文件到输出文件夹")
message("匹配的文件名示例: ", paste(head(basename(matched_files)), collapse = ", "))

# 5. 合并所有复制后的CSV文件（处理数据类型冲突）
output_csv_files <- dir_ls(
  path = output_dir,
  type = "file",
  regexp = "\\.csv$"
)

# 读取并合并所有CSV文件，统一处理数据类型
combined_data <- map_dfr(output_csv_files, function(file) {
  # 读取单个文件
  df <- read_csv(file, show_col_types = FALSE)

  # 处理可能的类型冲突列（提前转换）
  if ("effect_allele.exposure" %in% colnames(df)) {
    df <- df %>% mutate(effect_allele.exposure = as.character(effect_allele.exposure))
  }
  if ("effect_allele.outcome" %in% colnames(df)) {
    df <- df %>% mutate(effect_allele.outcome = as.character(effect_allele.outcome))
  }
  if ("SNP" %in% colnames(df)) {
    df <- df %>% mutate(SNP = as.character(SNP))
  }
  if ("other_allele.exposure" %in% colnames(df)) {
    df <- df %>% mutate(other_allele.exposure = as.character(other_allele.exposure))
  }
  if ("other_allele.outcome" %in% colnames(df)) {
    df <- df %>% mutate(other_allele.outcome = as.character(other_allele.outcome))
  }

  # 从文件名提取id.exposure并添加来源列
  id_exposure <- str_extract(basename(file), "ebi-a-GCST9000\\d+")
  df <- df %>% mutate(
    id.exposure = id_exposure,
    source_file = basename(file)
  )

  return(df)
})

# 6. 根据SNP和id.exposure组合去重
if (all(c("SNP", "id.exposure") %in% colnames(combined_data))) {
  combined_data_unique <- combined_data %>%
    distinct(SNP, id.exposure, .keep_all = TRUE)

  message("去重前数据行数: ", nrow(combined_data))
  message("去重后数据行数: ", nrow(combined_data_unique))
} else {
  stop("数据中缺少'SNP'列或'id.exposure'列，无法按组合去重！")
}

# 新增：处理第2-4列中的TRUE为T，再强制转换前5列为字符型
if (ncol(combined_data_unique) >= 5) {
  # 第一步：将第2-4列中的TRUE（字符型或逻辑型）替换为T
  combined_data_unique <- combined_data_unique %>%
    mutate(across(2:4, ~ ifelse(.x == TRUE | .x == "TRUE", "T", as.character(.x))))

  # 第二步：强制前5列为字符型（确保所有元素均为字符型）
  combined_data_unique <- combined_data_unique %>%
    mutate(across(1:5, as.character))

  message("已将第2-4列的TRUE替换为T，并将前5列强制转换为字符型")
} else {
  warning("数据列数不足5列，无法执行第2-4列替换及前5列转换")
}

# 定义target_cols为前5列的列名
if (ncol(combined_data_unique) >= 5) {
  target_cols <- colnames(combined_data_unique)[1:5]
  message("用于统计的前5列列名: ", paste(target_cols, collapse = ", "))
} else {
  target_cols <- colnames(combined_data_unique)
  warning("数据列数不足5列，将使用所有列进行统计")
}

# 统计前5列组合的唯一行数
unique_combination_count <- if (length(target_cols) > 0) {
  unique_combinations <- combined_data_unique %>%
    distinct(!!!syms(target_cols))  # 计算唯一组合数
  nrow(unique_combinations)
} else {
  NA
}

if (!is.na(unique_combination_count)) {
  message("前5列组合的唯一行数: ", unique_combination_count)
} else {
  warning("无法计算前5列组合的唯一行数")
}

# 7. 保存文件（文件名包含前5列唯一值数量，不追加统计行）
# 处理唯一值数量为NA的情况（用"unknown"代替）
suffix <- ifelse(is.na(unique_combination_count), "unknown", unique_combination_count)
combined_file <- here(output_dir, paste0("combined_sig_SNPs_", suffix, ".csv"))

# 仅保存去重后的数据（不添加统计行）
write_csv(combined_data_unique, combined_file, append = FALSE)

message("成功保存文件，前5列唯一值数量: ", suffix)
message("合并后文件路径: ", combined_file)

# 筛选pleiotropy和heterogeneity相关文件 ------------------------------------------
# 提取id.exposure的唯一值（用于后续匹配）
id_exposure_unique <- combined_data_unique %>%
  pull(id.exposure) %>%
  unique() %>%
  na.omit()
message("待筛选的id.exposure唯一值数量: ", length(id_exposure_unique))


# 一、处理pleiotropy文件（pval < 0.05）
pleio_target_dir <- here("3_outputs", paste0(today_date, "_存在pleiotropy的exposure"))
dir_create(pleio_target_dir)
message("pleiotropy文件目标文件夹: ", pleio_target_dir)

pleio_count <- 0  # 统计符合条件的文件数

for (id in id_exposure_unique) {
  pleio_pattern <- paste0(".*", id, ".*pleiotropy.*\\.csv$")
  pleio_files <- dir_ls(
    path = source_dir,
    type = "file",
    regexp = pleio_pattern,
    ignore.case = TRUE
  )

  if (length(pleio_files) == 0) next

  for (file in pleio_files) {
    df <- read_csv(file, show_col_types = FALSE)
    if (!"pval" %in% colnames(df)) next

    pvals <- df$pval
    if (any(is.na(pvals))) next
    if (any(pvals < 0.05)) {
      file_copy(file, pleio_target_dir, overwrite = TRUE)
      pleio_count <- pleio_count + 1
    }
  }
}
message("符合条件的pleiotropy文件数量: ", pleio_count)


# 二、处理heterogeneity文件（Q_pval < 0.05）
hetero_target_dir <- here("3_outputs", paste0(today_date, "_存在heterogeneity的exposure"))
dir_create(hetero_target_dir)
message("heterogeneity文件目标文件夹: ", hetero_target_dir)

hetero_count <- 0  # 统计符合条件的文件数

for (id in id_exposure_unique) {
  hetero_pattern <- paste0(".*", id, ".*heterogeneity.*\\.csv$")
  hetero_files <- dir_ls(
    path = source_dir,
    type = "file",
    regexp = hetero_pattern,
    ignore.case = TRUE
  )

  if (length(hetero_files) == 0) next

  for (file in hetero_files) {
    df <- read_csv(file, show_col_types = FALSE)
    if (!"Q_pval" %in% colnames(df)) next

    q_pvals <- df$Q_pval
    if (any(is.na(q_pvals))) next
    if (any(q_pvals < 0.05)) {
      file_copy(file, hetero_target_dir, overwrite = TRUE)
      hetero_count <- hetero_count + 1
    }
  }
}
message("符合条件的heterogeneity文件数量: ", hetero_count)


message("所有操作完成！")
