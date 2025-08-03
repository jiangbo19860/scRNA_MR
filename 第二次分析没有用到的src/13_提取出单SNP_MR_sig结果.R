rm(list = ls())
# 加载所需包
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, fs)

here()

# 定义主目录路径
main_dir <- here("3_outputs/25个sig_exposures/")
main_dir <- normalizePath(main_dir, winslash = "/", mustWork = FALSE)

# 检查主目录是否存在
if (!dir_exists(main_dir)) {
  stop("主目录不存在，请检查路径：", main_dir)
}

# 查找所有目标文件
input_files <- dir_ls(
  path = main_dir,
  recurse = TRUE,
  type = "file",
  regexp = "single_snp_results_default\\.csv$"
)

# 统计总文件夹数量
total_folders <- length(dir_ls(main_dir, type = "directory"))
message("总子文件夹数量：", total_folders)

if (length(input_files) == 0) {
  stop("未找到任何single_snp_results_default.csv文件")
}

# 初始化变量
success_folders <- 0
combined_data <- tibble()  # 用于存储合并后的所有符合条件的行（包含所有列）

# 处理每个文件
for (file in sort(input_files)) {
  # 提取暴露和结局ID（用于追溯来源）
  folder_name <- str_remove(file, "/single_snp_results_default.csv$") %>% basename()
  ids <- str_split(folder_name, "_")[[1]]
  exposure_id <- ids[1]
  outcome_id <- ids[2]

  tryCatch({
    # 读取完整数据（保留所有列）
    df <- read_csv(file, show_col_types = FALSE)

    # 检查必要列是否存在
    if (!all(c("SNP", "p") %in% colnames(df))) {
      message("⚠️ 跳过", file, "：缺少SNP或p列")
      next
    }

    # 筛选符合条件的行（SNP以rs开头且p<0.05），保留所有原始列
    filtered_rows <- df %>%
      filter(str_detect(SNP, "^rs"), p < 0.05) %>%
      # 添加来源标识列（不影响原始数据列）
      mutate(
        exposure_id = exposure_id,
        outcome_id = outcome_id,
        source_file = file
      )

    # 合并数据（按行追加）
    if (nrow(filtered_rows) > 0) {
      success_folders <- success_folders + 1
      combined_data <- bind_rows(combined_data, filtered_rows)
      message("✅ 处理", file, "：找到", nrow(filtered_rows), "行符合条件数据")
    } else {
      message("ℹ️ ", file, "：未找到符合条件的行")
    }

  }, error = function(e) {
    message("❌ 错误处理", file, "：", e$message)
  })
}

# 输出合并结果
if (nrow(combined_data) > 0) {
  # 保存包含所有列的合并数据
  output_file <- file.path(main_dir, "combined_单个SNP_MR分析结果有统计学差异的SNPs.csv")
  write_csv(combined_data, output_file)

  # 输出统计信息
  message("\n===== 合并结果 =====")
  message("总符合条件的行数：", nrow(combined_data))
  message("保留的列数（含原始列+来源标识）：", ncol(combined_data))
  message("合并文件保存至：", output_file)
} else {
  message("\n⚠️ 未找到任何符合条件的SNP行")
}

# 输出文件夹统计
message("\n===== 文件夹统计 =====")
message("总子文件夹数量：", total_folders)
message("含符合条件数据的文件夹数量：", success_folders)
