rm(list = ls())

# 加载所需包
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  here,
  tidyverse,
  readr,
  dplyr
)

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

# 创建输出文件夹
output_dir <- here("3_outputs", "汇总有显著意义的MR结果")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  cat("已创建输出文件夹：", output_dir, "\n\n")
} else {
  cat("输出文件夹已存在：", output_dir, "\n\n")
}

# 初始化存储结果的数据框
sig_heterogeneity <- tibble()  # 存储显著的异质性结果
all_ivw_results <- tibble()    # 存储所有IVW方法结果
sig_pleiotropy <- tibble()     # 存储显著的多效性结果
file_stats <- tibble(          # 存储文件统计信息
  folder = character(),
  hetero_file_exists = logical(),
  hetero_rows = integer(),
  mr_file_exists = logical(),
  ivw_rows = integer(),
  pleio_file_exists = logical(),
  pleio_rows = integer()
)

# 循环处理每个文件夹
for (folder in folder_paths) {
  folder_name <- basename(folder)
  cat("处理文件夹：", folder_name, "\n")

  # 初始化当前文件夹的统计信息
  current_stats <- list(
    folder = folder_name,
    hetero_file_exists = FALSE,
    hetero_rows = 0,
    mr_file_exists = FALSE,
    ivw_rows = 0,
    pleio_file_exists = FALSE,
    pleio_rows = 0
  )

  # 1. 处理heterogeneity_test_all_methods.csv
  hetero_path <- file.path(folder, "heterogeneity_test_all_methods.csv")
  if (file.exists(hetero_path)) {
    current_stats$hetero_file_exists <- TRUE
    hetero_df <- read_csv(hetero_path, show_col_types = FALSE)
    current_stats$hetero_rows <- nrow(hetero_df)

    # 筛选Q_pval <= 0.05的行，添加文件夹标识
    if ("Q_pval" %in% colnames(hetero_df)) {
      sig_hetero <- hetero_df %>%
        filter(Q_pval <= 0.05) %>%
        mutate(source_folder = folder_name)
      if (nrow(sig_hetero) > 0) {
        sig_heterogeneity <- bind_rows(sig_heterogeneity, sig_hetero)
        cat("  发现", nrow(sig_hetero), "行显著的异质性结果（Q_pval ≤ 0.05）\n")
      } else {
        cat("  无显著的异质性结果（Q_pval > 0.05）\n")
      }
    } else {
      cat("  警告：", hetero_path, "中未找到Q_pval列，跳过分析\n")
    }
  } else {
    cat("  警告：异质性文件不存在，跳过：", hetero_path, "\n")
  }

  # 2. 处理mr_results.csv（提取IVW方法结果）
  mr_path <- file.path(folder, "mr_results.csv")
  if (file.exists(mr_path)) {
    current_stats$mr_file_exists <- TRUE
    mr_df <- read_csv(mr_path, show_col_types = FALSE)

    # 提取Inverse variance weighted方法的结果
    if ("method" %in% colnames(mr_df)) {
      ivw_df <- mr_df %>%
        filter(method == "Inverse variance weighted") %>%
        mutate(source_folder = folder_name)
      current_stats$ivw_rows <- nrow(ivw_df)
      all_ivw_results <- bind_rows(all_ivw_results, ivw_df)
      cat("  提取到", nrow(ivw_df), "行IVW方法结果\n")
    } else {
      cat("  警告：", mr_path, "中未找到method列，跳过IVW提取\n")
    }
  } else {
    cat("  警告：MR结果文件不存在，跳过：", mr_path, "\n")
  }

  # 3. 处理pleiotropy_test.csv
  pleio_path <- file.path(folder, "pleiotropy_test.csv")
  if (file.exists(pleio_path)) {
    current_stats$pleio_file_exists <- TRUE
    pleio_df <- read_csv(pleio_path, show_col_types = FALSE)
    current_stats$pleio_rows <- nrow(pleio_df)

    # 筛选pval <= 0.05的行，添加文件夹标识
    if ("pval" %in% colnames(pleio_df)) {
      sig_plei <- pleio_df %>%
        filter(pval <= 0.05) %>%
        mutate(source_folder = folder_name)
      if (nrow(sig_plei) > 0) {
        sig_pleiotropy <- bind_rows(sig_pleiotropy, sig_plei)
        cat("  发现", nrow(sig_plei), "行显著的多效性结果（pval ≤ 0.05）\n")
      } else {
        cat("  无显著的多效性结果（pval > 0.05）\n")
      }
    } else {
      cat("  警告：", pleio_path, "中未找到pval列，跳过分析\n")
    }
  } else {
    cat("  警告：多效性文件不存在，跳过：", pleio_path, "\n")
  }

  # 更新统计信息
  file_stats <- bind_rows(file_stats, as_tibble(current_stats))
  cat("-------------------------\n")
}

# 汇总统计信息
cat("\n=== 总体统计结果 ===\n")
cat("总文件夹数量：", nrow(file_stats), "\n")
cat("存在异质性文件的文件夹：", sum(file_stats$hetero_file_exists), "\n")
cat("存在MR结果文件的文件夹：", sum(file_stats$mr_file_exists), "\n")
cat("存在多效性文件的文件夹：", sum(file_stats$pleio_file_exists), "\n")
cat("提取到的IVW方法总行数：", nrow(all_ivw_results), "\n")

# 筛选IVW结果中pval <= 0.05的显著结果
if (nrow(all_ivw_results) > 0 && "pval" %in% colnames(all_ivw_results)) {
  sig_ivw <- all_ivw_results %>%
    filter(pval <= 0.05)
  cat("IVW方法中显著的结果行数：", nrow(sig_ivw), "\n")
} else {
  sig_ivw <- tibble()
  cat("未找到有效的IVW结果或pval列\n")
}

# 保存所有结果
# 1. 显著的异质性结果
if (nrow(sig_heterogeneity) > 0) {
  hetero_out <- file.path(output_dir, "显著异质性结果.csv")
  write_csv(sig_heterogeneity, hetero_out)
  cat("已保存显著异质性结果至：", hetero_out, "\n")
} else {
  cat("无显著的异质性结果，不保存\n")
}

# 2. 所有IVW结果及显著IVW结果
if (nrow(all_ivw_results) > 0) {
  ivw_all_out <- file.path(output_dir, "所有IVW方法结果.csv")
  write_csv(all_ivw_results, ivw_all_out)
  cat("已保存所有IVW方法结果至：", ivw_all_out, "\n")

  if (nrow(sig_ivw) > 0) {
    ivw_sig_out <- file.path(output_dir, "显著IVW结果.csv")
    write_csv(sig_ivw, ivw_sig_out)
    cat("已保存显著IVW结果至：", ivw_sig_out, "\n")
  } else {
    cat("无显著的IVW结果，不保存\n")
  }
} else {
  cat("无有效的IVW结果，不保存\n")
}

# 3. 显著的多效性结果
if (nrow(sig_pleiotropy) > 0) {
  pleio_out <- file.path(output_dir, "显著多效性结果.csv")
  write_csv(sig_pleiotropy, pleio_out)
  cat("已保存显著多效性结果至：", pleio_out, "\n")
} else {
  cat("无显著的多效性结果，不保存\n")
}

# 4. 保存文件统计信息
stats_out <- file.path(output_dir, "文件统计信息.csv")
write_csv(file_stats, stats_out)
cat("已保存文件统计信息至：", stats_out, "\n")

cat("\n所有分析完成！\n")
