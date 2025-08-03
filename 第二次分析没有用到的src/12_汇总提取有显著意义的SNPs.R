rm(list = ls())
# 加载所需包
library(tidyverse)
library(fs)

# 定义关键路径
summary_file_path <- "/Users/lijiangbo/scRNA_MR/3_outputs/汇总有显著意义的MR结果/显著IVW结果.csv"
output_path <- "/Users/lijiangbo/scRNA_MR/3_outputs/汇总有显著意义的MR结果/合并的harmonised_data.csv"

# 待筛选的文件夹路径列表（完整包含用户提供的所有路径）
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

# 读取汇总文件中的id.exposure列
summary_data <- read_csv(summary_file_path, show_col_types = FALSE)
target_ids <- unique(summary_data$id.exposure)

# 筛选出前缀匹配id.exposure的文件夹
filtered_folders <- folder_paths %>%
  enframe(name = NULL, value = "folder_path") %>%
  mutate(
    prefix = sub("_.*", "", basename(folder_path))  # 提取文件夹名中"_"前的部分
  ) %>%
  filter(prefix %in% target_ids) %>%
  pull(folder_path)

if (length(filtered_folders) == 0) {
  stop("未找到匹配的文件夹，请检查id.exposure与文件夹前缀是否一致")
}

# 读取并合并harmonised_data.csv，确保碱基列为大写
combined_data <- filtered_folders %>%
  map(function(path) {
    csv_file <- path_join(c(path, "harmonised_data.csv"))
    if (!file_exists(csv_file)) {
      warning(paste("跳过不存在的文件:", csv_file))
      return(NULL)
    }
    read_csv(csv_file, show_col_types = FALSE) %>%
      mutate(
        across(
          c(effect_allele.exposure, other_allele.exposure,
            effect_allele.outcome, other_allele.outcome),
          toupper  # 转换为大写字母确保格式一致
        )
      )
  }) %>%
  compact() %>%  # 移除空值
  bind_rows()

# 保存合并结果
write_csv(combined_data, output_path)
cat("合并完成，结果已保存至:", output_path, "\n")
