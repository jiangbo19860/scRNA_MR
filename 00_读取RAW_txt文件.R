library(tidyverse)

# 1. 读取文件并查看基本信息
series_file <- "/Users/lijiangbo/scRNA_MR/1_data/pancreatic_GEO/GSE119794/GSE119794_series_matrix.txt"

# 检查文件存在性
if (!file.exists(series_file)) {
  stop("文件不存在，请确认路径正确！")
}

# 读取所有行（保留原始格式）
lines <- readLines(series_file, warn = FALSE)  # 避免换行符警告
cat("文件总共有", length(lines), "行\n\n")


# 2. 解析文件结构：分离不同类型的注释行
series_lines <- lines[str_detect(lines, "^!Series_")]    # 数据集整体信息
sample_lines <- lines[str_detect(lines, "^!Sample_")]    # 样本注释信息
platform_lines <- lines[str_detect(lines, "^!Platform_")]  # 测序平台信息

cat("文件结构分类：\n")
cat("- 系列信息行（!Series_）：", length(series_lines), "行\n")
cat("- 样本信息行（!Sample_）：", length(sample_lines), "行\n")
cat("- 平台信息行（!Platform_）：", length(platform_lines), "行\n\n")


# 3. 解析系列信息（!Series_开头）：修复管道操作语法
if (length(series_lines) > 0) {
  series_info <- tibble(raw = series_lines) %>%
    mutate(
      key = str_remove(raw, "^!Series_"),  # 提取键（如title、summary）
      # 修复值提取逻辑：先去除键部分，再处理剩余内容
      value = str_remove(raw, "^!Series_[^\\t]+\\t") %>%
        str_remove_all('"')  # 去除可能的引号
    ) %>%
    # 修复键名清理逻辑（不使用管道在select中直接操作）
    mutate(key = str_remove(key, "\\t.*")) %>%
    select(key, value)  # 选择最终的键和值

  cat("===== 数据集基本信息（前5条）：=====\n")
  print(head(series_info, 5))
  cat("\n")
}


# 4. 解析样本信息（!Sample_开头）：核心样本注释
if (length(sample_lines) > 0) {
  # 提取样本数量（第一个样本行的制表符数量=样本数）
  sample_count <- str_count(sample_lines[1], "\t")
  cat("样本总数：", sample_count, "个\n\n")

  # 解析为数据框（每行对应一个样本的所有属性）
  sample_tbl <- tibble(raw = sample_lines) %>%
    mutate(
      key = str_extract(raw, "(?<=^!Sample_)[^\\t]+"),  # 提取属性名
      # 提取值并处理引号和分割
      values = str_remove(raw, "^!Sample_[^\\t]+\\t") %>%
        str_remove_all('"') %>%  # 去除引号
        str_split("\\t")  # 按制表符分割为样本列表
    ) %>%
    select(-raw) %>%
    unnest(values) %>%  # 展开列表为行
    group_by(key) %>%
    mutate(sample_id = row_number()) %>%  # 样本编号（1到40）
    pivot_wider(names_from = key, values_from = values) %>%
    relocate(sample_id)  # 样本编号放第一列

  # 展示样本信息结构（前2个样本的前6个属性）
  cat("===== 样本信息结构（前2个样本，前6个属性）：=====\n")
  print(sample_tbl %>% slice(1:2) %>% select(1:6))

  # 列出所有样本属性
  cat("\n===== 样本包含的所有属性（共", ncol(sample_tbl)-1, "个）：=====\n")
  print(colnames(sample_tbl)[-1])  # 排除sample_id
}


# 5. 提取WGC编号并分组（解决之前的NA问题）
if (exists("sample_tbl") && nrow(sample_tbl) > 0) {
  # 从title中提取WGC编号并分组
  sample_groups <- sample_tbl %>%
    mutate(
      wgc_id = str_extract(title, "WGC\\d+"),  # 提取WGC编号
      # 根据source_name_ch1分组（Normal=正常，PC=肿瘤）
      group = case_when(
        source_name_ch1 == "Normal" ~ "normal",
        source_name_ch1 == "PC" ~ "tumor",
        TRUE ~ "unknown"
      )
    ) %>%
    filter(group != "unknown" & !is.na(wgc_id))  # 过滤有效分组和编号

  # 输出分组结果
  cat("\n===== 样本分组及WGC编号提取结果：=====\n")
  cat("正常组样本数：", sum(sample_groups$group == "normal"), "\n")
  cat("肿瘤组样本数：", sum(sample_groups$group == "tumor"), "\n")

  # 预览正常组和肿瘤组的WGC编号
  cat("\n正常组WGC编号：\n")
  print(sample_groups %>% filter(group == "normal") %>% pull(wgc_id))

  cat("\n肿瘤组WGC编号：\n")
  print(sample_groups %>% filter(group == "tumor") %>% pull(wgc_id))

  # 保存结果
  writeLines(
    sample_groups %>% filter(group == "normal") %>% pull(wgc_id),
    "normal_wgc.txt"
  )
  writeLines(
    sample_groups %>% filter(group == "tumor") %>% pull(wgc_id),
    "tumor_wgc.txt"
  )
  cat("\n分组WGC编号已分别保存至normal_wgc.txt和tumor_wgc.txt\n")
}
