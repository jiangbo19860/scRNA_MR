# 将输入的sig_SNP数据与Ensembl数据库关联，获取对应的基因信息（基因 ID、基因符号等），并对结果进行筛选、统计和输出，最终得到sig_genes.txt。
rm(list = ls())
pacman::p_load(
  here,
  biomaRt,
  dplyr
)

# 1. 设置文件路径
input_file <- here("3_outputs/20250802_sig_SNPs/combined_sig_SNPs_84.csv")
output_dir <- here("3_outputs")

# 创建输出目录（如果不存在）
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 2. 读取数据并检查SNP列
df <- read.csv(input_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

# 检查是否存在"SNP"列
if (!"SNP" %in% colnames(df)) {
  stop("输入文件中未找到名为'SNP'的列，请检查列名是否正确")
}

# 检查SNP是否唯一
snp_duplicates <- sum(duplicated(df$SNP))
if (snp_duplicates > 0) {
  warning("发现", snp_duplicates, "个重复的SNP，将保留唯一值进行处理")
  df <- df %>% distinct(SNP, .keep_all = TRUE)
  cat("去重后保留", nrow(df), "个唯一SNP\n")
} else {
  cat("所有SNP均为唯一值，共", nrow(df), "个\n")
}

# 筛选以rs开头的SNP
rs_snps_df <- df %>%
  filter(grepl("^rs", SNP, ignore.case = FALSE))

rs_count <- nrow(rs_snps_df)
cat("共筛选出", rs_count, "个以rs开头的SNP\n")

if (rs_count == 0) {
  stop("未找到以rs开头的SNP，请检查数据格式")
}

# 提取唯一的rs开头SNP列表
unique_snps <- unique(rs_snps_df$SNP)

# 3. 连接Ensembl SNP数据库（获取ensembl_gene_id）
tryCatch({
  snp_mart <- useMart(
    biomart = "ENSEMBL_MART_SNP",
    dataset = "hsapiens_snp",
    host = "https://asia.ensembl.org"
  )

  # 查询SNP对应的ensembl_gene_id
  snp_annot <- getBM(
    attributes = c(
      "refsnp_id",                  # SNP的rs编号
      "ensembl_gene_stable_id"      # 基因的Ensembl ID
    ),
    filters = "snp_filter",
    values = unique_snps,
    mart = snp_mart
  )

  # 处理查询结果（去重，保留每个SNP的第一个基因ID）
  snp_annot <- snp_annot %>%
    distinct(refsnp_id, .keep_all = TRUE)

}, error = function(e) {
  stop("SNP数据库连接或查询失败：", e$message)
})

# 4. 连接Ensembl基因数据库（将ensembl_gene_id转换为标准基因符号）
tryCatch({
  gene_mart <- useMart(
    biomart = "ENSEMBL_MART_ENSEMBL",
    dataset = "hsapiens_gene_ensembl",
    host = "https://asia.ensembl.org"
  )

  # 提取需要转换的基因ID（排除空值和空字符串）
  valid_gene_ids <- snp_annot$ensembl_gene_stable_id[
    !is.na(snp_annot$ensembl_gene_stable_id) & snp_annot$ensembl_gene_stable_id != ""
  ]

  if (length(valid_gene_ids) > 0) {
    # 查询基因ID对应的标准符号（external_gene_name）
    gene_symbols <- getBM(
      attributes = c(
        "ensembl_gene_id",           # 基因ID
        "external_gene_name",         # 标准基因符号（如TP53）
        "hgnc_symbol"           # 基因名（HGNC 符号，即目标 gene_symbol）
      ),
      filters = "ensembl_gene_id",
      values = valid_gene_ids,
      mart = gene_mart
    ) %>%
      distinct(ensembl_gene_id, .keep_all = TRUE)  # 去重
  } else {
    gene_symbols <- data.frame(ensembl_gene_id = character(), external_gene_name = character(), hgnc_symbol = character())
  }

}, error = function(e) {
  stop("基因数据库连接或查询失败：", e$message)
})

# 5. 合并基因符号到SNP注释中
snp_annot_with_symbol <- snp_annot %>%
  left_join(gene_symbols, by = c("ensembl_gene_stable_id" = "ensembl_gene_id")) %>%
  rename(
    SNP = refsnp_id,
    ensembl_gene_id = ensembl_gene_stable_id,
    gene_symbol = external_gene_name,  # 标准基因符号
    hgnc_symbol = hgnc_symbol          # 保留hgnc_symbol
  )

# 6. 合并到原始数据中（修改部分，新增hgnc_symbol的位置）
result_df <- df %>%
  left_join(snp_annot_with_symbol, by = "SNP") %>%
  relocate(ensembl_gene_id, .after = SNP) %>%        # SNP后第1列：基因ID
  relocate(gene_symbol, .after = ensembl_gene_id) %>% # 第2列：external_gene_name
  relocate(hgnc_symbol, .after = gene_symbol)         # 第3列：hgnc_symbol

# 7. 修正匹配统计逻辑（gene_symbol为空即视为未匹配）
total <- nrow(result_df)

# 未匹配条件：gene_symbol是NA或空字符串（""）
unmatched <- sum(
  is.na(result_df$gene_symbol) |
    result_df$gene_symbol == "" |
    nchar(trimws(result_df$gene_symbol)) == 0  # 排除仅含空格的情况
)

matched <- total - unmatched  # 匹配数 = 总数 - 未匹配数

# 8. 准备输出文件名（包含日期、匹配数和SNP总行数）
today_date <- format(Sys.Date(), "%Y%m%d")
# 文件名格式：日期_snp_with_genes_matched_匹配数_总SNP数
output_basename <- paste0(today_date, "_snps_", total, "_genes_matched_", matched)
output_csv <- file.path(output_dir, paste0(output_basename, ".csv"))
output_txt <- file.path(output_dir, paste0(output_basename, ".txt"))

# 9. 保存结果
write.csv(result_df, output_csv, row.names = FALSE, quote = FALSE)
write.table(result_df, output_txt, row.names = FALSE, sep = "\t", quote = FALSE)

cat("CSV结果已保存至：", output_csv, "\n")
cat("TXT结果已保存至：", output_txt, "\n")

# 10. 输出修正后的匹配统计
cat("\n匹配统计：\n")
cat("总SNP数：", total, "\n")
cat("成功匹配基因符号的SNP数：", matched, "(", round(matched/total*100, 2), "%)\n")
cat("未匹配到基因符号的SNP数：", unmatched, "(", round(unmatched/total*100, 2), "%)\n")

# 11. 结果预览
cat("\n结果预览（前5行）：\n")
print(head(result_df %>% dplyr::select(SNP, ensembl_gene_id, gene_symbol), 5))

# 12. 比对gene_symbol和hgnc_symbol是否一致（修正版）
comparison_df <- result_df %>%
  # 第一步：统一处理空值（将""、纯空格转为NA）
  mutate(
    gene_symbol_clean = ifelse(
      is.na(gene_symbol) | gene_symbol == "" | nchar(trimws(gene_symbol)) == 0,
      NA,
      trimws(gene_symbol)
    ),
    hgnc_symbol_clean = ifelse(
      is.na(hgnc_symbol) | hgnc_symbol == "" | nchar(trimws(hgnc_symbol)) == 0,
      NA,
      trimws(hgnc_symbol)
    )
  ) %>%
  # 第二步：判断一致性（核心修正）
  mutate(
    is_consistent = case_when(
      # 两列均为NA（已统一处理空值后）→ 视为一致
      is.na(gene_symbol_clean) & is.na(hgnc_symbol_clean) ~ TRUE,
      # 两列均为非NA且值相同 → 视为一致
      !is.na(gene_symbol_clean) & !is.na(hgnc_symbol_clean) &
        gene_symbol_clean == hgnc_symbol_clean ~ TRUE,
      # 其他情况（一列为NA另一列非NA，或均非NA但值不同）→ 视为不一致
      TRUE ~ FALSE
    )
  )

# 提取并打印不一致的行
inconsistent_rows <- comparison_df %>%
  filter(!is_consistent) %>%
  dplyr::select(SNP, gene_symbol, hgnc_symbol, gene_symbol_clean, hgnc_symbol_clean)

# 输出比对结果
cat("\n=== gene_symbol与hgnc_symbol比对结果 ===\n")
if (nrow(inconsistent_rows) == 0) {
  cat("所有行的gene_symbol与hgnc_symbol一致（包括均为空值的情况）\n")
} else {
  cat("发现", nrow(inconsistent_rows), "行不一致，详情如下：\n")
  print(inconsistent_rows, n = nrow(inconsistent_rows), na.print = "NA")  # 解决NA打印问题
}

# 13. 筛选出gene_symbol不为NA且非空字符串的行
filtered_result <- result_df %>%
  filter(
    !is.na(gene_symbol),
    gene_symbol != "",
    nchar(trimws(gene_symbol)) > 0  # 排除仅含空格的情况
  )

# 14. 保存筛选后的结果为新的CSV和TXT文件
filtered_output_basename <- paste0(today_date, "_snp_gene_valid_matched_", nrow(filtered_result))
filtered_output_csv <- file.path(output_dir, paste0(filtered_output_basename, ".csv"))
filtered_output_txt <- file.path(output_dir, paste0(filtered_output_basename, ".txt"))

write.csv(filtered_result, filtered_output_csv, row.names = FALSE, quote = FALSE)
write.table(filtered_result, filtered_output_txt, row.names = FALSE, sep = "\t", quote = FALSE)

cat("\n筛选后（gene_symbol有效）的结果已保存至：\n")
cat("CSV：", filtered_output_csv, "\n")
cat("TXT：", filtered_output_txt, "\n")

# 15. 提取gene_symbol列的唯一基因名，保存到包含日期和基因数量的文件中
if (nrow(filtered_result) > 0) {
  # 提取唯一基因名（去重）并排序
  sig_genes <- filtered_result %>%
    pull(gene_symbol) %>%
    unique() %>%
    sort()  # 按字母顺序排序（可选）

  # 获取基因数量和当前日期
  gene_count <- length(sig_genes)
  today_date <- format(Sys.Date(), "%Y%m%d")

  # 构建包含日期和基因数量的文件名
  sig_genes_filename <- paste0(today_date, "_sig_genes_", gene_count, ".txt")
  sig_genes_path <- file.path(output_dir, sig_genes_filename)

  # 保存到文本文件，每行一个基因名
  writeLines(sig_genes, sig_genes_path)

  cat("\n已提取", gene_count, "个唯一基因名至：", sig_genes_path, "\n")
} else {
  cat("\n未找到有效的gene_symbol，无法生成基因文件\n")
}

# 去重前的所有基因名（含重复）
all_genes_before_unique <- filtered_result %>% pull(gene_symbol)
# 每个基因出现的次数
gene_counts <- table(all_genes_before_unique)

# 筛选出现次数 > 1的重复基因
duplicated_genes <- gene_counts[gene_counts > 1]

cat("\n=== 基因名重复情况（导致52→44的原因） ===\n")
cat("去重前总基因数（含重复）：", length(all_genes_before_unique), "\n")
cat("去重后唯一基因数：", length(sig_genes), "\n")
cat("重复基因总数：", sum(duplicated_genes) - length(duplicated_genes), "\n")  # 总重复次数 = 总和 - 基因个数
cat("重复基因明细（基因名：出现次数）：\n")
print(duplicated_genes)
