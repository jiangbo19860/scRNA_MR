# 加载已安装的包
library(biomaRt)
library(dplyr)

# 1. 设置文件路径
input_file <- "/Users/lijiangbo/scRNA_MR/3_outputs/combined_harmonised_data.csv"
output_file <- "/Users/lijiangbo/scRNA_MR/3_outputs/snp_with_genes.csv"

# 2. 读取数据并筛选以rs开头的SNP
df <- read.csv(input_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

# 检查是否存在"SNP"列
if (!"SNP" %in% colnames(df)) {
  stop("输入文件中未找到名为'SNP'的列，请检查列名是否正确")
}

# 筛选以rs开头的SNP
df <- df %>%
  filter(grepl("^rs", SNP, ignore.case = FALSE))

# 去重SNP列表
unique_snps <- unique(df$SNP)
cat("共筛选出", length(unique_snps), "个以rs开头的唯一SNP\n")

if (length(unique_snps) == 0) {
  stop("未找到以rs开头的SNP，请检查数据格式")
}

# 3. 连接Ensembl数据库（使用HTTPS和亚洲镜像）
tryCatch({
  snp_mart <- useMart(
    biomart = "ENSEMBL_MART_SNP",
    dataset = "hsapiens_snp",
    host = "https://asia.ensembl.org"  # 关键修正：添加HTTPS并使用亚洲镜像
  )

  # 查询SNP对应的基因（使用之前筛选到的有效属性名）
  snp_annot <- getBM(
    attributes = c(
      "refsnp_id",                  # SNP的rs编号
      "ensembl_gene_stable_id",     # 基因的Ensembl ID
      "ensembl_gene_name"           # 修正：使用有效的基因名属性
    ),
    filters = "snp_filter",
    values = unique_snps,
    mart = snp_mart
  )

  # 处理查询结果
  snp_annot <- snp_annot %>%
    filter(!is.na(ensembl_gene_name)) %>%  # 过滤无基因名的记录
    distinct(refsnp_id, ensembl_gene_name, .keep_all = TRUE)  # 去重

  cat("成功匹配到", length(unique(snp_annot$refsnp_id)), "个SNP的基因信息\n")

}, error = function(e) {
  stop("数据库连接或查询失败：", e$message)
})

# 4. 合并基因信息到原始数据
colnames(snp_annot) <- c("SNP", "ensembl_gene_id", "gene_name")  # 重命名列名

result_df <- df %>%
  dplyr::left_join(snp_annot, by = "SNP") %>%
  dplyr::select(SNP, gene_name, ensembl_gene_id, everything())

# 5. 保存结果
write.csv(result_df, output_file, row.names = FALSE, quote = FALSE)
cat("结果已保存至：", output_file, "\n")

# 6. 匹配统计
cat("\n匹配统计：\n")
total <- nrow(df)
sum(is.na(result_df$gene_name))
matched <- sum(!is.na(result_df$gene_name))
cat("总SNP数：", total, "\n")
cat("成功匹配基因的SNP数：", matched, "(", round(matched/total*100, 2), "%)\n")
# 提取gene_name列中的空值样本（前5个）
empty_samples <- result_df$gene_name[result_df$gene_name == "" |
                                       result_df$gene_name == " " |  # 空格
                                       nchar(trimws(result_df$gene_name)) == 0]  # 去除空格后为空

head(empty_samples)  # 查看空值的实际存储形式
# 统计gene_name列中空字符串的行数
empty_count <- sum(result_df$gene_name == "")
cat("gene_name列中的空值（空字符串）数量：", empty_count, "\n")


# 提取未匹配的SNP
unmatched_snps <- result_df %>%
  filter(result_df$gene_name == "") %>%
  pull(SNP) %>%
  unique()

# 单独查询这些SNP，查看是否有隐藏注释
if (length(unmatched_snps) > 0) {
  snp_annot_unmatched <- getBM(
    attributes = c("refsnp_id", "ensembl_gene_stable_id", "ensembl_gene_name", "chrom_start", "chrom_end"),
    filters = "snp_filter",
    values = unmatched_snps,
    mart = snp_mart  # 之前连接的snp_mart
  )
  print(snp_annot_unmatched)  # 若仍无结果，则确认无注释
}


# 1. 从snp_annot_unmatched中提取有效基因注释（生成new_annot）
# （这一步是核心，之前可能漏跑了）
new_annot <- snp_annot_unmatched %>%
  dplyr::filter(!is.na(ensembl_gene_stable_id)) %>%  # 保留有基因ID的记录
  dplyr::select(
    refsnp_id,  # SNP编号（与result_df中的SNP列对应）
    ensembl_gene_stable_id,  # 基因的Ensembl ID
    ensembl_gene_name  # 基因名
  ) %>%
  distinct(refsnp_id, .keep_all = TRUE)  # 每个SNP只保留第一条有效注释

# 2. 重命名new_annot的列名，与result_df匹配
colnames(new_annot) <- c("SNP", "ensembl_gene_id_new", "gene_name_new")

# 3. 将新注释合并到result_df（而非原始df，因为result_df已有之前的匹配结果）
result_df_updated <- result_df %>%
  dplyr::left_join(new_annot, by = "SNP") %>%  # 通过SNP编号关联
  # 填充空值：如果原gene_name为空，则用新查询的gene_name_new补充
  mutate(
    gene_name = ifelse(
      gene_name == "" | is.na(gene_name),  # 原基因名为空或NA
      gene_name_new,                       # 用新基因名填充
      gene_name                            # 否则保留原有基因名
    ),
    # 同步更新ensembl_gene_id
    ensembl_gene_id = ifelse(
      gene_name == "" | is.na(gene_name),  # 原基因ID为空或NA
      ensembl_gene_id_new,                 # 用新基因ID填充
      ensembl_gene_id                      # 否则保留原有ID
    )
  ) %>%
  # 移除临时列（新基因名和新ID）
  dplyr::select(-gene_name_new, -ensembl_gene_id_new)

# 4. 查看更新效果
cat("更新前空值数量：", sum(result_df$gene_name == "" | is.na(result_df$gene_name)), "\n")
cat("更新后空值数量：", sum(result_df_updated$gene_name == "" | is.na(result_df_updated$gene_name)), "\n")

write.csv(result_df_updated, output_file, row.names = FALSE, quote = FALSE)
cat("结果已保存至：", output_file, "\n")
head(result_df_updated)

# 初始化基因相关的 BioMart（人类基因数据库示例）
ensembl_mart <- useMart(
  biomart = "ENSEMBL_MART_ENSEMBL",  # Ensembl 核心数据库
  dataset = "hsapiens_gene_ensembl"   # 人类基因数据集
)

# 从 result_df_updated 提取唯一的 ENSG 基因 ID
unique_ensg <- result_df_updated %>%
  distinct(ensembl_gene_id) %>%  # 避免重复查询
  pull(ensembl_gene_id) %>%
  na.omit()  # 忽略缺失值

# 查询符号
gene_annot <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),  # ENSG → 符号
  filters = "ensembl_gene_id",
  values = unique_ensg,
  mart = ensembl_mart
) %>%
  # 确保每个 ENSG 仅保留第一个匹配的符号（处理多映射）
  distinct(ensembl_gene_id, .keep_all = TRUE)

result_df_with_symbol <- result_df_updated %>%
  dplyr::left_join(gene_annot, by = c("ensembl_gene_id" = "ensembl_gene_id")) %>%  # 关联 ENSG
  # 替换或保留原有名称：若原 gene_name 是 ENSG 格式，则覆盖为符号；否则保留原值
  dplyr::mutate(
    gene_name = ifelse(
      grepl("^ENSG", gene_name),  # 判断是否为 ENSG 编号格式
      external_gene_name,          # 替换为符号
      gene_name                    # 保留非 ENSG 的原始名称（如手动注释）
    )
  ) %>%
  # 移除临时列
  dplyr::select(-external_gene_name) %>%
  filter(gene_name != "")  # 保留gene_name不为空的行

nrow(result_df_with_symbol) # 一共有85行数据

# 保存更新后的数据（包含符号）
write.csv(result_df_with_symbol, output_file, row.names = FALSE)

# 查看符号覆盖情况
table(grepl("^ENSG", result_df_with_symbol$gene_name))
# 输出应为 0，表示所有 ENSG 已被符号替换。

colnames(result_df_with_symbol)  # 查看列名，确认是否包含 gene_name 和 SNP)

# 提取未匹配的SNP
unmatched_snps <- result_df_with_symbol %>%
  filter(result_df_with_symbol$gene_name == "") %>%
  pull(SNP) %>%
  unique()
print(unmatched_snps)
length(unmatched_snps)
