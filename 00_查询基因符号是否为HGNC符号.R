rm(list = ls())  # 清空工作空间
library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")  # 人类基因数据库

# 提取表达数据中的基因名列表
current_genes <- rownames(filtered_matrix)

# 查询每个基因的HGNC符号及有效性
gene_info <- getBM(
  attributes = c("external_gene_name", "hgnc_symbol", "hgnc_id", "ensembl_gene_id"),
  filters = "external_gene_name",
  values = current_genes,
  mart = ensembl
)

# 输出统计结果：
cat("当前基因名的有效性验证结果：\n")
valid_hgnc <- sum(!is.na(gene_info$hgnc_symbol))
invalid_genes <- sum(is.na(gene_info$hgnc_symbol))
cat("- 有效HGNC符号数量：", valid_hgnc, "\n")
cat("- 无法匹配HGNC的基因数量：", invalid_genes, "\n")
if (invalid_genes > 0) {
  cat("示例无效基因名：", head(gene_info[is.na(gene_info$hgnc_symbol), "external_gene_name"]), "\n")
}


sig_genes <- read.table("/Users/lijiangbo/scRNA_MR/3_outputs/汇总有显著意义的MR结果/sig_genes.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
target_genes <- unique(as.vector(sig_genes[, 1]))  # 提取基因名并去重
# 若连接失败，尝试更换镜像（如useMart("ensembl", host = "uswest.ensembl.org")）
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# 3. 批量查询基因名对应的HGNC官方符号
# 注：external_gene_name包含别名，hgnc_symbol为官方符号
gene_info <- getBM(
  attributes = c("external_gene_name", "hgnc_symbol", "hgnc_id"),  # hgnc_id存在说明是官方符号
  filters = "external_gene_name",
  values = target_genes,
  mart = ensembl
)

# 4. 统计验证结果
# 有效HGNC符号：hgnc_symbol不为NA且与输入基因名一致
valid_hgnc <- gene_info[!is.na(gene_info$hgnc_symbol) & gene_info$external_gene_name == gene_info$hgnc_symbol, ]
# 别名或非官方符号：能匹配到HGNC符号但名称不一致
alias_genes <- gene_info[!is.na(gene_info$hgnc_symbol) & gene_info$external_gene_name != gene_info$hgnc_symbol, ]
# 无效符号：无法匹配到HGNC信息
invalid_genes <- setdiff(target_genes, gene_info$external_gene_name)
