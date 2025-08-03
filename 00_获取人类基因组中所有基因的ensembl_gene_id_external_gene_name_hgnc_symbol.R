## 连接到Ensembl数据库中获取人类基因的基本注释信息并保存为all_human_genes_basic.csv，48405行，5列。列名依次为ensembl_gene_id,	external_gene_name,hgnc_symbol,chromosome_name,start_position.

rm(list = ls())
pacman::p_load(
  here,
  biomaRt,
  data.table
)

# 连接Ensembl基因数据库
gene_mart <- useMart(
  biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "hsapiens_gene_ensembl",
  host = "https://asia.ensembl.org"
)

# 查看当前数据库中所有可用的属性（包含基因名称相关）
attributes <- listAttributes(gene_mart)
# 搜索与基因名称相关的属性
gene_name_attrs <- attributes[grep("gene_name|symbol", attributes$description, ignore.case = TRUE), ]
print(gene_name_attrs)
# 直接检查"external_gene_name"是否在属性名称中
attributes[attributes$name == "external_gene_name", ]

# 获取所有基因的基础注释（无SNP，数据量小）
all_genes <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "external_gene_name",
    "hgnc_symbol",
    "chromosome_name",  # 可选：染色体位置，用于后续筛选
    "start_position"    # 可选：基因起始位置
  ),
  mart = gene_mart
)

colnames(all_genes)
class(all_genes)
dim(all_genes)

# 筛选出external_gene_name不为NA且不为空字符串的行
filtered_genes <- all_genes %>%
  filter(
    !is.na(external_gene_name),  # 排除NA
    trimws(external_gene_name) != ""  # 排除空字符串（去除空格后）
  )



# 保存到本地（使用fread/fwrite更高效）
file_path <- here("all_human_genes_basic.csv")
fwrite(filtered_genes, file_path, row.names = FALSE)
cat("文件已保存至：", file_path, "\n")
# 验证文件是否存在
if (file.exists(file_path)) {
  cat("文件保存成功，文件大小：", file.size(file_path), "字节\n")
} else {
  cat("警告：文件保存失败！\n")
}
# 打印保存信息
cat("原始数据总行数：", nrow(all_genes), "\n")
cat("筛选后（external_gene_name非空）行数：", nrow(filtered_genes), "\n")
cat("文件已保存至：", file_path, "\n")
