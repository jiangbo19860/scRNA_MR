# 加载必需包（若未安装，先自动安装）
if (!require("readxl")) install.packages("readxl", quiet = TRUE)
if (!require("dplyr")) install.packages("dplyr", quiet = TRUE)
if (!require("purrr")) install.packages("purrr", quiet = TRUE)
if (!require("tidyr")) install.packages("tidyr", quiet = TRUE)

library(readxl)
library(dplyr)
library(purrr)
library(tidyr)

# ===================== 1. 定义文件路径 =====================
# Excel文件路径（含待比对基因）
excel_path <- "/Users/lijiangbo/1_Projects/Epi-eQTL-MR/results_old/1004coloc/coloc阳性_103行_54genes.xlsx"

# 8个细胞群DEGs CSV文件路径（按顺序排列）
csv_paths <- c(
  "/Users/lijiangbo/1_Projects/Epi1106/5_scRNA/1113重新分析/1113outputs/1115/1116_Astrocytes_DEGs.csv",
  "/Users/lijiangbo/1_Projects/Epi1106/5_scRNA/1113重新分析/1113outputs/1115/1116_Endothelial_Cells_DEGs.csv",
  "/Users/lijiangbo/1_Projects/Epi1106/5_scRNA/1113重新分析/1113outputs/1115/1116_Excitatory_Neurons_DEGs.csv",
  "/Users/lijiangbo/1_Projects/Epi1106/5_scRNA/1113重新分析/1113outputs/1115/1116_Inhibitory_Neurons_DEGs.csv",
  "/Users/lijiangbo/1_Projects/Epi1106/5_scRNA/1113重新分析/1113outputs/1115/1116_Microglia_DEGs.csv",
  "/Users/lijiangbo/1_Projects/Epi1106/5_scRNA/1113重新分析/1113outputs/1115/1116_Oligodendrocytes_DEGs.csv",
  "/Users/lijiangbo/1_Projects/Epi1106/5_scRNA/1113重新分析/1113outputs/1115/1116_OPCs_DEGs.csv",
  "/Users/lijiangbo/1_Projects/Epi1106/5_scRNA/1113重新分析/1113outputs/1115/1116_Pericytes_DEGs.csv"
)

# 输出报告路径
output_path <- "/Users/lijiangbo/1_Projects/Epi1106/5_scRNA/1113重新分析/1113outputs/1115/基因比对详细报告.csv"

# ===================== 2. 读取并整理待比对基因（Excel文件） =====================
# 读取Excel中的geneSymble列，去重、统一格式（大写，避免大小写差异）
excel_genes <- read_excel(excel_path) %>%
  select(geneSymble) %>%  # 选择目标列
  filter(!is.na(geneSymble)) %>%  # 去除NA值
  distinct() %>%  # 去重
  mutate(geneSymble = toupper(geneSymble)) %>%  # 统一转为大写
  pull(geneSymble)  # 转为向量

# 统计Excel中有效基因数
excel_gene_count <- length(excel_genes)
message(paste0("📊 Excel文件中有效基因数：", excel_gene_count, "个"))

# ===================== 3. 批量读取并整理DEGs文件（CSV） =====================
# 定义函数：读取单个CSV，提取gene列并整理
read_deg_genes <- function(csv_path) {
  # 读取CSV
  deg_df <- read.csv(csv_path, stringsAsFactors = FALSE)
  
  # 提取gene列，去重、统一格式（大写）
  deg_genes <- deg_df %>%
    select(gene) %>%
    filter(!is.na(gene)) %>%
    distinct() %>%
    mutate(gene = toupper(gene)) %>%
    pull(gene)
  
  # 返回结果：细胞群名称、总DEGs数、基因列表
  cell_type <- gsub("1116_|_DEGs.csv", "", basename(csv_path))  # 从文件名提取细胞群名称
  list(
    cell_type = cell_type,
    total_deg_count = length(deg_genes),
    deg_genes = deg_genes
  )
}

# 批量处理所有CSV文件
deg_list <- map(csv_paths, read_deg_genes)

# ===================== 4. 基因比对（统计交集） =====================
# 定义函数：比对单个细胞群的DEGs与Excel基因
compare_genes <- function(deg_data) {
  # 计算交集基因
  common_genes <- intersect(excel_genes, deg_data$deg_genes)
  common_count <- length(common_genes)
  
  # 生成结果行
  tibble(
    细胞群 = deg_data$cell_type,
    Excel中总基因数 = excel_gene_count,
    该细胞群总DEGs数 = deg_data$total_deg_count,
    交集基因数 = common_count,
    交集基因占Excel基因比例 = round(common_count / excel_gene_count * 100, 2),  # 百分比，保留2位小数
    交集基因占该细胞群DEGs比例 = round(common_count / deg_data$total_deg_count * 100, 2),  # 百分比，保留2位小数
    交集基因列表 = paste(common_genes, collapse = ", ")  # 交集基因用逗号分隔
  )
}

# 批量比对所有细胞群
comparison_result <- map_dfr(deg_list, compare_genes)

# ===================== 5. 添加全局统计行 =====================
# 计算所有细胞群中与Excel基因交集的总基因数（去重）
all_common_genes <- comparison_result %>%
  pull(交集基因列表) %>%
  strsplit(", ") %>%
  unlist() %>%
  unique() %>%
  .[. != ""]  # 去除空字符串（无交集时）
all_common_count <- length(all_common_genes)

# 添加全局统计行
global_stats <- tibble(
  细胞群 = "全局统计",
  Excel中总基因数 = excel_gene_count,
  该细胞群总DEGs数 = sum(map_dbl(deg_list, "total_deg_count")),  # 所有细胞群DEGs总数（去重前）
  交集基因数 = all_common_count,
  交集基因占Excel基因比例 = round(all_common_count / excel_gene_count * 100, 2),
  交集基因占该细胞群DEGs比例 = round(all_common_count / sum(map_dbl(deg_list, "total_deg_count")) * 100, 2),
  交集基因列表 = paste(all_common_genes, collapse = ", ")
)

# 合并结果（细胞群比对结果 + 全局统计）
final_report <- bind_rows(comparison_result, global_stats)

# ===================== 6. 保存报告 =====================
write.csv(
  final_report,
  file = output_path,
  row.names = FALSE,
  quote = FALSE,
  fileEncoding = "UTF-8"
)

# ===================== 7. 输出完成提示 =====================
message("✅ 基因比对报告生成完成！")
message(paste0("📁 报告保存路径：", output_path))
message(paste0("🔍 全局交集基因数：", all_common_count, "个（去重后）"))
message(paste0("📋 交集基因列表：", paste(all_common_genes, collapse = ", ")))