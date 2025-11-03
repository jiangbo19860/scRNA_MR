# 脚本2：GSEA（基因集富集分析）
# 功能：基于所有基因的log2FC排序，分析预定义基因集（如Hallmark）的富集趋势，并可视化关键通路

rm(list = ls()) # 清空工作空间
# 1. 加载依赖包
pacman::p_load(
  here, Seurat, tidyverse, fgsea, enrichplot, clusterProfiler,
  AnnotationDbi, org.Hs.eg.db, ggplot2, ggrepel
)


# 2. 数据加载与预处理
# 加载单细胞数据对象
load(here("1_data", "E-MTAB-6149", "scRNA_endo.Rds"))

# 设置细胞分组（与脚本1一致）
scRNA$group <- ifelse(scRNA$CellFromTumor == "1", "Tumor", "Non-malignant")
Idents(scRNA) <- "group"


# 3. 差异基因分析（获取所有基因的log2FC，用于排序）
# 计算所有基因的差异（不筛选，保留全基因用于GSEA）
deg_all <- FindMarkers(
  scRNA,
  ident.1 = "Tumor",
  ident.2 = "Non-malignant",
  logfc.threshold = 0 # 保留所有基因
)
cat("全基因差异分析结果维度：", dim(deg_all), "\n") # 行=基因，列=差异指标


# 4. GSEA输入准备：排序的基因列表（按log2FC降序，基因ID与基因集匹配）
# 4.1 基因ID转换：确保与Hallmark基因集的ID类型一致（此处为基因符号）
# 注：若Hallmark基因集为Entrez ID，需转换为Entrez ID（见注释）
# 提取基因符号（行名）和对应的log2FC
gene_df <- data.frame(
  symbol = rownames(deg_all),
  log2fc = deg_all$avg_log2FC,
  stringsAsFactors = FALSE
)

# 4.2 按log2FC降序排序（GSEA要求基因按差异程度排序）
sorted_genes <- gene_df %>%
  arrange(desc(log2fc)) %>% # 从高表达到低表达排序
  filter(!is.na(symbol)) # 过滤无基因名的条目

# 4.3 构建命名向量（基因符号→log2FC）：GSEA核心输入
ranked_genes <- setNames(sorted_genes$log2fc, sorted_genes$symbol)
cat("排序后基因数量：", length(ranked_genes), "\n")
head(ranked_genes) # 查看前几个基因的log2FC


# 5. 加载预定义基因集（如Hallmark）
# 从GMT文件读取Hallmark基因集（确保文件路径正确）
hallmark_gmt <- read.gmt(here("1_data", "GSEA", "h.all.v2025.1.Hs.symbols.gmt"))
# 转换为列表格式（通路名称→基因列表）
hallmark_list <- hallmark_gmt %>%
  split(.$term) %>% # 按通路名称分组
  lapply("[[", "gene") # 每个通路保留基因列表

# 查看基因集示例（如凋亡通路）
cat("凋亡通路基因数量：", length(hallmark_list[["HALLMARK_APOPTOSIS"]]), "\n")
head(hallmark_list[["HALLMARK_APOPTOSIS"]])


# 6. 用fgsea算法执行GSEA分析
fgsea_res <- fgsea(
  pathways = hallmark_list,
  stats = ranked_genes,
  minSize = 15, # 保留其他参数
  maxSize = 500 # 保留其他参数
  # 移除 nperm 参数，自动使用 fgseaMultilevel
)

# 转换为数据框查看结果
gsea_df <- data.frame(fgsea_res)
cat("GSEA结果前5行：\n")
head(gsea_df[, c("pathway", "NES", "padj", "size")]) # 关键列：通路、NES、校正后P值、基因数


# 7. GSEA结果筛选与可视化
# 确保加载必要的包（特别是 dplyr 用于数据处理）
library(clusterProfiler)
library(dplyr) # 提供 select 函数（用于筛选列）
library(enrichplot)
library(here)

# 1. 读取 Hallmark 基因集（GMT 文件）
# 注意：read.gmt 读取的列名通常是 "term"、"description"、"gene"（不同包可能有差异）
hallmark_gmt <- read.gmt(here("1_data", "GSEA", "h.all.v2025.1.Hs.symbols.gmt"))

# 查看读取的列名（关键步骤：确认实际列名）
colnames(hallmark_gmt) # 通常为 c("term", "description", "gene") 或类似

# 2. 转换为 TERM2GENE 格式（仅保留通路名和基因名两列）
# 若列名为 "term" 和 "gene"：
hallmark <- hallmark_gmt %>%
  dplyr::select(term, gene) # 用 dplyr::select 明确指定包，避免冲突

# 若列名不同（如 "pathway" 和 "gene"），根据实际列名修改：
# hallmark <- hallmark_gmt %>%
#   dplyr::select(pathway, gene) %>%  # 替换为实际列名
#   rename(term = pathway)  # 重命名为 "term"（TERM2GENE 要求第一列为通路名）

# 3. 确认 TERM2GENE 格式正确（两列：term 和 gene）
head(hallmark)
colnames(hallmark) # 必须是 c("term", "gene")

# 4. 执行 GSEA 分析（使用 clusterProfiler::GSEA）
# 确保 ranked_genes 是命名向量（基因名 -> log2FC，已排序）
gsea_result <- GSEA(
  geneList = ranked_genes,
  TERM2GENE = hallmark, # 正确的 TERM2GENE 数据框
  pvalueCutoff = 0.25
)

# 5. 绘制气泡图
gsea_bubble <- dotplot(
  gsea_result,
  showCategory = 15,
  x = "NES",
  color = "p.adjust"
) +
  labs(title = "GSEA Significant Enriched Pathways") +
  theme_minimal()

print(gsea_bubble)

# 7.3 单个通路富集曲线（如缺氧通路和凋亡通路）
# 缺氧通路
# 使用 clusterProfiler 的 GSEA 结果（而非 fgsea 的结果）
gseaplot2(
  gsea_result, # 关键修改：使用 clusterProfiler 的 GSEA 结果
  geneSetID = "HALLMARK_HYPOXIA",
  title = "Hypoxia Pathway GSEA Enrichment Plot",
  pvalue_table = TRUE,
  color = "red"
)

# 凋亡通路
gseaplot2(
  gsea_result, # 统一使用 clusterProfiler 的 GSEA 结果
  geneSetID = "HALLMARK_APOPTOSIS",
  title = "Apoptosis Pathway GSEA Enrichment Plot", # 建议统一为英文标题
  pvalue_table = TRUE,
  color = "blue"
)


# 8. 结果保存
if (!dir.exists(here("3_outputs", "GSEA"))) {
  dir.create(here("3_outputs", "GSEA"), recursive = TRUE)
}

# 检查 gsea_df 的列类型（确认哪些列是列表）
sapply(gsea_df, class) # 会显示 leadingEdge 为 "list"

# 将列表列（如 leadingEdge）转换为字符串（用分号分隔基因）
gsea_df$leadingEdge <- sapply(gsea_df$leadingEdge, function(x) paste(x, collapse = ";"))

# 再次保存 CSV（此时所有列均为字符串/数值类型，可正常写入）
write.csv(
  gsea_df,
  here("3_outputs", "GSEA", "gsea_all_results.csv"),
  row.names = FALSE
)
