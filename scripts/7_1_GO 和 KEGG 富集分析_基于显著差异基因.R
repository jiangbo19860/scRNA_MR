# 脚本1：GO和KEGG富集分析
# 功能：对差异表达基因进行GO（生物过程、细胞组分、分子功能）和KEGG通路富集分析，并可视化结果

rm(list = ls())
# 1. 加载依赖包
pacman::p_load(
  here, Seurat, tidyverse, clusterProfiler, AnnotationDbi, org.Hs.eg.db,
  ggplot2, ggpubr, patchwork
)


# 2. 数据加载与预处理
# 加载单细胞数据对象
load(here("1_data", "E-MTAB-6149", "scRNA_endo.Rds"))

# 查看数据结构
str(scRNA, max.level = 2) # 简化显示结构
dim(scRNA) # 查看细胞和基因数量

# 设置细胞分组（肿瘤vs非肿瘤）
scRNA$group <- ifelse(scRNA$CellFromTumor == "1", "Tumor", "Non-malignant")
Idents(scRNA) <- "group" # 设置分组为分析身份
table(Idents(scRNA)) # 确认分组正确


# 3. 差异基因分析（为富集分析提供输入）
# 计算肿瘤vs非肿瘤的差异基因
deg <- FindMarkers(
  scRNA,
  ident.1 = "Tumor",
  ident.2 = "Non-malignant",
  logfc.threshold = 0 # 保留所有差异基因（后续筛选）
)

# 筛选显著差异基因（FDR<0.01且|log2FC|>1）
sig_deg <- subset(deg, p_val_adj < 0.01 & abs(avg_log2FC) > 1)
cat("显著差异基因数量：", nrow(sig_deg), "\n")

# 可视化差异基因分布（火山图）
ggplot(deg, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = factor(ifelse(p_val_adj < 0.01 & abs(avg_log2FC) > 1, "Significant", "Non-significant"),
    levels = c("Significant", "Non-significant")
  ))) +
  scale_color_manual(values = c("red", "gray")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "blue") +
  labs(
    title = "Volcano Plot of Differential Expression",
    x = "Log2 Fold Change",
    y = "-Log10(Adjusted p-value)",
    color = "Significance"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(), # 移除主要网格线
    panel.grid.minor = element_blank(), # 移除次要网格线
    axis.line = element_line(color = "black", size = 0.5) # 添加x轴和y轴的黑线
  )


# 4. GO富集分析
# 基于显著差异基因进行GO富集（包含BP、CC、MF）
ego_ALL <- enrichGO(
  gene = rownames(sig_deg), # 输入显著差异基因（基因符号）
  OrgDb = org.Hs.eg.db, # 人类基因注释数据库
  keyType = "SYMBOL", # 输入基因ID类型为基因符号
  ont = "ALL", # 分析所有GO类别（BP+CC+MF）
  pAdjustMethod = "BH", # FDR校正方法
  pvalueCutoff = 0.01, # P值阈值
  qvalueCutoff = 0.05 # 校正后P值阈值
)

# 转换为数据框便于查看
go_results <- data.frame(ego_ALL)
cat("GO富集结果前5行：\n")
head(go_results[, 1:6]) # 显示关键列

# 可视化GO富集结果（分面展示BP、CC、MF）
# Visualize GO enrichment results (faceted by BP, CC, MF)
go_dotplot <- dotplot(ego_ALL, split = "ONTOLOGY") +
  facet_grid(ONTOLOGY ~ ., scales = "free") +
  labs(title = "Gene Ontology Enrichment Analysis") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12), # 增大x轴标签字体
    axis.text.y = element_text(size = 12), # 增大y轴标签字体
    axis.title = element_text(size = 14, face = "bold"), # 增大轴标题字体
    legend.text = element_text(size = 12), # 增大图例文本字体
    legend.title = element_text(size = 13, face = "bold"), # 增大图例标题字体
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), # 增大标题字体并居中
    strip.text = element_text(size = 26, face = "bold") # 增大分面标签字体
  )
print(go_dotplot)

# 5. KEGG富集分析
# 5.1 基因ID转换：KEGG需要Entrez ID，将基因符号转为Entrez ID
kegg_genes <- bitr(
  geneID = rownames(sig_deg),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)
# 提取转换后的Entrez ID（过滤无法转换的基因）
entrez_ids <- kegg_genes$ENTREZID
cat("可转换为Entrez ID的基因数量：", length(entrez_ids), "\n")

# 5.2 执行KEGG富集分析
ekegg <- enrichKEGG(
  gene = entrez_ids, # 输入Entrez ID
  organism = "hsa", # 人类（hsa = Homo sapiens）
  pAdjustMethod = "BH", # FDR校正
  pvalueCutoff = 0.05 # 阈值
)

# 转换为数据框
kegg_results <- data.frame(ekegg)
cat("KEGG富集结果前5行：\n")
head(kegg_results[, 1:6])

# 可视化KEGG富集结果
kegg_dotplot <- dotplot(ekegg, showCategory = 20) +
  labs(title = "KEGG通路富集分析结果", x = "基因比率", y = "通路名称") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(kegg_dotplot)


# 6. 结果保存
# 创建结果目录
if (!dir.exists(here("3_outputs", "GO_KEGG"))) {
  dir.create(here("3_outputs", "GO_KEGG"), recursive = TRUE)
}

# 保存富集结果表格
write.csv(go_results, here("3_outputs", "GO_KEGG", "go_enrichment_results.csv"), row.names = FALSE)
write.csv(kegg_results, here("3_outputs", "GO_KEGG", "kegg_enrichment_results.csv"), row.names = FALSE)

# 保存可视化图表
ggsave(here("3_outputs", "GO_KEGG", "go_dotplot.png"), go_dotplot, width = 10, height = 8, dpi = 300)
ggsave(here("3_outputs", "GO_KEGG", "kegg_dotplot.png"), kegg_dotplot, width = 10, height = 8, dpi = 300)

cat("GO和KEGG富集分析完成！结果已保存至3_outputs/GO_KEGG目录\n")
