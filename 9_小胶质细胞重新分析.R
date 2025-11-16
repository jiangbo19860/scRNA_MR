library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

setwd("~/1_Projects/Epi1106/5_scRNA/1113重新分析/1113outputs/1115/microglia")

mg4 <- readRDS("/Users/lijiangbo/1_Projects/Epi1106/5_scRNA/2_processed_data/1113/00_processed_data/Microglia_亚分4类.rds")
head(mg4)
# 聚类与UMAP降维
mg4 <- FindNeighbors(mg4, dims = 1:15)
mg4 <- FindClusters(mg4, resolution = 0.03)
mg4 <- RunUMAP(mg4, dims = 1:15)
unique(mg4$SCT_snn_res.0.03)

p1 <- DimPlot(mg4, label = TRUE)
p2 <- DimPlot(mg4, group.by = "orig.ident") + 
  theme(plot.title = element_blank())  # 隐藏右边图的标题
p1 | p2

ggsave("mg4_subclusters_UMAP.pdf", width = 10, height = 5)

table(mg4$SCT_snn_res.0.03)

DefaultAssay(mg4) <- "SCT"
mg_markers <- FindAllMarkers(
  mg4, 
  only.pos = TRUE,           # 只看上调基因
  min.pct = 0.25,            # 至少25%细胞表达
  logfc.threshold = 0.25,    # log2FC > 0.25
  test.use = "wilcox"        # 或用 "MAST"（更严格但慢）
)

top10_markers_full <- mg_markers %>%
  filter(p_val_adj < 0.05) %>%  # 仅保留显著差异基因
  group_by(cluster) %>%         # 按簇分组
  # 排序规则：显著性优先，再按表达差异降序
  arrange(p_val_adj, desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 10) %>%        # 取每个簇Top10
  ungroup() %>%
  mutate(rank = rep(1:10, length(unique(cluster))))  # 新增排名列（1-10）

# 基于完整数据，生成“汇总查看用”的表格（合并基因名、log2FC、P值）
top10_markers_per_cluster <- top10_markers_full %>%
  group_by(cluster) %>%
  summarise(
    top10_genes = paste(gene, collapse = ", "),
    top10_log2FC = paste(round(avg_log2FC, 2), collapse = ", "),  # 保留2位小数
    top10_padj = paste(format(p_val_adj, scientific = TRUE, digits = 2), collapse = ", ")
  )

print(top10_markers_per_cluster)
table(mg4$seurat_clusters)  # 此处mg替换为mg4
write.csv(
  top10_markers_full, 
  "mg_top10_cluster_markers_full.csv",  # 文件名可根据需要修改（如加mg4前缀）
  row.names = FALSE,
  quote = FALSE  # 避免基因名带引号
)