rm(list = ls())

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(scales)

# 设置输出目录
out_dir <- "/Users/lijiangbo/1_Projects/Epi1106/5_scRNA/1113重新分析/1113outputs/8类细胞类型注释后完整可视化"
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
setwd(out_dir)

all10 <- readRDS("/Users/lijiangbo/1_Projects/Epi1106/5_scRNA/2_processed_data/1113/00_processed_data/5_combined_10_去掉簇9和簇12后.rds")
head(all10)
unique(all10$seurat_clusters)
current_levels <- levels(all10$seurat_clusters)
current_levels

# # 定义替换规则：将"10"改为"9"，"11"改为"10"
# new_levels <- ifelse(current_levels == "10", "9", 
#                      ifelse(current_levels == "11", "10", current_levels))
# 
# # 重新设置因子水平
# levels(all10$seurat_clusters) <- new_levels
# 
# # 验证修改结果
# unique(all10$seurat_clusters)
# 
# cluster_order <- as.character(0:10)
# 
# # 直接修改Seurat对象中seurat_clusters的因子水平
# all10$seurat_clusters <- factor(all10$seurat_clusters, levels = cluster_order)
# 
# # 验证修改结果（可选）
# unique(all10$seurat_clusters)
# head(all10)
# unique(all10$Cell.type)

# 第一步：先给Seurat对象添加 cell_type_8class 列（基于之前的cluster映射规则）
# 映射关系：seurat_clusters（0-10）→ 8类细胞类型
all10$cell_type_8class <- dplyr::case_when(
  all10$seurat_clusters == "0" ~ "Oligodendrocytes",
  all10$seurat_clusters == "1" ~ "Astrocytes",
  all10$seurat_clusters == "2" ~ "Microglia",
  all10$seurat_clusters == "3" ~ "Excitatory_Neurons",
  all10$seurat_clusters == "4" ~ "OPCs",
  all10$seurat_clusters == "5" ~ "Excitatory_Neurons",
  all10$seurat_clusters == "6" ~ "Inhibitory_Neurons",
  all10$seurat_clusters == "7" ~ "Inhibitory_Neurons",
  all10$seurat_clusters == "8" ~ "Pericytes",
  all10$seurat_clusters == "9" ~ "Endothelial_Cells",  # 原10重命名后
  all10$seurat_clusters == "10" ~ "Oligodendrocytes", # 原11重命名后
  TRUE ~ "Unknown"  # 兜底（避免遗漏）
)

table(all10$cell_type_8class)  # 查看8类细胞的数量分布

# 1. 统计分析
cat("=== 细胞类型统计分析 ===\n")
stats_overall <- all10@meta.data %>%
  group_by(cell_type_8class) %>%
  summarise(Cell_Count = n(), .groups = 'drop') %>%
  mutate(Percentage = round(Cell_Count / sum(Cell_Count) * 100, 2)) %>%
  arrange(desc(Cell_Count))

stats_detailed <- all10@meta.data %>%
  group_by(cell_type_8class, seurat_clusters) %>%  # 按8类细胞+聚类序号分组
  summarise(Cell_Count = n(), .groups = 'drop') %>%
  arrange(cell_type_8class, seurat_clusters)  # 按细胞类型+聚类序号排序

stats_qc <- all10@meta.data %>%
  group_by(cell_type_8class) %>%
  summarise(
    Cell_Count = n(),
    Mean_nFeature = round(mean(nFeature_RNA), 0),
    Median_nFeature = round(median(nFeature_RNA), 0),
    Mean_nCount = round(mean(nCount_RNA), 0),
    Mean_pct_mt = round(mean(percent.mt), 2),
    .groups = 'drop'
  ) %>%
  arrange(desc(Cell_Count))

# 第四步：保存统计表格
write.csv(stats_overall, "Table_S1_CellType_Overall_Statistics.csv", row.names = FALSE)
write.csv(stats_detailed, "Table_S2_CellType_Detailed_Statistics.csv", row.names = FALSE)
write.csv(stats_qc, "Table_S3_CellType_QC_Statistics.csv", row.names = FALSE)

cat("=== 统计分析完成，表格已保存 ===\n")

# saveRDS(all10, "/Users/lijiangbo/1_Projects/Epi1106/5_scRNA/2_processed_data/1113/00_processed_data/5_combined_10_去掉簇9和簇12后.rds")

# 2. 颜色配置
celltype_colors <- c(
  "Astrocytes" = "#E64B35FF", "Microglia" = "#4DBBD5FF", 
  "Oligodendrocytes" = "#00A087FF", "OPCs" = "#3C5488FF",
  "Excitatory_Neurons" = "#F39B7FFF", "Inhibitory_Neurons" = "#8491B4FF",
  "Endothelial_Cells" = "#91D1C2FF", "Pericytes" = "#DC0000FF"
)
cluster_colors <- colorRampPalette(brewer.pal(11, "Spectral"))(11)
cell_type_order <- c("Astrocytes", "Microglia", "Oligodendrocytes", "OPCs",
                     "Excitatory_Neurons", "Inhibitory_Neurons", 
                     "Endothelial_Cells", "Pericytes")

# 3. Marker基因验证
markers_list <- list(
  Astrocytes = c("AQP4", "ALDH1A1", "GJA1", "GFAP", "SLC1A2", "SLCO1C1"),
  Microglia = c("P2RY12", "CX3CR1", "TMEM119", "C1QA", "C1QB", "CSF1R"),
  Oligodendrocytes = c("MOBP", "OPALIN", "MOG", "PLP1", "MBP", "MAG"),
  OPCs = c("PDGFRA", "CSPG4", "VCAN", "COL9A1", "GPR17", "MYT1"),
  Excitatory_Neurons = c("SLC17A7", "CAMK2A", "CUX2", "SATB2", "TBR1", "NWD2"),
  Inhibitory_Neurons = c("GAD1", "GAD2", "SLC32A1", "KCNC2", "SST", "VIP", "PVALB"),
  Endothelial_Cells = c("CLDN5", "FLT1", "PECAM1", "VWF", "CDH5", "MECOM"),
  Pericytes = c("PDGFRB", "RGS5", "ANPEP", "DCN", "ITIH5", "ABCC9")
)
markers_available <- lapply(markers_list, function(ms) ms[ms %in% rownames(all10)])
all_markers <- unique(unlist(markers_available))

# DotPlot
Idents(all10) <- "cell_type_8class"
p_dotplot <- DotPlot(all10, features = all_markers, cols = c("lightgrey", "red"), 
                     dot.scale = 8, cluster.idents = FALSE) +
  RotatedAxis() + theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, face = "italic"),
        axis.text.y = element_text(face = "bold"),
        panel.border = element_rect(color = "black", fill = NA)) +
  labs(title = "Canonical Marker Expression", x = "Marker Genes", y = "Cell Types") +
  scale_y_discrete(limits = rev(cell_type_order))
ggsave("Fig1_DotPlot_Marker_Genes.pdf", p_dotplot, width = 20, height = 7)
ggsave("Fig1_DotPlot_Marker_Genes.png", p_dotplot, width = 20, height = 7)

head(all10)
# FeaturePlot
# 定义x轴刻度的有序顺序（0到12）
cluster_order <- as.character(0:10)
all10$seurat_clusters <- factor(all10$seurat_clusters, levels = cluster_order)

# 关键marker列表
key_markers <- c("MOBP", "AQP4", "P2RY12", "SLC17A7", "PDGFRA",  "GAD2", "DCN", "FLT1")
# 替换SLC17A7的逻辑保持不变
key_markers <- if(!"SLC17A7" %in% rownames(all10) && "CAMK2A" %in% rownames(all10)) 
  replace(key_markers, key_markers == "SLC17A7", "CAMK2A") else key_markers

# 绘制FeaturePlot并设置x轴顺序
p_feature <- FeaturePlot(
  all10, 
  features = key_markers, 
  ncol = 4, 
  pt.size = 0.1, 
  order = TRUE,
  split.by = "seurat_clusters",  # 按自定义的cluster顺序拆分
) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # 调整x轴标签角度，避免重叠

# 保存图片
ggsave("Fig2_FeaturePlot_Key_Markers_Ordered.pdf", p_feature, width = 20, height = 8)
ggsave("Fig2_FeaturePlot_Key_Markers_Ordered.png", p_feature, width = 20, height = 8)


# ViolinPlot
# 定义x轴刻度的有序顺序（0到10）
cluster_order <- as.character(0:10)
all10$seurat_clusters <- factor(all10$seurat_clusters, levels = cluster_order)

# 关键marker列表
key_markers <- c("MOBP", "AQP4", "P2RY12", "SLC17A7", "PDGFRA",  "GAD2", "DCN", "FLT1")
# 绘制小提琴图并设置x轴样式
Idents(all10) <- "seurat_clusters"
p_violin <- VlnPlot(
  all10, 
  features = key_markers, 
  ncol = 4, 
  pt.size = 0, 
  cols = cluster_colors  # 若需保留自定义配色，否则可移除
) &
  theme_classic() &
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10)  # 标签水平居中，不倾斜
  )

# 保存图片
ggsave("Fig3_ViolinPlot_Markers_by_Cluster_0-10.pdf", p_violin, width = 17, height = 8)
ggsave("Fig3_ViolinPlot_Markers_by_Cluster_0-10.png", p_violin, width = 17, height = 8)

# 热图
# 保持细胞类型顺序（如需默认顺序可移除此行，但建议保留以确保一致性）
all10$cell_type_8class <- factor(all10$cell_type_8class, levels = cell_type_order)
Idents(all10) <- "cell_type_8class"

# 绘制热图并调整样式
p_heatmap <- DoHeatmap(
  all10, 
  features = all_markers, 
  group.by = "cell_type_8class", 
  raster = FALSE  # 非栅格化保证清晰度
) +
  theme(
    strip.text.x = element_blank(),  # 完全隐藏顶部细胞类型标签
    axis.text.y = element_text(size = 15, face = "plain"),  # 进一步放大基因名字体
    legend.position = "right",       # 保留右侧图例
    legend.text = element_text(size = 14)  # 图例字体放大（可根据需求调整size值）
  )

# 保存调整后的热图
ggsave("Fig4_Heatmap_Marker_Genes.pdf", p_heatmap, width = 20, height = 20, device = cairo_pdf)
ggsave("Fig4_Heatmap_Marker_Genes.png", p_heatmap, width = 20, height = 20, dpi = 600)

# 4. UMAP可视化
# 原始cluster
Idents(all10) <- "seurat_clusters"
p_umap_original <- DimPlot(all10, reduction = "umap", label = TRUE, cols = cluster_colors) +
  theme_minimal() + theme(panel.border = element_rect(color = "black", fill = NA)) +
  labs(title = "Original Clusters")
ggsave("Fig5A_UMAP_seurat_clusterss.pdf", p_umap_original, width = 12, height = 9)
ggsave("Fig5A_UMAP_seurat_clusterss.png", p_umap_original, width = 12, height = 9)

# 8类细胞
Idents(all10) <- "cell_type_8class"
p_umap_8class <- DimPlot(
  all10, 
  reduction = "umap", 
  cols = celltype_colors, 
  label = TRUE, 
  label.size = 6  # 增大label字体（数值可调整，如7、8等）
) +
  theme_minimal() + theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # 边框加粗（size值调大）
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  ) +
  labs(title = "Eight Major Cell Types", x = "UMAP_1", y = "UMAP_2")
ggsave("Fig5B_UMAP_8_Major_CellTypes.pdf", p_umap_8class, width = 12, height = 9)
ggsave("Fig5B_UMAP_8_Major_CellTypes.png", p_umap_8class, width = 12, height = 9)

# 对比图
p_comparison <- p_umap_original + p_umap_8class + 
  plot_annotation(title = "Original Clusters vs Annotated Types")
ggsave("Fig5C_UMAP_Comparison.pdf", p_comparison, width = 24, height = 9)
ggsave("Fig5C_UMAP_Comparison.png", p_comparison, width = 24, height = 9)

# 分面图
p_split <- DimPlot(all10, split.by = "cell_type_8class", ncol = 4, cols = celltype_colors) +
  theme_minimal() + theme(legend.position = "none") + labs(title = "Cell Type Distribution")
ggsave("Fig5D_UMAP_Split_by_CellType.pdf", p_split, width = 16, height = 8)
ggsave("Fig5D_UMAP_Split_by_CellType.png", p_split, width = 16, height = 8)

# 发表级UMAP
p_umap_pub <- DimPlot(all10, cols = celltype_colors, label = FALSE) +
  theme_classic() + theme(panel.border = element_rect(color = "black", fill = NA)) +
  labs(title = "Single Cell Transcriptomic Atlas")
ggsave("Fig5E_UMAP_Publication.pdf", p_umap_pub, width = 11, height = 9)
ggsave("Fig5E_UMAP_Publication.png", p_umap_pub, width = 11, height = 9)

# 5. 统计图表
# 柱状图
p_bar <- ggplot(stats_overall, aes(x = reorder(cell_type_8class, -Cell_Count), 
                                   y = Cell_Count, fill = cell_type_8class)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = paste0(Cell_Count, "\n(", Percentage, "%)")), vjust = -0.3) +
  scale_fill_manual(values = celltype_colors) + theme_classic() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Cell Type Composition", y = "Number of Cells", x = "")
ggsave("Fig6A_Barplot_CellComposition.pdf", p_bar, width = 12, height = 8)
ggsave("Fig6A_Barplot_CellComposition.png", p_bar, width = 12, height = 8)

# 饼图
p_pie <- ggplot(stats_overall, aes(x = "", y = Percentage, fill = cell_type_8class)) +
  geom_bar(stat = "identity", width = 1, color = "white") + coord_polar("y") +
  scale_fill_manual(values = celltype_colors) + theme_void() +
  geom_text(aes(label = paste0(cell_type_8class, "\n", Percentage, "%")), 
            position = position_stack(vjust = 0.5), color = "white") +
  labs(title = "Cell Type Distribution")
ggsave("Fig6B_Piechart_CellDistribution.pdf", p_pie, width = 11, height = 10)
ggsave("Fig6B_Piechart_CellDistribution.png", p_pie, width = 11, height = 10)

# 堆叠柱状图
p_stacked <- ggplot(stats_detailed, aes(x = cell_type_8class, y = Cell_Count, 
                                        fill = factor(seurat_clusters))) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = paste0("C", seurat_clusters)), position = position_stack(vjust = 0.5)) +
  scale_fill_brewer(palette = "Set3") + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Cluster Composition Details", y = "Cell Count", x = "")
ggsave("Fig6C_StackedBar_ClusterDetails.pdf", p_stacked, width = 14, height = 8)
ggsave("Fig6C_StackedBar_ClusterDetails.png", p_stacked, width = 14, height = 8)

# QC小提琴图
p_qc1 <- VlnPlot(all10, features = "nFeature_RNA", cols = celltype_colors, pt.size = 0) +
  theme_classic() + labs(title = "Gene Count Distribution")
p_qc2 <- VlnPlot(all10, features = "nCount_RNA", cols = celltype_colors, pt.size = 0) +
  theme_classic() + labs(title = "UMI Count Distribution")
p_qc3 <- VlnPlot(all10, features = "percent.mt", cols = celltype_colors, pt.size = 0) +
  theme_classic() + labs(title = "Mitochondrial Percentage")
p_qc_combined <- p_qc1 / p_qc2 / p_qc3
ggsave("Fig7_QC_Metrics.pdf", p_qc_combined, width = 14, height = 15)
ggsave("Fig7_QC_Metrics.png", p_qc_combined, width = 14, height = 15)

# 6. 验证cluster合并
# 兴奋性神经元
excitatory_markers <- c("SLC17A7", "CAMK2A", "CUX2")
if(any(excitatory_markers %in% rownames(all10))) {
  Idents(all10) <- "seurat_clusters"
  p_ex_vln <- VlnPlot(all10, features = excitatory_markers[excitatory_markers %in% rownames(all10)],
                      idents = c("3", "5"), cols = rep("#F39B7FFF", 2)) + theme_classic()
  ggsave("Fig8A_Excitatory_MergeValidation.pdf", p_ex_vln, width = 12, height = 6)
}

# 抑制性神经元
inhibitory_markers <- c("GAD1", "GAD2", "SST")
if(any(inhibitory_markers %in% rownames(all10))) {
  p_in_vln <- VlnPlot(all10, features = inhibitory_markers[inhibitory_markers %in% rownames(all10)],
                      idents = c("6", "7"), cols = rep("#8491B4FF", 2)) + theme_classic()
  ggsave("Fig8B_Inhibitory_MergeValidation.pdf", p_in_vln, width = 12, height = 6)
}

# 少突胶质细胞
oligo_markers <- c("MOBP", "MOG", "PLP1")
if(any(oligo_markers %in% rownames(all10))) {
  p_ol_vln <- VlnPlot(all10, features = oligo_markers[oligo_markers %in% rownames(all10)],
                      idents = c("0", "11"), cols = rep("#00A087FF", 2)) + theme_classic()
  ggsave("Fig8C_Oligo_MergeValidation.pdf", p_ol_vln, width = 12, height = 6)
}

# 7. 导出数据表格
write.csv(data.frame(
  Cell_ID = colnames(all10),
  seurat_clusters = all10$seurat_clusters,
  Cell_Type = all10$cell_type_8class,
  UMAP_1 = Embeddings(all10, "umap")[,1],
  UMAP_2 = Embeddings(all10, "umap")[,2],
  nFeature_RNA = all10$nFeature_RNA,
  nCount_RNA = all10$nCount_RNA,
  Percent_MT = all10$percent.mt
), "Table_S4_Complete_CellAnnotations.csv", row.names = FALSE)

write.csv(data.frame(
  Cell_Type = rep(names(markers_list), sapply(markers_list, length)),
  Marker_Gene = unlist(markers_list),
  Available = unlist(markers_list) %in% rownames(all10)
), "Table_S5_MarkerGeneList.csv", row.names = FALSE)

write.csv(data.frame(
  seurat_clusters = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "10", "11"),
  Cell_Type = c("Oligodendrocytes", "Astrocytes", "Microglia", "Excitatory_Neurons",
                "OPCs", "Excitatory_Neurons", "Inhibitory_Neurons", "Inhibitory_Neurons",
                "Pericytes", "Endothelial_Cells", "Oligodendrocytes")
), "Table_S6_ClusterMapping.csv", row.names = FALSE)

# 8. 保存Seurat对象和报告
saveRDS(all10, "all10_Annotated_8CellTypes.rds")

sink("Analysis_Report.txt")
cat("单细胞RNA测序细胞类型注释报告\n")
cat("分析日期:", Sys.Date(), "\n")
cat("总细胞数:", ncol(all10), "\n")
cat("细胞类型组成:\n")
print(stats_overall)
sink()

cat("=== 分析完成，结果保存至", out_dir, "===\n")