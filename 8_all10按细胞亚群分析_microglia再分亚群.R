library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
all10 <- readRDS("/Users/lijiangbo/1_Projects/Epi1106/5_scRNA/2_processed_data/1113/00_processed_data/5_combined_10_去掉簇9和簇12后.rds")
head(all10)
DefaultAssay(all10)
head(all10@meta.data)
unique(all10$cell_type_8class)

# 拆分注释后的8个细胞群 -------------------------------------------------------------
# 设置保存路径（确保路径存在，若不存在则创建）
save_dir <- "/Users/lijiangbo/1_Projects/Epi1106/5_scRNA/2_processed_data/1113/00_processed_data"
if (!dir.exists(save_dir)) {
  dir.create(save_dir, recursive = TRUE)
}

# 获取所有细胞类型
cell_types <- unique(all10$cell_type_8class)

# 循环拆分并保存
for (ct in cell_types) {
  # 按细胞类型子集化Seurat对象
  subset_obj <- subset(all10, subset = cell_type_8class == ct)
  
  # 计算细胞数（列数）和基因数（行数，基于当前assay）
  n_cells <- ncol(subset_obj)
  n_genes <- nrow(subset_obj)  # 基因数为对象的行数（所有assay共享基因集）
  
  # 构建文件名（替换空格/特殊字符，确保合规）
  ct_clean <- gsub(" ", "_", ct)  # 这里细胞类型已用下划线，可省略但保留以防万一
  filename <- sprintf(
    "%s_%s_cells_%d_genes_%d.rds",
    ct_clean,
    current_date,
    n_cells,
    n_genes
  )
  
  # 完整保存路径
  save_path <- file.path(save_dir, filename)
  
  # 保存rds文件
  saveRDS(subset_obj, file = save_path)
  
  # 打印进度信息
  cat(sprintf("已保存：%s（细胞数：%d，基因数：%d）\n", filename, n_cells, n_genes))
}


# 分析小胶质细胞亚群：核心功能：吞噬、抗炎、突触调控 -----------------------------------------------------------------
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)


mg <- readRDS("/Users/lijiangbo/1_Projects/Epi1106/5_scRNA/2_processed_data/1113/00_processed_data/Microglia_20251115_cells_5662_genes_14851.rds")
head(mg)
nrow(mg)
ncol(mg)

# 执行SCTransform，使用标准名称"SCT"（避免后续识别问题）
mg <- SCTransform(mg, assay = "RNA", new.assay.name = "SCT", verbose = FALSE)
DefaultAssay(mg) <- "SCT"  # 仅设置一次默认assay

# 高变基因筛选（基于"SCT" assay）
mg <- FindVariableFeatures(mg, selection.method = "vst", nfeatures = 2000)
# 绘制高变基因图
top20_var_genes <- head(VariableFeatures(mg), 20)
VariableFeaturePlot(mg) %>% 
  LabelPoints(points = top20_var_genes, repel = TRUE, size = 3) +
  labs(title = "Top 20 Variable Features + All Variable Genes")
ggsave("mg_variable_features.png", width = 10, height = 6)
# 归一化
mg <- ScaleData(mg, verbose = FALSE)
# PCA及elbow图
mg <- RunPCA(mg, features = VariableFeatures(mg), verbose = FALSE)
ElbowPlot(mg, ndims = 30)

# 聚类与UMAP降维
mg <- FindNeighbors(mg, dims = 1:15)
mg <- FindClusters(mg, resolution = 0.1)
mg <- RunUMAP(mg, dims = 1:15)

unique(mg$orig.ident)

# UMAP聚类图（带亚群标签）
DimPlot(mg, label = TRUE, repel = TRUE) + 
  labs(title = "Microglia Subclusters (UMAP)")
# 保存图片
ggsave("mg_subclusters_UMAP.png", width = 10, height = 8)

DimPlot(
  mg, 
  group.by = "orig.ident",  # 核心：按健康/癫痫分组
  label = FALSE,  # 若不需要亚群标签可设为FALSE（避免与分组颜色混淆）
  repel = TRUE,
  cols = c("Healthy" = "#2E86AB", "Epilepsy" = "#A23B72")  # 自定义颜色（健康蓝、癫痫紫）
)

cluster_counts <- table(mg$seurat_clusters)
print(cluster_counts)
# 0    1    2    3    4    5    6    7    8 
# 1811 1738  842  500  240  201  128  127   75 

# 按Sample_ID拆分绘图，直接在DimPlot中使用split.by
p <- DimPlot(mg, split.by = "Sample_ID", label = TRUE, repel = TRUE) +
  facet_wrap(~Sample_ID, nrow = 2, ncol = 5) +  # ★核心修改
  labs(title = "Microglia Subclusters by Sample (10 Samples)") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    strip.text = element_text(size = 10),
    axis.line.x = element_line(color = "black", size = 0.5),  # ★添加X轴线
    axis.line.y = element_line(color = "black", size = 0.5)   # Y轴线（可选）
  )

print(p)
# 保存
ggsave("sample_subclusters_2x5.pdf", p, width = 20, height = 8)
ggsave("sample_subclusters_2x5.png", p, width = 20, height = 8, dpi = 300)


# 标记基因鉴定 ------------------------------------------------------------------
# 找每个cluster的marker基因
DefaultAssay(mg) <- "SCT"
mg_markers <- FindAllMarkers(
  mg, 
  only.pos = TRUE,           # 只看上调基因
  min.pct = 0.25,            # 至少25%细胞表达
  logfc.threshold = 0.25,    # log2FC > 0.25
  test.use = "wilcox"        # 或用 "MAST"（更严格但慢）
)

# 保存结果
write.csv(mg_markers, "mg_cluster_markers.csv", row.names = FALSE)

# 根据文献，小胶质细胞的主要状态包括
#      状态	        代表性 Marker
# 稳态/Homeostatic	P2RY12, TMEM119, CX3CR1, SALL1
# 激活/Activated	CD68, CD74, APOE, SPP1, GPNMB
# 炎症/Inflammatory	IL1B, TNF, CCL3, CCL4, CXCL10
# 增殖/Proliferative	MKI67, TOP2A, PCNA
# 疾病相关/DAM	APOE, TREM2, LPL, TYROBP, CST7
# 抗炎/Anti-inflammatory	IL10, MERTK, CD163

# 提取每个簇的Top10显著标志物（核心优化：先过滤显著基因，再排序）
top10_markers_full <- mg_markers %>%
  filter(p_val_adj < 0.05) %>%  # 仅保留显著差异基因（仅需改一次）
  group_by(cluster) %>%         # 按簇分组
  # 排序规则（仅需定义一次）：显著性优先，再按表达差异降序
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
table(mg$seurat_clusters)
write.csv(
  top10_markers_full, 
  "mg_top10_cluster_markers_full.csv", 
  row.names = FALSE,
  quote = FALSE  # 避免基因名带引号
)


# 删除簇6-是少突胶质细胞污染，高表达MOBP, OPALIN ------------------------------------------
mg <- subset(mg, seurat_clusters != 6)
table(mg$seurat_clusters)
cluster8_sample <- mg@meta.data %>%
  filter(seurat_clusters == 8) %>%
  group_by(orig.ident) %>%
  summarise(cell_count = n())
print(cluster8_sample)  # 预期结果：仅Healthy组有数据，Epilepsy组无


# 删除后重新SCT标准化 -------------------------------------------------------------
# 1. 重新运行SCT标准化（自带高变基因筛选和标准化，无需额外步骤）
mg <- SCTransform(
  mg,
  assay = "RNA",
  new.assay.name = "SCT",
  vars.to.regress = "percent.mt",  # 回归线粒体基因
  verbose = FALSE
)

# 2. 查看SCT自动筛选的高变基因（无需用vst重复筛选）
cat("SCT自动筛选的高变基因数量：", length(VariableFeatures(mg, assay = "SCT")), "\n")
top20_var_genes <- head(VariableFeatures(mg, assay = "SCT"), 20)  # 直接用SCT的高变基因

# 3. 绘制高变基因图（基于SCT的高变基因，确保一致性）
VariableFeaturePlot(mg, assay = "SCT") %>%  # 指定assay为SCT
  LabelPoints(points = top20_var_genes, repel = TRUE, size = 3) +
  labs(title = "Top 20 Variable Features (SCT) + All Variable Genes")
ggsave("删除少突污染细胞后的小胶高变基因图mg_variable_features_SCT.png", width = 10, height = 6)

# 4. PCA及elbow图（直接用SCT的高变基因和标准化数据，无需额外ScaleData）
mg <- RunPCA(mg, assay = "SCT", features = VariableFeatures(mg, assay = "SCT"), verbose = FALSE)
ElbowPlot(mg, ndims = 30)

# 5. 聚类与UMAP降维（流程不变，基于SCT的PCA结果）
mg <- FindNeighbors(mg, dims = 1:12)
mg <- FindClusters(mg, resolution = 0.05)
mg <- RunUMAP(mg, dims = 1:12)

# 验证样本来源
unique(mg$orig.ident)
table(mg$seurat_clusters)

DimPlot(mg, label = TRUE, repel = TRUE)

saveRDS(mg, "小胶质细胞去少突后6群1115晚.rds")

# 1. 提取每个簇的Top5标志物（验证功能特异性）
mg_markers <- FindAllMarkers(  # ✅ 改为 mg_markers
  mg, 
  only.pos = TRUE, 
  min.pct = 0.25, 
  logfc.threshold = 0.25, 
  test.use = "wilcox"
)

# 2. 检查结果
cat("找到的marker基因数:", nrow(mg_markers), "\n")
cat("涉及的cluster数:", length(unique(mg_markers$cluster)), "\n")
head(mg_markers)

# 3. 提取Top10标志物（无需修改，变量名已统一）
top10_markers_full <- mg_markers %>%
  filter(p_val_adj < 0.05) %>%  
  group_by(cluster) %>%
  arrange(p_val_adj, desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  ungroup() %>%
  mutate(rank = rep(1:10, length(unique(cluster))))

# 4. 生成汇总表（无需修改）
top10_markers_per_cluster <- top10_markers_full %>%
  group_by(cluster) %>%
  summarise(
    top10_genes = paste(gene, collapse = ", "),
    top10_log2FC = paste(round(avg_log2FC, 2), collapse = ", "),
    top10_padj = paste(format(p_val_adj, scientific = TRUE, digits = 2), collapse = ", ")
  )

print(top10_markers_per_cluster)
# cluster top10_genes                                                                    top10_log2FC top10_padj
# <fct>   <chr>                                                                          <chr>        <chr>     
#   1 0       SYNDIG1, IPCEF1, PLCXD3, MICOS10, RASGEF1C, DLEU7, ZNF710, ADAM28, MAP4K4, DO… 0.86, 0.64,… "1.2e-147…
# 2 1       FKBP5, SLC1A3, SMAP2, B2M, ITM2B, RPS27A, RPL21, IFNGR1, RPS27, RPS23          2.37, 1.62,… " 0.0e+00…
# 3 2       PLPPR1, RUNX1T1, CSMD2, MLLT3, SLC44A5, CDH2, ENOX1, FRMD5, GRIP1, BACH2       3.73, 3.68,… "0e+00, 0…
# 4 3       F13A1, MS4A4E, CR1, COLEC12, CD163, SIGLEC1, IQGAP2, MS4A6A, VAV3, IL15        9.6, 9.45, … " 0.0e+00…
# 5 4       CACHD1, SPON1, COL5A3, LRRC3B, MMD2, ACSBG1, RFX4, SLCO1C1, BMPR1B, HIF3A      7.22, 6.93,… "0e+00, 0…
# 6 5       RGS1, DUSP1, RGS2, SRGN, RHOB, SPATS2L, HSPH1, FOS, HIF1A, ZC3HAV1             5.28, 3.79,… "3.8e-137…

# 5. 保存结果
write.csv(
  top10_markers_full, 
  "去除少突后的小胶6簇_top10_markers.csv", 
  row.names = FALSE,
  quote = FALSE
)


# 按Sample_ID拆分绘图，直接在DimPlot中使用split.by
p <- DimPlot(mg, split.by = "Sample_ID", label = TRUE, repel = TRUE) +
  facet_wrap(~Sample_ID, nrow = 2, ncol = 5) +  # ★核心修改
  labs(title = "Microglia Subclusters by Sample (10 Samples)") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    strip.text = element_text(size = 10),
    axis.line.x = element_line(color = "black", linewidth = 0.5),  # ★添加X轴线
    axis.line.y = element_line(color = "black", linewidth = 0.5)   # Y轴线（可选）
  )

print(p)

# 注释：“健康稳态（簇 0、2、5）→应激启动（簇 1）→抗炎消退（簇 3）→慢性修复（簇 4）”，完整覆盖癫痫 “健康 - 急性  --------
# 1. 生物学原理：基于癫痫 "健康稳态→急性应激→抗炎消退→慢性修复" 的病理进程，将 6 个原始亚群按功能关联性整合，避免亚群过多导致的解读复杂
# 2. Marker 依据：
# - 健康稳态维持群（0/2/5）：SYNDIG1（突触调控）、PLPPR1（血管保护）、RGS1（应激刹车），均为健康脑内小胶质细胞稳态功能标志物
# - 急性应激启动群（1）：FKBP5（应激响应）、SLC1A3（谷氨酸清除）、IFNGR1（免疫预备），对应癫痫急性发作期应激损伤响应
# - 抗炎清除修复群（3）：CD163（抗炎金标）、CR1（炎症清除）、IL15（免疫调节），是炎症消退期核心功能标志物
# - 慢性病理修复群（4）：SPON1（基质重构）、COL5A3（胶质瘢痕）、BMPR1B（慢性修复），仅在癫痫慢性期高表达

# 稳态群（0+2）：健康状态占优  
# 急性激活群（1+5）：发作期应激响应  
# 抗炎修复群（3）：炎症消退期  
# 慢性病理群（4）：慢性期瘢痕形成  

# 检查经典小胶质细胞marker
稳态marker <- c("P2RY12", "TMEM119", "CX3CR1")
激活marker <- c("CD68", "CD74", "APOE", "SPP1")
抗炎marker <- c("CD163", "F13A1", "MS4A4E", "MERTK")

FeaturePlot(mg, features = c(稳态marker, 激活marker, 抗炎marker))

# 检查您提供的稳态marker在这些cluster中的表达  
FeaturePlot(mg, features = c("SYNDIG1", "PLPPR1", "RGS1"))  

# 2. 给Seurat对象添加4类功能群标签（按癫痫病理进程整合）
mg$integrated_cluster <- case_when(
  mg$seurat_clusters %in% c(0, 2, 5) ~ "1_Homeostatic Maintenance Group",
  mg$seurat_clusters == 1 ~ "2_Acute Stress Initiation Group",
  mg$seurat_clusters == 3 ~ "3_Anti-inflammatory Clearance Group",
  mg$seurat_clusters == 4 ~ "4_Chronic Pathological Repair Group",
  TRUE ~ "Others"
)
# 固定标签顺序（避免绘图乱序）
mg$integrated_cluster <- factor(
  mg$integrated_cluster,
  levels = c("1_Homeostatic Maintenance Group", "2_Acute Stress Initiation Group", 
             "3_Anti-inflammatory Clearance Group", "4_Chronic Pathological Repair Group")
)

# 3. 统计各功能群在健康/癫痫组的细胞数量和占比
cluster_sample_stats <- mg@meta.data %>%
  group_by(integrated_cluster, orig.ident) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  mutate(cell_pct = round(cell_count / sum(cell_count) * 100, 2)) %>%
  ungroup()
# 打印统计结果（验证数据）
print("4 functional groups count & percentage (Healthy/Epilepsy):")
print(cluster_sample_stats)

# 4. 绘制健康vs癫痫比例堆叠图
stacked_plot <- ggplot(
  data = cluster_sample_stats,
  aes(x = orig.ident, y = cell_pct, fill = integrated_cluster)
) +
  geom_col(width = 0.6, color = "white", size = 0.5) +
  geom_text(
    aes(label = paste0(cell_pct, "%")),
    position = position_stack(vjust = 0.5),
    size = 4, color = "black", fontface = "bold"
  ) +
  scale_fill_manual(
    values = c(
      "1_Homeostatic Maintenance Group" = "#2E86AB",
      "2_Acute Stress Initiation Group" = "#F18F01",
      "3_Anti-inflammatory Clearance Group" = "#C73E1D",
      "4_Chronic Pathological Repair Group" = "#6A994E"
    ),
    name = "4 Functional Groups of Microglia"
  ) +
  labs(
    x = "Sample Group",
    y = "Cell Percentage (%)",
    title = "Proportion of 4 Microglial Functional Groups (Epilepsy vs Healthy)",
    subtitle = "Integrated from 6 Original Clusters"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 11, color = "#666666"),
    axis.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 13, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid = element_blank(),
    panel.border = element_rect(size = 0.8, color = "black")
  )
# 显示图表
print(stacked_plot)
table(mg$integrated_cluster)
# 1_Homeostatic Maintenance Group     2_Acute Stress Initiation Group 3_Anti-inflammatory Clearance Group 
# 3579                                1700                                 130 
# 4_Chronic Pathological Repair Group 
# 125 
# 5. 保存图表（默认保存在当前工作目录，可修改filename路径）
ggsave(
  filename = "小胶质细胞分6亚簇4亚群.png",
  plot = stacked_plot,
  width = 10,
  height = 7,
  dpi = 300,
  bg = "white"
)

# 6. 保存带功能群标签的Seurat对象（默认保存在当前工作目录）
saveRDS(mg, "Microglia_亚分4类.rds")

head(mg)


# 加载数据与定义marker
mg <- readRDS("Microglia_亚分4类.rds")
markers <- c("SYNDIG1", "FKBP5", "CD163", "SPON1")

# 仅生成一次子图列表
plot_list <- lapply(markers, function(g) {
  FeaturePlot(mg, features = g, pt.size = 0.4, cols = c("lightgrey", "#C73E1D")) +
    labs(title = g) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
          axis.text = element_blank(), axis.title = element_blank(),
          axis.ticks = element_blank(), panel.grid = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, size = 0.6))
})

# 组合图表并保存（无重复代码）
p <- wrap_plots(plot_list, ncol = 2, nrow = 2) +
  plot_annotation(title = "Marker Expression in Microglial Functional Subgroups",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14)))
p

ggsave("Microglia_Markers_2x2_FeaturePlot.png", width = 8, height = 8, dpi = 300, bg = "white")

VlnPlot(
  mg, 
  features = c("P2RY12", "DUSP1", "CD163", "COL5A3"), 
  group.by = "integrated_cluster", 
  ncol = 2
) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
