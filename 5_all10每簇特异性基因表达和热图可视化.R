rm(list = ls())
gc()

# 加载必需包
pacman::p_load(Seurat, dplyr, ggplot2, cowplot, SingleR, BiocManager)
library(Seurat)  
library(dplyr)  
library(ggplot2)  
library(patchwork)  
library(RColorBrewer)  
library(scales)  

setwd("/Users/lijiangbo/1_Projects/Epi1106/5_scRNA/1113重新分析/1113outputs/")  


# 1. 加载处理后Seurat对象
all10 <- readRDS("/Users/lijiangbo/1_Projects/Epi1106/5_scRNA/2_processed_data/1113/00_processed_data/4_combined_10.rds")
head(all10)


# 解决基因索引不一致导致的越界访问 --------------------------------------------------------
# 问题链条：
# 整合7个样本 → 每个样本的基因数不同 → merge后基因是并集(16874) 
# → 每个SCT模型只有子集基因(15381-16673) → PrepSCTFindMarkers访问时索引越界
cat("对象总基因数:", nrow(all10), " | 总细胞数:", ncol(all10), "\n")
cat("SCT模型数量:", length(all10@assays$SCT@SCTModel.list), "\n\n")
model_stats <- data.frame(
  Model = 1:length(all10@assays$SCT@SCTModel.list),
  Genes = sapply(all10@assays$SCT@SCTModel.list, function(m) 
    nrow(m@feature.attributes)),
  Cells = sapply(all10@assays$SCT@SCTModel.list, function(m) 
    nrow(m@cell.attributes))
)
print(model_stats)
#           Model Genes Cells
# model1       1 14851 16889
# model1.6     2 14851  3783
# model1.1     3 14851  3189
# model1.2     4 14851  3136
# model1.3     5 14851  4693
# model1.4     6 14851  5433
# model1.5     7 14851  3638
cat("细胞总数验证:", sum(model_stats$Cells), "\n\n")  # 40761 


cat("=== SCT Assay组件 ===\n")
cat("counts:", dim(all10@assays$SCT@counts), "\n")
cat("data:", dim(all10@assays$SCT@data), "\n")
cat("meta.features:", nrow(all10@assays$SCT@meta.features), "行\n")
cat("var.features:", length(all10@assays$SCT@var.features), "个\n\n")

counts_genes <- rownames(all10@assays$SCT@counts)
model_genes_list <- lapply(all10@assays$SCT@SCTModel.list, function(m) 
  rownames(m@feature.attributes))

mismatch_summary <- sapply(1:length(model_genes_list), function(i) {
  model_genes <- model_genes_list[[i]]
  data.frame(
    Model = i,
    ModelGenes = length(model_genes),
    MissingInCounts = sum(!model_genes %in% counts_genes),
    ExtraInCounts = sum(!counts_genes %in% model_genes)
  )
})
cat("=== 基因不匹配统计 ===\n")
print(t(mismatch_summary))

if (any(mismatch_summary["MissingInCounts",] > 0) || 
    any(mismatch_summary["ExtraInCounts",] > 0)) {
  
  cat("\n⚠️  检测到基因不一致，开始修复...\n")
  
  # 2.1 计算所有模型的共有基因
  common_genes <- Reduce(intersect, model_genes_list)
  cat("所有模型共有基因数:", length(common_genes), "\n")
  
  # 2.2 同步SCT assay的所有矩阵
  all10@assays$SCT@counts <- all10@assays$SCT@counts[common_genes, ]
  all10@assays$SCT@data <- all10@assays$SCT@data[common_genes, ]
  all10@assays$SCT@meta.features <- all10@assays$SCT@meta.features[common_genes, , drop = FALSE]
  
  # 2.3 同步var.features
  all10@assays$SCT@var.features <- intersect(
    all10@assays$SCT@var.features, 
    common_genes
  )
  
  # 2.4 同步scale.data（如果存在）
  if (nrow(all10@assays$SCT@scale.data) > 0) {
    scale_genes <- intersect(rownames(all10@assays$SCT@scale.data), common_genes)
    all10@assays$SCT@scale.data <- all10@assays$SCT@scale.data[scale_genes, ]
  }
  
  # 2.5 同步所有SCT模型
  all10@assays$SCT@SCTModel.list <- lapply(
    all10@assays$SCT@SCTModel.list,
    function(model) {
      model@feature.attributes <- model@feature.attributes[common_genes, , drop = FALSE]
      
      # 同步细胞（通常已对齐，但保险起见）
      keep_cells <- intersect(rownames(model@cell.attributes), colnames(all10))
      model@cell.attributes <- model@cell.attributes[keep_cells, , drop = FALSE]
      
      return(model)
    }
  )
  
  cat("✅ 修复完成！\n\n")
  
} else {
  cat("\n✅ 未检测到基因不一致，无需修复\n\n")
}


# ===== 步骤3：验证修复结果 =====
cat("=== 修复后验证 ===\n")
cat("SCT counts:", nrow(all10@assays$SCT@counts), "基因\n")
cat("SCT data:", nrow(all10@assays$SCT@data), "基因\n")
cat("meta.features:", nrow(all10@assays$SCT@meta.features), "行\n")
cat("var.features:", length(all10@assays$SCT@var.features), "个\n")

for (i in seq_along(all10@assays$SCT@SCTModel.list)) {
  cat(sprintf("模型%d: %d基因\n", i, 
              nrow(all10@assays$SCT@SCTModel.list[[i]]@feature.attributes)))
}

cluster_counts <- table(Idents(all10))
cluster_counts
#     0     1     2     3     4     5     6     7     8     9    10    11    12 
# 10539  6127  5662  4646  4339  3733  1921  1631   599   584   478   440    62 

cluster_prop <- prop.table(cluster_counts) * 100
print("\n各Cluster占比（%）：")
print(round(cluster_prop, 2))
#    0     1     2     3     4     5     6     7     8     9    10    11    12 
# 25.86 15.03 13.89 11.40 10.64  9.16  4.71  4.00  1.47  1.43  1.17  1.08  0.15 

cluster_stats <- data.frame(
  Cluster = names(cluster_counts),
  CellCount = as.vector(cluster_counts),
  Percentage = round(as.vector(cluster_prop), 2)
)
write.csv(cluster_stats, "/Users/lijiangbo/1_Projects/Epi1106/5_scRNA/1113重新分析/1113outputs/5_cluster_cell_counts.csv", row.names = FALSE)

# ===== 步骤4：执行PrepSCTFindMarkers =====
cat("\n=== 执行PrepSCTFindMarkers ===\n")
all10 <- PrepSCTFindMarkers(all10)
cat("✅ 成功！对象已准备好进行差异分析\n")


DefaultAssay(all10) <- "SCT"
cat("当前默认assay:", DefaultAssay(all10), "\n")


cat("开始鉴定差异基因...\n")
diff.wilcox <- FindAllMarkers(
  object = all10,
  assay = "SCT",
  test.use = "wilcox",
  only.pos = TRUE, # 只保留目标簇中上调的基因（默认，符合marker基因特征）
  logfc.threshold = 1, # 差异倍数阈值（可根据需求提高，如0.5）
  min.pct = 0.1, # 基因表达比例阈值（过滤低表达基因）
  verbose = TRUE
)
# 保存所有差异基因
write.csv(
  diff.wilcox, 
  file = "/Users/lijiangbo/1_Projects/Epi1106/5_scRNA/1113重新分析/1113outputs/5_diff.wilcox.csv",
  row.names = FALSE
)
cat("总共发现", nrow(diff.wilcox), "个差异基因\n")
cat("涉及", length(unique(diff.wilcox$cluster)), "个cluster\n\n")


# 第一轮：严格筛选
specific_markers_strict <- diff.wilcox %>%
  dplyr::filter(
    p_val_adj < 0.01, # 使用校正后的p值更严谨
    pct.1 > 0.6,   # 目标簇表达比例>60%
    pct.2 < 0.1,   # 其他簇表达比例<10%
    avg_log2FC > 1.5   # log2FC>1 即2倍差异（若需2.83倍则用1.5）
  )

# 统计各cluster通过数量
cluster_counts_strict <- table(specific_markers_strict$cluster)
cat("严格标准下各cluster基因数：\n")
print(cluster_counts_strict)
# 0  1  2  3  4  5  6  7  8  9 10 11 12 
# 38 35 60  4 13  0  4  3 16  0 37  7  3 
# 第二轮：对基因数<5的cluster放宽标准
clusters_need_relaxed <- names(cluster_counts_strict[cluster_counts_strict < 5])
clusters_need_relaxed <- c(
  clusters_need_relaxed,
  setdiff(as.character(0:12), names(cluster_counts_strict))  # 包括0个基因的cluster
)

cat("\n需要放宽标准的cluster:", clusters_need_relaxed, "\n")
# 需要放宽标准的cluster: 3 5 6 7 9 12 
specific_markers_relaxed <- diff.wilcox %>%
  dplyr::filter(
    cluster %in% clusters_need_relaxed,
    p_val_adj < 0.05,      # 放宽p值
    pct.1 > 0.5,           # 降低表达比例
    pct.2 < 0.15,           # 放宽背景表达
    avg_log2FC > 0.8       # 降低差异倍数（1.4倍）
  )

# 合并两轮结果
specific_markers_combined <- rbind(
  specific_markers_strict %>% 
    dplyr::filter(!cluster %in% clusters_need_relaxed),  # 保留严格标准的cluster
  specific_markers_relaxed  # 添加放宽标准的cluster
)

cat("\n放宽后各cluster基因数：\n")
print(table(specific_markers_combined$cluster))

# 保存特异性marker
write.csv(
  specific_markers_combined, 
  file = "/Users/lijiangbo/1_Projects/Epi1106/5_scRNA/1113重新分析/1113outputs/5_cluster_specific_markers.csv",
  row.names = FALSE
)


p1 <- DimPlot(all10,   
              reduction = "umap",  
              label = TRUE)
p1
ggsave("/Users/lijiangbo/1_Projects/Epi1106/5_scRNA/1113重新分析/1113outputs/5_UMAP_plot.png", plot = p1, width = 10, height = 8, dpi = 300)  
           
# 4.1 先从严格筛选的结果中取top5
cluster_specific_top5 <- specific_markers_combined %>%
  dplyr::group_by(cluster) %>%
  dplyr::arrange(desc(avg_log2FC)) %>%
  dplyr::slice_head(n = 5) %>%
  dplyr::ungroup() %>%
  split(.$cluster)
cluster_specific_top5_genes <- lapply(cluster_specific_top5, function(x) x$gene)
sorted_specific_genes <- unlist(
  cluster_specific_top5_genes[order(as.numeric(names(cluster_specific_top5_genes)))]
)
cat("\n最终用于可视化的基因数:", length(sorted_specific_genes), "\n")
cat("各cluster基因数：\n")
print(sapply(cluster_specific_top5_genes, length))

cluster_specific_top5_genes
# `0` = c("OPALIN", "FOLH1", "MOBP", "CNDP1", "ST18"),        
# # 成熟少突胶质细胞 (Mature Oligodendrocytes)
# # 髓鞘化轴突的主要细胞，提高神经传导速度
# 
# `1` = c("ALDH1A1", "AQP4", "SLCO1C1", "GJA1", "RANBP3L"),   
# # 星形胶质细胞 (Astrocytes)
# # 代谢支持、神经递质循环、血脑屏障维持
# 
# `2` = c("CX3CR1", "ADAM28", "P2RY12", "APBB1IP", "CSF2RA"), 
# # 小胶质细胞 (Microglia)
# # 脑内常驻免疫细胞，免疫监视和突触修剪
# 
# `3` = c("NWD2", "CCBE1", "SV2B", "PCSK2", "ZNF804B"),       
# # 深层兴奋性神经元 (Deep Layer Excitatory Neurons)
# # 谷氨酸能神经元，可能位于皮层第5/6层，长程投射
# 
# `4` = c("PDGFRA", "MYT1", "COL9A1", "STK32A", "MEGF11"),    
# # 少突胶质前体细胞 (Oligodendrocyte Precursor Cells, OPCs)
# # PDGFRA+，成体脑中唯一保持增殖能力的胶质细胞
# 
# `5` = c("SOX11", "EPHA3", "CUX2", "DPYSL3"),                
# # 上层兴奋性神经元 (Upper Layer Excitatory Neurons, Layer 2/3)
# # CUX2+标志上层皮层神经元，参与皮层内信息处理
# 
# `6` = c("GRIN3A", "KCNC2", "ANK1", "CDH9", "GRIP2"),        
# # PV+快速放电抑制性中间神经元 (PV+ Fast-spiking Inhibitory Interneurons)
# # KCNC2 (Kv3.2)+，快速无适应放电(>100Hz)，释放GABA，篮状细胞或枝晶细胞
# 
# `7` = c("CXCL14", "SLC35F4", "SYNPR", "GAD2", "CNR1"),      
# # SST+/VIP+ GABA能抑制性神经元 (SST+/VIP+ GABAergic Inhibitory Interneurons)
# # GAD2+确认GABA能，可能包含SST+或VIP+亚型，树突抑制或去抑制功能
# 
# `8` = c("ITIH5", "EBF1", "DCN", "COLEC12", "RBPMS"),        
# # 周细胞/血管周围成纤维细胞 (Pericytes/Perivascular Fibroblasts)
# # 包裹血管，维持血脑屏障完整性，调节脑血流
# 
# `9` = c("RIPOR2", "ARHGAP15", "RPL18A", "PTPRC", "RPS4X"),  
# # 外周免疫细胞/双细胞 (Peripheral Immune Cells/Doublets) ⚠️
# # PTPRC (CD45)+，高核糖体蛋白，可能是污染或双细胞，建议过滤
# 
# `10` = c("MECOM", "FLT1", "CLDN5", "ATP10A", "ABCB1"),      
# # 血管内皮细胞 (Vascular Endothelial Cells)
# # CLDN5+紧密连接，FLT1+血管生成，构成血脑屏障
# 
# `11` = c("GPR17", "KCNS3", "SLC16A10", "SEMA5B", "PLAAT1"), 
# # 新生/分化中少突胶质细胞 (Newly Formed Oligodendrocytes)
# # GPR17+是OPC向成熟少突胶质细胞过渡的分子开关
# 
# `12` = c("FERMT1", "NTN1", "PDGFRA", "MYT1", "COL9A1")      
# # 少突胶质前体细胞亚群 (OPCs Subpopulation)
# # 与Cluster 4共享PDGFRA，可能是不同脑区或分化阶段的OPC亚群

VlnPlot(all10, features = c("PTPRC", "percent.mt", "nFeature_RNA"), 
        idents = c(1, 2, 9), ncol = 3)


# 删除cluster9 --------------------------------------------------------------

cluster9_cells <- sum(Idents(all10) == 9)
total_cells <- ncol(all10)
cat("\n=== Cluster 9 删除统计 ===\n")
cat("Cluster 9 细胞数:", cluster9_cells, "\n")
cat("占总细胞比例:", round(cluster9_cells/total_cells*100, 2), "%\n")
cat("删除原因：\n")
cat("  1. PTPRC高表达 - 外周免疫细胞污染\n")
cat("  2. 线粒体比例偏高 - 细胞质量差\n")
cat("  3. 检测基因数偏低 - 信息含量不足\n")
cat("  4. 无明确生物学意义\n\n")

# 3. 执行删除
all10 <- subset(all10, idents = setdiff(0:12, 9))
cat("删除后细胞数:", ncol(all10), "\n")
cat("保留的cluster:", paste(sort(unique(Idents(all10))), collapse = ", "), "\n\n")
rm(all10_cleaned)
table(Idents(all10))
# 0     1     2     3     4     5     6     7     8    10    11    12 
# 10539  6127  5662  4646  4339  3733  1921  1631   599   478   440    62
# 删除Cluster 12
all10 <- subset(all10, idents = setdiff(c(0:8, 10:12), 12))
table(Idents(all10))

saveRDS(all10, "/Users/lijiangbo/1_Projects/Epi1106/5_scRNA/2_processed_data/1113/00_processed_data/5_combined_10_去掉簇9和簇12后.rds")


# 细胞类型注释 ------------------------------------------------------------------
# 定义8大类细胞的marker基因
markers_8class <- list(
  Astrocytes = c("AQP4", "ALDH1A1", "GJA1", "GFAP", "SLC1A2", "SLCO1C1"),
  Microglia = c("P2RY12", "CX3CR1", "TMEM119", "C1QA", "C1QB", "ADAM28"),
  Oligodendrocytes = c("MOBP", "OPALIN", "MOG", "PLP1", "MBP", "MAG"),
  OPCs = c("PDGFRA", "CSPG4", "VCAN", "COL9A1", "GPR17", "MYT1"),
  Excitatory_Neurons = c("SLC17A7", "CAMK2A", "CUX2", "SATB2", "TBR1", "NWD2"),
  Inhibitory_Neurons = c("GAD1", "GAD2", "SLC32A1", "KCNC2", "SST", "VIP", "PVALB"),
  Endothelial_Cells = c("CLDN5", "FLT1", "PECAM1", "VWF", "CDH5", "MECOM"),
  Pericytes = c("PDGFRB", "RGS5", "ANPEP", "DCN", "ITIH5", "COLEC12")
)

cat("\nMarker基因可用性检查：\n")
markers_available <- lapply(markers_8class, function(x) {
  present <- x[x %in% rownames(all10)]
  cat("  ", names(markers_8class)[which(sapply(markers_8class, identical, x))], ":", 
      length(present), "/", length(x), "\n")
  return(present)
})
all_markers <- unique(unlist(markers_available))

library(Seurat)
library(dplyr)

# 1. 确认原始cluster已保存
all10$original_cluster <- as.character(Idents(all10))

cat("原始cluster分布：\n")
table(all10$original_cluster)

# 2. 创建映射向量
celltype_8class_map <- c(
  "0" = "Oligodendrocytes",
  "1" = "Astrocytes",
  "2" = "Microglia",
  "3" = "Excitatory_Neurons",
  "4" = "OPCs",
  "5" = "Excitatory_Neurons",
  "6" = "Inhibitory_Neurons",
  "7" = "Inhibitory_Neurons",
  "8" = "Pericytes",
  "10" = "Endothelial_Cells",
  "11" = "Oligodendrocytes"
)

celltype_8class_cn_map <- c(
  "0" = "少突胶质细胞",
  "1" = "星形胶质细胞",
  "2" = "小胶质细胞",
  "3" = "兴奋性神经元",
  "4" = "少突胶质前体细胞",
  "5" = "兴奋性神经元",
  "6" = "抑制性神经元",
  "7" = "抑制性神经元",
  "8" = "周细胞",
  "10" = "内皮细胞",
  "11" = "少突胶质细胞"
)

# 3. 正确的映射方式：使用 unname() 去掉多余的names
all10$cell_type_8class <- unname(celltype_8class_map[all10$original_cluster])
all10$cell_type_8class_cn <- unname(celltype_8class_cn_map[all10$original_cluster])

# 4. 验证
cat("\n✅ 映射成功！\n\n")

cat("8大类细胞分布（英文）：\n")
table(all10$cell_type_8class)
# Astrocytes  Endothelial_Cells Excitatory_Neurons Inhibitory_Neurons 
# 6127                478               8379               3552 
# Microglia   Oligodendrocytes               OPCs          Pericytes 
# 5662              10979               4339                599 

cat("\n8大类细胞分布（中文）：\n")
table(all10$cell_type_8class_cn)

# 5. 查看前几行
cat("\n前10个细胞的映射结果：\n")
head(all10@meta.data[, c("original_cluster", "cell_type_8class", "cell_type_8class_cn")], 10)

# 兴奋性神经元         内皮细胞           周细胞       小胶质细胞 少突胶质前体细胞 
# 8379              478              599             5662             4339 
# 少突胶质细胞     抑制性神经元     星形胶质细胞 
# 10979             3552             6127 

Idents(all10) <- "cell_type_8class" 
unique(Idents(all10))
key_markers <- c(  
  "AQP4",      # Astrocytes  
  "P2RY12",    # Microglia  
  "MOBP",      # Oligodendrocytes  
  "PDGFRA",    # OPCs  
  "SLC17A7",   # Excitatory Neurons (如果不可用则用CAMK2A)  
  "GAD2",      # Inhibitory Neurons  
  "CLDN5",     # Endothelial Cells  
  "PDGFRB"     # Pericytes  
)  

# 检查可用性并在需要时替换  
key_markers_present <- key_markers[key_markers %in% rownames(all10)]  

head(all10)

# 基于您的DotPlot结果，挑选最特异的marker
top_markers_list <- list(
  Astrocytes = c("AQP4", "GFAP", "SLCO1C1"),           # 星形胶质细胞
  Microglia = c("P2RY12", "CX3CR1", "C1QB"),        # 小胶质细胞
  Oligodendrocytes = c("MOBP", "MOG", "PLP1"),         # 少突胶质细胞
  OPCs = c("PDGFRA", "COL9A1", "VCAN"),                 # 少突胶质前体细胞
  Excitatory_Neurons = c("SLC17A7", "CUX2", "SATB2"), # 兴奋性神经元
  Inhibitory_Neurons = c("GAD1", "GAD2", "KCNC2"),       # 抑制性神经元
  Endothelial_Cells = c("CLDN5", "FLT1", "MECOM"),    # 内皮细胞
  Pericytes = c("PDGFRB", "DCN", "COLEC12")               # 周细胞
)
# 显示选择的marker
for(celltype in names(top_markers_list)) {
  cat(sprintf("%-25s: %s\n", celltype, 
              paste(top_markers_list[[celltype]], collapse = ", ")))
}
cat("\n")
top_markers_available <- lapply(names(top_markers_list), function(celltype) {
  markers <- top_markers_list[[celltype]]
  present <- markers[markers %in% rownames(all10)]
  missing <- markers[!markers %in% rownames(all10)]
  
  cat(sprintf("%-25s: %d/%d 可用", 
              celltype, length(present), length(markers)))
  
  if(length(missing) > 0) {
    cat(sprintf("  (缺失: %s)", paste(missing, collapse = ", ")))
  }
  cat("\n")
  
  return(present)
})
names(top_markers_available) <- names(top_markers_list)

top_markers_final <- unique(unlist(top_markers_available))
# 设置Identity
Idents(all10) <- "cell_type_8class"

# 细胞类型顺序
cell_type_order <- c("Astrocytes", "Microglia", "Oligodendrocytes", "OPCs",
                     "Excitatory_Neurons", "Inhibitory_Neurons", 
                     "Endothelial_Cells", "Pericytes")
p_dotplot_top3 <- DotPlot(all10, 
                          features = top_markers_final,
                          cols = c("lightgrey", "red"),
                          dot.scale = 6,  # 增大点的尺寸
                          cluster.idents = FALSE) +
  RotatedAxis() +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.5, 
                               face = "plain", color = "black"),
    axis.text.y = element_text(size = 8, face = "bold", color = "black"),
    axis.title = element_text(size = 10, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 8),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.line = element_blank(),
    legend.position = "right"
  ) +
  labs(title = "Top 3 Specific Markers for Each Cell Type",
       x = "Marker Genes",
       y = "Cell Types") +
  scale_y_discrete(limits = rev(cell_type_order))

# 显示图形
print(p_dotplot_top3)

# 保存图形
ggsave("5_Fig1_DotPlot_Top3_Markers.pdf", plot = p_dotplot_top3, 
       width = 10, height = 6, dpi = 300)
ggsave("5_Fig1_DotPlot_Top3_Markers.png", plot = p_dotplot_top3, 
       width = 10, height = 6, dpi = 300)

cat("\n✓ 精简版 DotPlot 已保存\n")
cat("  - Fig1_DotPlot_Top3_Markers.pdf\n")
cat("  - Fig1_DotPlot_Top3_Markers.png\n\n")







