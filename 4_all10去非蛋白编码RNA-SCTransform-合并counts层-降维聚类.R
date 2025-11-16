## 筛选人类蛋白编码基因 → 合并对象 → 修复 RNA assay 的 counts 层（唯一结构修复步骤） → 所有分析步骤（高变基因、PCA、批次校正、聚类）
rm(list = ls())

library(Seurat)
library(dplyr)
library(harmony)
library(ggplot2)
library(patchwork)
library(biomaRt)

# 1. 读取数据并处理样本标识 ------------------------------------------------
health4 <- readRDS("/Users/lijiangbo/1_Projects/Epi1106/5_scRNA/2_processed_data/1113/00_processed_data/3_healthy4_hq.rds")
epi6 <- readRDS("/Users/lijiangbo/1_Projects/Epi1106/5_scRNA/2_processed_data/1113/00_processed_data/3_epi6_hq.rds")
head(health4)
head(epi6)
epi6$Sample_Name <- epi6$orig.ident  # 新增一列存储原始样本名（如P010_frontal）
epi6$orig.ident <- "Epilepsy" # 将epi6的orig.ident统一改为"Epilepsy"（与健康组的"Healthy"对应）
unique(health4$orig.ident)
unique(epi6$orig.ident)
nrow(health4)
head(rownames(health4), 10)
ncol(health4)
nrow(epi6)
head(rownames(epi6), 10)
ncol(epi6)

# 2. 筛选人类蛋白编码基因 ---------------------------------------------------
common_genes <- intersect(rownames(health4), rownames(epi6)) # 29543
length(common_genes)
round(length(common_genes)/nrow(health4)*100, 2)
round(length(common_genes)/nrow(epi6)*100, 2)

# 获取人类蛋白编码基因列表（用biomaRt，确保准确性）
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
protein_coding_genes <- getBM(
  attributes = "external_gene_name",
  filters = "biotype",
  values = "protein_coding",
  mart = ensembl
)$external_gene_name

# （2）健康组筛选蛋白编码基因
health4 <- subset(health4, features = intersect(protein_coding_genes, rownames(health4)))

# （3）癫痫组筛选蛋白编码基因
epi6 <- subset(epi6, features = intersect(protein_coding_genes, rownames(epi6)))

# 验证筛选结果（确保两组基因集基本一致）
cat("健康组筛选后基因数：", nrow(health4), "\n") # 16910
cat("癫痫组筛选后基因数：", nrow(epi6), "\n") # 17573
cat("共同基因比例：", round(length(intersect(rownames(health4), rownames(epi6)))/nrow(health4)*100, 2), "%\n")  # 99.89 %

# 3. SCT标准化 + 合并对象 ---------------------------------------------------
# 对筛选后的样本分别进行SCTransform标准化
health4 <- SCTransform(health4, verbose = FALSE)  # 基于蛋白编码基因标准化
epi6 <- SCTransform(epi6, verbose = FALSE)       # 基于蛋白编码基因标准化

# 合并对象
combined_obj <- merge(
  x = health4,
  y = epi6,
  add.cell.ids = c("Healthy", "Epilepsy"),
  project = "Healthy_vs_Epilepsy"  
)

# 验证合并后的默认assay是否为SCT（关键！）
cat("合并后默认assay：", DefaultAssay(combined_obj), "\n")  # "SCT"

# 整合样本信息到Sample_ID列
combined_obj@meta.data$Sample_ID <- ifelse(
  !is.na(combined_obj$Donor.ID),  # 若Donor.ID非NA（健康组），用Donor.ID
  combined_obj$Donor.ID,
  combined_obj$Sample_Name        # 若Donor.ID为NA（癫痫组），用Sample_Name
)
unique(combined_obj$Sample_ID)


# 合并RNA assay的counts层 

# 4. 合并RNA assay的counts层（彻底修复数据结构） ------------------------
# 4.1 查看RNA assay的拆分layer
cat("combined_obj 的 RNA assay 所有 layer：\n")
print(Layers(combined_obj[["RNA"]]))
# [1] "counts.Healthy4" 等7个拆分层

# 4.2 查看各layer详情
if (length(Layers(combined_obj[["RNA"]])) > 1) {
  layer_info <- data.frame()
  for (layer in Layers(combined_obj[["RNA"]])) {
    mat <- GetAssayData(combined_obj, assay = "RNA", layer = layer)
    info <- data.frame(
      Layer = layer,
      Gene_Count = nrow(mat),
      Cell_Count = ncol(mat),
      First_Cell = colnames(mat)[1]
    )
    layer_info <- rbind(layer_info, info)
  }
  cat("\n各 layer 详细信息：\n")
  print(layer_info)
}

# 4.3 合并为统一counts层
all_unique_cells <- colnames(combined_obj)  # 全局唯一细胞（40761个）
counts_layers <- Layers(combined_obj[["RNA"]])[grep("counts", Layers(combined_obj[["RNA"]]))]

# 提取矩阵→过滤重复→合并
mat_list <- lapply(counts_layers, function(layer) {
  mat <- GetAssayData(combined_obj, assay = "RNA", layer = layer)
  mat[, colnames(mat) %in% all_unique_cells, drop = FALSE]
})
common_genes <- Reduce(intersect, lapply(mat_list, rownames))  # 16305个共同基因
mat_list_filtered <- lapply(mat_list, function(mat) mat[common_genes, , drop = FALSE])
merged_counts_clean <- do.call(cbind, mat_list_filtered)

# 替换RNA assay（彻底清除旧层）
combined_obj[["RNA"]] <- CreateAssayObject(counts = merged_counts_clean)

# 验证合并结果（确保结构正确后，再往下走）
cat("===== counts层合并验证（必须通过！）=====\n")
cat("RNA assay 的所有 layer：\n")
print(Layers(combined_obj[["RNA"]]))  # 应输出 [1] "counts"
cat("RNA assay counts 层细胞数：", ncol(combined_obj[["RNA"]]$counts), "\n")  # 40761
cat("RNA assay counts 层基因数：", nrow(combined_obj[["RNA"]]$counts), "\n")  # 16305
cat("细胞匹配：", all(colnames(combined_obj[["RNA"]]$counts) == colnames(combined_obj)), "\n")  # TRUE


# 【数据结构修复完成，再开展所有分析步骤】

# 5. 高变基因筛选→PCA→批次校正→降维→聚类 ----------------------------------
# 5.1 筛选高变基因（基于SCT assay，不受RNA合并影响）
combined_obj <- FindVariableFeatures(combined_obj, assay = "SCT", verbose = FALSE)  # 筛选SCT后的高可变基因
hvf_count <- length(VariableFeatures(combined_obj, assay = "SCT"))
cat("基于SCT assay筛选出的高变基因数量：", hvf_count, "\n")  # 2000 
head(VariableFeatures(combined_obj, assay = "SCT"), 10)  # 前10个高变基因：LHFPL3等

# 5.2 PCA降维（基于SCT高变基因）
combined_obj <- RunPCA(
  combined_obj,
  assay = "SCT",               # 基于SCT标准化数据
  features = VariableFeatures(combined_obj),  # 用高可变基因做PCA
  npcs = 30,                   # 保留30个主成分
  verbose = FALSE
)
ElbowPlot(combined_obj, ndims = 30)  # 查看PCA效果

# 5.3 Harmony批次校正
combined_obj <- RunHarmony(
  object = combined_obj,
  group.by.vars = "Sample_ID",  # 按最细的技术批次（sample_ID）校正
  reduction.use = "pca",
  dims = 1:30,
  verbose = FALSE
)

# 验证校正效果
p_before <- DimPlot(combined_obj, group.by = "Sample_ID", reduction = "pca") + 
  ggtitle("Sample Distribution (Before Harmony, PCA)") + theme(plot.title = element_text(hjust = 0.5))
p_after <- DimPlot(combined_obj, group.by = "Sample_ID", reduction = "harmony") + 
  ggtitle("Sample Distribution (After Harmony)") + theme(plot.title = element_text(hjust = 0.5))
p_before | p_after

# 5.4 UMAP降维（基于校正后的Harmony数据）
combined_obj <- RunUMAP(
  combined_obj,
  reduction = "harmony",
  dims = 1:30,
  verbose = FALSE
)

# 5.5 细胞聚类
combined_obj <- FindNeighbors(
  combined_obj,
  reduction = "harmony",
  dims = 1:30,
  verbose = FALSE
)
combined_obj <- FindClusters(
  combined_obj,
  resolution = 0.1,
  verbose = FALSE
)

# 查看聚类结果
unique(combined_obj$seurat_clusters)
# [1] 3  6  5  7  4  8  10 1  2  9  11 0  12
# Levels: 0 1 2 3 4 5 6 7 8 9 10 11 12

# 可视化聚类结果
p_cluster <- DimPlot(combined_obj, group.by = "seurat_clusters", reduction = "umap", label = TRUE) +
  ggtitle("Cell Clusters (UMAP)") +
  theme(plot.title = element_text(hjust = 0.5))
print(p_cluster)


p4 <- DimPlot(combined_obj, group.by = "Sample_ID", reduction = "umap", 
              split.by = "orig.ident") +  # 按健康/癫痫分组显示
  ggtitle("Sample Distribution in UMAP (Split by Group)")
print(p4)

# 6. 确认默认assay（后续分析用SCT） -----------------------------------------
DefaultAssay(combined_obj) <- "SCT"
cat("最终默认assay：", DefaultAssay(combined_obj), "\n")  # 应返回 "SCT"

saveRDS(combined_obj, "/Users/lijiangbo/1_Projects/Epi1106/5_scRNA/2_processed_data/1113/00_processed_data/4_combined_10.rds")

# 7. 另存一份去除P025的对象------------------------------------------------------------------
unique(combined_obj$Sample_ID)
combined_obj_filtered <- subset(
  x = combined_obj,
  subset = Sample_ID != "P025"  # 核心条件：排除P025样本
)
# 1. 检查剩余样本ID中是否已无P025
cat("删除后剩余样本ID：\n")
print(unique(combined_obj_filtered$Sample_ID))  # 应无"P025"

# 2. 检查细胞数是否减少（P025原始细胞数为3638，总细胞数应减少3638）
cat("删除前总细胞数：", ncol(combined_obj), "\n")          # 原40761
cat("删除后总细胞数：", ncol(combined_obj_filtered), "\n")  # 应为40761 - 3638 = 37123

# 3. 检查数据结构是否正常（RNA和SCT assay的layer未被破坏）
cat("删除后RNA assay的layer：\n")
print(Layers(combined_obj_filtered[["RNA"]]))  # 仍为["counts"]
cat("删除后SCT assay的layer：\n")
print(Layers(combined_obj_filtered[["SCT"]]))  # 仍为["counts", "data", "scale.data"]
DefaultAssay(combined_obj_filtered)  # [1] "SCT"

saveRDS(
  combined_obj_filtered,
  "/Users/lijiangbo/1_Projects/Epi1106/5_scRNA/2_processed_data/1113/00_processed_data/4_combined_noP025.rds"
)


# 8. 存能后续用python处理的文件格式Seurat → h5ad（AnnData 格式） -----------------------------------------------
# 1. 安装并加载SeuratDisk包（若未安装）
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("mojaveazure/seurat-disk")
library(SeuratDisk)

# 2. 加载已处理的combined_obj（确保前文代码已运行并保存对象）
# 若已保存为RDS，可先读取：
# combined_obj <- readRDS("/Users/lijiangbo/1_Projects/Epi1106/5_scRNA/2_processed_data/1113/00_processed_data/4_combined_10.rds")

# 3. 保存为h5Seurat格式，再转换为h5ad
SaveH5Seurat(
  object = combined_obj, 
  filename = "/Users/lijiangbo/1_Projects/Epi1106/5_scRNA/2_processed_data/1113/00_processed_data/combined_obj.h5Seurat"
)
Convert(
  source = "/Users/lijiangbo/1_Projects/Epi1106/5_scRNA/2_processed_data/1113/00_processed_data/combined_obj.h5Seurat", 
  dest = "h5ad"  # 转换为Python兼容的h5ad格式
)




