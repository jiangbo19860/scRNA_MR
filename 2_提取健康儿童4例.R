library(Seurat)
library(dplyr)

health12 <- readRDS("/Users/lijiangbo/1_Projects/Epi1106/5_scRNA/2_processed_data/1113/healthy12_seurat.rds")
head(health12)
unique(health12$Donor.ID)

# 1. 验证匹配的细胞数量（确认有可筛选的细胞）
cat("匹配的细胞总数：", sum(health12$Donor.ID %in% c("Child2", "Adol2", "Child1", "Adol1")), "\n")

# 2. 直接按Donor.ID筛选（无需手动提取细胞名，Seurat自动匹配）
health4 <- subset(health12, subset = Donor.ID %in% c("Child2", "Adol2", "Child1", "Adol1"))

# 3. 验证结果
cat("筛选后细胞数：", ncol(health4), "\n")
cat("筛选后Donor.ID分布：\n")
print(table(health4$Donor.ID))
# Adol1  Adol2 Child1 Child2 
# 3079   4972   5266   3572 

# 保存筛选后的对象
saveRDS(health4, "/Users/lijiangbo/1_Projects/Epi1106/5_scRNA/2_processed_data/1113/healthy4_seurat.rds")


# 加载包
pacman::p_load(Seurat, dplyr)

# 读取原始健康4样本（未经过初步过滤的对象，或从原对象提取counts）
health4 <- readRDS("/Users/lijiangbo/1_Projects/Epi1106/5_scRNA/2_processed_data/1113/healthy4_seurat.rds")

# 1. 提取原始表达矩阵（counts）
counts <- GetAssayData(health4, assay = "RNA", slot = "counts")  # 获取原始计数矩阵

# 2. 重新创建Seurat对象，应用min.cells和min.features（与癫痫样本一致）
health4_initial <- CreateSeuratObject(
  counts = counts,
  project = "Healthy4",
  min.cells = 3,  # 保留至少在3个细胞中表达的基因
  min.features = 200  # 保留至少表达200个基因的细胞
)

# 3. 恢复原元数据（Donor.ID、Cell.type等关键信息）
# 仅保留与初步过滤后细胞匹配的元数据
health4_initial@meta.data <- health4@meta.data[colnames(health4_initial), ]

# 4. 后续步骤：计算线粒体比例（与癫痫流程一致）
health4_initial <- PercentageFeatureSet(health4_initial, pattern = "^MT-", col.name = "percent.mt")

# 5. 按癫痫标准筛选高质量细胞（完整质控）
high_quality <- health4_initial$nFeature_RNA > 200 & 
  health4_initial$nFeature_RNA < 5000 & 
  health4_initial$percent.mt < 10
health4_hq <- subset(health4_initial, cells = rownames(health4_initial@meta.data[high_quality, ]))

# 验证：输出初步过滤和质控后的统计（对比癫痫样本）
cat("健康4初步过滤后（应用min参数）- 细胞数：", ncol(health4_initial), "\n")
cat("健康4初步过滤后（应用min参数）- 基因数：", nrow(health4_initial), "\n")
cat("健康4最终质控后 - 细胞数：", ncol(health4_hq), "\n")
cat("健康4最终质控后 - 基因数：", nrow(health4_hq), "\n")

# 保存处理后对象
saveRDS(health4_hq, "/Users/lijiangbo/1_Projects/Epi1106/5_scRNA/2_processed_data/1113/healthy4_hq_final.rds")

head(health4)
