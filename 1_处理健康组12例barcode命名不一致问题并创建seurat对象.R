# 加载必要的包
library(Seurat)
library(dplyr)
library(tidyr)
library(stringr)
library(SeuratDisk)

health_matrix <- readRDS("/Users/lijiangbo/1_Projects/Epi1106/5_scRNA/1_raw_data/GSE204683_healthy/GSE204683_count_matrix.RDS")
health_barcodes <- read.table("/Users/lijiangbo/1_Projects/Epi1106/5_scRNA/1_raw_data/GSE204683_healthy/GSE204683_barcodes.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
head(health_barcodes)
head(health_matrix)

# 表达矩阵列名barcode前有序号，与barcode文件中的barcode不一致，前缀出现次数与Donor.ID 细胞数一直
matrix_cell_ids <- colnames(health_matrix)
extracted_prefixes <- str_extract(matrix_cell_ids, "^[0-9]+")

# 统计前缀数字的出现次数
prefix_counts <- table(extracted_prefixes)
# 过滤掉前缀数字为空的情况
prefix_counts <- prefix_counts[names(prefix_counts) != ""]
# 按出现次数降序排序
prefix_counts <- sort(prefix_counts, decreasing = TRUE)

print("矩阵中的前缀数字统计:")
print(prefix_counts)

# 检查列名并使用正确的列名
donor_col <- if("Donor.ID" %in% colnames(health_barcodes)) "Donor.ID" else if("Donor ID" %in% colnames(health_barcodes)) "Donor ID" else stop("错误：找不到Donor.ID或Donor ID列")
donor_counts <- table(health_barcodes[[donor_col]])
donor_counts <- sort(donor_counts, decreasing = TRUE)
print(donor_counts)


# 用前缀关联矩阵和元数据 -----------------------------------------------------------
# 1. 去除表达矩阵列名前的序号（保留纯Barcode）
colnames(health_matrix) <- sub("^[0-9]+_", "", colnames(health_matrix))

head(health_matrix)
cat("元数据行名与Barcode列是否一致：", all(rownames(health_barcodes) == health_barcodes$Barcode), "\n")
cat("矩阵重复Barcode数：", sum(duplicated(colnames(health_matrix))), "\n")
cat("元数据Barcode列重复数：", sum(duplicated(health_barcodes$Barcode)), "\n")
cat("元数据行名重复数：", sum(duplicated(rownames(health_barcodes))), "\n")

# 2. 元数据行名设为Barcode（Seurat匹配要求）
health_barcodes <- health_barcodes %>%
  group_by(Barcode) %>%  # 按Barcode分组
  slice_max(UMI.counts, n = 1) %>%  # 每组保留UMI最高的1个细胞
  ungroup()  # 取消分组

# 2. 再对齐矩阵与元数据（此时元数据Barcode已无重复）
common_bc <- intersect(colnames(health_matrix), health_barcodes$Barcode)  # 取共同Barcode
health_matrix <- health_matrix[, common_bc]  # 筛选矩阵列
health_barcodes <- health_barcodes[health_barcodes$Barcode %in% common_bc, ]  # 筛选元数据行（用Barcode列匹配，避免行名问题）
head(health_barcodes)
health_barcodes <- as.data.frame(health_barcodes)

# 4. 创建Seurat对象
seurat_obj <- CreateSeuratObject(counts = health_matrix, meta.data = health_barcodes, project = "Healthy")

# 5. 保存Seurat对象（.RDS）
saveRDS(seurat_obj, file = "/Users/lijiangbo/1_Projects/Epi1106/5_scRNA/2_processed_data/1113/healthy12_seurat.rds")


str(seurat_obj) # Seurat 对象-Assay(RNA)-槽（Slot）
head(seurat_obj@meta.data)
# assays：这是一个包含一个或多个 Assay 对象的列表。在你的情况下，只有一个名称为 "RNA" 的 Assay，其中包含了表达数据、元数据等。
# 
# RNA：
# layers：包含多个层次的表达数据（如 counts, data 等）。
# cells：表示每个细胞的元信息。
# features：包含基因或其他特征的信息。
# meta.data：这是一个数据框，包含与细胞相关的元信息。在你的对象中，有 44222 个细胞和 10 个变量（如 orig.ident, nCount_RNA, Cell.type 等）。
# 
# active.assay：指示 seurat_obj 当前活跃的 Assay。你的对象中当前活跃的 Assay 是 "RNA"。
# 
# active.ident：指示活跃细胞身份的因子变量，通常用于群体分析。
# 
# graphs, neighbors, reductions, images：这些槽通常用于存储细胞之间的关系图、降维结果（如 UMAP 或 PCA）以及图像数据。
# 
# project.name：存储项目的名称（在你的例子中为 "Healthy"）。
# 
# commands, tools, misc：这些槽通常包含分析过程中使用的工具或额外信息。
# Seurat 对象和其组件的层级结构示意图：
# Seurat对象
# │
# ├── Assays (List)
# │   └── RNA (Assay)
# │       ├── counts (dgCMatrix)
# │       ├── data (dgCMatrix)
# │       ├── scale.data (dgCMatrix)
# │       ├── meta.data (DataFrame)
# │       ├── layers (List)
# │       └── misc (List)
# │
# ├── meta.data (DataFrame)
# ├── active.assay (字符)
# ├── active.ident (Factor)
# ├── project.name (字符)
# ├── graphs (List)
# ├── neighbors (List)
# ├── reductions (List)
# ├── images (List)
# ├── commands (List)
# └── tools (List)


# 获取原始计数数据，使用 layer 参数
counts_data <- GetAssayData(seurat_obj[["RNA"]], layer = "counts")

# 获取元数据
meta_data <- seurat_obj@meta.data

# 将计数数据保存为 CSV 文件
write.csv(as.data.frame(counts_data), file = "/Users/lijiangbo/1_Projects/Epi1106/5_scRNA/2_processed_data/1113/healthy12_rna_counts.csv", row.names = TRUE)



