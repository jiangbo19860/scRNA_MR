rm(list = ls())

# 单细胞RNA测序（scRNA-seq）数据预处理、聚类分析及可视化完整流程
# 核心目标：通过对单细胞基因表达数据的系统处理，揭示细胞异质性并可视化细胞聚类特征
# 主要流程包括：
# 1. 数据加载与整合
#   - 读取基因表达矩阵（log2cpm，对数转换的每百万计数）
#   - 整合基因注释信息（featuredata，包含基因ID、名称等元数据）
#   - 导入细胞层面元数据（如细胞来源CellFromTumor、患者编号、样本信息等）
#   - 加载预计算的t-SNE降维坐标，最终整合成Seurat对象用于统一分析
#
# 2. 质量控制与预处理
#   - 数据归一化：采用LogNormalize方法，缩放因子设为10000
#   - 高变基因筛选：使用VST（方差稳定变换）方法，筛选2000个高变基因
#   - 数据中心化：对高变基因进行标准化（减均值、除标准差），消除技术偏差
#
# 3. 降维与聚类分析
#   - 主成分分析（PCA）：基于高变基因提取主要变异信息
#   - 聚类分析：使用前15个主成分进行细胞邻居查找，以resolution=0.5进行聚类
#   - 降维可视化：通过t-SNE和UMAP算法将高维数据映射到二维空间，展示细胞聚类分布
#
# 4. 特征基因分析与可视化
#   - 高变基因展示：绘制前10个高变基因的散点图和分聚类热图
#   - 标记基因鉴定：通过差异表达分析（FindAllMarkers）获取每个聚类的特征标记基因
#   - 聚类特征可视化：结合PCA、t-SNE、UMAP图展示细胞异质性及标记基因表达模式

# # 安装 remotes 包（如果尚未安装）
# if (!requireNamespace("remotes", quietly = TRUE)) {
#   install.packages("remotes")
# }
#
# # 从 GitHub 安装 SCopeLoomR
# remotes::install_github("aertslab/SCopeLoomR")
pacman::p_load(
  here, Seurat, SCopeLoomR, loomR, dplyr, ggplot2, cowplot, patchwork, ggrepel, SingleR, celldex, BiocParallel, scater, scran, SingleCellExperiment, data.table, hdf5r, Matrix, tidyr, stringr, forcats, pheatmap, RColorBrewer, viridis, ggridges, tibble
)
here()


# 0. 定义路径，加载数据 ------------------------------------------------------------
rds_file <- here("1_data", "E-MTAB-6149", "Allsamples.Cellview.Rds")
loom_path <- here("1_data", "E-MTAB-6149", "Thienpont_Tumors_52k_v4_R_fixed.loom")

# 定义输出路径函数，所有输出文件保存在3_outputs目录下的E-MTAB-6149_outputs子目录中。
output_path <- function(...) {
  here("3_outputs", "E-MTAB-6149_outputs", ...)
}

# 创建输出目录（如果不存在）
if (!dir.exists(output_path())) {
  dir.create(output_path(), recursive = TRUE)
  message("已创建输出目录:", output_path())
}

# 读取Rds文件
load(rds_file) # 读取Rds文件，加载数据到当前工作空间
ls() # 显示当前工作空间中的所有对象："featuredata" "log2cpm"     "loom_path"   "output_path" "rds_file"    "tsne.data"。 featuredata存储基因 / 特征的元数据，如基因名称、染色体位置、生物类型等，用于注释表达矩阵中的行。log2cpm存储基因表达矩阵（行 = 基因，列 = 细胞），通常是对数转换后的CPM值（Counts Per Million），即每百万reads的计数。log2cpm_sparse是log2cpm的稀疏矩阵版本，适用于大规模单细胞数据分析。loom_path是loom文件的路径，包含单细胞数据的元数据和降维坐标。meta存储细胞元数据，如细胞类型、批次信息等。nsclc存储从loom文件中读取的单细胞数据对象。nsclc.meta存储从loom文件中提取的细胞元数据。nsclc.tsne存储从loom文件中提取的t-SNE坐标。output_path是一个函数，用于生成输出文件的路径。rds_file是Rds文件的路径，包含单细胞数据。scRNA是一个Seurat对象，包含处理后的单细胞数据。temp_env是一个临时环境，用于加载Rds文件而不污染全局变量。

exists("tsne.data") # 检查对象tsne.data是否存在，tsne.data存储了 t-SNE 降维的坐标和细胞类型注释。

class(tsne.data) # "data.frame"，tsne.data是一个数据框，是一个有4列的DataFrame，包含t-SNE坐标（V3都是NA）和细胞类型注释。
str(tsne.data) # 'data.frame':	52698 obs. of  4 variables: V1, V2, V3, dbCluster（细胞类型注释）
unique(tsne.data$dbCluster) # 查看细胞类型注释的唯一值，表示不同的细胞类型或聚类结果, 有8个：Alveolar, B_cell, EC, Epi, Fibro, Myeloid, T_cell, tumor。
unique(tsne.data$V3) # "NA"。


# 1. 对齐基因表达矩阵log2cpm和基因注释信息featuredata，确保两者的基因标识符一致。 -------------------------------------------
# 基因表达矩阵log2cpm（行 = 基因，列 = 细胞），行名通常是基因 ID（如 Ensembl ID）。
# featuredata：基因注释信息（行 = 基因，列 = 各种注释字段），行名也是基因 ID，且包含更多注释（如基因名、染色体位置）。
colnames(featuredata)
featuredata[1:3, 1:6] # 基因注释信息的前3行和前6列
log2cpm[1:3, 1:6] # 基因表达矩阵（通常是对数转换后的CPM值（Counts Per Million），即每百万reads的计数）的前3行和前6列
featuredata <- featuredata[rownames(log2cpm), ] # 按log2cpm的行名筛选featuredata

rownames(log2cpm) <- featuredata[, 5] # 将log2cpm的行名设置为featuredata的第5列（通常是基因名或ID）
log2cpm[1:3, 1:6] # 查看更新后的基因表达矩阵的前3行和前6列

library(Matrix)

# 转换为稀疏矩阵
log2cpm_sparse <- as(as.matrix(log2cpm), "dgCMatrix")

# 用稀疏矩阵创建 Seurat 对象
scRNA <- CreateSeuratObject(
  counts = log2cpm_sparse,
  project = "scRNA"
)
scRNA # 查看创建的 Seurat 对象信息：An object of class Seurat, 22180 features (基因) across 52698 samples(单细胞) within 1 assay (分析模块), 1个assay，说明当前对象只存储了一种类型的数据（默认是 scRNA-seq 的基因表达数据）。Active assay: RNA (22180 features, 0 variable features) , “Active assay” 指当前正在使用的 assay, 该assay包含22180个基因。“RNA” 是 assay 的名称，是 CreateSeuratObject 函数的默认命名，对应 scRNA-seq 的基因表达数据。若后续添加其他 assay（如 “ATAC”），可通过 DefaultAssay(scRNA) <- "ATAC" 切换活跃 assay。0 variable features：表示当前未标记任何 “可变基因”。可变基因（variable features）是指在不同细胞中表达差异较大的基因，是后续降维（如PCA)和聚类分析的核心依据（因为它们更能反映细胞间的异质性）。执行 FindVariableFeatures 函数（筛选可变基因）后这里会显示筛选出的可变基因数量（如 “2000 variable features”）。 1 layer present: counts. “layer” 指 assay 中存储的数据层，即同一类数据的不同处理状态（如原始计数、归一化后的数据等）。这里显示有1个数据层：counts，表示该 assay 目前只存储了原始表达矩阵（即 log2cpm_sparse 传入的数据，虽然变量名含 “log2cpm”，但 Seurat 会将其作为 “原始计数” 层存储）。
summary(scRNA@meta.data) # 查看 Seurat 对象的元数据摘要：nFeature_RNA（每个细胞检测到的基因数），nCount_RNA（每个细胞的UMI总数），percent.mt（线粒体基因比例）等。

str(scRNA) # 查看 Seurat 对象的结构：Formal class 'Seurat' [package "SeuratObject"] with 13 slots


# Idents(scRNA)  # 返回当前聚类标识（cluster labels），如果没有进行聚类相关的分析步骤（如 FindClusters 函数），所有细胞暂时被标记为同一个 “聚类”（标签为 "scRNA"）；
table(Idents(scRNA)) # 统计每个聚类的细胞数

# 查看替换前的基因名（取前5个）
head(rownames(log2cpm_sparse), 5)

# 查看替换后的基因名（Seurat对象中的基因名）
head(rownames(scRNA), 5)

# 计算线粒体基因百分比（假设是人类数据，线粒体基因以"MT-"开头）
scRNA <- PercentageFeatureSet(scRNA, pattern = "^MT-", col.name = "percent.mt")

VlnPlot(scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0, ncol = 3) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("QC Metrics Violin Plot") +
  xlab("Metrics") +
  ylab("Value")

### loom文件读取，loom文件是一种 HDF5 结构的单细胞数据存储格式，通常包含基因表达矩阵、细胞元数据（如细胞类型、聚类结果）和降维坐标（如 t-SNE/UMAP）。
nsclc <- connect(filename = loom_path, mode = "r") # mode = "r"表示以只读模式打开loom文件
nsclc # loom文件, 包含单细胞数据的元数据和降维坐标。3个顶层属性Attributes: title, MetaData, LOOM_SPEC_VERSION

nsclc.meta <- get_cell_annotation(nsclc) # 细胞元数据如细胞类型标签、批次信息、实验条件等。
str(nsclc.meta)
class(nsclc.meta) # "data.frame"
colnames(nsclc.meta) # 查看元数据列名:"CellFromTumor" "ClusterName"   "PatientNumber" "Sample"        "nGene"         "nUMI"

nsclc.tsne <- get_embeddings(nsclc) # 提取降维后的坐标（如 t-SNE 或 UMAP），用于可视化细胞聚类结果。
class(nsclc.tsne) # list，是无序的 “集合”，本质是 “元素的集合”，不具备二维表格的 “行 / 列” 概念。元素通过索引位置（如 [[1]]）或名称（如 $name）访问。
str(nsclc.tsne) # list 的整体结构
names(nsclc.tsne) # 元素名称："Seurat t-SNE", "SCENIC 50PC, 50perp", SCENIC（Single-Cell Regulatory Network Inference and Clustering） 算法计算的 t-SNE 降维结果。SCENIC 是一种专门用于推断单细胞中转录调控网络的工具，它分析的不是原始表达量，而是转录因子（TF）对靶基因的调控活性。50PC：使用了 50 个主成分进行降维（PCA 的结果）。50perp：t-SNE 的 perplexity 参数设为 50（控制局部与全局结构的平衡，较大的值强调全局结构）。

save(nsclc.meta, file = output_path("meta.Rds")) # 保存元数据
save(nsclc.tsne, file = output_path("tsne.Rds")) # 保存t-SNE坐标

### 非常重要！！！把nsclc.meta中的细胞来源细胞类型病人编号样本编号基因数UMI数等信息与tsne.data中的t-SNE坐标和细胞类型注释合并，形成一个新的元数据框meta，包括完整的元数据框：实验设计信息（来自nsclc.meta）、可视化坐标（来自tsne.data）、细胞类型标签（来自tsne.data）。
meta <- cbind(nsclc.meta, tsne.data) # nsclc.meta是52698 obs. of 6 variables的data.frame, 有1细胞来源，2细胞类型，3病人编号，4样本编号，5基因数，6UMI数。tsne.data是一个data.frame，有4个变量：V1, V2, V3（t-SNE坐标），dbCluster（细胞类型注释）。cbind函数将这两个数据框按列合并，要求行数必须相同。

colnames(meta) # 查看合并后的元数据列名:"CellFromTumor", "ClusterName", "PatientNumber" "Sample", "nGene"         "nUMI", "V1", "V2", "V3", "dbCluster".
head(meta) # 查看合并后的元数据前几行
str(meta)
unique(meta$CellFromTumor) # "CellFromTumor"列：细胞来源（如肿瘤或正常组织）0和1。
unique(meta$ClusterName) # "ClusterName"列：细胞类型或聚类结果，有40种不同类型的细胞。
unique(meta$PatientNumber) # "PatientNumber"列的唯一值，表示患者编号，1-5有5个。
length(unique(meta$Sample)) # "Sample"列的唯一值，表示样本编号，有19个。


scRNA <- AddMetaData(scRNA, metadata = meta) # 将元数据添加到 Seurat 对象中，便于后续分析和可视化。
summary(scRNA@meta.data)
colnames(scRNA@meta.data)
scRNA@meta.data[1:5, 1:10]
dim(scRNA@assays$RNA) # 查看基因表达矩阵的维度：22180 features (基因) across 52698 samples(单细胞)。

# 2. 对 Seurat 对象进行预处理和质量控制 -------------------------------------------
# 2.1 对数据进行归一化处理，使用对数归一化方法。
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)

# 2.2 选择高变基因FindVariableFeatures
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000) # 选择2000个高变基因，使用方差稳定变换（vst）方法。

# 2.3 提取前10个高变基因
top10 <- head(VariableFeatures(scRNA), 10)
print(top10)

# 2.4 数据中心化与回归
scRNA <- ScaleData(scRNA, features = VariableFeatures(scRNA)) # 对高变基因进行中心化处理，减去均值并除以标准差。

# 2.5 运行主成分分析（PCA），提取主要的变异信息。
scRNA <- RunPCA(scRNA, features = VariableFeatures(scRNA))

# 2.6 可视化PCA结果
DimPlot(scRNA, reduction = "pca", group.by = "orig.ident") +
  ggtitle("PCA Plot") +
  xlab("PC1") +
  ylab("PC2")

# 2.7 绘制肘部图，帮助选择合适的主成分数量。
ElbowPlot(scRNA, ndims = 20, reduction = "pca") +
  ggtitle("Elbow Plot") +
  xlab("Principal Components") +
  ylab("Variance Explained")
# 选择前15个主成分
pc.num <- 1:15

# 2.8 使用选定的主成分进行邻居查找
scRNA <- FindNeighbors(scRNA, dims = pc.num)

# 2.9 使用选定的主成分进行聚类分析
scRNA <- FindClusters(scRNA, resolution = 0.5)
# resolution参数控制聚类的粒度，0.5是一个常用的值，可以根据数据集的特点进行调整。

# 2.10 查看每个聚类的细胞数量
table(scRNA@meta.data$seurat_clusters)

# 2.11 保存细胞的注释信息：将Seurat对象中的元数据（metadata，包括各种注释信息如聚类标签、样本来源、细胞类型等）提取出来，保存到一个独立的数据框（metadata）中
metadata <- scRNA@meta.data

# 2.12 两种降维：RunTSNE() 和 RunUMAP()。在 Seurat 中，RunTSNE() 和 RunUMAP() 两个降维函数不会相互覆盖或冲突，因为它们会将结果存储在 Seurat 对象的不同位置（reductions 槽中），并使用不同的名称标识。
scRNA <- RunTSNE(scRNA, dims = pc.num) # 运行t-SNE降维，RunTSNE() 的结果会以 tsne 为名称存储（可通过 reduction.name 参数自定义，默认是 "tsne"）。
scRNA <- RunUMAP(scRNA, dims = pc.num) # 运行UMAP降维，RunUMAP() 的结果会以 umap 为名称存储（默认名称是 "umap"）。

# 查看所有降维方法的名称
names(scRNA@reductions)

scRNA

# 查看元数据列（确认聚类标签、样本来源等是否存在）
colnames(scRNA@meta.data)

# 保存处理后的 Seurat 对象
save(scRNA, file = output_path("scRNA_processed.Rds")) # 把经过预处理后的scRNA对象保存到指定路径的scRNA_processed.Rds 文件中，方便后续复用（无需重新运行数据预处理步骤）。


# 3. 可视化 ------------------------------------------------------------------
# 1. PCA 聚类图（按聚类标签着色）
DimPlot(scRNA, reduction = "pca", label = TRUE) + ggtitle("PCA Clusters")

# 2. t-SNE 聚类图
DimPlot(scRNA, reduction = "tsne", label = TRUE) + ggtitle("t-SNE Clusters")

# 3. UMAP 聚类图（推荐，效果通常更好）
DimPlot(scRNA, reduction = "umap", label = TRUE) + ggtitle("UMAP Clusters")

# 4. 按其他元数据分组（如样本来源、病人编号）
# DimPlot(scRNA, reduction = "umap", group.by = "orig.ident") + ggtitle("By Sample")

### 查看高变基因和标记基因
# 1. 高变基因散点图（显示平均表达量 vs. 变异程度）
# 可视化高变异基因
# 可视化高变异基因
top10 <- head(VariableFeatures(scRNA), 10)
plot1 <- VariableFeaturePlot(scRNA)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
plot2

top10 <- head(VariableFeatures(scRNA), 10)

# 1.1 高变基因热图
# 创建一个存储每个聚类高变基因的列表
cluster_hvgs <- list()

# 对每个聚类分别计算高变基因
for (cluster in unique(Idents(scRNA))) {
  # 提取当前聚类的细胞
  cluster_cells <- WhichCells(scRNA, idents = cluster)

  # 对当前聚类的细胞计算高变基因
  cluster_subset <- subset(scRNA, cells = cluster_cells)
  cluster_subset <- FindVariableFeatures(cluster_subset, nfeatures = 2000)

  # 存储当前聚类的高变基因
  cluster_hvgs[[as.character(cluster)]] <- VariableFeatures(cluster_subset)
}

# 为每个聚类绘制高变基因热图
for (cluster in names(cluster_hvgs)) {
  top10 <- head(cluster_hvgs[[cluster]], 10) # 获取当前聚类的前10个高变基因

  # 绘制热图（只显示当前聚类的细胞）
  p <- DoHeatmap(
    subset(scRNA, idents = cluster),
    features = top10
  ) + ggtitle(paste("Cluster", cluster, "Top 10 HVGs"))

  print(p)
}

# 2. 查看每个聚类的标记基因（差异表达基因）
markers <- FindAllMarkers(scRNA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top_markers <- markers %>%
  group_by(cluster) %>%
  top_n(n = 2, wt = avg_log2FC)
print(top_markers)

# 3. 热图可视化标记基因表达
DoHeatmap(scRNA, features = top_markers$gene) + NoLegend()



# Fig1_t-SNE降维图，按CellFromTumor着色 -----------------------------------------------------------
unique(scRNA$CellFromTumor)
DimPlot(scRNA, reduction = "tsne", group.by = "CellFromTumor") +
  # 根据unique(scRNA$CellFromTumor)返回的实际值 "1" 和 "0" 调整颜色映射
  scale_color_manual(
    values = c("1" = "#2D5474", "0" = "#72C667"), # "1" 对应肿瘤，"0" 对应非肿瘤
    labels = c("Tumor", "Non-malignant"), # 自定义图例标签（顺序与 values 一致）
    limits = c("1", "0") # 控制图例顺序（Tumor 在上）
  ) +
  # 将图例形状改为方块（通过guides设置）
  guides(color = guide_legend(override.aes = list(shape = 15, size = 5))) + # 15=方块，16=圆，17=三角.
  theme(
    panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid"),
    legend.position = c(.01, .05), # 图例位置（左下角），(0,0)为左下角，(1,1)为右上角
    legend.key = element_rect(fill = NA) # 去除图例背景
  ) +
  labs(title = "Sample Origin")

ggsave(filename = output_path("Fig1_tSNE_CellFromTumor.pdf"), height = 8, width = 8)
